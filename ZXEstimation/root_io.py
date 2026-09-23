"""Stream scalar skim observables with PyROOT TTreeFormula, file by file."""
import math
import uuid
from collections import Counter


def formula(ROOT, tree, expression):
    result = ROOT.TTreeFormula("expr_" + uuid.uuid4().hex, str(expression), tree)
    if result.GetNdim() <= 0:
        raise ValueError(f"Invalid expression or missing branch: {expression}")
    return result


def scalar(compiled):
    if compiled.GetNdata() != 1:
        raise ValueError(f"Expected scalar expression: {compiled.GetTitle()}")
    value = float(compiled.EvalInstance())
    if not math.isfinite(value):
        raise ValueError(f"Nonfinite value in {compiled.GetTitle()}")
    return value


class EventReader:
    def __init__(self, ROOT, cfg):
        self.ROOT, self.cfg = ROOT, cfg
        self.audit = Counter()

    def events(self, regions, sample_filter=None, variables=()):
        cfg, ROOT = self.cfg, self.ROOT
        seen = set()
        selected = {name: cfg["regions"][name] for name in regions}
        for sample in cfg["samples"]:
            if sample_filter and not sample_filter(sample):
                continue
            for filename in sample["_files"]:
                print(f"Reading {sample['name']}: {filename}", flush=True)
                handle = ROOT.TFile.Open(filename, "READ")
                if not handle or handle.IsZombie():
                    raise OSError(f"Cannot open ROOT file {filename}")
                cuts, channels, probes, observables = {}, {}, [], {}
                common = sample_cut = jets = weight = None
                try:
                    tree = handle.Get(sample.get("tree", cfg["tree"]))
                    if not tree or not tree.InheritsFrom("TTree"):
                        raise ValueError(f"Missing TTree in {filename}")
                    branches = {b.GetName() for b in tree.GetListOfBranches()}
                    for spec in selected.values():
                        if spec["branch"] not in branches:
                            raise ValueError(f"{filename}: missing {spec['branch']}. For 3P0F, regenerate the skim or use an external fake-rate map.")
                    cuts = {name: formula(ROOT, tree, spec.get("cut", "1")) for name, spec in selected.items()}
                    common = formula(ROOT, tree, cfg["selection"]["cut"])
                    sample_cut = formula(ROOT, tree, sample.get("selection", "1"))
                    jets = formula(ROOT, tree, cfg["selection"]["jet_count"])
                    weight = formula(ROOT, tree, sample.get("weight", "1") if sample["kind"] == "mc" else "1")
                    channels = {name: formula(ROOT, tree, expr) for name, expr in cfg["channels"].items()}
                    probes = [{k: formula(ROOT, tree, expr) for k, expr in probe.items()} for probe in cfg["probes"]]
                    observables = {name: formula(ROOT, tree, cfg["variables"][name]["expression"]) for name in variables}
                    is_data = sample["kind"] == "data"
                    identifiers = cfg["selection"]["event_id"]
                    if is_data and cfg["selection"]["deduplicate_data"]:
                        if set(identifiers) - branches:
                            raise ValueError(f"{filename}: missing event identity branches {identifiers}")
                    for entry in range(tree.GetEntries()):
                        if tree.GetEntry(entry) < 0:
                            raise OSError(f"Read error at {filename}:{entry}")
                        self.audit["entries_read"] += 1
                        active = [name for name, spec in selected.items() if bool(getattr(tree, spec["branch"]))]
                        if len(active) > 1:
                            raise ValueError(f"Overlapping CR flags at {filename}:{entry}: {active}")
                        if not active:
                            continue
                        region = active[0]
                        if not scalar(common) or not scalar(sample_cut) or not scalar(cuts[region]):
                            continue
                        if scalar(jets) < int(cfg["selection"]["min_jets"]):
                            continue
                        if is_data and cfg["selection"]["deduplicate_data"]:
                            # Read identifiers as integers, never through double-valued formulas.
                            key = (sample["era"], *(int(getattr(tree, name)) for name in identifiers))
                            if key in seen:
                                self.audit["duplicate_data"] += 1
                                continue
                            seen.add(key)
                        role = selected[region]["role"]
                        count = 1 if role.startswith("rate_") else 2
                        legs = [{k: scalar(expr) for k, expr in probe.items()} for probe in probes[:count]]
                        n_fail = sum(not bool(p["tight"]) for p in legs)
                        expected = {"rate_fail": 1, "rate_pass": 0, "single_fail": 1, "double_fail": 2}[role]
                        if n_fail != expected or any(abs(int(p["pdgid"])) not in (11, 13) for p in legs):
                            raise ValueError(f"Inconsistent probe metadata at {filename}:{entry} ({region})")
                        value_weight = scalar(weight) * sample["_scale"]
                        if not math.isfinite(value_weight):
                            raise ValueError(f"Nonfinite event weight at {filename}:{entry}")
                        self.audit[f"selected/{sample['name']}/{region}"] += 1
                        yield {"sample": sample, "region": region, "role": role,
                               "weight": value_weight, "probes": legs,
                               "channels": [name for name, expr in channels.items() if scalar(expr)],
                               "values": {name: scalar(expr) for name, expr in observables.items()
                                          if region in cfg["variables"][name]["regions"]}}
                finally:
                    # Destroy formulas while the file-owned tree is still alive.
                    cuts.clear()
                    channels.clear()
                    probes.clear()
                    observables.clear()
                    common = sample_cut = jets = weight = None
                    handle.Close()


def mkdir(root_file, path):
    directory = root_file
    for part in path.split("/"):
        directory = directory.GetDirectory(part) or directory.mkdir(part)
    return directory


def write_object(root_file, path, obj, name):
    mkdir(root_file, path).WriteObject(obj, name)

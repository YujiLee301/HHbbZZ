"""Prompt-subtracted T/(T+F) rates, separately for each era and flavor."""
import math
import os

from configuration import local_path
from histograms import h2, clone
from root_io import EventReader, write_object


class FakeRates:
    def __init__(self, cfg):
        self.cfg = cfg
        self.maps = {}
        self.errors = {}

    def locate(self, era, probe):
        flavor = abs(int(probe["pdgid"]))
        hist, valid = self.maps[(era, flavor)]
        pt, eta = probe["pt"], abs(probe["eta"])
        ix, iy = hist.GetXaxis().FindFixBin(pt), hist.GetYaxis().FindFixBin(eta)
        if ix > hist.GetNbinsX() and self.cfg["fake_rate"]["pt_overflow"] == "last_bin":
            ix = hist.GetNbinsX()
        if not (1 <= ix <= hist.GetNbinsX() and 1 <= iy <= hist.GetNbinsY()):
            raise ValueError(f"Probe outside FR acceptance: era={era}, flavor={flavor}, pt={pt}, abs_eta={eta}")
        if valid.GetBinContent(ix, iy) != 1:
            raise ValueError(f"Invalid/unpopulated FR bin {era}/{flavor}/{ix}/{iy}; merge bins in YAML or supply a validated external rate map")
        return hist.GetBinContent(ix, iy), (era, flavor, ix, iy)

    def register(self, era, flavor, hist, valid):
        spec = self.cfg["fake_rate"]["binning"][flavor]
        for axis, expected in ((hist.GetXaxis(), spec["pt"]), (hist.GetYaxis(), spec["abs_eta"]),
                               (valid.GetXaxis(), spec["pt"]), (valid.GetYaxis(), spec["abs_eta"])):
            actual = [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 2)]
            if len(actual) != len(expected) or any(abs(a - b) > 1e-9 for a, b in zip(actual, expected)):
                raise ValueError(f"Rate map axes differ from YAML for {era}/{flavor}")
        for ix in range(1, hist.GetNbinsX() + 1):
            for iy in range(1, hist.GetNbinsY() + 1):
                f, error = hist.GetBinContent(ix, iy), hist.GetBinError(ix, iy)
                if valid.GetBinContent(ix, iy) == 1:
                    if not math.isfinite(f) or not 0 <= f < 1 or not math.isfinite(error) or error < 0:
                        raise ValueError(f"Invalid fake rate/error in {era}/{flavor}/{ix}/{iy}")
                self.errors[(era, flavor, ix, iy)] = error
        self.maps[(era, flavor)] = (hist, valid)


def measure(ROOT, cfg, output):
    regions = [k for k, v in cfg["regions"].items() if v["role"].startswith("rate_")]
    if {cfg["regions"][r]["role"] for r in regions} != {"rate_pass", "rate_fail"}:
        raise ValueError("Measuring fake rates needs both rate_pass (3P0F) and rate_fail (2P1F)")
    subtract = cfg["fake_rate"]["subtract_processes"]
    counts = {}
    for era in cfg["eras"]:
        for flavor, bins in cfg["fake_rate"]["binning"].items():
            for process in ["data"] + subtract:
                for state in ("pass", "fail"):
                    counts[era, flavor, process, state] = h2(ROOT, bins["pt"], bins["abs_eta"])
    reader = EventReader(ROOT, cfg)
    def needed(sample):
        return sample["kind"] == "data" or sample["process"] in subtract
    for event in reader.events(regions, needed):
        sample, probe = event["sample"], event["probes"][0]
        era, flavor = sample["era"], abs(int(probe["pdgid"]))
        process = "data" if sample["kind"] == "data" else sample["process"]
        state = "pass" if event["role"] == "rate_pass" else "fail"
        hist = counts[era, flavor, process, state]
        pt, eta = probe["pt"], abs(probe["eta"])
        if pt >= hist.GetXaxis().GetXmax() and cfg["fake_rate"]["pt_overflow"] == "last_bin":
            pt = hist.GetXaxis().GetBinCenter(hist.GetNbinsX())
        if not (hist.GetXaxis().GetXmin() <= pt < hist.GetXaxis().GetXmax() and
                hist.GetYaxis().GetXmin() <= eta < hist.GetYaxis().GetXmax()):
            raise ValueError(f"Rate measurement probe outside configured bins: {era}/{flavor}: {pt}, {eta}")
        hist.Fill(pt, eta, event["weight"])
    rates, audit = FakeRates(cfg), {"input": dict(reader.audit), "bins": []}
    handle = ROOT.TFile.Open(str(output), "RECREATE")
    if not handle or handle.IsZombie():
        raise OSError(f"Cannot create {output}")
    try:
        for era in cfg["eras"]:
            for flavor, bins in cfg["fake_rate"]["binning"].items():
                passed = clone(counts[era, flavor, "data", "pass"])
                failed = clone(counts[era, flavor, "data", "fail"])
                for process in subtract:
                    passed.Add(counts[era, flavor, process, "pass"], -1.)
                    failed.Add(counts[era, flavor, process, "fail"], -1.)
                rate, valid = h2(ROOT, bins["pt"], bins["abs_eta"]), h2(ROOT, bins["pt"], bins["abs_eta"])
                for ix in range(1, rate.GetNbinsX() + 1):
                    for iy in range(1, rate.GetNbinsY() + 1):
                        t, f = passed.GetBinContent(ix, iy), failed.GetBinContent(ix, iy)
                        vt, vf = passed.GetBinError(ix, iy) ** 2, failed.GetBinError(ix, iy) ** 2
                        good = t > 0 and f > 0 and t + f > 0
                        value = t / (t + f) if good else None
                        # T and F are disjoint samples. This accounts for
                        # numerator/denominator correlation in T/(T+F).
                        variance = (f*f*vt + t*t*vf) / (t+f)**4 if good else None
                        if good:
                            rate.SetBinContent(ix, iy, value)
                            rate.SetBinError(ix, iy, math.sqrt(variance))
                            valid.SetBinContent(ix, iy, 1)
                        audit["bins"].append({"era": era, "flavor": flavor, "pt_bin": ix, "eta_bin": iy,
                                              "tight_subtracted": t, "fail_subtracted": f,
                                              "rate": value, "variance": variance, "valid": good})
                path = f"{era}/{flavor}"
                for name, hist in (("rate", rate), ("valid", valid), ("tight_subtracted", passed), ("fail_subtracted", failed)):
                    write_object(handle, path, hist, name)
                for process in ["data"] + subtract:
                    for state in ("pass", "fail"):
                        write_object(handle, path, counts[era, flavor, process, state], f"{process}_{state}")
                rates.register(era, flavor, rate, valid)
        ROOT.TObjString("OS Z1+l; >=2 cleaned jets; f=T/(T+F); independent pass/fail sumw2").Write("definition")
    finally:
        handle.Close()
    return rates, audit


def load(ROOT, cfg, output):
    external = cfg["fake_rate"]["source"] == "external"
    options = cfg["fake_rate"].get("external", {}) if external else {}
    filename = output
    if external:
        configured = os.path.expandvars(os.path.expanduser(str(options["file"])))
        if "$" in configured:
            raise ValueError(f"Unresolved environment variable in {configured}")
        filename = configured if configured.startswith(("root://", "https://")) else local_path(configured, cfg["_base"])
    handle = ROOT.TFile.Open(str(filename), "READ")
    if not handle or handle.IsZombie():
        raise OSError(f"Cannot open fake-rate file {filename}; run --stage fake-rate first")
    rates = FakeRates(cfg)
    try:
        for era in cfg["eras"]:
            for flavor in cfg["fake_rate"]["binning"]:
                objects = []
                for kind in ("rate", "valid"):
                    pattern = options.get(kind, "{era}/{flavor}/" + kind)
                    name = pattern.format(era=era, flavor=flavor)
                    hist = handle.Get(name)
                    if not hist or not hist.InheritsFrom("TH2"):
                        raise ValueError(f"Missing TH2 {name} in {filename}")
                    objects.append(clone(hist))
                rates.register(era, flavor, *objects)
    finally:
        handle.Close()
    return rates

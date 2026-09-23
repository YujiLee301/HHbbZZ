"""OS inclusion-exclusion transfer and correlated FR statistical errors."""
import math
from collections import defaultdict

from histograms import h1, h2, clone, fill, total
from root_io import EventReader, write_object


def transfer(event, rates):
    failed = [p for p in event["probes"] if not bool(p["tight"])]
    values = [rates.locate(event["sample"]["era"], p) for p in failed]
    odds = [f / (1. - f) for f, _ in values]
    sign = 1. if event["role"] == "single_fail" else -1.
    weight = sign * math.prod(odds)
    derivatives = defaultdict(float)
    for i, (f, key) in enumerate(values):
        other = math.prod(r for j, r in enumerate(odds) if i != j)
        derivatives[key] += sign * other / (1. - f)**2
    return weight, derivatives


def estimate(ROOT, cfg, rates, output):
    regions = [k for k, v in cfg["regions"].items() if v["role"] in ("single_fail", "double_fail")]
    if {cfg["regions"][r]["role"] for r in regions} != {"single_fail", "double_fail"}:
        raise ValueError("OS estimation needs both single_fail and double_fail regions")
    variables = {name: spec for name, spec in cfg["variables"].items() if spec.get("estimate", False)}
    for name, spec in variables.items():
        if not set(regions) <= set(spec["regions"]):
            raise ValueError(f"Estimator variable {name} must be defined in both application regions")
    specs = dict(variables)
    specs["yield"] = {"edges": [0., 1.], "fold": True}
    templates, components, gradients = {}, {}, {}
    eras = list(cfg["eras"]) + ["combined"]
    if "combined" in cfg["eras"]:
        raise ValueError("Era name 'combined' is reserved for summed templates")
    for era in eras:
        for channel in cfg["channels"]:
            for variable, spec in specs.items():
                key = (era, channel, variable)
                templates[key] = h1(ROOT, spec["edges"])
                components[key] = {part: h1(ROOT, spec["edges"]) for part in
                                   ("data_single", "prompt_single", "data_double", "prompt_double")}
                gradients[key] = {}
    reader = EventReader(ROOT, cfg)
    processes = set(cfg["estimation"]["subtract_single"] + cfg["estimation"]["subtract_double"])
    def needed(sample):
        return sample["kind"] == "data" or sample["process"] in processes
    for event in reader.events(regions, needed, variables):
        sample = event["sample"]
        single = event["role"] == "single_fail"
        category = "single" if single else "double"
        is_data = sample["kind"] == "data"
        if not is_data and sample["process"] not in cfg["estimation"]["subtract_" + category]:
            continue
        base = event["weight"] * (1. if is_data else -1.)
        weight, derivative = transfer(event, rates)
        weight *= base
        values = dict(event["values"], yield_value=0.5)
        for era in (sample["era"], "combined"):
            for channel in event["channels"]:
                for variable, spec in specs.items():
                    value = values["yield_value"] if variable == "yield" else values[variable]
                    key = era, channel, variable
                    fill(templates[key], value, weight, spec.get("fold", True))
                    fill(components[key][("data_" if is_data else "prompt_") + category],
                         value, weight, spec.get("fold", True))
                    for nuisance, slope in derivative.items():
                        if nuisance not in gradients[key]:
                            gradients[key][nuisance] = h1(ROOT, spec["edges"])
                        fill(gradients[key][nuisance], value, base * slope, spec.get("fold", True))
    handle = ROOT.TFile.Open(str(output), "RECREATE")
    if not handle or handle.IsZombie():
        raise OSError(f"Cannot create {output}")
    audit = {"input": dict(reader.audit), "yields": {}, "negative_bins": [],
             "uncertainties": "application sumw2 plus first-order independent era/flavor/FR-bin errors; same-bin probes and events are correlated"}
    try:
        for key, nominal in templates.items():
            era, channel, variable = key
            path = "/".join(key)
            spec = specs[variable]
            combined = clone(nominal)
            covariance = h2(ROOT, spec["edges"], spec["edges"])
            # Include under/overflow in covariance if folding is disabled.
            for i in range(nominal.GetNbinsX() + 2):
                for j in range(nominal.GetNbinsX() + 2):
                    fr_cov = sum(hist.GetBinContent(i) * hist.GetBinContent(j) * rates.errors[nuisance]**2
                                 for nuisance, hist in gradients[key].items())
                    value = fr_cov + (nominal.GetBinError(i)**2 if i == j else 0.)
                    covariance.SetBinContent(i, j, value)
                    if i == j:
                        combined.SetBinError(i, math.sqrt(max(0., value)))
                if nominal.GetBinContent(i) < 0:
                    audit["negative_bins"].append({"path": path, "bin": i, "yield": nominal.GetBinContent(i)})
            write_object(handle, path, nominal, "zx_stat")
            write_object(handle, path, combined, "zx_total")
            write_object(handle, path, covariance, "covariance_stat_plus_fr")
            for name, hist in components[key].items():
                write_object(handle, path, hist, name)
            for nuisance, gradient in gradients[key].items():
                label = "fr_" + "_".join(map(str, nuisance))
                for sign, suffix in ((1., "Up"), (-1., "Down")):
                    shifted = clone(nominal)
                    shifted.Add(gradient, sign * rates.errors[nuisance])
                    # Variation bin errors describe the nominal application statistics.
                    for i in range(nominal.GetNbinsX() + 2):
                        shifted.SetBinError(i, nominal.GetBinError(i))
                    write_object(handle, path + "/variations", shifted, label + suffix)
            if variable == "yield":
                audit["yields"].setdefault(era, {})[channel] = {
                    **total(nominal), "stat_plus_fr_error": combined.GetBinError(1),
                    "components": {name: total(hist) for name, hist in components[key].items()}}
        ROOT.TObjString("ZX = data(single)*r - prompt(single)*r - data(double)*r1*r2 + optional prompt(double)*r1*r2; r=f/(1-f)").Write("definition")
    finally:
        handle.Close()
    return audit

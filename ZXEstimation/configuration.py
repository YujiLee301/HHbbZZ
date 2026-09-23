"""YAML input contract; all sample paths and normalization live outside code."""
import glob
import math
import os
from pathlib import Path
import re

import yaml


def edges(values, label):
    result = [float(x) for x in values]
    if len(result) < 2 or not all(math.isfinite(x) for x in result):
        raise ValueError(f"{label}: need at least two finite bin edges")
    if any(a >= b for a, b in zip(result, result[1:])):
        raise ValueError(f"{label}: bin edges must increase strictly")
    return result


def local_path(value, base):
    value = os.path.expandvars(os.path.expanduser(str(value)))
    if "$" in value:
        raise ValueError(f"Unresolved environment variable in {value}")
    path = Path(value)
    return path if path.is_absolute() else base / path


def load_config(path):
    path = Path(path).resolve()
    with path.open(encoding="utf-8-sig") as stream:
        cfg = yaml.safe_load(stream)
    if not isinstance(cfg, dict) or cfg.get("version") != 1:
        raise ValueError("Expected a YAML mapping with version: 1")
    cfg["_base"] = path.parent
    cfg["_path"] = str(path)
    cfg["_output"] = local_path(cfg["output"], path.parent).resolve()
    for section in ("eras", "processes", "regions", "channels", "variables"):
        if not isinstance(cfg.get(section), dict) or not cfg[section]:
            raise ValueError(f"{section}: expected a nonempty mapping")
        for name in cfg[section]:
            if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9_-]+", name):
                raise ValueError(f"Unsafe ROOT/output name: {section}/{name}")
    if "combined" in cfg["eras"] or "data" in cfg["processes"] or "yield" in cfg["variables"]:
        raise ValueError("Reserved names: era 'combined', process 'data', variable 'yield'")
    for era, settings in cfg["eras"].items():
        if (settings.get("lumi_fb") is None or not math.isfinite(float(settings["lumi_fb"]))
                or float(settings["lumi_fb"]) <= 0):
            raise ValueError(f"{era}: supply lumi_fb in YAML")
    roles = [r["role"] for r in cfg["regions"].values()]
    if len(roles) != len(set(roles)) or not set(roles) <= {"rate_fail", "rate_pass", "single_fail", "double_fail"}:
        raise ValueError("Region roles must be unique: rate_fail/rate_pass/single_fail/double_fail")
    for variable, spec in cfg["variables"].items():
        spec["edges"] = edges(spec["edges"], variable)
        unknown = set(spec["regions"]) - set(cfg["regions"])
        if unknown:
            raise ValueError(f"{variable}: unknown regions {unknown}")
    if int(cfg["selection"]["min_jets"]) < 2:
        raise ValueError("Every CR must require at least two selected jets")
    if len(cfg["probes"]) != 2 or any(set(p) != {"pt", "eta", "pdgid", "tight"} for p in cfg["probes"]):
        raise ValueError("Configure two probes with pt, eta, pdgid and tight expressions")
    if not cfg["plotting"]["regions"] or set(cfg["plotting"]["regions"]) - set(cfg["regions"]):
        raise ValueError("plotting.regions must name configured regions")
    if len(cfg["plotting"]["regions"]) != len(set(cfg["plotting"]["regions"])):
        raise ValueError("Duplicate plotting regions")
    if not cfg["plotting"]["formats"] or set(cfg["plotting"]["formats"]) - {"png", "pdf", "svg"}:
        raise ValueError("Plot formats must be png, pdf or svg")
    ratio = cfg["plotting"]["ratio_range"]
    if len(ratio) != 2 or not all(math.isfinite(float(v)) for v in ratio) or float(ratio[0]) >= float(ratio[1]):
        raise ValueError("plotting.ratio_range must be [low, high]")
    fr = cfg["fake_rate"]
    if fr["source"] not in ("measure", "external"):
        raise ValueError("fake_rate.source must be measure or external")
    for flavor, spec in fr["binning"].items():
        spec["pt"] = edges(spec["pt"], f"{flavor}/pt")
        spec["abs_eta"] = edges(spec["abs_eta"], f"{flavor}/eta")
    if {int(x) for x in fr["binning"]} != {11, 13}:
        raise ValueError("fake_rate.binning requires PDG flavors 11 and 13")
    # Normalize YAML numeric keys to integers.
    fr["binning"] = {int(k): v for k, v in fr["binning"].items()}
    for group in (fr["subtract_processes"], cfg["estimation"]["subtract_single"], cfg["estimation"]["subtract_double"]):
        if len(group) != len(set(group)) or set(group) - set(cfg["processes"]):
            raise ValueError(f"Unknown or duplicate subtraction process in {group}")
    if fr["pt_overflow"] not in ("last_bin", "error"):
        raise ValueError("fake_rate.pt_overflow must be last_bin or error")
    if not cfg.get("samples"):
        raise ValueError("No samples configured")
    names, used_files = set(), set()
    for sample in cfg["samples"]:
        name = sample["name"]
        if name in names:
            raise ValueError(f"Duplicate sample name {name}")
        names.add(name)
        if sample["era"] not in cfg["eras"] or sample["kind"] not in ("data", "mc"):
            raise ValueError(f"{name}: invalid era/kind")
        if sample["kind"] == "mc" and sample["process"] not in cfg["processes"]:
            raise ValueError(f"{name}: unknown MC process")
        if sample["kind"] == "mc" and "weight" not in sample:
            raise ValueError(f"{name}: explicitly configure an MC weight expression")
        if not isinstance(sample["files"], list):
            raise ValueError(f"{name}: files must be a YAML list")
        files = []
        for pattern in sample["files"]:
            pattern = os.path.expandvars(os.path.expanduser(str(pattern)))
            if "$" in pattern:
                raise ValueError(f"Unresolved environment variable in {pattern}")
            if str(pattern).startswith(("root://", "https://")):
                if any(char in str(pattern) for char in "*?["):
                    raise ValueError("Remote files must be explicit URLs, not globs")
                matched = [str(pattern)]
            else:
                matched = sorted(glob.glob(str(local_path(pattern, path.parent))))
                matched = [str(Path(p).resolve()) for p in matched]
            if not matched:
                raise ValueError(f"{name}: no files match {pattern}")
            for filename in matched:
                if filename in used_files:
                    raise ValueError(f"Input file listed more than once: {filename}")
                used_files.add(filename)
                files.append(filename)
        if not files:
            raise ValueError(f"{name}: fill the files list in YAML")
        sample["_files"] = files
        sample["_scale"] = normalization(sample, cfg)
    for era in cfg["eras"]:
        samples = [s for s in cfg["samples"] if s["era"] == era]
        if not any(s["kind"] == "data" for s in samples):
            raise ValueError(f"{era}: missing data")
        missing = set(cfg["processes"]) - {s.get("process") for s in samples if s["kind"] == "mc"}
        if missing:
            raise ValueError(f"{era}: missing MC processes {sorted(missing)}")
    return cfg


def normalization(sample, cfg):
    if sample["kind"] == "data":
        return 1.0
    norm = sample["normalization"]
    if norm["mode"] == "preweighted":
        scale = float(norm.get("scale", 1.0))
    elif norm["mode"] == "xsec":
        if norm.get("xsec_pb") is None or norm.get("sum_gen_weights") is None:
            raise ValueError(f"{sample['name']}: supply xsec_pb and full-production sum_gen_weights")
        xsec, denominator = float(norm["xsec_pb"]), float(norm["sum_gen_weights"])
        if not math.isfinite(xsec) or not math.isfinite(denominator) or xsec <= 0 or denominator == 0:
            raise ValueError(f"{sample['name']}: invalid cross section or normalization denominator")
        scale = (1000.0 * float(cfg["eras"][sample["era"]]["lumi_fb"]) * xsec /
                 denominator * float(norm.get("k_factor", 1.0)) * float(norm.get("filter_efficiency", 1.0)))
    else:
        raise ValueError("normalization.mode must be xsec or preweighted")
    if not math.isfinite(scale):
        raise ValueError(f"{sample['name']}: nonfinite normalization")
    return scale

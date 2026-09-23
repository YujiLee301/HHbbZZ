#!/usr/bin/env python3
"""Run fake-rate measurement, OS Z+X estimation, and CR data/MC plots."""
import argparse
import json
from pathlib import Path
import sys

from configuration import load_config


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True, help="External YAML sample/analysis configuration")
    parser.add_argument("--stage", choices=("all", "plots", "fake-rate", "estimate"), default="all")
    args = parser.parse_args()
    cfg = load_config(args.config)
    # Import ROOT only after argument/config parsing; --help works without ROOT.
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptStat(0)
    from fake_rates import measure, load
    from estimation import estimate
    from plotting import control_plots, rate_plots

    directory = cfg["_output"]
    directory.mkdir(parents=True, exist_ok=True)
    manifest = {"config": cfg["_path"], "stage": args.stage,
                "root_version": ROOT.gROOT.GetVersion(), "status": "running",
                "samples": [{"name": s["name"], "era": s["era"], "kind": s["kind"],
                             "process": s.get("process"), "files": s["_files"], "scale": s["_scale"],
                             "weight": s.get("weight", "1")} for s in cfg["samples"]]}
    manifest_path = directory / f"manifest_{args.stage}.json"
    write_json(manifest_path, manifest)
    (directory / f"config_{args.stage}.yaml").write_text(args.config.read_text(encoding="utf-8-sig"), encoding="utf-8")
    try:
        # Draw the CRs first so diagnostic outputs survive a sparse FR map.
        if args.stage in ("all", "plots"):
            print("Producing CR data/MC plots", flush=True)
            audit = control_plots(ROOT, cfg, directory / "control_regions.root")
            write_json(directory / "control_regions.json", audit)
        if args.stage in ("all", "fake-rate", "estimate"):
            rate_file = directory / "fake_rates.root"
            if cfg["fake_rate"]["source"] == "measure" and args.stage in ("all", "fake-rate"):
                print("Measuring prompt-subtracted fake rates", flush=True)
                rates, audit = measure(ROOT, cfg, rate_file)
                write_json(directory / "fake_rates.json", audit)
            else:
                rates = load(ROOT, cfg, rate_file)
            rate_plots(ROOT, cfg, rates)
            if args.stage in ("all", "estimate"):
                print("Building OS Z+X templates", flush=True)
                audit = estimate(ROOT, cfg, rates, directory / "zx_estimation.root")
                write_json(directory / "zx_estimation.json", audit)
        manifest["status"] = "complete"
    except Exception as error:
        manifest["status"] = "failed"
        manifest["error"] = str(error)
        raise
    finally:
        write_json(manifest_path, manifest)
    print(f"Outputs: {directory}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, ImportError, KeyError) as error:
        print(f"ZXEstimation: {error}", file=sys.stderr)
        sys.exit(1)

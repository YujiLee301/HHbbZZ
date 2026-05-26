#!/usr/bin/env python3
"""
Build b-tag efficiency JSON files (correctionlib v2 schema) for use with
PhysicsTools.NATModules.modules.jetBtag.

The jetBtag module needs, for MC, a "btag_efficiency" correction binned in
(flavor, pt, abseta). These efficiencies are analysis-specific (they depend
on the jet pT spectrum after the analysis preselection), so each analysis
must derive them from its own MC.

Usage
-----
Run separately per era. Pass MC files that already have your analysis
preselection applied (e.g. the post-skim trees), so the eff is representative
of the phase space where the b-tag weights will be used.

    python3 scripts/build_btag_eff.py \\
        --era 2022_Summer22EE \\
        --out data/btag/2022_Summer22EE_RPT_M_eff.json.gz \\
        --tree Events \\
        /path/to/skim/Run3Summer22EE/*.root

Output files are placed in data/btag/ with names matching what
corrections_config.py expects.

Requirements
------------
- correctionlib (already in CMSSW). Only used to look up the WP value from the
  central BTV JSON.
- uproot + numpy. Available in CMSSW; otherwise:  pip install --user uproot numpy
"""

import argparse
import gzip
import json
import os
import sys
import numpy as np

# ---- Configuration -----------------------------------------------------------

# Binning matches the reference example file shipped with NATModules
# (PhysicsTools/NATModules/test/btag_2022EE.json.gz). Adjust if you need finer
# resolution and have enough MC stats per bin.
PT_EDGES = [30.0, 40.0, 50.0, 70.0, 100.0, 150.0, 200.0, 300.0, 600.0]
ABSETA_EDGES = [0.0, 0.6, 1.2, 2.1, 2.5]
FLAVORS = ["light", "c", "b"]

# Tagger to evaluate. Must match _RPT_TAGGER / _RPT_TAGGER_NAME / _RPT_WP in
# corrections_config.py.
TAGGER_JSON_KEY = "robustParticleTransformer"  # key in BTV correctionlib JSON
TAGGER_BRANCH = "Jet_btagRobustParTAK4B"       # NanoAOD branch
WP_NAME = "M"

# Mapping era -> central BTV SF JSON. Used only to look up the WP value, so the
# script and corrections_config.py stay consistent if BTV updates the WPs.
BTV_SF_JSON = {
    "2022_Summer22":     "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2022_Summer22/btagging.json.gz",
    "2022_Summer22EE":   "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2022_Summer22EE/btagging.json.gz",
    "2023_Summer23":     "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2023_Summer23/btagging.json.gz",
    "2023_Summer23BPix": "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2023_Summer23BPix/btagging.json.gz",
}


def get_wp_value(era):
    """Look up the numeric WP cut from the central BTV correctionlib JSON."""
    import correctionlib
    cs = correctionlib.CorrectionSet.from_file(BTV_SF_JSON[era])
    return float(cs["%s_wp_values" % TAGGER_JSON_KEY].evaluate(WP_NAME))


def flavor_label(hf):
    """Map Jet_hadronFlavour to the string categories used in the JSON."""
    hf = np.abs(hf)
    out = np.full(hf.shape, "light", dtype=object)
    out[hf == 4] = "c"
    out[hf == 5] = "b"
    return out


def fill_from_file(path, tree_name, wp, n_tot, n_tag):
    """Fill numerator/denominator histograms from one input file."""
    import uproot
    branches = [
        "Jet_pt",
        "Jet_eta",
        "Jet_hadronFlavour",
        TAGGER_BRANCH,
    ]
    try:
        with uproot.open(path) as f:
            if tree_name not in f:
                print("  [warn] tree '%s' not found in %s, skipping" % (tree_name, path))
                return 0
            arrays = f[tree_name].arrays(branches, library="np")
    except Exception as e:
        print("  [warn] could not read %s: %s" % (path, e))
        return 0

    # Flatten across events
    pt = np.concatenate(arrays["Jet_pt"]) if len(arrays["Jet_pt"]) else np.array([])
    if pt.size == 0:
        return 0
    aeta = np.abs(np.concatenate(arrays["Jet_eta"]))
    hf = np.concatenate(arrays["Jet_hadronFlavour"])
    disc = np.concatenate(arrays[TAGGER_BRANCH])

    # Restrict to the histogram support; clamp behavior is set in the JSON.
    in_range = (
        (pt >= PT_EDGES[0]) & (pt < PT_EDGES[-1]) & (aeta < ABSETA_EDGES[-1])
    )
    pt, aeta, hf, disc = pt[in_range], aeta[in_range], hf[in_range], disc[in_range]
    if pt.size == 0:
        return 0

    lab = flavor_label(hf)
    for fl in FLAVORS:
        sel = lab == fl
        if not sel.any():
            continue
        tot, _, _ = np.histogram2d(pt[sel], aeta[sel], bins=[PT_EDGES, ABSETA_EDGES])
        sel_tag = sel & (disc > wp)
        tag, _, _ = np.histogram2d(
            pt[sel_tag], aeta[sel_tag], bins=[PT_EDGES, ABSETA_EDGES]
        )
        n_tot[fl] += tot
        n_tag[fl] += tag

    return int(pt.size)


def make_correctionlib_json(era, wp, n_tot, n_tag):
    """Assemble the correctionlib v2 JSON dictionary."""

    def make_pt_block(fl):
        # Zero-stats bins get eff=0; jetBtag.py guards against eff<=0 and falls
        # back to the unweighted decision. Widen bins or merge samples if many
        # bins are empty.
        with np.errstate(divide="ignore", invalid="ignore"):
            eff = np.divide(
                n_tag[fl],
                n_tot[fl],
                out=np.zeros_like(n_tot[fl]),
                where=n_tot[fl] > 0,
            )
        return {
            "nodetype": "binning",
            "input": "pt",
            "edges": list(PT_EDGES),
            "flow": "clamp",
            "content": [
                {
                    "nodetype": "binning",
                    "input": "abseta",
                    "edges": list(ABSETA_EDGES),
                    "flow": "clamp",
                    "content": eff[ipt].tolist(),
                }
                for ipt in range(len(PT_EDGES) - 1)
            ],
        }

    return {
        "schema_version": 2,
        "description": (
            "b-tag MC efficiency map for %s WP=%s, %s "
            "(derived from analysis-preselected MC)" % (TAGGER_JSON_KEY, WP_NAME, era)
        ),
        "corrections": [
            {
                "name": "btag_efficiency",
                "description": (
                    "MC efficiency for %s WP=%s (WP cut value = %.6f)"
                    % (TAGGER_JSON_KEY, WP_NAME, wp)
                ),
                "version": 1,
                "inputs": [
                    {"name": "flavor", "type": "string",
                     "description": "Jet hadron flavour label: 'b', 'c', or 'light'"},
                    {"name": "pt", "type": "real", "description": "Jet pT [GeV]"},
                    {"name": "abseta", "type": "real",
                     "description": "Jet |eta|"},
                ],
                "output": {
                    "name": "efficiency",
                    "type": "real",
                    "description": "MC b-tag efficiency in [0, 1]",
                },
                "data": {
                    "nodetype": "category",
                    "input": "flavor",
                    "content": [
                        {"key": fl, "value": make_pt_block(fl)} for fl in FLAVORS
                    ],
                },
            }
        ],
    }


def print_summary(n_tot, n_tag):
    """Print a small per-flavor coverage report."""
    print("\nBin coverage summary:")
    for fl in FLAVORS:
        total = float(n_tot[fl].sum())
        tagged = float(n_tag[fl].sum())
        empty = int(np.sum(n_tot[fl] == 0))
        nbins = n_tot[fl].size
        avg = (tagged / total) if total > 0 else 0.0
        print(
            "  flavor=%-5s  jets=%-10d  tagged=%-10d  <eff>=%.4f  empty_bins=%d/%d"
            % (fl, int(total), int(tagged), avg, empty, nbins)
        )


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--era", required=True, choices=sorted(BTV_SF_JSON),
                   help="Era label; used both to pick the WP and to name the output.")
    p.add_argument("--out", required=True,
                   help="Output path (.json.gz).")
    p.add_argument("--tree", default="Events",
                   help="Input TTree name (default: Events).")
    p.add_argument("files", nargs="+",
                   help="Input ROOT files (post-selection MC).")
    args = p.parse_args()

    wp = get_wp_value(args.era)
    print("Era:        %s" % args.era)
    print("Tagger key: %s  WP=%s  cut=%.6f" % (TAGGER_JSON_KEY, WP_NAME, wp))
    print("Branch:     %s" % TAGGER_BRANCH)
    print("Inputs:     %d file(s)" % len(args.files))

    shape = (len(PT_EDGES) - 1, len(ABSETA_EDGES) - 1)
    n_tot = {fl: np.zeros(shape, dtype=np.float64) for fl in FLAVORS}
    n_tag = {fl: np.zeros(shape, dtype=np.float64) for fl in FLAVORS}

    n_jets_total = 0
    for i, path in enumerate(args.files, 1):
        n = fill_from_file(path, args.tree, wp, n_tot, n_tag)
        n_jets_total += n
        if i % 10 == 0 or i == len(args.files):
            print("  processed %d/%d files, %d jets so far"
                  % (i, len(args.files), n_jets_total))

    if n_jets_total == 0:
        print("ERROR: no jets read from any input file. Check --tree and inputs.",
              file=sys.stderr)
        return 1

    print_summary(n_tot, n_tag)

    out_dict = make_correctionlib_json(args.era, wp, n_tot, n_tag)
    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    with gzip.open(args.out, "wt") as f:
        json.dump(out_dict, f)
    print("\nWrote %s" % args.out)

    # Quick self-test: re-open with correctionlib and probe one point.
    try:
        import correctionlib
        cs = correctionlib.CorrectionSet.from_file(args.out)
        e = cs["btag_efficiency"].evaluate("b", 50.0, 1.0)
        print("Self-test: eff(b, pT=50 GeV, |eta|=1.0) = %.4f" % e)
    except Exception as e:
        print("Self-test failed (file written but not loadable): %s" % e,
              file=sys.stderr)
        return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())

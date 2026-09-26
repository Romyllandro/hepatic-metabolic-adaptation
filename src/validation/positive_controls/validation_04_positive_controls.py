#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validation_04_positive_controls.py
==================================
PRONG 4 of the validation plan: does the pipeline recover KNOWN diet physiology?

This is a "recovery of ground truth" check where the ground truth is textbook
hepatic biochemistry rather than a benchmark dataset. It requires NO solver and
NO model file -- it runs entirely on the pipeline's existing output CSVs, so it
is fully reproducible by anyone with the repo outputs.

SOLVER-STEP NOTE (revised)
--------------------------
No change was required to this script. It reads RQ1_* output CSVs and performs no
optimisation of its own, so the FBA/pFBA distinction does not arise here. It inherits
whatever solver step produced those CSVs -- plain FBA, for RQ1/RQ2. State this when
citing these results, so the positive-control evidence is not implicitly attributed
to a pFBA pipeline.

WHY THIS APPROACH
-----------------
A composed pipeline can be internally consistent yet biologically wrong. Positive
controls catch that: if dietary-fat loading does not raise predicted beta-oxidation
flux, the customisation has broken something the engine should capture. Encoding the
expectations in an explicit, auditable table (diet x pathway x expected direction)
makes the check falsifiable and lets a reviewer see exactly what we claimed and
whether it held.

WHAT IT PRODUCES
----------------
1. Pathway-level fold-change vs SCD baseline (mean |flux| over each control set).
2. Per-marker signed flux + effect size + FDR from the reaction-stats files.
3. A pass/fail report against the expectations table -> positive_control_report.csv

USAGE
-----
    python validation_04_positive_controls.py --data_dir /path/to/RQ_outputs_files
"""
import argparse, os, sys
import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Control reaction sets (iMM1415 IDs). Pathway sets are matched by subsystem;
# marker reactions are the cleanest single-reaction readouts.
# --------------------------------------------------------------------------- #
PATHWAY_SUBSYSTEMS = {
    "Fatty acid oxidation": ["Fatty acid oxidation"],          # beta-oxidation
    "Cholesterol synthesis": ["Cholesterol Metabolism"],       # filtered to synth markers below
}
CHOLESTEROL_SYNTH_IDS = ["HMGCOARr", "HMGCOARx", "SQLSr", "SQLEr", "LNSTLSr",
                         "DHCR241r", "DHCR242r", "DHCR243r", "DHCR71r", "DHCR72r",
                         "C14STRr", "C3STKR2r"]
KETOGENESIS_IDS = ["HMGCOASm", "BDHm", "ACACT1m", "ADCim", "ACACt2m", "EX_acac_e"]

# Expectations: (pathway, diet, expected_direction vs SCD). 'up'/'down'.
# Grounded in standard hepatic physiology; see VALIDATION_STRATEGY.md.
EXPECTATIONS = [
    ("Fatty acid oxidation", "KD",  "up"),
    ("Fatty acid oxidation", "HFD", "up"),
    ("Fatty acid oxidation", "WD",  "up"),
    ("Cholesterol synthesis", "WD", "up"),   # fructose/carbohydrate-driven DNL & sterol synthesis
]
# Marker-level expectations (single rate-limiting reactions), checked with stats.
MARKER_EXPECTATIONS = [
    # (marker_id, diet, expected_direction, note)
    ("EX_acac_e", "WD", "up", "ketone-body (acetoacetate) export under Western diet"),
    ("HMGCOASm",  "KD", "up", "Hmgcs2 = rate-limiting ketogenic enzyme (textbook KD response)"),
]


def load_inputs(data_dir):
    ext = pd.read_csv(os.path.join(data_dir, "RQ1_reaction_flux_comparison_extended.csv"))
    ann = pd.read_csv(os.path.join(data_dir, "RQ4_reaction_annotations.csv"))
    stats = {}
    for d in ["HFD", "KD", "WD"]:
        f = os.path.join(data_dir, f"RQ1_reaction_stats_{d}_vs_SCD.csv")
        if os.path.exists(f):
            stats[d] = pd.read_csv(f).set_index("ReactionID")
    return ext.set_index("ReactionID"), ann.set_index("reaction_id"), stats


def pathway_fold_changes(ext, ann):
    """Mean |flux| per condition for each control pathway, expressed as fold vs SCD."""
    meancols = {c.split("_MeanFlux")[0]: c for c in ext.columns if c.endswith("_MeanFlux")}
    conds = [c for c in ["SCD", "HFD", "KD", "WD"] if c in meancols]

    sets = {}
    fao = ann.index[ann["subsystem"].astype(str).str.contains("Fatty acid oxidation", case=False)]
    sets["Fatty acid oxidation"] = [r for r in fao if r in ext.index]
    sets["Cholesterol synthesis"] = [r for r in CHOLESTEROL_SYNTH_IDS if r in ext.index]
    sets["Ketogenesis (mito)"] = [r for r in KETOGENESIS_IDS if r in ext.index]

    rows = []
    for label, rxns in sets.items():
        sub = ext.loc[rxns]
        activity = {c: float(np.nanmean(np.abs(sub[meancols[c]].values))) for c in conds}
        base = activity.get("SCD", np.nan)
        row = {"pathway": label, "n_reactions": len(rxns)}
        for c in conds:
            row[f"{c}_mean_abs_flux"] = round(activity[c], 5)
        for c in conds:
            if c != "SCD":
                row[f"{c}_fold_vs_SCD"] = round(activity[c] / base, 3) if base else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def evaluate(fold_df, stats):
    """Score pathway- and marker-level expectations -> pass/fail rows."""
    report = []
    # pathway-level
    fold_idx = fold_df.set_index("pathway")
    for pathway, diet, expected in EXPECTATIONS:
        if pathway not in fold_idx.index:
            continue
        col = f"{diet}_fold_vs_SCD"
        fold = fold_idx.loc[pathway, col] if col in fold_idx.columns else np.nan
        observed = "up" if (pd.notna(fold) and fold > 1.0) else "down"
        report.append({
            "level": "pathway", "target": pathway, "diet": diet,
            "expected": expected, "observed": observed,
            "metric": f"fold={fold}", "pass": bool(observed == expected),
        })
    # marker-level (signed direction + FDR)
    for marker, diet, expected, note in MARKER_EXPECTATIONS:
        if diet not in stats or marker not in stats[diet].index:
            continue
        r = stats[diet].loc[marker]
        observed = "up" if r.get("log2FC", 0) > 0 else "down"
        report.append({
            "level": "marker", "target": f"{marker} ({note})", "diet": diet,
            "expected": expected, "observed": observed,
            "metric": f"log2FC={r.get('log2FC'):+.3f}; Cohen_d={r.get('Cohen_d'):+.3f}; "
                      f"q={r.get('q_value'):.3g}; sig={r.get('Significant')}",
            "pass": bool(observed == expected),
        })
    return pd.DataFrame(report)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_dir", default=".", help="Dir with RQ1_*/RQ4_* output CSVs")
    ap.add_argument("--out", default="positive_control_report.csv")
    args = ap.parse_args()

    ext, ann, stats = load_inputs(args.data_dir)
    fold_df = pathway_fold_changes(ext, ann)

    print("=" * 70)
    print("PRONG 4: POSITIVE-CONTROL PHYSIOLOGY RECOVERY")
    print("=" * 70)
    print("\nPathway activity (mean |flux|) and fold-change vs SCD:\n")
    print(fold_df.to_string(index=False))

    report = evaluate(fold_df, stats)
    print("\nExpectation scorecard:\n")
    print(report.to_string(index=False))
    n_pass = int(report["pass"].sum()); n_tot = len(report)
    print(f"\nPASSED {n_pass}/{n_tot} expectations.")
    print("NOTE: marker-level KD ketogenesis is a known caveat at the pooled-bulk "
          "level (see VALIDATION_STRATEGY.md). Report transparently.")

    report.to_csv(args.out, index=False)
    fold_df.to_csv(args.out.replace(".csv", "_pathway_folds.csv"), index=False)
    print(f"\n[OK] Wrote {args.out} and {args.out.replace('.csv','_pathway_folds.csv')}")


if __name__ == "__main__":
    sys.exit(main())

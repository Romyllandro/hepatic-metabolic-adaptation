#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_measured_modules_21.py
=============================
Restricts the external 13C-MFA ground truth used by benchmark_ecoli_13c_methods.py
to the SAME 21 curated central-carbon reactions already used to produce the
manuscript's existing Table SL1 (mean r = 0.817 DC / 0.789 AC for the reduced
E-Flux/pFBA engine).

WHY
---
benchmark_ecoli_13c_methods.py's --measured_modules input (ecoli_inputs/measured_flux_modules.csv)
evaluates ~35 Gerosa reference modules per condition via their rnx_relationship
(sums/medians of one or more model reactions). 8 of those 35 modules never resolve
to an iJO1366 reaction ID at all (see ecoli_inputs/gerosa_module_mapping_audit.csv:
Suc_Ex, Fum_Ex, PFK-FBP, PYK-PPS, SUCDH3, MDH+MQO, ME1+ME2 are fully unresolved, and
Growth_rate has "missing relationship" entirely). That is a DIFFERENT, broader ground
truth than the 21-reaction panel Table SL1 was built from, which is why running eflux
through it does not reproduce 0.817/0.789 (it gives ~0.66 DC / ~0.51 AC instead).

The 21-reaction panel itself already exists on disk, built earlier by
prepare_benchmark_inputs.py, at prepared_def/meas_Ecoli_<condition>.csv (one row per
reaction: reaction_id,measured_flux; identical 21 reaction IDs across all 8 conditions).
This is the same input run_benchmark.py --prepared prepared_def feeds into
validation_01_benchmark_recovery.py -- the original engine-fidelity script that produced
Table SL1.

This script converts those 8 flat files into the tidy module-relationship format
benchmark_ecoli_13c_methods.py's --measured_modules expects (organism,condition,
module_short,module_name,measured_flux,sd,rnx_relationship), with each of the 21
reactions treated as its own trivial 1:1 module (module_short == rnx_relationship ==
the exact reaction ID already used in prepared_def). No aggregation, no ambiguity:
every row should resolve to exactly one model reaction with zero "unresolved" entries
in the mapping audit that benchmark_ecoli_13c_methods.py writes back out.

USAGE
-----
    python build_measured_modules_21.py --prepared_def prepared_def \
        --out ecoli_inputs/measured_flux_modules_21reaction.csv

Then run the existing multi-method benchmark against this restricted, previously-verified
ground truth in place of the full ecoli_inputs/measured_flux_modules.csv:

    python benchmarking\\benchmark_ecoli_13c_methods.py ^
      --model iJO1366.json ^
      --objective BIOMASS_Ec_iJO1366_core_53p95M ^
      --transcriptomes ecoli_inputs\\transcriptomes.csv ^
      --measured_modules ecoli_inputs\\measured_flux_modules_21reaction.csv ^
      --uptakes ecoli_inputs\\uptakes.csv ^
      --regime both ^
      --solver gurobi ^
      --methods eflux,gimme,pfba ^
      --output reviewer_results\\benchmark_ecoli_13c_21reaction

Inspect the new gerosa_module_mapping_audit.csv from that run first: with this input,
every one of the 21 reactions x 8 conditions = 168 rows should show a fully resolved
single reaction ID and NO issues. If eflux's mean_uncentered_pearson under this run is
close to 0.817 (DC) / 0.789 (AC), that verifies the new script reproduces the original
reduced-engine result on a matched ground truth, and GIMME/pFBA (and RIPTiDe once
implemented) can be trusted on the same restricted panel for the Major Comment 2
closure table. If it does NOT land near 0.817/0.789, that is a second-level
discrepancy (e.g. E-Flux settings/solver differences between the two scripts) that
needs to be resolved before quoting any of these numbers to the reviewer -- do not
average over this in the writeup either way; report exactly what this run gives.
"""
import argparse
import csv
import glob
import os
import sys

ORGANISM = "E. coli"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepared_def", default="prepared_def",
                     help="Directory with meas_Ecoli_<condition>.csv (reaction_id,measured_flux)")
    ap.add_argument("--out", default="ecoli_inputs/measured_flux_modules_21reaction.csv")
    args = ap.parse_args()

    meas_files = sorted(glob.glob(os.path.join(args.prepared_def, "meas_Ecoli_*.csv")))
    if not meas_files:
        sys.exit(f"[ERROR] no meas_Ecoli_*.csv found in {args.prepared_def}")

    rows = []
    reaction_sets = {}
    for path in meas_files:
        condition = os.path.basename(path)[len("meas_Ecoli_"):-len(".csv")]
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            rxns = []
            for r in reader:
                rxn_id = r["reaction_id"].strip()
                flux = r["measured_flux"].strip()
                rxns.append(rxn_id)
                rows.append({
                    "organism": ORGANISM,
                    "condition": condition,
                    "module_short": rxn_id,
                    "module_name": rxn_id,
                    "measured_flux": flux,
                    "sd": "",
                    "rnx_relationship": rxn_id,
                })
            reaction_sets[condition] = tuple(rxns)

    # Sanity check: every condition must use the identical 21-reaction panel,
    # exactly like Table SL1 (otherwise the comparison across conditions is not apples-to-apples).
    distinct = set(reaction_sets.values())
    if len(distinct) != 1:
        print("[WARN] conditions do not share an identical reaction panel:")
        for cond, rxns in reaction_sets.items():
            print(f"  {cond}: {len(rxns)} reactions")
    else:
        n = len(next(iter(distinct)))
        print(f"[OK] all {len(reaction_sets)} conditions share the same {n}-reaction panel.")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "organism", "condition", "module_short", "module_name",
            "measured_flux", "sd", "rnx_relationship"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] wrote {len(rows)} rows ({len(meas_files)} conditions x "
          f"{len(rows)//len(meas_files)} reactions) to {args.out}")


if __name__ == "__main__":
    main()

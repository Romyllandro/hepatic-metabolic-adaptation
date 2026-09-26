#!/usr/bin/env python3
"""
flatten_outputs.py
==================
Copy pipeline outputs into one flat folder with the file names used by the
manuscript tables and figure scripts (the production folder was called
`RQs_Output_Files_final`): RQ1_<table>.csv and RQ2_<STRAIN>_<table>.csv.

This is bookkeeping only (files are copied byte-for-byte, nothing is recomputed).

Usage
-----
python flatten_outputs.py --out results/RQs_Output_Files --rq1 results/RQ1_primary_FBA
python flatten_outputs.py --out results/RQs_Output_Files --rq2 results/RQ2_primary_FBA_aggregate \
       [--annotations inputs/bulk/RQ4_reaction_annotations.csv]
"""
import argparse, shutil
from pathlib import Path

STRAINS = ["129S1SvImJ", "AJ", "C57BL6J", "CASTEiJ", "DBA2J", "NODShiLtJ", "NZOHlLtJ", "PWKPhJ", "WSBEiJ"]


def cp(src, dst):
    if Path(src).exists():
        Path(dst).parent.mkdir(parents=True, exist_ok=True); shutil.copy2(src, dst); return 1
    print(f"[skip] missing {src}"); return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rq1"); ap.add_argument("--rq2"); ap.add_argument("--annotations")
    a = ap.parse_args(); o = Path(a.out); o.mkdir(parents=True, exist_ok=True); n = 0
    if a.rq1:
        r = Path(a.rq1)
        n += cp(r / "flux_analysis/reaction_flux_comparison_extended.csv", o / "RQ1_reaction_flux_comparison_extended.csv")
        n += cp(r / "flux_analysis/reaction_flux_comparison_extended_batch_corrected.csv",
                o / "RQ1_reaction_flux_comparison_extended_batch_corrected.csv")
        n += cp(r / "stats_comparison/flux_pairwise_stats.csv", o / "RQ1_flux_pairwise_stats.csv")
        n += cp(r / "stats_comparison/flux_global_distance_matrix.csv", o / "RQ1_flux_global_distance_matrix.csv")
        n += cp(r / "rank_product.csv", o / "RQ1_rank_product.csv")
        for f in (r / "cytoscape_edges").glob("edges_*.csv"):
            n += cp(f, o / f"RQ1_{f.name}")
        for f in (r / "comprehensive_analysis/csv_outputs").glob("*.csv"):
            n += cp(f, o / f"RQ1_{f.name}")
    if a.rq2:
        r = Path(a.rq2)
        for s in STRAINS:
            d = r / f"results_{s}_GSE182668"
            n += cp(d / "flux_analysis/reaction_flux_comparison_extended.csv", o / f"RQ2_{s}_reaction_flux_comparison_extended.csv")
            n += cp(d / "stats_comparison/flux_pairwise_stats.csv", o / f"RQ2_{s}_flux_pairwise_stats.csv")
            n += cp(d / "rank_product.csv", o / f"RQ2_{s}_rank_product.csv")
            n += cp(d / "cytoscape_edges/edges_HFD_vs_SCD.csv", o / f"RQ2_{s}_edges_HFD_vs_SCD.csv")
    if a.annotations:
        n += cp(a.annotations, o / "RQ4_reaction_annotations.csv")
    print(f"[flatten] {n} files -> {o}")


if __name__ == "__main__":
    main()

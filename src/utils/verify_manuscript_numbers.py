#!/usr/bin/env python3
"""
verify_manuscript_numbers.py
============================
Recompute the manuscript headline numbers from pipeline outputs (or from the frozen
tables in reference_outputs/) and compare them with EXPECTED_MANUSCRIPT_RESULTS.json.
No optimisation solver is needed. Exit code 0 = every checked value matches.

Usage
-----
# frozen reference tables shipped with the repository (runs in seconds)
python src/utils/verify_manuscript_numbers.py --expected reference_outputs/EXPECTED_MANUSCRIPT_RESULTS.json \
   --rq1_stats reference_outputs/RQ1_primary_FBA --rq1_prefix RQ1_reaction_stats_ \
   --rq2_flat reference_outputs/RQ2_primary_FBA_aggregate --rq2_replicate_flat reference_outputs/RQ2_FBA_replicate \
   --benchmark reference_outputs/benchmarks/benchmark_5method_summary.csv
# a fresh pipeline run
python src/utils/verify_manuscript_numbers.py --expected ... --rq1_stats results/RQ1_primary_FBA/comprehensive_analysis/csv_outputs \
   --rq2_root results/RQ2_primary_FBA_aggregate --rq2_replicate_root results/RQ2_FBA_replicate
"""
import argparse, json, subprocess, sys, tempfile
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
SRC = HERE.parent


def run_json(cmd, outdir, name):
    subprocess.run([sys.executable] + [str(c) for c in cmd], check=True, stdout=subprocess.DEVNULL)
    return json.loads((Path(outdir) / name).read_text())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--expected", required=True)
    ap.add_argument("--rq1_stats"); ap.add_argument("--rq1_prefix", default="reaction_stats_")
    ap.add_argument("--rq2_root"); ap.add_argument("--rq2_flat")
    ap.add_argument("--rq2_replicate_root"); ap.add_argument("--rq2_replicate_flat")
    ap.add_argument("--benchmark")
    a = ap.parse_args()
    exp = json.loads(Path(a.expected).read_text())
    tmp = Path(tempfile.mkdtemp(prefix="verify_"))
    fails, checks = [], 0

    def check(label, got, want):
        nonlocal checks; checks += 1
        ok = (abs(float(got) - float(want)) < 1e-6) if isinstance(want, (int, float)) and not isinstance(want, bool) else (got == want)
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}: got {got}  expected {want}")
        if not ok: fails.append(label)

    if a.rq1_stats:
        print("RQ1 primary (standard FBA)")
        r = run_json([SRC / "rq1_bulk/primary_fba/rq1_signature_summary.py", "--stats_dir", a.rq1_stats,
                      "--prefix", a.rq1_prefix, "--output", tmp / "rq1"], tmp / "rq1", "RQ1_signature_summary.json")
        e = exp["rq1_primary_fba"]
        for c, v in e["significant_per_contrast"].items():
            check(f"significant {c}", r["significant_per_contrast"].get(c), v)
        for k in ["union_SCD_referenced", "union_all_six_contrasts_including_KD_vs_WD", "cross_diet_signature",
                  "signature_directionally_concordant", "signature_reactions"]:
            check(k, r[k], e[k])

    def rq2(root, flat, key, tag):
        args = ["--rq2_root", root] if root else ["--flat_dir", flat]
        return run_json([SRC / "rq2_strains/primary_fba/rq2_conservation_tiers.py", *args, "--output", tmp / tag],
                        tmp / tag, "RQ2_tiers_summary.json")
    agg = None
    if a.rq2_root or a.rq2_flat:
        print("RQ2 primary (standard FBA, aggregate design)")
        agg = rq2(a.rq2_root, a.rq2_flat, "rq2_primary_fba_aggregate", "rq2a")
        for k, v in exp["rq2_primary_fba_aggregate"].items():
            if k == "theta": continue
            check(k, agg[k], v)
    if a.rq2_replicate_root or a.rq2_replicate_flat:
        print("RQ2 replicate-preserving design")
        rep = rq2(a.rq2_replicate_root, a.rq2_replicate_flat, "rq2_fba_replicate_preserving", "rq2r")
        e = exp["rq2_fba_replicate_preserving"]
        check("union", rep["union"], e["union"]); check("universal", rep["universal"], e["universal"])
        if agg is not None:
            ua = set(pd.read_csv(tmp / "rq2a/RQ2_union_membership.csv", index_col=0).index)
            ur = set(pd.read_csv(tmp / "rq2r/RQ2_union_membership.csv", index_col=0).index)
            check("overlap_with_primary_union", len(ua & ur), e["overlap_with_primary_union"])
    if a.benchmark:
        print("Five-method E. coli 13C-MFA benchmark")
        b = pd.read_csv(a.benchmark)
        for reg, d in exp["benchmark_ecoli_5method_mean_r"].items():
            for m, v in d.items():
                got = b.loc[(b.method == m) & (b.regime == reg), "mean_r"]
                check(f"{reg} {m}", round(float(got.iloc[0]), 4) if len(got) else None, v)
    print(f"\n{checks - len(fails)}/{checks} checks passed")
    if fails:
        print("FAILED:", ", ".join(fails)); sys.exit(1)


if __name__ == "__main__":
    main()

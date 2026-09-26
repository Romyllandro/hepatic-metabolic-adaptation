#!/usr/bin/env python3
"""
merge_5method_benchmark.py
==========================
Rebuild the five-method E. coli 13C-MFA benchmark table (Supplementary Table SL3)
from the two long-format outputs produced by the pipeline:

  1. benchmark_ecoli_13c_methods.py  (conventional E-Flux, GIMME, unconstrained pFBA,
     RIPTiDe; zero-filled so every method is scored on the identical 21-reaction panel)
  2. run_benchmark.py -> validation_01_benchmark_recovery.py  (reduced three-layer engine)

Score: uncentered Pearson r between |measured| and |predicted| flux per condition,
then the mean across the 8 growth conditions, per regime
(AC = uptake-unconstrained, DC = uptake-constrained).

Usage
-----
python merge_5method_benchmark.py \
    --methods_long  results/benchmark_ecoli_13c/benchmark_long_format_flux.csv \
    --reduced_long  results/benchmark_ecoli_reduced_engine/benchmark_long_format_flux_reduced_engine.csv \
    --output        results/benchmark_ecoli_13c_5method_final
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def uncentered_r(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    d = np.sqrt((x * x).sum() * (y * y).sum())
    return float((x * y).sum() / d) if d > 0 else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--methods_long", required=True)
    ap.add_argument("--reduced_long", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    m = pd.read_csv(a.methods_long)
    r = pd.read_csv(a.reduced_long)
    r["condition"] = r["condition"].astype(str).str.replace(r"^Ecoli_", "", regex=True)
    keep = ["condition", "regime", "method", "reaction_id", "measured_flux", "predicted_flux", "status"]
    long = pd.concat([m[keep], r[keep]], ignore_index=True)
    long = long.rename(columns={"measured_flux": "measured_flux_raw", "predicted_flux": "predicted_flux_raw"})
    long["measured_flux_abs"] = long["measured_flux_raw"].abs()
    long["predicted_flux_abs"] = long["predicted_flux_raw"].abs()

    rows = []
    for (cond, reg, meth), g in long.groupby(["condition", "regime", "method"]):
        g = g.dropna(subset=["measured_flux_abs", "predicted_flux_abs"])
        rows.append({"condition": cond, "regime": reg, "method": meth,
                     "n_reactions": len(g), "uncentered_r": uncentered_r(g.measured_flux_abs, g.predicted_flux_abs)})
    per_cond = pd.DataFrame(rows)
    summ = (per_cond.groupby(["method", "regime"])
            .agg(mean_r=("uncentered_r", "mean"), n_conditions=("uncentered_r", "size"))
            .reset_index().sort_values(["regime", "method"]))
    summ["mean_r"] = summ["mean_r"].round(4)

    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)
    long.to_csv(out / "benchmark_long_format_flux_5method_merged.csv", index=False)
    per_cond.to_csv(out / "benchmark_5method_by_condition.csv", index=False)
    summ.to_csv(out / "benchmark_5method_summary.csv", index=False)
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()

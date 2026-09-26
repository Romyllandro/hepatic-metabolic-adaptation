#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_benchmark.py
================
Portable replacement for run_benchmark.sh. Iterates the per-condition inputs produced by
prepare_benchmark_inputs.py, runs validation_01_benchmark_recovery.py for each (AC and/or DC),
parses the uncentered-Pearson r from its output, and reports per-condition values plus the
MEAN per regime against the Bhadra-Lobo 2020 targets (DC 0.725 / AC 0.385).

WHY A PYTHON RUNNER INSTEAD OF A .sh
------------------------------------
The generated .sh breaks on Windows: text-mode writing adds CRLF (Git Bash then chokes on
'\r'), and a bash shell often can't see the conda `python`. This runner invokes the SAME
interpreter that runs it (sys.executable), needs no shell, and works identically on Windows,
macOS and Linux -- so the benchmark is reproducible regardless of environment.

USAGE
-----
    python run_benchmark.py --prepared prepared_def --model iJO1366.json \
        --objective BIOMASS_Ec_iJO1366_core_53p95M --regime both
    # --regime {AC,DC,both}   --validation path/to/validation_01_benchmark_recovery.py
"""
import argparse, glob, os, re, subprocess, sys
from pathlib import Path
import pandas as pd

R_RE = re.compile(r"uncentered Pearson\s*:\s*([-\d.]+)")
N_RE = re.compile(r"reactions compared\s*:\s*(\d+)")


def run_one(validation, model, objective, solver, expr, meas, uptake=None):
    # Read validation_01's own per-reaction recovery file (if present from a
    # PREVIOUS call) out of the way first, so a run that errors before writing
    # its own copy can never be mistaken for this call's detail.
    recovery_path = Path("benchmark_flux_recovery.csv")
    if recovery_path.exists():
        recovery_path.unlink()

    cmd = [sys.executable, validation, "--mode", "flux", "--model", model,
           "--objective", objective, "--solver", solver,
           "--expression", expr, "--measured", meas]
    if uptake:
        cmd += ["--uptake", uptake]
    p = subprocess.run(cmd, capture_output=True, text=True)
    out = p.stdout + p.stderr
    r = R_RE.search(out)
    n = N_RE.search(out)

    # DESIGN 2026-09-23: capture this call's own "benchmark_flux_recovery.csv"
    # (written by validation_01_benchmark_recovery.py itself, untouched by this
    # patch) before the NEXT call's subprocess overwrites the same fixed
    # filename. This is the only way to get all 8 conditions' per-reaction
    # detail out of a script that was written to report one condition at a
    # time. validation_01_benchmark_recovery.py's own FBA/pFBA solve logic is
    # not modified by this patch in any way.
    detail_rows = []
    if recovery_path.exists():
        try:
            det = pd.read_csv(recovery_path)
            detail_rows = det.to_dict("records")
        except Exception:
            detail_rows = []

    return (float(r.group(1)) if r else None,
            int(n.group(1)) if n else None,
            p.returncode, out, detail_rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepared", required=True, help="dir from prepare_benchmark_inputs.py")
    ap.add_argument("--model", required=True)
    ap.add_argument("--objective", required=True)
    ap.add_argument("--validation", default="validation_01_benchmark_recovery.py")
    ap.add_argument("--solver", default="gurobi")
    ap.add_argument("--regime", choices=["AC", "DC", "both"], default="both")
    ap.add_argument("--out", default="benchmark_results.csv")
    ap.add_argument("--verbose", action="store_true", help="print each validation_01 run's output")
    args = ap.parse_args()

    meas_files = sorted(glob.glob(os.path.join(args.prepared, "meas_*.csv")))
    if not meas_files:
        sys.exit(f"[ERROR] no meas_*.csv in {args.prepared}; run prepare_benchmark_inputs.py first.")

    regimes = ["AC", "DC"] if args.regime == "both" else [args.regime]
    rows = []
    long_rows = []
    for meas in meas_files:
        cond = os.path.basename(meas)[len("meas_"):-len(".csv")]
        expr = os.path.join(args.prepared, f"expr_{cond}.csv")
        uptake = os.path.join(args.prepared, f"uptake_{cond}.csv")
        if not os.path.exists(expr):
            print(f"[WARN] no expr file for {cond}; skipping")
            continue
        for regime in regimes:
            up = uptake if (regime == "DC") else None
            if regime == "DC" and not os.path.exists(uptake):
                print(f"[WARN] DC requested but no uptake file for {cond}; skipping DC")
                continue
            r, n, rc, out, detail_rows = run_one(args.validation, args.model, args.objective,
                                    args.solver, expr, meas, up)
            if args.verbose or rc != 0 or r is None:
                print(f"\n----- {cond} [{regime}] (rc={rc}) -----\n{out}")
            print(f"  {cond:<22} {regime}  r={r if r is not None else 'NA':<7} (n={n})")
            rows.append({"condition": cond, "regime": regime, "n_reactions": n,
                         "uncentered_pearson": r, "returncode": rc})
            for det in detail_rows:
                long_rows.append({
                    "condition": cond, "regime": regime, "method": "reduced_engine",
                    "reaction_id": det.get("reaction_id"),
                    "measured_flux": det.get("measured_abs"),
                    "predicted_flux": det.get("predicted_abs"),
                    "status": "OK", "solve_mode": det.get("solve_mode"),
                })

    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    print("\n" + "=" * 60)
    print("BENCHMARK SUMMARY (uncentered Pearson)")
    print("=" * 60)
    targets = {"AC": 0.385, "DC": 0.725}
    for regime in regimes:
        sub = df[(df["regime"] == regime) & df["uncentered_pearson"].notna()]
        if len(sub):
            mean_r = sub["uncentered_pearson"].mean()
            tgt = targets[regime]
            verdict = ">=" if mean_r >= tgt else "<"
            print(f"  {regime}: mean r = {mean_r:.3f} over {len(sub)} conditions "
                  f"(Bhadra-Lobo 2020 target {tgt}; {verdict} target)")
        else:
            print(f"  {regime}: no successful runs")
    if long_rows:
        long_out = "benchmark_long_format_flux_reduced_engine.csv"
        pd.DataFrame(long_rows).to_csv(long_out, index=False)
        print(f"[OK] wrote {long_out} ({len(long_rows)} rows, "
              f"{pd.DataFrame(long_rows).condition.nunique()} conditions)")
    else:
        print("[WARN] no per-reaction detail captured -- benchmark_flux_recovery.csv "
              "was never found after any run_one() call; long-format file not written.")

    print(f"\n[OK] wrote {args.out}")
    print("Report the per-regime MEAN against the target; cite Bhadra-Lobo, Kim & Lun 2020 "
          "and note Gerosa-native mapping (band comparison, not bit-identical).")


if __name__ == "__main__":
    sys.exit(main())

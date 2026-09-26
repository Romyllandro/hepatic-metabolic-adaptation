#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prepare_benchmark_inputs.py
===========================
Convert the Bhadra-Lobo, Kim & Lun (2020) PLOS ONE 11:e0238689 Supporting-Information
tables into the per-condition CSVs that validation_01_benchmark_recovery.py consumes:

    expr_<organism>_<condition>.csv   ->  columns: gene_id, expression
    meas_<organism>_<condition>.csv   ->  columns: reaction_id, measured_flux
    uptake_<organism>_<condition>.csv ->  columns: exchange_id, lower_bound, upper_bound  (DC regime)

WHY A SEPARATE ADAPTER
----------------------
The benchmark ships its data as supplementary spreadsheets, not as model-ready CSVs.
Keeping the format conversion in one auditable, parameterised place (rather than hard-coded
inside the validation script) means the validation engine stays generic and the messy,
source-specific parsing is isolated and reviewable -- the same separation-of-concerns
principle you already use between the flux engine and the batch-correction step.

INPUT (recommended "tidy" intermediate)
---------------------------------------
After downloading the SI (see DOWNLOAD_BHADRA_LOBO_2020.md), export to two tidy CSVs:
    transcriptomes.csv : organism, condition, gene_id, expression
    measured_fluxes.csv: organism, condition, reaction_id, measured_flux
    (optional) uptakes.csv: organism, condition, exchange_id, lower_bound, upper_bound

This script then emits one input pair per (organism, condition) and prints the exact
validation_01 commands for both DC (with --uptake) and AC (without) regimes.

USAGE
-----
  python prepare_benchmark_inputs.py \
      --transcriptomes transcriptomes.csv \
      --measured_fluxes measured_fluxes.csv \
      --uptakes uptakes.csv \
      --organism "E. coli" --model iJO1366.json \
      --objective BIOMASS_Ec_iJO1366_core_53p95M \
      --outdir prepared/
"""
import argparse, os, sys
import pandas as pd


def _norm_cols(df):
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    return df


def emit_for_condition(sub, key_col, val_col, out_path, out_cols):
    t = sub[[key_col, val_col]].dropna()
    t.columns = out_cols
    t.to_csv(out_path, index=False)
    return len(t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--transcriptomes", required=True)
    ap.add_argument("--measured_fluxes", required=True)
    ap.add_argument("--uptakes", default=None)
    ap.add_argument("--organism", required=True, help="filter value matching the 'organism' column")
    ap.add_argument("--model", required=True)
    ap.add_argument("--objective", required=True)
    ap.add_argument("--outdir", default="prepared")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    tx = _norm_cols(pd.read_csv(args.transcriptomes))
    fx = _norm_cols(pd.read_csv(args.measured_fluxes))
    up = _norm_cols(pd.read_csv(args.uptakes)) if args.uptakes else None

    for df, need in [(tx, {"organism", "condition", "gene_id", "expression"}),
                     (fx, {"organism", "condition", "reaction_id", "measured_flux"})]:
        missing = need - set(df.columns)
        if missing:
            sys.exit(f"[ERROR] input missing columns {missing}. Got {list(df.columns)}. "
                     f"See header of this script for the expected tidy schema.")

    org = args.organism
    tx_o = tx[tx["organism"].astype(str).str.strip().str.lower() == org.strip().lower()]
    fx_o = fx[fx["organism"].astype(str).str.strip().str.lower() == org.strip().lower()]
    conditions = sorted(set(tx_o["condition"]) & set(fx_o["condition"]))
    if not conditions:
        sys.exit(f"[ERROR] no shared (organism={org}) conditions between transcriptome and flux tables.")

    print(f"[INFO] organism={org}: {len(conditions)} conditions -> {conditions}")
    cmds_dc, cmds_ac = [], []
    safe_org = org.replace(" ", "").replace(".", "")
    for cond in conditions:
        safe = f"{safe_org}_{str(cond).replace(' ', '')}"
        e_path = os.path.join(args.outdir, f"expr_{safe}.csv")
        m_path = os.path.join(args.outdir, f"meas_{safe}.csv")
        n_e = emit_for_condition(tx_o[tx_o["condition"] == cond], "gene_id", "expression",
                                 e_path, ["gene_id", "expression"])
        n_m = emit_for_condition(fx_o[fx_o["condition"] == cond], "reaction_id", "measured_flux",
                                 m_path, ["reaction_id", "measured_flux"])
        u_arg = ""
        if up is not None:
            up_o = up[(up["organism"].astype(str).str.strip().str.lower() == org.strip().lower())
                      & (up["condition"] == cond)]
            if len(up_o):
                u_path = os.path.join(args.outdir, f"uptake_{safe}.csv")
                up_o[["exchange_id", "lower_bound", "upper_bound"]].to_csv(u_path, index=False)
                u_arg = f" --uptake {u_path}"
        print(f"  {cond}: {n_e} genes, {n_m} measured reactions"
              + (" (+uptake)" if u_arg else ""))
        base = (f"python validation_01_benchmark_recovery.py --mode flux "
                f"--model {args.model} --objective {args.objective} "
                f"--expression {e_path} --measured {m_path}")
        cmds_dc.append(base + u_arg)        # DC regime (target ~0.725)
        cmds_ac.append(base)                # AC regime (target ~0.385)

    runner = os.path.join(args.outdir, "run_benchmark.sh")
    with open(runner, "w") as f:
        f.write("#!/usr/bin/env bash\nset -e\n")
        f.write("\n# === DC regime (uptake known; Bhadra-Lobo 2020 E-Flux2 mean ~0.725) ===\n")
        f.write("\n".join(cmds_dc) + "\n")
        f.write("\n# === AC regime (uptake unknown; E-Flux2 mean ~0.385) ===\n")
        f.write("\n".join(cmds_ac) + "\n")
    os.chmod(runner, 0o755)
    print(f"\n[OK] wrote per-condition inputs to {args.outdir}/ and runner {runner}")
    print("    Aggregate the per-condition uncentered-Pearson values and report the MEAN "
          "against the regime target (DC 0.725 / AC 0.385).")


if __name__ == "__main__":
    sys.exit(main())

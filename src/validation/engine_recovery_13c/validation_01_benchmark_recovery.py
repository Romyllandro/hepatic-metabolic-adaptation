#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validation_01_benchmark_recovery.py
===================================
PRONG 1: prove the ENGINE our pipeline is built on is implemented faithfully, by
reproducing published behaviour on a standard benchmark. This isolates the engine from
our diet-aware customisation.

SOLVER-STEP NOTE (revised)
--------------------------
The production pipeline does NOT use one uniform solver step, and this script must
mirror that rather than assume it:

    RQ1 / RQ2 (map_fixv5_multigroupsv8_layered_manuscript_run.py) : plain FBA, model.optimize()
    RQ3       (rq3_FINAL_COMPLETE_REVISED.py)                     : pFBA, pfba(model)
    RQ4 hep.  (rq4_hepatic_integration_CORRECTED_v13.py)          : pFBA, fraction_of_optimum=0.99
    RQ4 MICOM (rq4_microbiome_community_modeling_v2026.py)        : cooperative_tradeoff(pfba=False)

An earlier version of this script validated the engine under pFBA only. Because the
RQ1/RQ2 results are produced with plain FBA, a pFBA-only benchmark does not strictly
cover them. --solve_mode now defaults to "both", so the recovery statistic is reported
under each configuration and the manuscript can cite the one matching each analysis arm.

TWO MODES
---------
  --mode flux         Uncentered Pearson between predicted and 13C-MFA MEASURED fluxes.
                      Regime-aware:
                        * DC (carbon-source uptake KNOWN)  -> pass --uptake; Bhadra-Lobo 2020
                          E-Flux2 mean r ~ 0.725
                        * AC (uptake UNKNOWN)              -> omit --uptake; mean r ~ 0.385
                      Kim et al. 2016 reference band: r in [0.59, 0.87] (glucose, E. coli/yeast).

  --mode essentiality iMM1415 single-gene deletion vs an experimental ground-truth list.
                      Reports sensitivity / specificity / precision / MCC.

USAGE
-----
  # AC regime (uptake unknown)
  python validation_01_benchmark_recovery.py --mode flux \
      --model iJO1366.json --objective BIOMASS_Ec_iJO1366_core_53p95M \
      --expression prepared/expr_Ecoli_glucose.csv --measured prepared/meas_Ecoli_glucose.csv

  # DC regime (uptake known)
  python validation_01_benchmark_recovery.py --mode flux \
      --model iJO1366.json --objective BIOMASS_Ec_iJO1366_core_53p95M \
      --expression prepared/expr_Ecoli_glucose.csv --measured prepared/meas_Ecoli_glucose.csv \
      --uptake prepared/uptake_Ecoli_glucose.csv

INPUT FORMATS
-------------
  --expression : CSV [gene_id, expression]                       (single condition)
  --measured   : CSV [reaction_id, measured_flux]                (13C-MFA fluxes)
  --uptake     : CSV [exchange_id, lower_bound, upper_bound]     (DC regime)
  --ground_truth : CSV [gene_id, essential]                      (essentiality mode)
"""
import argparse, os, re, sys
import numpy as np
import pandas as pd

import cobra
from cobra.io import load_json_model
from cobra.flux_analysis import pfba
from cobra.flux_analysis import single_gene_deletion


# --------------------------------------------------------------------------- #
# Vanilla E-Flux2 engine -- the REDUCED form of the manuscript pipeline:
# expression -> reaction capacity via GPR (AND=min, OR=max), normalise, then FBA
# followed by a parsimonious (minimal-flux) selection step. Exchanges are left to
# the medium/diet layer, exactly as in the full pipeline.
# --------------------------------------------------------------------------- #
def gpr_reaction_expression(rxn, expr):
    rule = rxn.gene_reaction_rule
    if not rule:
        return None
    vals = []
    for or_part in re.sub(r"[()]", "", rule).split(" or "):
        and_vals = [expr.get(g.strip()) for g in or_part.split(" and ") if g.strip()]
        and_vals = [v for v in and_vals if v is not None]
        if and_vals:
            vals.append(min(and_vals))      # AND (complex) = min subunit
    return max(vals) if vals else None       # OR (isozyme) = max


def eflux2_solve(model, expr, q=0.95, floor=0.1, cap=1e6, capacity=1000.0,
                 constrain_exchanges=False, solve_mode="pfba"):
    """Return a flux solution under E-Flux2 capacity constraints.

    solve_mode:
      "pfba" -- maximise objective then minimise total flux (matches RQ3 / RQ4 hepatic)
      "fba"  -- plain model.optimize() (matches the RQ1 / RQ2 production pipeline)
    """
    pos = np.array([v for v in expr.values() if v > 0])
    denom = np.quantile(pos, q) if len(pos) else 1.0
    norm = {g: float(np.clip(v / max(denom, 1e-9), floor, cap)) * capacity
            for g, v in expr.items()}
    with model:
        for rxn in model.reactions:
            if (not constrain_exchanges) and (rxn.id.startswith("EX_") or rxn.boundary):
                continue
            e = gpr_reaction_expression(rxn, norm)
            if e is None:
                continue
            if rxn.lower_bound < 0:
                rxn.lower_bound = max(rxn.lower_bound, -e)
            rxn.upper_bound = min(rxn.upper_bound, e)
        if solve_mode == "pfba":
            model.slim_optimize()
            sol = pfba(model)
        elif solve_mode == "fba":
            sol = model.optimize()
        else:
            raise ValueError(f"unknown solve_mode: {solve_mode}")
    return sol


def uncentered_pearson(pred, meas):
    pred = np.asarray(pred, float)
    meas = np.asarray(meas, float)
    return float(np.dot(pred, meas) / (np.linalg.norm(pred) * np.linalg.norm(meas) + 1e-12))


def _resolve_exchange_id(model, ex_id):
    """Find an exchange id in the model, tolerating SBML<->BiGG naming differences:
    'EX_ac_LPAREN_e_RPAREN_' (libSBML) <-> 'EX_ac_e' (BiGG json)."""
    candidates = [ex_id]
    if "_LPAREN_e_RPAREN_" in ex_id:
        candidates.append(ex_id.replace("_LPAREN_e_RPAREN_", "_e"))
    if ex_id.endswith("_e"):
        candidates.append(ex_id[:-2] + "_LPAREN_e_RPAREN_")
    if ex_id.startswith("R_"):
        candidates.append(ex_id[2:])
    # BiGG often suffixes stereoisomers: EX_glc_e -> EX_glc__D_e ; also handle the _e base
    base = None
    for c in list(candidates):
        if c.endswith("_e"):
            base = c[:-2]
            candidates += [base + "__D_e", base + "__L_e"]
    for c in candidates:
        try:
            return model.reactions.get_by_id(c)
        except KeyError:
            continue
    return None


def apply_uptake(model, uptake_csv):
    """DC regime: set measured carbon-source / key-metabolite uptake (exchange) bounds.
    CSV columns: [exchange_id, lower_bound, upper_bound]. Uptake is a negative lower_bound."""
    up = pd.read_csv(uptake_csv)
    up.columns = [c.strip().lower() for c in up.columns]
    n, miss = 0, []
    for _, r in up.iterrows():
        rxn = _resolve_exchange_id(model, str(r["exchange_id"]))
        if rxn is None:
            miss.append(str(r["exchange_id"]))
            continue
        rxn.lower_bound = float(r["lower_bound"])
        rxn.upper_bound = float(r["upper_bound"])
        n += 1
    print(f"[INFO] DC regime: applied uptake bounds to {n} exchanges"
          + (f"; {len(miss)} not found ({miss})" if miss else ""))
    return model


def run_flux_mode(args):
    model = load_json_model(args.model)
    if args.objective:
        model.objective = args.objective
    regime = "AC (uptake unknown)"
    if args.uptake:                       # DC regime: known carbon-source uptake
        model = apply_uptake(model, args.uptake)
        regime = "DC (uptake known)"
    expr_df = pd.read_csv(args.expression)
    expr = dict(zip(expr_df.iloc[:, 0].astype(str), expr_df.iloc[:, 1].astype(float)))

    modes = ["fba", "pfba"] if args.solve_mode == "both" else [args.solve_mode]
    meas_df = pd.read_csv(args.measured)
    meas_df.columns = ["reaction_id", "measured_flux"]
    meas_df = meas_df.set_index("reaction_id")

    print("=" * 66)
    print("PRONG 1 (flux mode): engine vs 13C-MFA measured fluxes")
    print("=" * 66)
    print(f"  regime             : {regime}")

    summary, frames = [], []
    for mode in modes:
        sol = eflux2_solve(model, expr, q=args.eflux_quantile,
                           floor=args.eflux_floor, capacity=args.capacity,
                           solve_mode=mode)
        common = [r for r in meas_df.index if r in sol.fluxes.index]
        if len(common) < 3:
            print(f"[ERROR] Only {len(common)} reactions overlap measured set; check IDs.")
            sys.exit(1)
        pred = np.abs(sol.fluxes.loc[common].values)   # magnitudes (E-Flux2 convention)
        meas = np.abs(meas_df.loc[common, "measured_flux"].values)
        r_unc = uncentered_pearson(pred, meas)
        r_pear = float(np.corrcoef(pred, meas)[0, 1])
        ref = 0.725 if args.uptake else 0.385
        covers = {"fba": "RQ1 / RQ2 (production bulk + per-strain)",
                  "pfba": "RQ3 / RQ4 hepatic"}[mode]
        print(f"\n  --- solve_mode = {mode.upper()}  [covers: {covers}] ---")
        print(f"  reactions compared : {len(common)}")
        print(f"  uncentered Pearson : {r_unc:.3f}")
        print(f"  centered  Pearson  : {r_pear:.3f}")
        print(f"  Bhadra-Lobo 2020 E-Flux2 mean for this regime: {ref}")
        print(f"  Kim et al. 2016 reference band: [0.59, 0.87]")
        print(f"  => engine is {'>=' if r_unc >= ref else '<'} the published E-Flux2 mean")
        summary.append({"solve_mode": mode, "covers_analysis_arms": covers,
                        "n_reactions": len(common),
                        "uncentered_pearson": round(r_unc, 4),
                        "centered_pearson": round(r_pear, 4),
                        "reference_mean": ref, "regime": regime})
        frames.append(pd.DataFrame({"solve_mode": mode, "reaction_id": common,
                                    "predicted_abs": pred, "measured_abs": meas}))

    pd.concat(frames).to_csv("benchmark_flux_recovery.csv", index=False)
    pd.DataFrame(summary).to_csv("benchmark_flux_recovery_summary.csv", index=False)
    print("\n" + pd.DataFrame(summary).to_string(index=False))
    print("\n  [OK] wrote benchmark_flux_recovery.csv and benchmark_flux_recovery_summary.csv")
    if len(modes) > 1:
        print("  NOTE: cite the FBA row when describing RQ1/RQ2 and the pFBA row for RQ3/RQ4.")


def run_essentiality_mode(args):
    model = load_json_model(args.model)
    if args.objective:
        model.objective = args.objective
    wt = model.slim_optimize()
    thr = args.essential_fraction * wt
    print(f"[INFO] WT objective={wt:.4f}; essential if KO objective < {thr:.4f}")

    res = single_gene_deletion(model)
    res = res.reset_index(drop=True)
    res["gene"] = res["ids"].apply(lambda s: list(s)[0] if len(s) else None)
    res["pred_essential"] = (res["growth"].fillna(0) < thr).astype(int)
    pred = dict(zip(res["gene"].astype(str), res["pred_essential"]))

    gt = pd.read_csv(args.ground_truth)
    gt.columns = ["gene_id", "essential"]
    gt["gene_id"] = gt["gene_id"].astype(str)
    common = [g for g in gt["gene_id"] if g in pred]
    if not common:
        print("[ERROR] No gene-ID overlap between model and ground truth. Check ID namespace "
              "(model uses Entrez; convert ground truth accordingly).")
        sys.exit(1)
    y_true = gt.set_index("gene_id").loc[common, "essential"].astype(int).values
    y_pred = np.array([pred[g] for g in common])

    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())
    sens = tp / (tp + fn) if (tp + fn) else float("nan")
    spec = tn / (tn + fp) if (tn + fp) else float("nan")
    prec = tp / (tp + fp) if (tp + fp) else float("nan")
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / denom if denom else float("nan")

    print("=" * 66)
    print("PRONG 1 (essentiality mode): iMM1415 single-gene deletion vs in-vivo")
    print("=" * 66)
    print(f"  genes compared : {len(common)}   (TP={tp} TN={tn} FP={fp} FN={fn})")
    print(f"  sensitivity    : {sens:.3f}")
    print(f"  specificity    : {spec:.3f}")
    print(f"  precision      : {prec:.3f}")
    print(f"  MCC            : {mcc:.3f}")
    res[["gene", "growth", "pred_essential"]].to_csv("benchmark_essentiality.csv", index=False)
    print("  [OK] wrote benchmark_essentiality.csv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["flux", "essentiality"], required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--objective", default=None)
    ap.add_argument("--solver", default="gurobi")
    # flux mode
    ap.add_argument("--expression")
    ap.add_argument("--measured")
    ap.add_argument("--uptake", default=None,
                    help="DC regime: CSV [exchange_id,lower_bound,upper_bound] of measured uptake")
    ap.add_argument("--eflux_quantile", type=float, default=0.95)
    ap.add_argument("--eflux_floor", type=float, default=0.1)
    ap.add_argument("--capacity", type=float, default=1000.0)
    ap.add_argument("--solve_mode", choices=["fba", "pfba", "both"], default="both",
                    help="Solver step for the engine benchmark. 'fba' matches RQ1/RQ2, "
                         "'pfba' matches RQ3/RQ4 hepatic, 'both' (default) reports each.")
    # essentiality mode
    ap.add_argument("--ground_truth")
    ap.add_argument("--essential_fraction", type=float, default=0.10)
    args = ap.parse_args()

    try:
        cobra.Configuration().solver = args.solver
    except Exception as e:
        print(f"[WARN] solver '{args.solver}' unavailable ({e}); using default.")

    if args.mode == "flux":
        if not (args.expression and args.measured):
            ap.error("--mode flux requires --expression and --measured")
        run_flux_mode(args)
    else:
        if not args.ground_truth:
            ap.error("--mode essentiality requires --ground_truth")
        run_essentiality_mode(args)


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from cobra.io import load_json_model

from fva_refinement_common import (
    fva_difference,
    load_modeling_helpers,
    nanmean_expression,
    run_parsimonious_fva,
    run_standard_fva,
    safe_name,
)

DEFAULT_CONTRASTS = [
    ("HFD", "SCD"),
    ("KD", "SCD"),
    ("WD", "SCD"),
]


def parse_dataset_and_group(column: str):
    m = re.match(r"^(SCD|HFD|KD|WD|LFD)_(GSE\d+)_", str(column), flags=re.IGNORECASE)
    if not m:
        return None
    group = m.group(1).upper()
    if group == "LFD":
        group = "SCD"
    return m.group(2).upper(), group


def load_primary_stats(rq1_results: Path, treatment: str, control: str) -> pd.DataFrame:
    p = rq1_results / "comprehensive_analysis" / "csv_outputs" / f"reaction_stats_{treatment}_vs_{control}.csv"
    if not p.exists():
        raise FileNotFoundError(f"Primary reaction-statistics file not found: {p}")
    df = pd.read_csv(p)
    if "ReactionID" not in df.columns or "Significant" not in df.columns:
        raise ValueError(f"Unexpected reaction-statistics format: {p}")
    sig = df[df["Significant"].astype(str).str.lower().eq("true")].copy()
    return sig


def parse_eps(text: str) -> list[float]:
    vals = []
    for x in re.split(r"[,;\s]+", str(text).strip()):
        if not x:
            continue
        v = float(x)
        if v < 0:
            raise ValueError("Parsimony epsilon values must be >= 0")
        vals.append(v)
    return vals


def main():
    ap = argparse.ArgumentParser(
        description="Cohort-specific representative FVA for the multi-cohort RQ1 rerun."
    )
    ap.add_argument("--expression", default="GSEMERGED_SCD_HFD_KD_WD_gene_expression.csv")
    ap.add_argument("--model", default="iMM1415.json")
    ap.add_argument("--diet_bounds", default="expanded_diet_bounds_flat.json")
    ap.add_argument("--mapping", default="mouse_entrez_to_symbol.csv")
    ap.add_argument("--modeling_script", default="map_fixv5_multigroupsv8_layered_manuscript_rerun.py")
    ap.add_argument("--rq1_results", default="RERUN_RQ1_PFBA")
    ap.add_argument("--output", default="RERUN_RQ1_PFBA/cohort_representative_FVA")
    ap.add_argument("--objective_id", default="BIOMASS_mm_1_no_glygln")
    ap.add_argument("--objective_sense", default="max", choices=["max", "min"])
    ap.add_argument("--eflux_quantile", type=float, default=0.95)
    ap.add_argument("--eflux_floor", type=float, default=0.1)
    ap.add_argument("--eflux_cap", type=float, default=1000.0)
    ap.add_argument("--transporter_strategy", default="either", choices=["e_to_non_e", "regex", "either"])
    ap.add_argument("--solver", default="gurobi")
    ap.add_argument("--fva_fraction", type=float, default=1.0)
    ap.add_argument("--parsimony_eps", default="0.01,0.05",
                    help="Comma-separated near-pFBA total-flux tolerances, e.g. 0.01,0.05.")
    ap.add_argument("--tolerance", type=float, default=1e-7)
    ap.add_argument("--processes", type=int, default=1)
    ap.add_argument("--skip_standard", action="store_true")
    ap.add_argument("--skip_parsimonious", action="store_true")
    args = ap.parse_args()

    expr_path = Path(args.expression)
    model_path = Path(args.model)
    diet_path = Path(args.diet_bounds)
    mapping_path = Path(args.mapping)
    rq1_results = Path(args.rq1_results)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    (out / "per_model").mkdir(exist_ok=True)
    (out / "per_cohort_contrast").mkdir(exist_ok=True)

    helpers = load_modeling_helpers(args.modeling_script)

    # Primary significant targets for each diet-vs-SCD contrast.
    primary_stats = {}
    targets_by_contrast = {}
    for treatment, control in DEFAULT_CONTRASTS:
        contrast = f"{treatment}_vs_{control}"
        stats = load_primary_stats(rq1_results, treatment, control)
        primary_stats[contrast] = stats.set_index("ReactionID")
        targets_by_contrast[contrast] = set(stats["ReactionID"].astype(str))
        print(f"[TARGETS] {contrast}: {len(targets_by_contrast[contrast])} primary significant reactions")

    shared_core = set.intersection(*(targets_by_contrast[c] for c in targets_by_contrast))
    print(f"[TARGETS] Shared HFD/KD/WD-vs-SCD core: {len(shared_core)} reactions")
    pd.DataFrame({"ReactionID": sorted(shared_core)}).to_csv(out / "RQ1_shared_core_targets.csv", index=False)

    # Load merged expression and derive cohort/group membership from headers.
    rna = pd.read_csv(expr_path)
    if "Gene_Symbol" not in rna.columns:
        raise ValueError("Expression file must contain Gene_Symbol.")
    sample_meta = []
    for col in rna.columns:
        parsed = parse_dataset_and_group(col)
        if parsed:
            dataset, group = parsed
            sample_meta.append({"Sample": col, "Dataset": dataset, "Group": group})
    meta = pd.DataFrame(sample_meta)
    if meta.empty:
        raise ValueError("No sample columns matched the expected DIET_GSE####_* naming convention.")
    meta.to_csv(out / "sample_metadata_cohort_fva.csv", index=False)

    count_table = meta.groupby(["Dataset", "Group"]).size().unstack(fill_value=0)
    print("\n[DESIGN] Sample counts by dataset/group:")
    print(count_table.to_string())

    eligible = defaultdict(list)
    for treatment, control in DEFAULT_CONTRASTS:
        contrast = f"{treatment}_vs_{control}"
        for dataset, g in meta.groupby("Dataset"):
            groups = set(g["Group"])
            if treatment in groups and control in groups:
                eligible[contrast].append(dataset)
        print(f"[DESIGN] {contrast}: {eligible[contrast]}")

    expected = {
        "HFD_vs_SCD": {"GSE101657", "GSE159090", "GSE160646", "GSE188344", "GSE246221", "GSE248297"},
        "KD_vs_SCD": {"GSE101657", "GSE248297"},
        "WD_vs_SCD": {"GSE159090"},
    }
    for contrast, exp in expected.items():
        got = set(eligible[contrast])
        if got != exp:
            raise ValueError(f"RQ1 cohort preflight failed for {contrast}: expected {sorted(exp)}, got {sorted(got)}")

    # Gene IDs and base model exactly as in the primary rerun.
    symbol_map = helpers.load_symbol_to_entrez_mapping(str(mapping_path))
    gene_symbols = rna["Gene_Symbol"].astype(str).str.lower().tolist()
    gene_ids = [symbol_map.get(sym.strip().lower(), sym) for sym in gene_symbols]

    with diet_path.open("r", encoding="utf-8") as fh:
        raw = json.load(fh)
    diet_bounds = {
        helpers.canonical_code(k): {rid: [float(v[0]), float(v[1])] for rid, v in d.items()}
        for k, d in raw.items()
    }

    base = load_json_model(str(model_path))
    base = helpers.configure_solver(base, solver=args.solver, deterministic=True)
    base, obj_id = helpers.set_objective_reaction(
        base, objective_id=args.objective_id, objective_regex=None, sense=args.objective_sense
    )
    ex_set, trans_set, int_set = helpers.classify_reactions(
        base,
        transporter_strategy=args.transporter_strategy,
        transporter_regex=None,
        transporter_subsystem_regex=None,
        compartments_for_transport=("e",),
    )
    print(f"[MODEL] reactions={len(base.reactions)} EX={len(ex_set)} transporters={len(trans_set)} internal={len(int_set)} objective={obj_id}")

    def build_model(dataset: str, group: str):
        cols = meta.loc[(meta.Dataset == dataset) & (meta.Group == group), "Sample"].tolist()
        if not cols:
            raise ValueError(f"No samples for {dataset}/{group}")
        vec = nanmean_expression(rna, cols)
        mdl = base.copy()
        mdl = helpers.configure_solver(mdl, solver=args.solver, deterministic=True, logger=lambda *_: None)
        mdl, l1 = helpers.apply_diet_bounds_layer1(
            mdl, code=group, diet_bounds=diet_bounds, diet_units="model", mw_map={}, gDW=1.0, hours_per_day=24.0
        )
        mdl, l2, *_ = helpers.apply_expression_constraints_scoped(
            mdl, gene_ids, vec, trans_set,
            eflux_quantile=args.eflux_quantile, eflux_floor=args.eflux_floor,
            eflux_cap=args.eflux_cap, label="L2", symbol_to_entrez=symbol_map
        )
        mdl, l3, *_ = helpers.apply_expression_constraints_scoped(
            mdl, gene_ids, vec, int_set,
            eflux_quantile=args.eflux_quantile, eflux_floor=args.eflux_floor,
            eflux_cap=args.eflux_cap, label="L3", symbol_to_entrez=symbol_map
        )
        mdl = helpers.validate_model(mdl)
        psol = helpers.solve_flux(mdl, solve_mode="pfba", pfba_fraction=1.0)
        return mdl, {
            "Dataset": dataset,
            "Group": group,
            "N_replicates": len(cols),
            "BiologicalObjective": float(psol.objective_value),
            "pFBATotalFlux": float(getattr(psol, "total_flux", np.abs(psol.fluxes).sum())),
            "L1": int(l1), "L2": int(l2), "L3": int(l3),
        }

    # Determine which target reactions each unique representative model needs.
    model_targets = defaultdict(set)
    for treatment, control in DEFAULT_CONTRASTS:
        contrast = f"{treatment}_vs_{control}"
        for dataset in eligible[contrast]:
            model_targets[(dataset, treatment)] |= targets_by_contrast[contrast]
            model_targets[(dataset, control)] |= targets_by_contrast[contrast]

    eps_values = parse_eps(args.parsimony_eps)
    model_results = {}
    model_qc = []
    model_manifest = []

    for (dataset, group), target_set in sorted(model_targets.items()):
        targets = sorted(target_set)
        print(f"\n[MODEL] {dataset}/{group}: building representative model from {len(meta[(meta.Dataset==dataset)&(meta.Group==group)])} replicates; targets={len(targets)}")
        mdl, info = build_model(dataset, group)
        model_manifest.append(info)
        key = (dataset, group)
        model_results[key] = {"standard": None, "pfva": {}}

        if not args.skip_standard:
            fva, qc = run_standard_fva(
                mdl, targets, fraction=args.fva_fraction,
                processes=args.processes, tolerance=args.tolerance
            )
            fva.reset_index().to_csv(out / "per_model" / f"{dataset}_{group}_standard_FVA.csv", index=False)
            qc["Dataset"] = dataset; qc["Group"] = group; qc["Mode"] = "standard"
            model_qc.append(qc)
            model_results[key]["standard"] = fva

        if not args.skip_parsimonious:
            pfva_targets = sorted(set(targets) & shared_core)
            if pfva_targets:
                for eps in eps_values:
                    pfva, qc, pfmeta = run_parsimonious_fva(
                        mdl, pfva_targets, fraction=args.fva_fraction,
                        parsimony_epsilon=eps, processes=args.processes,
                        tolerance=args.tolerance
                    )
                    label = f"pfva_eps{eps:g}".replace(".", "p")
                    pfva.reset_index().to_csv(out / "per_model" / f"{dataset}_{group}_{label}.csv", index=False)
                    qc["Dataset"] = dataset; qc["Group"] = group; qc["Mode"] = label
                    model_qc.append(qc)
                    model_results[key]["pfva"][eps] = pfva
                    mrow = dict(info)
                    mrow.update({"Mode": label, **pfmeta})
                    pd.DataFrame([mrow]).to_csv(out / "per_model" / f"{dataset}_{group}_{label}_metadata.csv", index=False)

    pd.DataFrame(model_manifest).to_csv(out / "representative_model_manifest.csv", index=False)
    if model_qc:
        pd.concat(model_qc, ignore_index=True).to_csv(out / "FVA_numerical_QC.csv", index=False)

    # Cohort-level difference intervals and cross-cohort summaries.
    all_summaries = []
    for treatment, control in DEFAULT_CONTRASTS:
        contrast = f"{treatment}_vs_{control}"
        pstats = primary_stats[contrast]
        cohort_frames = []
        for dataset in eligible[contrast]:
            if not args.skip_standard:
                d = fva_difference(
                    model_results[(dataset, treatment)]["standard"],
                    model_results[(dataset, control)]["standard"],
                    tolerance=args.tolerance,
                ).reset_index()
                d["Dataset"] = dataset
                d["Mode"] = "standard"
                d["Contrast"] = contrast
                d["PrimaryDirection"] = d["ReactionID"].map(pstats["Direction"])
                d["MatchesPrimary"] = d["Robust"] & (d["RobustDirection"] == d["PrimaryDirection"])
                d.to_csv(out / "per_cohort_contrast" / f"{contrast}_{dataset}_standard.csv", index=False)
                cohort_frames.append(d)

            if not args.skip_parsimonious:
                for eps in eps_values:
                    tdf = model_results[(dataset, treatment)]["pfva"].get(eps)
                    cdf = model_results[(dataset, control)]["pfva"].get(eps)
                    if tdf is None or cdf is None:
                        continue
                    d = fva_difference(tdf, cdf, tolerance=args.tolerance).reset_index()
                    label = f"pfva_eps{eps:g}".replace(".", "p")
                    d["Dataset"] = dataset
                    d["Mode"] = label
                    d["Contrast"] = contrast
                    d["PrimaryDirection"] = d["ReactionID"].map(pstats["Direction"])
                    d["MatchesPrimary"] = d["Robust"] & (d["RobustDirection"] == d["PrimaryDirection"])
                    d.to_csv(out / "per_cohort_contrast" / f"{contrast}_{dataset}_{label}.csv", index=False)
                    cohort_frames.append(d)

        if not cohort_frames:
            continue
        long = pd.concat(cohort_frames, ignore_index=True)
        long.to_csv(out / f"{contrast}_cohort_FVA_long.csv", index=False)

        for mode, mdf in long.groupby("Mode"):
            rows = []
            for rid, rdf in mdf.groupby("ReactionID"):
                primary_dir = pstats.loc[rid, "Direction"] if rid in pstats.index else ""
                n = rdf["Dataset"].nunique()
                nrob = int(rdf["Robust"].sum())
                nmatch = int(rdf["MatchesPrimary"].sum())
                robust_dirs = sorted(set(rdf.loc[rdf["Robust"], "RobustDirection"].astype(str)))
                rows.append({
                    "Contrast": contrast,
                    "Mode": mode,
                    "ReactionID": rid,
                    "PrimaryDirection": primary_dir,
                    "N_cohorts": n,
                    "N_robust": nrob,
                    "N_match_primary": nmatch,
                    "AllCohortsRobustSameAsPrimary": bool(n > 0 and nmatch == n),
                    "AnyFVARobust": bool(nrob > 0),
                    "RobustDirectionsObserved": ";".join(robust_dirs),
                })
            sdf = pd.DataFrame(rows)
            sdf.to_csv(out / f"{contrast}_{mode}_cross_cohort_summary.csv", index=False)
            all_summaries.append(sdf)

    combined = pd.concat(all_summaries, ignore_index=True) if all_summaries else pd.DataFrame()
    if not combined.empty:
        combined.to_csv(out / "RQ1_cross_cohort_FVA_summary_all_modes.csv", index=False)
        high = combined.groupby(["Contrast", "Mode"]).agg(
            N_reactions=("ReactionID", "nunique"),
            N_all_cohorts_robust=("AllCohortsRobustSameAsPrimary", "sum"),
            N_any_robust=("AnyFVARobust", "sum"),
        ).reset_index()
        high.to_csv(out / "RQ1_cross_cohort_FVA_high_level_summary.csv", index=False)
        print("\n[SUMMARY]")
        print(high.to_string(index=False))

    run_report = {
        "expression": str(expr_path),
        "model": str(model_path),
        "objective": args.objective_id,
        "fva_fraction": args.fva_fraction,
        "parsimony_eps": eps_values,
        "numerical_tolerance": args.tolerance,
        "eligible_cohorts": dict(eligible),
        "targets_per_contrast": {k: len(v) for k, v in targets_by_contrast.items()},
        "shared_core_targets": len(shared_core),
        "interpretation": (
            "Standard FVA tests alternative primary-objective-optimal solutions. "
            "Parsimonious FVA additionally restricts total flux to within epsilon of the pFBA minimum. "
            "Neither replaces replicate-level RQ1 statistical inference."
        ),
    }
    (out / "run_report.json").write_text(json.dumps(run_report, indent=2), encoding="utf-8")
    print(f"\n[DONE] Cohort-specific RQ1 FVA written to: {out.resolve()}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
rq2_layer3_strain_specificity.py
=================================
RQ2 "Layer 3" statistical test: is one genetic-background strain's HFD-vs-SCD
flux response statistically different from the pooled response of the other
eight strains?

Background
----------
RQ2's union/core/strain-unique reaction classification (e.g. union=243,
core=11, strain-unique=80) is a DETERMINISTIC MAGNITUDE THRESHOLD
(|Diff(HFD-SCD)| >= --edge_threshold), not a statistical test. Two further,
complementary statistical layers were requested (SBBI lab meeting,
2026-09-22, Dr. Juan Cui):

  Layer 2 (within-strain diet effect; NOT computed here, see
  flux_pairwise_stats.csv from the replicate-level pipeline run): does HFD
  differ from SCD *within one strain*, given that strain's own replicates?

  Layer 3 (between-strain specificity; THIS SCRIPT): when a reaction is
  called strain-specific / strain-unique (i.e. "one of nine numbers differs
  from the other eight"), is that difference actually supported by a
  statistical test, at both levels of data granularity Dr. Cui described:

    3a. Aggregate-level test ("you may have one summarized number for each
        strain... that comparison still requires statistical evidence
        appropriate to those group-level values"): treat each strain's
        Diff(HFD-SCD) as a single point estimate (9 numbers total) and test
        whether the target strain's number is an outlier relative to the
        other eight, via a one-sample t-test (t-distribution, df = n_other-1).

    3b. Replicate-level test ("at the individual level... assess ...
        between-group differences using the replicate-level distributions"):
        use the actual per-replicate flux values (3 HFD + 3 SCD per strain)
        for the target strain and the pooled other eight strains, and fit
        flux ~ C(diet) * C(is_target_strain) by OLS with HC3 (heteroscedasticity-
        robust) standard errors, matching the HC3-robust-OLS approach already
        used elsewhere in this manuscript's RQ1 statistics. The interaction
        term's p-value is the between-strain specificity test: it asks
        whether the diet effect is different for the target strain than for
        the pooled rest, using genuine replicate-level variance rather than
        collapsing each strain to one number.

Both tests are run for three reaction-set scenarios:
  - union243      : every reaction flagged in >=1 strain; tested against
                     every strain in which it was flagged.
  - strainunique80 : every reaction flagged in EXACTLY one strain; tested
                     against that one strain.
  - core11         : every reaction flagged in ALL nine strains; tested
                     against every one of the nine strains (a sanity check —
                     true conservation should show few/no significant
                     interactions).

Union/core/strain-unique are recomputed directly from the aggregate-run flux
files inside this script (not read from a cached list), so this script is
fully self-contained and reproducible from the raw pipeline outputs alone.

Inputs
------
--aggregate_dir   Directory containing results_<strain>_GSE182668/flux_analysis/
                   reaction_flux_comparison_extended.csv for the AGGREGATE
                   (mean-collapsed, one flux value per group) pipeline run,
                   e.g. 972026_Step_2_RQ2_aggregated/
--replicate_dir   Same layout, for the REPLICATE-PRESERVING pipeline run
                   (3 HFD + 3 SCD per-sample flux columns per reaction),
                   e.g. 972026_Step_2_RQ2/
--strains         Comma-separated strain list (default: the 9 strains used
                   throughout RQ2).
--edge_threshold  Magnitude threshold for the union/core/strain-unique
                   classification (default 0.2, matching Table S3/S4).
--output_dir      Where to write the per-scenario result CSVs.

Output (per scenario, written to --output_dir)
-----------------------------------------------
<scenario>_layer3_results.csv with columns:
  reaction, target_strain, n_target_strains, n_other_strains,
  agg_target_value, agg_other_mean, agg_other_sd, agg_tstat, agg_pvalue, agg_fdr,
  rep_interaction_coef, rep_tstat, rep_pvalue, rep_fdr, rep_n_target_obs, rep_n_other_obs

layer3_scenario_summary.csv: one row per scenario x test-level, with counts of
raw-p<0.05 and FDR<0.05 tests, for a quick top-line read.

Example
-------
python rq2_layer3_strain_specificity.py \\
    --aggregate_dir 972026_Step_2_RQ2_aggregated \\
    --replicate_dir 972026_Step_2_RQ2 \\
    --output_dir reviewer_results/rq2_layer3_strain_specificity
"""
from __future__ import annotations
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

DEFAULT_STRAINS = [
    "129S1SvImJ", "AJ", "C57BL6J", "CASTEiJ", "DBA2J",
    "NODShiLtJ", "NZOHlLtJ", "PWKPhJ", "WSBEiJ",
]

FLUX_COMPARISON_REL = "results_{strain}_GSE182668/flux_analysis/reaction_flux_comparison_extended.csv"


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------

def load_aggregate_table(aggregate_dir: Path, strains: list[str]) -> pd.DataFrame:
    """Long table: reaction, strain, diff (aggregate Diff(HFD-SCD), one value/strain)."""
    rows = []
    for strain in strains:
        f = aggregate_dir / FLUX_COMPARISON_REL.format(strain=strain)
        df = pd.read_csv(f)
        diff_col = [c for c in df.columns if c.startswith("Diff")][0]
        for _, r in df.iterrows():
            rows.append({"reaction": r["ReactionID"], "strain": strain, "diff": r[diff_col]})
    return pd.DataFrame(rows)


def load_replicate_long_table(replicate_dir: Path, strains: list[str]) -> pd.DataFrame:
    """Long table: reaction, strain, diet, flux (one row per replicate observation)."""
    rows = []
    for strain in strains:
        f = replicate_dir / FLUX_COMPARISON_REL.format(strain=strain)
        df = pd.read_csv(f)
        hfd_cols = [c for c in df.columns if c.startswith("HFD_") and c.endswith("_Flux")
                    and c not in ("HFD_Flux",)]
        scd_cols = [c for c in df.columns if c.startswith("SCD_") and c.endswith("_Flux")
                    and c not in ("SCD_Flux",)]
        if not hfd_cols or not scd_cols:
            raise ValueError(
                f"{f}: could not find per-replicate HFD_*_Flux / SCD_*_Flux columns "
                f"(found columns: {list(df.columns)}). Is this an aggregate-run file "
                f"passed as --replicate_dir by mistake?"
            )
        for _, r in df.iterrows():
            rxn = r["ReactionID"]
            for c in hfd_cols:
                rows.append({"reaction": rxn, "strain": strain, "diet": "HFD", "flux": r[c]})
            for c in scd_cols:
                rows.append({"reaction": rxn, "strain": strain, "diet": "SCD", "flux": r[c]})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Reaction-set scenarios (recomputed from the aggregate table; not cached)
# --------------------------------------------------------------------------

def build_flag_matrix(agg_long: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """reaction x strain boolean matrix: |Diff(HFD-SCD)| >= threshold."""
    wide = agg_long.pivot(index="reaction", columns="strain", values="diff")
    return wide.abs() >= threshold


def scenario_union243(flags: pd.DataFrame) -> list[tuple[str, str]]:
    """(reaction, target_strain) for every reaction, every strain where it was flagged."""
    pairs = []
    for rxn, row in flags.iterrows():
        for strain in flags.columns:
            if bool(row[strain]):
                pairs.append((rxn, strain))
    return pairs


def scenario_strainunique(flags: pd.DataFrame) -> list[tuple[str, str]]:
    """(reaction, its single flagged strain) for reactions flagged in exactly one strain."""
    pairs = []
    for rxn, row in flags.iterrows():
        flagged = [s for s in flags.columns if bool(row[s])]
        if len(flagged) == 1:
            pairs.append((rxn, flagged[0]))
    return pairs


def scenario_core11(flags: pd.DataFrame) -> list[tuple[str, str]]:
    """(reaction, strain) for every strain, for reactions flagged in ALL strains (sanity check)."""
    pairs = []
    n_strains = len(flags.columns)
    for rxn, row in flags.iterrows():
        if row.sum() == n_strains:
            for strain in flags.columns:
                pairs.append((rxn, strain))
    return pairs


# --------------------------------------------------------------------------
# Layer 3a: aggregate-level "9 numbers, is one an outlier" test
# --------------------------------------------------------------------------

def aggregate_level_test(agg_wide: pd.Series, target_strain: str) -> dict:
    """agg_wide: Series indexed by strain, one Diff(HFD-SCD) value per strain, for one reaction."""
    target_val = agg_wide[target_strain]
    other_vals = agg_wide.drop(index=target_strain).values.astype(float)
    n_other = len(other_vals)
    other_mean = float(np.mean(other_vals))
    other_sd = float(np.std(other_vals, ddof=1)) if n_other > 1 else np.nan

    if other_sd is None or other_sd == 0 or np.isnan(other_sd):
        tstat, pval = np.nan, np.nan
    else:
        # One-sample t-test: is target_val consistent with the distribution
        # of the other n_other strain-level values?
        se = other_sd * np.sqrt(1 + 1.0 / n_other)
        tstat = (target_val - other_mean) / se
        pval = 2 * stats.t.sf(abs(tstat), df=n_other - 1)

    return {
        "agg_target_value": target_val,
        "agg_other_mean": other_mean,
        "agg_other_sd": other_sd,
        "agg_tstat": tstat,
        "agg_pvalue": pval,
        "n_other_strains": n_other,
    }


# --------------------------------------------------------------------------
# Layer 3b: replicate-level diet x strain-group interaction test (HC3 OLS)
# --------------------------------------------------------------------------

def replicate_level_test(rep_long_reaction: pd.DataFrame, target_strain: str) -> dict:
    """rep_long_reaction: rows for ONE reaction, columns strain/diet/flux, all strains."""
    d = rep_long_reaction.copy()
    d["is_target"] = (d["strain"] == target_strain).astype(int)
    d["diet_bin"] = (d["diet"] == "HFD").astype(int)

    n_target = int((d["is_target"] == 1).sum())
    n_other = int((d["is_target"] == 0).sum())

    if d["flux"].std(ddof=0) == 0:
        # No variance at all (e.g. reaction carries zero flux everywhere) -> undefined test
        return {
            "rep_interaction_coef": 0.0, "rep_tstat": np.nan, "rep_pvalue": np.nan,
            "rep_n_target_obs": n_target, "rep_n_other_obs": n_other,
        }

    try:
        model = smf.ols("flux ~ diet_bin * is_target", data=d).fit(cov_type="HC3")
        term = "diet_bin:is_target"
        coef = float(model.params.get(term, np.nan))
        tstat = float(model.tvalues.get(term, np.nan))
        pval = float(model.pvalues.get(term, np.nan))
    except Exception:
        coef, tstat, pval = np.nan, np.nan, np.nan

    return {
        "rep_interaction_coef": coef, "rep_tstat": tstat, "rep_pvalue": pval,
        "rep_n_target_obs": n_target, "rep_n_other_obs": n_other,
    }


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def run_scenario(name: str, pairs: list[tuple[str, str]], agg_wide_df: pd.DataFrame,
                  rep_long: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    rep_by_reaction = {rxn: g for rxn, g in rep_long.groupby("reaction")}
    rows = []
    for rxn, target_strain in pairs:
        row = {"reaction": rxn, "target_strain": target_strain}
        if rxn in agg_wide_df.index:
            row.update(aggregate_level_test(agg_wide_df.loc[rxn], target_strain))
        rep_g = rep_by_reaction.get(rxn)
        if rep_g is not None:
            row.update(replicate_level_test(rep_g, target_strain))
        rows.append(row)

    out = pd.DataFrame(rows)
    if len(out):
        for col, fdr_col in [("agg_pvalue", "agg_fdr"), ("rep_pvalue", "rep_fdr")]:
            valid = out[col].notna()
            out[fdr_col] = np.nan
            if valid.sum() > 0:
                out.loc[valid, fdr_col] = multipletests(out.loc[valid, col], method="fdr_bh")[1]

    out_path = output_dir / f"{name}_layer3_results.csv"
    out.to_csv(out_path, index=False)
    print(f"[{name}] n_tests={len(out)}  -> {out_path}")
    return out


def summarize(scenario_frames: dict[str, pd.DataFrame], output_dir: Path) -> None:
    rows = []
    for name, df in scenario_frames.items():
        for level, pcol, fdrcol in [("aggregate", "agg_pvalue", "agg_fdr"),
                                     ("replicate", "rep_pvalue", "rep_fdr")]:
            valid = df[pcol].notna()
            n = int(valid.sum())
            n_p05 = int((df.loc[valid, pcol] < 0.05).sum())
            n_fdr05 = int((df.loc[valid, fdrcol] < 0.05).sum())
            rows.append({
                "scenario": name, "test_level": level, "n_tests_valid": n,
                "n_raw_p_lt_0.05": n_p05, "n_fdr_lt_0.05": n_fdr05,
                "pct_raw_p_lt_0.05": round(100 * n_p05 / n, 1) if n else np.nan,
                "pct_fdr_lt_0.05": round(100 * n_fdr05 / n, 1) if n else np.nan,
            })
    summary = pd.DataFrame(rows)
    summary.to_csv(output_dir / "layer3_scenario_summary.csv", index=False)
    print()
    print(summary.to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aggregate_dir", required=True, type=Path)
    ap.add_argument("--replicate_dir", required=True, type=Path)
    ap.add_argument("--strains", default=",".join(DEFAULT_STRAINS))
    ap.add_argument("--edge_threshold", type=float, default=0.2)
    ap.add_argument("--output_dir", required=True, type=Path)
    args = ap.parse_args()

    strains = [s.strip() for s in args.strains.split(",") if s.strip()]
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading aggregate-run flux comparisons for {len(strains)} strains from {args.aggregate_dir} ...")
    agg_long = load_aggregate_table(args.aggregate_dir, strains)
    agg_wide = agg_long.pivot(index="reaction", columns="strain", values="diff")

    print(f"Loading replicate-run flux comparisons for {len(strains)} strains from {args.replicate_dir} ...")
    rep_long = load_replicate_long_table(args.replicate_dir, strains)

    flags = build_flag_matrix(agg_long, args.edge_threshold)
    union = flags.index[flags.any(axis=1)]
    core = flags.index[flags.sum(axis=1) == len(strains)]
    unique_pairs = scenario_strainunique(flags)
    print(f"Recomputed from aggregate run (|Diff|>={args.edge_threshold}): "
          f"union={len(union)}, core={len(core)}, strain-unique={len(unique_pairs)}")

    scenarios = {
        "union243": scenario_union243(flags),
        "strainunique80": unique_pairs,
        "core11": scenario_core11(flags),
    }

    frames = {}
    for name, pairs in scenarios.items():
        frames[name] = run_scenario(name, pairs, agg_wide, rep_long, args.output_dir)

    summarize(frames, args.output_dir)
    print(f"\n[OK] all outputs written to {args.output_dir}")


if __name__ == "__main__":
    main()

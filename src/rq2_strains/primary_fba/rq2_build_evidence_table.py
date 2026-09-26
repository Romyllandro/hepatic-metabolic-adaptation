#!/usr/bin/env python3
"""
rq2_build_evidence_table.py
============================
Builds the "Table S3/S4 addition" Dr. Cui asked for (SBBI lab meeting,
2026-09-22): keep the existing aggregate union/core/strain-unique
classification (magnitude-threshold table) AS-IS, and add columns carrying
the individual-replicate statistical evidence *alongside* it, rather than
replacing it.

For every one of the 243 union reactions this adds:
  - the magnitude-threshold classification (how many of 9 strains flag it;
    Core-11 / Strain-unique / Shared-partial; which strain(s) flag it)
  - Layer 2 (within-strain HFD-vs-SCD, from flux_pairwise_stats.csv):
    among the 8 estimable strains (NZOHlLtJ handled separately - see caveat
    in rq2-replicate-vs-aggregate-effect-2026-09-23.md), how many reach
    FDR<0.05 on their OWN replicates, and the strongest (min) FDR among them
  - Layer 3 (between-strain specificity, from rq2_layer3_strain_specificity.py
    output): the aggregate-level (9-number outlier) and replicate-level
    (diet x is_target_strain interaction, HC3-robust OLS) test results for
    the strain(s) that flag the reaction, plus a compact
    "supported_both_levels" flag

Nothing here recomputes Layer 3; it merges the already-produced
union243_layer3_results.csv / strainunique80_layer3_results.csv /
core11_layer3_results.csv with a freshly-computed Layer 2 summary and the
magnitude-threshold classification.

Usage
-----
python rq2_build_evidence_table.py \
    --aggregate_dir 972026_Step_2_RQ2_aggregated \
    --replicate_dir 972026_Step_2_RQ2 \
    --layer3_dir reviewer_results/rq2_layer3_strain_specificity \
    --output_dir reviewer_results/rq2_evidence_table
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rq2_layer3_strain_specificity import (
    DEFAULT_STRAINS, load_aggregate_table, build_flag_matrix,
)

L2_STATS_REL = "results_{strain}_GSE182668/stats_comparison/flux_pairwise_stats.csv"
NZO = "NZOHlLtJ"


def load_layer2(replicate_dir: Path, strains: list[str]) -> pd.DataFrame:
    """Long table: reaction, strain, FDR_BH (HFD vs SCD, within-strain, own replicates)."""
    rows = []
    rxn_names = {}
    for strain in strains:
        f = replicate_dir / L2_STATS_REL.format(strain=strain)
        df = pd.read_csv(f)
        df = df[(df["GroupA"] == "HFD") & (df["GroupB"] == "SCD")]
        for _, r in df.iterrows():
            rows.append({"reaction": r["ReactionID"], "strain": strain, "fdr": r["FDR_BH"]})
            if r["ReactionID"] not in rxn_names:
                rxn_names[r["ReactionID"]] = r["ReactionName"]
    return pd.DataFrame(rows), rxn_names


def summarize_layer2(l2_long: pd.DataFrame, reactions: list[str]) -> pd.DataFrame:
    out = []
    for rxn in reactions:
        sub = l2_long[l2_long["reaction"] == rxn]
        core8 = sub[sub["strain"] != NZO]
        n_tested8 = int(core8["fdr"].notna().sum())
        n_sig8 = int((core8["fdr"] < 0.05).sum())
        min_fdr8 = core8["fdr"].min() if n_tested8 else np.nan
        nzo_row = sub[sub["strain"] == NZO]
        nzo_fdr = nzo_row["fdr"].values[0] if len(nzo_row) else np.nan
        out.append({
            "reaction": rxn,
            "L2_n_strains_tested_of8": n_tested8,
            "L2_n_strains_FDR_lt_0.05_of8": n_sig8,
            "L2_min_FDR_of8": min_fdr8,
            "L2_NZOHlLtJ_FDR_excluded_nonestimable": nzo_fdr,
        })
    return pd.DataFrame(out).set_index("reaction")


def classify(flags: pd.DataFrame) -> pd.DataFrame:
    n_strains = len(flags.columns)
    rows = []
    for rxn, row in flags.iterrows():
        flagged = [s for s in flags.columns if bool(row[s])]
        n = len(flagged)
        if n == 0:
            continue
        if n == n_strains:
            cls = "Core-11 (all 9 strains)"
        elif n == 1:
            cls = "Strain-unique (1 of 9)"
        else:
            cls = f"Shared-partial ({n} of 9)"
        rows.append({
            "reaction": rxn,
            "N_strains_flagged_of9": n,
            "Classification": cls,
            "Flagged_strains": ";".join(flagged),
        })
    return pd.DataFrame(rows).set_index("reaction")


def attach_layer3_unique(class_df: pd.DataFrame, unique_l3: pd.DataFrame) -> pd.DataFrame:
    """For strain-unique reactions: single target-strain Layer 3 row, direct merge."""
    u = unique_l3.set_index("reaction")[
        ["target_strain", "agg_pvalue", "agg_fdr", "rep_pvalue", "rep_fdr"]
    ].rename(columns={
        "target_strain": "L3_target_strain",
        "agg_pvalue": "L3_agg_pvalue", "agg_fdr": "L3_agg_FDR",
        "rep_pvalue": "L3_rep_pvalue", "rep_fdr": "L3_rep_FDR",
    })
    return u


def attach_layer3_shared(class_df: pd.DataFrame, union_l3: pd.DataFrame, shared_reactions: list[str]) -> pd.DataFrame:
    """For shared-partial reactions (2-8 strains flagged): report the flagged strain with
    the STRONGEST (min) replicate-level FDR as the best-supported 'stands out most' case."""
    sub = union_l3[union_l3["reaction"].isin(shared_reactions)].copy()
    sub = sub.sort_values("rep_fdr", na_position="last")
    best = sub.groupby("reaction", as_index=True).first()
    return best[["target_strain", "agg_pvalue", "agg_fdr", "rep_pvalue", "rep_fdr"]].rename(columns={
        "target_strain": "L3_target_strain",
        "agg_pvalue": "L3_agg_pvalue", "agg_fdr": "L3_agg_FDR",
        "rep_pvalue": "L3_rep_pvalue", "rep_fdr": "L3_rep_FDR",
    })


def attach_layer3_core(class_df: pd.DataFrame, core_l3: pd.DataFrame) -> pd.DataFrame:
    """For core-11 reactions (sanity check): report the MOST significant of the 9 per-strain
    outlier tests. True conservation should show none/few reaching FDR<0.05."""
    sub = core_l3.copy().sort_values("rep_fdr", na_position="last")
    best = sub.groupby("reaction", as_index=True).first()
    out = best[["target_strain", "agg_pvalue", "agg_fdr", "rep_pvalue", "rep_fdr"]].rename(columns={
        "target_strain": "L3_target_strain",
        "agg_pvalue": "L3_agg_pvalue", "agg_fdr": "L3_agg_FDR",
        "rep_pvalue": "L3_rep_pvalue", "rep_fdr": "L3_rep_FDR",
    })
    out["L3_note"] = "core-11 sanity check: most-significant of 9 per-strain outlier tests shown; low significance expected"
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aggregate_dir", required=True, type=Path)
    ap.add_argument("--replicate_dir", required=True, type=Path)
    ap.add_argument("--layer3_dir", required=True, type=Path)
    ap.add_argument("--strains", default=",".join(DEFAULT_STRAINS))
    ap.add_argument("--edge_threshold", type=float, default=0.2)
    ap.add_argument("--output_dir", required=True, type=Path)
    args = ap.parse_args()

    strains = [s.strip() for s in args.strains.split(",") if s.strip()]
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading aggregate classification (magnitude threshold)...")
    agg_long = load_aggregate_table(args.aggregate_dir, strains)
    flags = build_flag_matrix(agg_long, args.edge_threshold)
    class_df = classify(flags)
    print(f"  union={len(class_df)}, core={sum(class_df.Classification.str.startswith('Core'))}, "
          f"unique={sum(class_df.Classification.str.startswith('Strain-unique'))}")

    print("Loading Layer 2 (within-strain, own replicates)...")
    l2_long, rxn_names = load_layer2(args.replicate_dir, strains)
    l2_summary = summarize_layer2(l2_long, list(class_df.index))

    print("Loading Layer 3 (between-strain specificity, pre-computed)...")
    union_l3 = pd.read_csv(args.layer3_dir / "union243_layer3_results.csv")
    unique_l3 = pd.read_csv(args.layer3_dir / "strainunique80_layer3_results.csv")
    core_l3 = pd.read_csv(args.layer3_dir / "core11_layer3_results.csv")

    is_unique = class_df["Classification"] == "Strain-unique (1 of 9)"
    is_core = class_df["Classification"] == "Core-11 (all 9 strains)"
    is_shared = ~is_unique & ~is_core

    l3_unique = attach_layer3_unique(class_df, unique_l3)
    l3_shared = attach_layer3_shared(class_df, union_l3, list(class_df.index[is_shared]))
    l3_core = attach_layer3_core(class_df, core_l3)

    l3_all = pd.concat([l3_unique, l3_shared, l3_core], axis=0)

    final = class_df.join(l2_summary, how="left").join(l3_all, how="left")
    final.insert(0, "ReactionName", [rxn_names.get(r, "") for r in final.index])
    final = final.reset_index().rename(columns={"reaction": "ReactionID"})

    final["L3_supported_both_levels_FDR_lt_0.05"] = (
        (final["L3_agg_FDR"] < 0.05) & (final["L3_rep_FDR"] < 0.05)
    )

    col_order = [
        "ReactionID", "ReactionName", "N_strains_flagged_of9", "Classification", "Flagged_strains",
        "L2_n_strains_tested_of8", "L2_n_strains_FDR_lt_0.05_of8", "L2_min_FDR_of8",
        "L2_NZOHlLtJ_FDR_excluded_nonestimable",
        "L3_target_strain", "L3_agg_pvalue", "L3_agg_FDR", "L3_rep_pvalue", "L3_rep_FDR",
        "L3_supported_both_levels_FDR_lt_0.05", "L3_note",
    ]
    for c in col_order:
        if c not in final.columns:
            final[c] = np.nan
    final = final[col_order].sort_values(
        ["Classification", "ReactionID"], ascending=[True, True]
    )

    out_path = args.output_dir / "TableS3_S4_addition_RQ2_evidence_layers.csv"
    final.to_csv(out_path, index=False)
    print(f"\n[OK] wrote {out_path}  ({len(final)} rows)")

    # Top-line counts for the write-up
    n_unique = int(is_unique.sum())
    n_unique_both = int(final.loc[final.Classification == "Strain-unique (1 of 9)",
                                   "L3_supported_both_levels_FDR_lt_0.05"].sum())
    n_unique_rep_only = int((
        (final.Classification == "Strain-unique (1 of 9)") &
        (final["L3_rep_FDR"] < 0.05)
    ).sum())
    n_core_flagged = int((
        (final.Classification == "Core-11 (all 9 strains)") &
        (final["L3_rep_FDR"] < 0.05)
    ).sum())
    print(f"\nStrain-unique (n={n_unique}): L3 replicate-level FDR<0.05 -> {n_unique_rep_only}; "
          f"supported at BOTH agg+rep FDR<0.05 -> {n_unique_both}")
    print(f"Core-11 sanity check: reactions where the best of 9 outlier tests still reaches "
          f"replicate-level FDR<0.05 -> {n_core_flagged} of 11")


if __name__ == "__main__":
    main()

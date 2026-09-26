#!/usr/bin/env python3
"""
run_cross_rq_final.py
=====================
Reproduce the final cross-RQ integration locally from the frozen RQ1-RQ4 outputs.

Inputs
------
RQ1:
  RERUN_RQ1_PFBA.zip (or extracted directory)
  must contain:
    reaction_stats_HFD_vs_SCD.csv
    reaction_stats_KD_vs_SCD.csv
    reaction_stats_WD_vs_SCD.csv

RQ2:
  RERUN_RQ2_PFBA_FVA.zip (or extracted directory)
  must contain:
    RQ2_cross_strain_reaction_summary.csv

RQ3:
  selected/frozen RQ3 setting directory, e.g.
    reviewer_results/RQ3_sensitivity/q99_BIOMASS_mm_1_no_glygln
  must contain:
    statistics/statistical_tests.csv

RQ4:
  frozen RQ4 v3 directory, e.g.
    reviewer_results/RQ4_final_full_v3
  must contain:
    host_scenarios/primary/condition_DD_HFD/DD_HFD_flux_comparison.csv
    host_scenarios/primary/attribution/flux_attribution_analysis.csv

Outputs
-------
cross_RQ_reaction_tracing_final.csv
cross_RQ_reaction_membership_matrix.csv
cross_RQ_pathway_matrix_final.csv
cross_RQ_summary_final.json
RQ4_primary_threshold_attributable_reactions.csv
RQ4_primary_microbiome_dominant_reactions.csv
RQ4_primary_pathway_attribution.csv
run_manifest.json

Scientific conventions
----------------------
- RQ1 HFD anchor:
    HFD-vs-SCD Significant == True AND |AdjustedMeanDiff| >= --anchor_abs_diff
  Default anchor threshold = 0.20.

- RQ1 shared all-diet signature:
    intersection of significant HFD-vs-SCD, KD-vs-SCD, WD-vs-SCD reactions.

- RQ2 universal core:
    MagnitudeUniversal == True.

- RQ3:
    threshold-responsive reactions from the frozen P99+biomass setting.
    RQ3 is treated as an ORTHOGONAL WesternDiet-vs-Chow cellular-context result,
    not as condition-matched HFD validation.

- RQ4 thresholded microbiome attribution:
    DD_HFD_flux_comparison.csv microbiome_attributable == True.
    In the final v3 pipeline this corresponds to |microbiome-induced delta flux|
    >= 0.01 by default.

- RQ4 microbiome-dominant:
    dominant_driver starts with "Microbiome" in flux_attribution_analysis.csv.
    This classification can include reactions below the |delta v| threshold, so
    threshold-attributable and microbiome-dominant are kept as separate concepts.

No COBRApy/Gurobi optimization is performed here. This is a downstream integration
script only.
"""
from __future__ import annotations

import argparse
import io
import json
import zipfile
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().isin(
        ["true", "1", "yes", "y", "t"]
    )


class CsvSource:
    """Read CSVs from either a directory tree or ZIP archive."""
    def __init__(self, path: str):
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(self.path)
        self._zip: Optional[zipfile.ZipFile] = None
        if self.path.is_file() and self.path.suffix.lower() == ".zip":
            self._zip = zipfile.ZipFile(self.path)

    def read_csv(self, basename: str, **kwargs) -> pd.DataFrame:
        if self._zip is not None:
            matches = [n for n in self._zip.namelist() if n.endswith(basename)]
            if not matches:
                raise FileNotFoundError(
                    f"{basename} not found inside {self.path}"
                )
            # Prefer shortest path if duplicate archived copies exist.
            name = sorted(matches, key=lambda x: (len(x), x))[0]
            return pd.read_csv(io.BytesIO(self._zip.read(name)), **kwargs)

        matches = list(self.path.rglob(basename))
        if not matches:
            raise FileNotFoundError(
                f"{basename} not found below directory {self.path}"
            )
        match = sorted(matches, key=lambda x: (len(str(x)), str(x)))[0]
        return pd.read_csv(match, **kwargs)

    def close(self):
        if self._zip is not None:
            self._zip.close()


def require_columns(df: pd.DataFrame, cols, label: str):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{label} missing required columns: {missing}")


def safe_get_annotation_table(hfd: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in ["ReactionID", "ReactionName", "Subsystem"] if c in hfd.columns]
    ann = hfd[cols].drop_duplicates("ReactionID").copy()
    if "ReactionName" not in ann.columns:
        ann["ReactionName"] = ""
    if "Subsystem" not in ann.columns:
        ann["Subsystem"] = "Unassigned"
    return ann.set_index("ReactionID")


def ensure_outdir(path: str) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def mapped_value(series: pd.Series, mapping: pd.Series, default=np.nan):
    x = series.map(mapping)
    if isinstance(default, str):
        return x.fillna(default)
    return x


# ---------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Final frozen RQ1-RQ4 cross-scale integration."
    )
    ap.add_argument("--rq1", required=True,
                    help="RQ1 rerun ZIP or extracted directory.")
    ap.add_argument("--rq2", required=True,
                    help="RQ2 rerun ZIP or extracted directory.")
    ap.add_argument("--rq3_results", required=True,
                    help="Frozen RQ3 setting directory.")
    ap.add_argument("--rq4_results", required=True,
                    help="Frozen RQ4 v3 results directory.")
    ap.add_argument("--rq3_comparison", default="WesternDiet_vs_Chow")
    ap.add_argument("--anchor_abs_diff", type=float, default=0.20)
    ap.add_argument("--rq4_flux_threshold", type=float, default=0.01,
                    help="Documentary threshold used by the frozen RQ4 pipeline.")
    ap.add_argument("--output", default="reviewer_results/final_cross_rq")
    ap.add_argument("--strict_expected_counts", action="store_true",
                    help="Fail unless frozen-study counts match expected final counts.")
    args = ap.parse_args()

    out = ensure_outdir(args.output)
    rq1src = CsvSource(args.rq1)
    rq2src = CsvSource(args.rq2)

    try:
        # -------------------------------------------------------------
        # RQ1
        # -------------------------------------------------------------
        hfd = rq1src.read_csv("reaction_stats_HFD_vs_SCD.csv")
        kd  = rq1src.read_csv("reaction_stats_KD_vs_SCD.csv")
        wd  = rq1src.read_csv("reaction_stats_WD_vs_SCD.csv")

        for label, df in [("RQ1 HFD", hfd), ("RQ1 KD", kd), ("RQ1 WD", wd)]:
            require_columns(
                df,
                ["ReactionID", "AdjustedMeanDiff", "q_value", "Cohen_d",
                 "Direction", "Significant"],
                label,
            )
            df["Significant"] = as_bool(df["Significant"])

        anchors = hfd[
            hfd["Significant"]
            & (pd.to_numeric(hfd["AdjustedMeanDiff"], errors="coerce").abs()
               >= args.anchor_abs_diff)
        ].copy()

        hfd_sig = set(hfd.loc[hfd["Significant"], "ReactionID"])
        kd_sig  = set(kd.loc[kd["Significant"], "ReactionID"])
        wd_sig  = set(wd.loc[wd["Significant"], "ReactionID"])
        rq1_shared = hfd_sig & kd_sig & wd_sig

        annotations = safe_get_annotation_table(hfd)

        # -------------------------------------------------------------
        # RQ2
        # -------------------------------------------------------------
        rq2 = rq2src.read_csv("RQ2_cross_strain_reaction_summary.csv")
        require_columns(
            rq2,
            ["ReactionID", "N_magnitude_responsive", "MagnitudeUniversal",
             "ConservationTier"],
            "RQ2",
        )
        rq2["MagnitudeUniversal"] = as_bool(rq2["MagnitudeUniversal"])
        rq2_universal = set(
            rq2.loc[rq2["MagnitudeUniversal"], "ReactionID"]
        )
        rq2i = rq2.set_index("ReactionID")

        # -------------------------------------------------------------
        # RQ3
        # -------------------------------------------------------------
        rq3_path = Path(args.rq3_results) / "statistics" / "statistical_tests.csv"
        if not rq3_path.exists():
            raise FileNotFoundError(rq3_path)
        rq3 = pd.read_csv(rq3_path)

        require_columns(
            rq3,
            ["reaction_id", "cell_type", "comparison", "abs_flux_change"],
            "RQ3",
        )

        if "significant" in rq3.columns:
            rq3["responsive"] = as_bool(rq3["significant"])
        else:
            # Fallback to the finalized deterministic threshold rule.
            require_columns(
                rq3,
                ["abs_fold_change", "relative_magnitude"],
                "RQ3 threshold reconstruction",
            )
            n = (
                (pd.to_numeric(rq3["abs_fold_change"], errors="coerce") > 1.5).astype(int)
                + (pd.to_numeric(rq3["abs_flux_change"], errors="coerce") > 0.10).astype(int)
                + (pd.to_numeric(rq3["relative_magnitude"], errors="coerce").abs() > 0.5).astype(int)
            )
            rq3["responsive"] = n >= 2

        rq3 = rq3[rq3["comparison"].astype(str) == args.rq3_comparison].copy()
        rq3_sig = rq3[rq3["responsive"]].copy()
        rq3_set = set(rq3_sig["reaction_id"])

        rq3_top = (
            rq3_sig.sort_values("abs_flux_change", ascending=False)
            .drop_duplicates("reaction_id")
            .set_index("reaction_id")
        )
        rq3_n_celltypes = rq3_sig.groupby("reaction_id")["cell_type"].nunique()

        # -------------------------------------------------------------
        # RQ4
        # -------------------------------------------------------------
        rq4root = Path(args.rq4_results)
        rq4_fc_path = (
            rq4root
            / "host_scenarios"
            / "primary"
            / "condition_DD_HFD"
            / "DD_HFD_flux_comparison.csv"
        )
        rq4_attr_path = (
            rq4root
            / "host_scenarios"
            / "primary"
            / "attribution"
            / "flux_attribution_analysis.csv"
        )
        if not rq4_fc_path.exists():
            raise FileNotFoundError(rq4_fc_path)
        if not rq4_attr_path.exists():
            raise FileNotFoundError(rq4_attr_path)

        rq4_fc = pd.read_csv(rq4_fc_path)
        rq4_attr = pd.read_csv(rq4_attr_path)

        require_columns(
            rq4_fc,
            ["reaction_id", "flux_delta", "microbiome_attributable"],
            "RQ4 primary HFD flux comparison",
        )
        require_columns(
            rq4_attr,
            ["reaction_id", "delta_diet", "delta_microbiome",
             "variance_explained_diet", "variance_explained_microbiome",
             "dominant_driver"],
            "RQ4 primary attribution",
        )

        rq4_fc["microbiome_attributable"] = as_bool(
            rq4_fc["microbiome_attributable"]
        )
        rq4_threshold_set = set(
            rq4_fc.loc[rq4_fc["microbiome_attributable"], "reaction_id"]
        )

        rq4_dom_mask = rq4_attr["dominant_driver"].astype(str).str.startswith(
            "Microbiome"
        )
        rq4_dominant_set = set(
            rq4_attr.loc[rq4_dom_mask, "reaction_id"]
        )

        # Add annotations to RQ4.
        rq4_ann = rq4_attr.copy()
        rq4_ann["ReactionName"] = rq4_ann["reaction_id"].map(
            annotations["ReactionName"]
        )
        rq4_ann["Subsystem"] = rq4_ann["reaction_id"].map(
            annotations["Subsystem"]
        ).fillna("Unassigned")

        rq4_fc_i = rq4_fc.set_index("reaction_id")
        rq4_ann["threshold_attributable_abs_delta_ge_threshold"] = (
            rq4_ann["reaction_id"].isin(rq4_threshold_set)
        )
        rq4_ann["HFD_microbiome_flux_delta"] = rq4_ann["reaction_id"].map(
            rq4_fc_i["flux_delta"]
        )

        rq4_threshold = rq4_ann[
            rq4_ann["threshold_attributable_abs_delta_ge_threshold"]
        ].copy()
        rq4_threshold = rq4_threshold.sort_values(
            "variance_explained_microbiome", ascending=False
        )
        rq4_threshold.to_csv(
            out / "RQ4_primary_threshold_attributable_reactions.csv",
            index=False,
        )

        rq4_dominant = rq4_ann[rq4_dom_mask].copy()
        rq4_dominant = rq4_dominant.sort_values(
            "variance_explained_microbiome", ascending=False
        )
        rq4_dominant.to_csv(
            out / "RQ4_primary_microbiome_dominant_reactions.csv",
            index=False,
        )

        rq4_pathway = (
            rq4_threshold.groupby("Subsystem", dropna=False)
            .agg(
                n_threshold_attributable=("reaction_id", "nunique"),
                n_microbiome_dominant=(
                    "dominant_driver",
                    lambda s: s.astype(str).str.startswith("Microbiome").sum(),
                ),
                mean_VE_microbiome=("variance_explained_microbiome", "mean"),
                median_abs_microbiome_delta=(
                    "delta_microbiome",
                    lambda s: float(np.nanmedian(np.abs(s))),
                ),
            )
            .reset_index()
            .sort_values(
                ["n_threshold_attributable", "n_microbiome_dominant"],
                ascending=[False, False],
            )
        )
        rq4_pathway.to_csv(
            out / "RQ4_primary_pathway_attribution.csv", index=False
        )

        # -------------------------------------------------------------
        # Cross-RQ reaction tracing of RQ1 HFD anchors
        # -------------------------------------------------------------
        trace_cols = [
            c for c in
            ["ReactionID", "ReactionName", "Subsystem", "AdjustedMeanDiff",
             "q_value", "Cohen_d", "Direction"]
            if c in anchors.columns
        ]
        trace = anchors[trace_cols].copy()

        trace["RQ2_N_magnitude_responsive"] = trace["ReactionID"].map(
            rq2i["N_magnitude_responsive"]
        )
        trace["RQ2_conservation_tier"] = trace["ReactionID"].map(
            rq2i["ConservationTier"]
        )
        trace["RQ2_universal_magnitude"] = trace["ReactionID"].isin(
            rq2_universal
        )

        trace["RQ3_WDChow_responsive"] = trace["ReactionID"].isin(rq3_set)
        trace["RQ3_N_responsive_celltypes"] = (
            trace["ReactionID"].map(rq3_n_celltypes).fillna(0).astype(int)
        )
        if len(rq3_top):
            trace["RQ3_primary_celltype"] = trace["ReactionID"].map(
                rq3_top["cell_type"]
            ).fillna("")
        else:
            trace["RQ3_primary_celltype"] = ""

        trace["RQ4_HFD_threshold_attributable"] = trace["ReactionID"].isin(
            rq4_threshold_set
        )
        trace["RQ4_HFD_microbiome_dominant"] = trace["ReactionID"].isin(
            rq4_dominant_set
        )

        rq4i = rq4_ann.set_index("reaction_id")
        for c in [
            "delta_diet", "delta_microbiome", "variance_explained_diet",
            "variance_explained_microbiome", "dominant_driver",
        ]:
            trace[f"RQ4_{c}"] = trace["ReactionID"].map(rq4i[c])

        trace.to_csv(
            out / "cross_RQ_reaction_tracing_final.csv", index=False
        )

        # -------------------------------------------------------------
        # Union membership matrix across all reaction-level sets
        # -------------------------------------------------------------
        union_reactions = sorted(
            hfd_sig
            | kd_sig
            | wd_sig
            | rq2_universal
            | rq3_set
            | rq4_threshold_set
            | rq4_dominant_set
        )
        membership = pd.DataFrame({"ReactionID": union_reactions})
        membership["ReactionName"] = membership["ReactionID"].map(
            annotations["ReactionName"]
        ).fillna("")
        membership["Subsystem"] = membership["ReactionID"].map(
            annotations["Subsystem"]
        ).fillna("Unassigned")
        membership["RQ1_HFD_significant"] = membership["ReactionID"].isin(hfd_sig)
        membership["RQ1_KD_significant"] = membership["ReactionID"].isin(kd_sig)
        membership["RQ1_WD_significant"] = membership["ReactionID"].isin(wd_sig)
        membership["RQ1_all3_shared"] = membership["ReactionID"].isin(rq1_shared)
        membership["RQ1_HFD_anchor"] = membership["ReactionID"].isin(
            set(anchors["ReactionID"])
        )
        membership["RQ2_universal_magnitude"] = membership["ReactionID"].isin(
            rq2_universal
        )
        membership["RQ3_WDChow_responsive"] = membership["ReactionID"].isin(
            rq3_set
        )
        membership["RQ4_HFD_threshold_attributable"] = membership[
            "ReactionID"
        ].isin(rq4_threshold_set)
        membership["RQ4_microbiome_dominant"] = membership["ReactionID"].isin(
            rq4_dominant_set
        )
        membership.to_csv(
            out / "cross_RQ_reaction_membership_matrix.csv", index=False
        )

        # -------------------------------------------------------------
        # Pathway-level integration
        # -------------------------------------------------------------
        p1 = (
            hfd[hfd["Significant"]]
            .groupby("Subsystem", dropna=False)
            .agg(
                RQ1_HFD_n_significant=("ReactionID", "nunique"),
                RQ1_HFD_median_abs_adjusted_diff=(
                    "AdjustedMeanDiff",
                    lambda s: float(np.nanmedian(np.abs(s))),
                ),
            )
            .reset_index()
        )

        p2 = (
            rq2[rq2["N_magnitude_responsive"] > 0]
            .groupby("Subsystem", dropna=False)
            .agg(
                RQ2_n_magnitude_responsive_union=("ReactionID", "nunique"),
                RQ2_n_universal=("MagnitudeUniversal", "sum"),
                RQ2_mean_N_strains=("N_magnitude_responsive", "mean"),
            )
            .reset_index()
        )

        p3 = (
            rq3_sig.groupby("subsystem", dropna=False)
            .agg(
                RQ3_WDChow_n_responsive=("reaction_id", "nunique"),
                RQ3_WDChow_n_celltypes=("cell_type", "nunique"),
            )
            .reset_index()
            .rename(columns={"subsystem": "Subsystem"})
        )

        p4 = (
            rq4_threshold.groupby("Subsystem", dropna=False)
            .agg(
                RQ4_HFD_n_threshold_attributable=("reaction_id", "nunique"),
                RQ4_HFD_n_microbiome_dominant=(
                    "dominant_driver",
                    lambda s: s.astype(str).str.startswith("Microbiome").sum(),
                ),
                RQ4_HFD_mean_VE_microbiome=(
                    "variance_explained_microbiome", "mean"
                ),
            )
            .reset_index()
        )

        pathways = (
            p1.merge(p2, on="Subsystem", how="outer")
              .merge(p3, on="Subsystem", how="outer")
              .merge(p4, on="Subsystem", how="outer")
        )

        for c in pathways.columns:
            if c != "Subsystem":
                pathways[c] = pd.to_numeric(
                    pathways[c], errors="coerce"
                ).fillna(0)

        pathways["N_RQ_layers_with_signal"] = (
            pathways["RQ1_HFD_n_significant"].gt(0).astype(int)
            + pathways["RQ2_n_magnitude_responsive_union"].gt(0).astype(int)
            + pathways["RQ3_WDChow_n_responsive"].gt(0).astype(int)
            + pathways["RQ4_HFD_n_threshold_attributable"].gt(0).astype(int)
        )

        pathways = pathways.sort_values(
            [
                "N_RQ_layers_with_signal",
                "RQ4_HFD_n_threshold_attributable",
                "RQ1_HFD_n_significant",
            ],
            ascending=[False, False, False],
        )
        pathways.to_csv(
            out / "cross_RQ_pathway_matrix_final.csv", index=False
        )

        # -------------------------------------------------------------
        # Final summary
        # -------------------------------------------------------------
        anchor_set = set(anchors["ReactionID"])
        all4_pathways = pathways.loc[
            pathways["N_RQ_layers_with_signal"] == 4, "Subsystem"
        ].dropna().astype(str).tolist()

        summary = {
            "RQ1_HFD_anchor_definition":
                f"q<0.05 and |dataset-adjusted mean difference|>={args.anchor_abs_diff}",
            "RQ1_HFD_anchor_n": int(len(anchor_set)),
            "RQ1_shared_all_three_SCD_referenced_n": int(len(rq1_shared)),
            "RQ2_universal_magnitude_core_n": int(len(rq2_universal)),
            "RQ3_P99_biomass_WDChow_responsive_unique_reactions":
                int(len(rq3_set)),
            "RQ4_primary_HFD_threshold_documentary_abs_delta":
                float(args.rq4_flux_threshold),
            "RQ4_primary_HFD_threshold_attributable_n":
                int(len(rq4_threshold_set)),
            "RQ4_primary_microbiome_dominant_n":
                int(len(rq4_dominant_set)),
            "RQ1_anchor_overlap_RQ2_universal_n":
                int(len(anchor_set & rq2_universal)),
            "RQ1_anchor_overlap_RQ3_WDChow_n":
                int(len(anchor_set & rq3_set)),
            "RQ1_anchor_overlap_RQ4_threshold_n":
                int(len(anchor_set & rq4_threshold_set)),
            "RQ1_anchor_overlap_RQ4_microbiome_dominant_n":
                int(len(anchor_set & rq4_dominant_set)),
            "RQ1_shared27_overlap_RQ2_universal_n":
                int(len(rq1_shared & rq2_universal)),
            "RQ1_shared27_overlap_RQ3_WDChow_n":
                int(len(rq1_shared & rq3_set)),
            "RQ1_shared27_overlap_RQ4_threshold_n":
                int(len(rq1_shared & rq4_threshold_set)),
            "pathways_with_signal_in_all_4_RQs":
                int(len(all4_pathways)),
            "all4_pathways": all4_pathways,
            "interpretation":
                "No single exact reaction-level conserved backbone is supported "
                "across RQ1-RQ4. Cross-scale convergence is stronger at pathway "
                "level. RQ3 is an orthogonal WesternDiet-vs-Chow cellular context, "
                "not condition-matched HFD validation.",
        }

        (out / "cross_RQ_summary_final.json").write_text(
            json.dumps(summary, indent=2), encoding="utf-8"
        )

        manifest = {
            "inputs": {
                "rq1": str(Path(args.rq1).resolve()),
                "rq2": str(Path(args.rq2).resolve()),
                "rq3_results": str(Path(args.rq3_results).resolve()),
                "rq4_results": str(Path(args.rq4_results).resolve()),
            },
            "parameters": {
                "rq3_comparison": args.rq3_comparison,
                "anchor_abs_diff": args.anchor_abs_diff,
                "rq4_flux_threshold_documentary": args.rq4_flux_threshold,
            },
            "counts": summary,
        }
        (out / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2), encoding="utf-8"
        )

        # -------------------------------------------------------------
        # Frozen-study count QC
        # -------------------------------------------------------------
        expected = {
            "RQ1_HFD_anchor_n": 20,
            "RQ1_shared_all_three_SCD_referenced_n": 27,
            "RQ2_universal_magnitude_core_n": 16,
            "RQ3_P99_biomass_WDChow_responsive_unique_reactions": 562,
            "RQ4_primary_HFD_threshold_attributable_n": 47,
        }

        mismatches = []
        for k, exp in expected.items():
            obs = summary.get(k)
            if obs != exp:
                mismatches.append(f"{k}: observed={obs}, expected={exp}")

        print("\n=== FINAL CROSS-RQ SUMMARY ===")
        for k in expected:
            print(f"{k}: {summary[k]}")
        print(
            "RQ1 anchor overlaps: "
            f"RQ2={summary['RQ1_anchor_overlap_RQ2_universal_n']}, "
            f"RQ3={summary['RQ1_anchor_overlap_RQ3_WDChow_n']}, "
            f"RQ4_threshold={summary['RQ1_anchor_overlap_RQ4_threshold_n']}, "
            f"RQ4_dominant={summary['RQ1_anchor_overlap_RQ4_microbiome_dominant_n']}"
        )
        print(
            "All-4-RQ pathways:",
            ", ".join(all4_pathways) if all4_pathways else "none",
        )

        if mismatches:
            msg = (
                "\n[WARNING] Frozen-study counts differ from the expected final "
                "revision counts:\n  " + "\n  ".join(mismatches)
            )
            if args.strict_expected_counts:
                raise RuntimeError(msg)
            print(msg)
        else:
            print("[QC PASS] Frozen-study headline counts match expected values.")

        print(f"\n[DONE] Outputs written to: {out.resolve()}")

    finally:
        rq1src.close()
        rq2src.close()


if __name__ == "__main__":
    main()

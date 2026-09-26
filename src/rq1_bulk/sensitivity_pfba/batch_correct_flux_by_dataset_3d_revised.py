#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch-adjust flux profiles for VISUALIZATION/QC only.

For each reaction:
    flux ~ diet_group + dataset
and subtract only the fitted dataset component.

Important:
- The corrected vectors are statistical feature vectors; they are NOT guaranteed
  to satisfy S v = 0 or the original reaction bounds.
- Use them for PCA, clustering, heatmaps, and batch diagnostics.
- Do NOT use them as the primary reaction-level inferential input. The revised
  comprehensive_flux_analysis.py tests diet effects on the original feasible
  pFBA flux vectors while adjusting for dataset in the statistical model.
"""

from __future__ import annotations
import argparse, json, os, re
from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


GROUP_ALIASES = {
    "SCD": {"SCD","SC","CD","CHOW","ND","CONTROL","STDCHOW","STANDARDCHOW","LFD","NCD","CTRL","CONTROL"},
    "HFD": {"HFD","HF","HIGHFAT","HIGH_FAT","HIGH-FAT"},
    "KD": {"KD","KETO","KETOGENIC","KETO_DIET","KETO-DIET"},
    "WD": {"WD","WESTERN","WESTERNDIET"},
}


def canonical_group(text: str) -> str:
    up = str(text).upper()
    tokens = re.split(r"[^A-Z0-9]+", up)
    for tok in tokens:
        for canon, aliases in GROUP_ALIASES.items():
            if tok == canon or tok in aliases:
                return canon
    # Fallback to prefix, but make it visible in QC.
    return str(text).split("_", 1)[0].upper()


def parse_dataset(text: str, pattern=r"(?:GSE|GSM)\d+") -> str:
    m = re.search(pattern, str(text), flags=re.I)
    return m.group(0).upper() if m else "UNKNOWN"


def detect_sample_columns(columns, suffix="_Flux"):
    bad = re.compile(r"(MeanFlux|StdFlux|SEMFlux|_N)$", re.I)
    return [str(c) for c in columns if str(c).endswith(suffix) and not bad.search(str(c))]


def design_matrix(groups, datasets):
    g = pd.get_dummies(pd.Categorical(groups), drop_first=True, dtype=float)
    b = pd.get_dummies(pd.Categorical(datasets), drop_first=True, dtype=float)
    n = len(groups)
    intercept = np.ones((n, 1), dtype=float)
    G = g.to_numpy(dtype=float)
    B = b.to_numpy(dtype=float)
    D = np.concatenate([intercept, G, B], axis=1)
    Db = np.concatenate([np.zeros_like(intercept), np.zeros_like(G), B], axis=1)
    return D, Db, list(g.columns), list(b.columns)


def remove_dataset_effect(X, groups, datasets):
    """Vectorized OLS across all reactions: X is samples x reactions."""
    D, Db, group_terms, batch_terms = design_matrix(groups, datasets)
    rank = int(np.linalg.matrix_rank(D))
    if rank < D.shape[1]:
        raise ValueError(
            f"Design matrix is rank deficient (rank={rank}, columns={D.shape[1]}). "
            "Diet and dataset are not sufficiently connected for this adjustment."
        )
    beta = np.linalg.pinv(D) @ X
    batch_component = Db @ beta
    return X - batch_component, {
        "design_rank": rank,
        "design_columns": int(D.shape[1]),
        "group_terms": [str(x) for x in group_terms],
        "dataset_terms": [str(x) for x in batch_terms],
    }


def pca_qc(X, meta, outdir, tag):
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)
    ncomp = min(10, Xs.shape[0]-1, Xs.shape[1])
    if ncomp < 2:
        return {}
    pca = PCA(n_components=ncomp, random_state=42)
    Z = pca.fit_transform(Xs)
    for key in ["Dataset", "Group"]:
        fig, ax = plt.subplots(figsize=(7, 6))
        for lvl in meta[key].unique():
            mask = (meta[key] == lvl).to_numpy()
            ax.scatter(Z[mask,0], Z[mask,1], s=45, alpha=.8, label=str(lvl))
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        ax.set_title(f"PCA {tag}: by {key.lower()}")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(Path(outdir)/f"PCA_{tag}_by_{key.lower()}.png", dpi=300)
        plt.close(fig)
    return {
        "PC1_variance": float(pca.explained_variance_ratio_[0]),
        "PC2_variance": float(pca.explained_variance_ratio_[1]),
    }


def recompute_derived_columns(df, sample_cols, meta):
    """Keep the wide table internally consistent after sample-column adjustment."""
    out = df.copy()
    group_to_cols = {
        g: meta.loc[meta["Group"] == g, "SampleID"].tolist()
        for g in sorted(meta["Group"].unique())
    }
    for g, cols in group_to_cols.items():
        vals = out[cols].apply(pd.to_numeric, errors="coerce")
        mean_col, std_col, n_col = f"{g}_MeanFlux", f"{g}_StdFlux", f"{g}_N"
        if mean_col in out.columns:
            out[mean_col] = vals.mean(axis=1)
        if std_col in out.columns:
            out[std_col] = vals.std(axis=1, ddof=0)
        if n_col in out.columns:
            out[n_col] = len(cols)

    diff_re = re.compile(r"^Diff\(([^)]+)-([^)]+)\)$")
    ratio_re = re.compile(r"^Ratio\(([^)]+)/([^)]+)\)$")
    for col in list(out.columns):
        md = diff_re.match(str(col))
        if md:
            a,b = md.groups()
            if a in group_to_cols and b in group_to_cols:
                ma = out[group_to_cols[a]].mean(axis=1)
                mb = out[group_to_cols[b]].mean(axis=1)
                out[col] = ma - mb
        mr = ratio_re.match(str(col))
        if mr:
            a,b = mr.groups()
            if a in group_to_cols and b in group_to_cols:
                aa = out[group_to_cols[a]].abs().mean(axis=1)
                bb = out[group_to_cols[b]].abs().mean(axis=1) + 1e-12
                out[col] = aa / bb
    return out


def main():
    ap = argparse.ArgumentParser(description="Dataset-effect adjustment for flux visualization/QC.")
    ap.add_argument("--input","-i",required=True)
    ap.add_argument("--output_csv","-o",required=True)
    ap.add_argument("--outdir",default="batch_correction_qc")
    ap.add_argument("--sample_suffix",default="_Flux")
    ap.add_argument("--dataset_regex",default=r"(?:GSE|GSM)\d+")
    ap.add_argument("--no_plots",action="store_true")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input)
    sample_cols = detect_sample_columns(df.columns, args.sample_suffix)
    if not sample_cols:
        raise SystemExit("No sample columns detected.")
    meta = pd.DataFrame({
        "SampleID": sample_cols,
        "Group": [canonical_group(c) for c in sample_cols],
        "Dataset": [parse_dataset(c, args.dataset_regex) for c in sample_cols],
    })
    if (meta["Dataset"] == "UNKNOWN").any():
        bad = meta.loc[meta["Dataset"]=="UNKNOWN","SampleID"].tolist()[:5]
        raise ValueError(f"Dataset ID missing from one or more sample names; first examples: {bad}")
    meta.to_csv(outdir/"sample_metadata.csv", index=False)

    X = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float).T
    pre = {} if args.no_plots else pca_qc(X, meta, outdir, "pre")
    Xc, design_qc = remove_dataset_effect(X, meta["Group"].values, meta["Dataset"].values)
    post = {} if args.no_plots else pca_qc(Xc, meta, outdir, "post")

    out = df.copy()
    out.loc[:, sample_cols] = Xc.T
    out = recompute_derived_columns(out, sample_cols, meta)
    Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output_csv, index=False)

    qc = {
        "purpose": "visualization_qc_only",
        "input": str(args.input),
        "output": str(args.output_csv),
        "n_samples": len(sample_cols),
        "n_reactions": int(X.shape[1]),
        "groups": meta["Group"].value_counts().to_dict(),
        "datasets": meta["Dataset"].value_counts().to_dict(),
        **design_qc,
        "pca_pre": pre,
        "pca_post": post,
        "warning": "Batch-adjusted vectors are not guaranteed to satisfy S*v=0 or model bounds."
    }
    (outdir/"batch_correction_qc.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")
    print(f"[OK] Corrected visualization matrix: {args.output_csv}")
    print(f"[OK] QC: {outdir/'batch_correction_qc.json'}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

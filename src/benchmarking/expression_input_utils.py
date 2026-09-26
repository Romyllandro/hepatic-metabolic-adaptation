#!/usr/bin/env python3
from __future__ import annotations
import numpy as np
import pandas as pd

def detect_gene_symbol_column(df: pd.DataFrame) -> str:
    preferred = ["Gene_Symbol", "gene_symbol", "GeneSymbol", "symbol", "Symbol"]
    for c in preferred:
        if c in df.columns:
            return c
    raise ValueError(
        "Could not find a gene-symbol column. Expected one of "
        f"{preferred}. Available first columns: {list(df.columns[:12])}"
    )

def sample_columns_for_dataset(df, dataset: str, conditions):
    annotation = {"Gene_Symbol","Gene_ID","gene_symbol","gene_id","GeneSymbol","symbol","Symbol"}
    cols = [c for c in df.columns if c not in annotation]
    groups = {}
    selected = []
    for cond in conditions:
        hits = [c for c in cols if c.startswith(f"{cond}_{dataset}_")]
        if not hits:
            hits = [c for c in cols if c.startswith(cond + "_") and dataset in c]
        for c in hits:
            groups[c] = cond
        selected.extend(hits)
    return selected, groups

def available_case_mean(df: pd.DataFrame, cols):
    X = df[list(cols)].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    valid = np.sum(np.isfinite(X), axis=1)
    sums = np.nansum(X, axis=1)
    out = np.full(X.shape[0], np.nan, dtype=float)
    ok = valid > 0
    out[ok] = sums[ok] / valid[ok]
    return out

def preflight_symbol_to_model(gene_symbols, symbol_to_entrez, model, min_matches=100):
    model_genes = {str(g.id) for g in model.genes}
    mapped = []
    for g in gene_symbols:
        key = str(g).lower().strip()
        gid = symbol_to_entrez.get(key)
        if gid is not None and str(gid) in model_genes:
            mapped.append(str(gid))
    n = len(set(mapped))
    print(f"[PRECHECK] Gene symbols mapping into iMM1415 genes: {n}")
    if n < min_matches:
        raise ValueError(
            f"Only {n} gene symbols map into model genes. The production RQ1/RQ2 "
            "runs matched ~1,100 genes. Stop: this benchmark would not apply "
            "expression constraints. Check Gene_Symbol and mouse_entrez_to_symbol.csv."
        )
    return n

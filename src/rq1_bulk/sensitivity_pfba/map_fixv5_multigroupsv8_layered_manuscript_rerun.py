#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import json
import argparse
import traceback
from typing import Dict, Set, Tuple, List, Iterable

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import cobra
from cobra.io import load_json_model
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis import pfba as _cobra_pfba
from scipy.stats import ttest_ind, ranksums, mannwhitneyu
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA

# =====================================
# Group handling
# =====================================
GROUP_ALIASES = {
    "SCD": {"SCD", "SC", "CD", "CHOW", "ND", "CONTROL", "STDCHOW", "STANDARDCHOW"},
    "WD":  {"WD", "WESTERN", "WESTERNDIET"},
    "HFD": {"HFD", "HF", "HIGHFAT", "HIGH_FAT", "HIGH-FAT"},
    "KD":  {"KD", "KETO", "KETOGENIC", "KETO_DIET", "KETO-DIET"},
    "LFD": {"LFD", "LOWFAT", "LOW_FAT", "LOW-FAT"},
}

def canonical_code(tag: str) -> str:
    t = re.sub(r"[^A-Za-z0-9]+", "", str(tag).upper())
    for canon, aliases in GROUP_ALIASES.items():
        if t == canon or t in aliases:
            return canon
    return t

def parse_groups_from_filename_multi(path: str) -> List[str]:
    base = os.path.basename(path)
    m = re.search(r"_([A-Za-z0-9_]+)_gene_expression\.csv$", base, re.IGNORECASE)
    if not m:
        raise ValueError(f"Filename must contain *_<GROUPS>_gene_expression.csv: {base}")
    chunk = m.group(1)
    raw_tags = [t for t in chunk.split("_") if t]
    tags = [canonical_code(t) for t in raw_tags]
    if len(tags) < 2:
        raise ValueError(f"Need at least 2 group tags in filename; got: {tags}")
    return tags

def choose_default_baseline(tags: List[str]) -> str:
    pref = ["SCD", "WD"]
    for p in pref:
        if p in tags:
            return p
    return tags[0]

# =====================================
# Gene ID mapping
# =====================================
def load_symbol_to_entrez_mapping(mapping_file: str) -> Dict[str, str]:
    """
    Load mapping from gene symbol to entrez ID.
    Returns dict: {lowercase_symbol: entrez_id_as_string}
    """
    df = pd.read_csv(mapping_file)
    mapping = {}
    for _, row in df.iterrows():
        symbol = str(row.get('symbol', '')).strip().lower()
        
        # Handle Entrez ID - convert from float to int to string to remove .0
        entrez_raw = row.get('entrez', '')
        if pd.isna(entrez_raw) or entrez_raw == '':
            continue
        try:
            # Convert to int first (handles floats like 239559.0)
            entrez = str(int(float(entrez_raw)))
        except (ValueError, TypeError):
            continue
            
        if symbol and entrez:
            mapping[symbol] = entrez
            
        # Also add aliases if present
        aliases_str = str(row.get('alias', '')).strip()
        if aliases_str and aliases_str != 'nan':
            for alias in aliases_str.split(','):
                alias = alias.strip().lower()
                if alias and alias not in mapping:
                    mapping[alias] = entrez
    return mapping

# =====================================
# Objective selection
# =====================================
def set_objective_reaction(model, objective_id=None, objective_regex=None, sense="max", logger=print):
    """
    Select objective by id or regex; fallback to biomass/ATPM/first.
    sense: 'max' or 'min'
    """
    chosen = None
    if objective_id:
        try:
            chosen = model.reactions.get_by_id(objective_id)
        except KeyError:
            logger(f"[WARN] Objective id '{objective_id}' not found; will try regex/fallbacks.")
    if (chosen is None) and objective_regex:
        rx = re.compile(objective_regex, re.IGNORECASE)
        for rxn in model.reactions:
            if rx.search(rxn.id) or rx.search(rxn.name or "") or (rxn.subsystem and rx.search(rxn.subsystem)):
                chosen = rxn; break
        if chosen is None:
            logger(f"[WARN] No reaction matched objective_regex '{objective_regex}'.")
    if chosen is None:
        for rxn in model.reactions:
            nm = (rxn.name or "").lower()
            if "biomass" in nm or (rxn.subsystem and "biomass" in rxn.subsystem.lower()):
                chosen = rxn; break
    if chosen is None:
        try:
            chosen = model.reactions.get_by_id("ATPM")
        except KeyError:
            chosen = model.reactions[0]
    # Apply sense
    try:
        model.objective = chosen
        if sense.lower().startswith("min"):
            model.objective_direction = "min"
        else:
            model.objective_direction = "max"
    except Exception:
        pass
    logger(f"[INFO] Objective set to: {chosen.id} ({chosen.name}) | sense={sense}")
    return model, chosen.id

# =====================================
# Reaction classification (EX / transporter / internal)
# =====================================
def classify_reactions(
    model,
    transporter_strategy="e_to_non_e",
    transporter_regex=None,
    transporter_subsystem_regex=None,
    compartments_for_transport=("e",),
) -> Tuple[Set[str], Set[str], Set[str]]:
    """
    transporter_strategy:
      - 'e_to_non_e' (default): any reaction with an 'e' metabolite and any non-'e' metabolite
      - 'regex': reactions whose id matches transporter_regex OR subsystem matches transporter_subsystem_regex
      - 'either': union of e_to_non_e and regex strategies
    """
    ex_rxns, transporter_rxns, internal_rxns = set(), set(), set()

    rx_id_re = re.compile(transporter_regex, re.IGNORECASE) if transporter_regex else None
    rx_subsys_re = re.compile(transporter_subsystem_regex, re.IGNORECASE) if transporter_subsystem_regex else None
    comp_set = set(compartments_for_transport or [])

    def is_transporter_by_e_non_e(rxn) -> bool:
        comps = set([m.compartment for m in rxn.metabolites])
        return (len(comp_set & comps) > 0) and any((c not in comp_set) for c in comps)

    def is_transporter_by_regex(rxn) -> bool:
        if rx_id_re and rx_id_re.search(rxn.id):
            return True
        if rx_subsys_re and rxn.subsystem and rx_subsys_re.search(rxn.subsystem):
            return True
        return False

    for rxn in model.reactions:
        rid = rxn.id
        if rid.startswith("EX_"):
            ex_rxns.add(rid); continue
        by_e_non_e = is_transporter_by_e_non_e(rxn)
        by_regex = is_transporter_by_regex(rxn)
        pick = False
        if transporter_strategy == "e_to_non_e":
            pick = by_e_non_e
        elif transporter_strategy == "regex":
            pick = by_regex
        elif transporter_strategy == "either":
            pick = by_e_non_e or by_regex
        else:
            pick = by_e_non_e
        if pick:
            transporter_rxns.add(rid)
        else:
            internal_rxns.add(rid)
    return ex_rxns, transporter_rxns, internal_rxns

# =====================================
# Expression mapping (E-flux style), scoped - FIXED VERSION
# =====================================
def _standardize_gene_id(gene_id: str) -> str:
    """Standardize gene ID by removing special characters and converting to lowercase"""
    return re.sub(r"[()\-'\",;]", "", str(gene_id).lower().strip())

def _build_gene_index(gene_names):
    """Build index from gene names/IDs to their position in the expression vector"""
    d = {}
    for idx, g in enumerate(gene_names):
        g_str = str(g).strip()
        
        # Store the original
        d[g_str] = idx
        
        # Store standardized version (lowercase, no special chars)
        standardized = _standardize_gene_id(g)
        if standardized != g_str:
            d[standardized] = idx
        
        # Store lowercase version
        g_lower = g_str.lower()
        if g_lower != g_str and g_lower != standardized:
            d[g_lower] = idx
            
        # If it looks like a number, also store without leading zeros and .0
        try:
            g_num = str(int(float(g_str)))
            if g_num != g_str:
                d[g_num] = idx
        except (ValueError, TypeError):
            pass
    
    return d

def _map_expression_to_reactions(model, gene_index, expr_vector, symbol_to_entrez=None):
    """
    Map gene expression to reactions using GPR rules.
    
    Args:
        model: COBRA model
        gene_index: dict mapping gene identifiers to expression vector indices
        expr_vector: numpy array of expression values
        symbol_to_entrez: dict mapping gene symbols to entrez IDs (optional)
    """
    rxn_expr = np.zeros(len(model.reactions))
    matched_genes = set()
    
    for i, rxn in enumerate(model.reactions):
        rule = rxn.gene_reaction_rule
        if not rule:
            continue
        
        rule_clean = re.sub(r"[()]", "", rule)
        or_parts = rule_clean.split(" or ")
        values = []
        
        for part in or_parts:
            and_parts = part.split(" and ")
            vals = []
            
            for gene in and_parts:
                gene = gene.strip()
                if not gene:
                    continue
                
                # Try multiple matching strategies
                expr_value = None
                matched_key = None
                
                # Strategy 1: Direct match with original gene ID
                if gene in gene_index:
                    expr_value = expr_vector[gene_index[gene]]
                    matched_key = gene
                
                # Strategy 2: Standardized match
                if expr_value is None:
                    key = _standardize_gene_id(gene)
                    if key in gene_index:
                        expr_value = expr_vector[gene_index[key]]
                        matched_key = key
                
                # Strategy 3: Lowercase match
                if expr_value is None:
                    key_lower = gene.lower()
                    if key_lower in gene_index:
                        expr_value = expr_vector[gene_index[key_lower]]
                        matched_key = key_lower
                
                # Strategy 4: Try as integer (remove .0 if present)
                if expr_value is None:
                    try:
                        gene_as_int = str(int(float(gene)))
                        if gene_as_int in gene_index:
                            expr_value = expr_vector[gene_index[gene_as_int]]
                            matched_key = gene_as_int
                    except (ValueError, TypeError):
                        pass
                
                if expr_value is not None and expr_value > 0:
                    vals.append(expr_value)
                    if matched_key:
                        matched_genes.add(matched_key)
            
            pos_vals = [v for v in vals if v > 0]
            if pos_vals:
                values.append(min(pos_vals))  # AND = min
        
        if values:
            rxn_expr[i] = max(values)  # OR = max
    
    return rxn_expr, matched_genes

def apply_expression_constraints_scoped(
    model, gene_names, expr_vector, scope_rxn_ids: Set[str],
    *, eflux_quantile=0.95, eflux_floor=0.1, eflux_cap=1000.0, logger=print, label="",
    symbol_to_entrez=None
):
    """
    Apply E-Flux style expression constraints to reactions in a specified scope.
    
    Args:
        symbol_to_entrez: Optional dict mapping gene symbols to Entrez IDs
    """
    # Normalize expression
    pos = expr_vector[expr_vector > 0]
    denom = np.quantile(pos, eflux_quantile) if len(pos) else 1.0
    x = expr_vector / max(denom, 1e-9)
    x = np.clip(x, eflux_floor, eflux_cap)

    # Build gene index
    gene_index = _build_gene_index(gene_names)
    
    # Map expression to reactions
    rxn_expr, matched_genes = _map_expression_to_reactions(
        model, gene_index, x, symbol_to_entrez=symbol_to_entrez
    )

    # Apply constraints
    changed = 0
    changed_rxn_ids = set()
    for rxn, val in zip(model.reactions, rxn_expr):
        if rxn.id not in scope_rxn_ids:
            continue
        if val > 0:
            old_lb, old_ub = rxn.lower_bound, rxn.upper_bound
            if rxn.lower_bound < 0:
                rxn.lower_bound = max(rxn.lower_bound, -float(val))
            rxn.upper_bound = min(rxn.upper_bound, float(val))
            if (rxn.lower_bound != old_lb) or (rxn.upper_bound != old_ub):
                changed += 1
                changed_rxn_ids.add(rxn.id)
    
    logger(f"[INFO] Layer {label}: expression constraints applied to {changed} reactions in scope={len(scope_rxn_ids)}. Genes matched: {len(matched_genes)}")
    return model, changed, len(matched_genes), len(gene_index), matched_genes, changed_rxn_ids

# =====================================
# Layer 1: Diet bounds (with optional unit conversion)
# =====================================
def load_mw_table(path: str) -> Dict[str, float]:
    tbl = pd.read_csv(path)
    out = {}
    for _, row in tbl.iterrows():
        rid = str(row["exchange_id"]).strip()
        mw = float(row["mw_g_per_mol"])
        out[rid] = mw
    return out

def convert_bound_value(lb_model_units: float, *, units: str, ex_id: str, mw_map: Dict[str, float], gDW: float, hours_per_day: float):
    if units == "model":
        return float(lb_model_units)
    if units == "mmol_per_day":
        return float(lb_model_units) / max(gDW, 1e-9) / max(hours_per_day, 1e-9)
    if units == "g_per_day":
        mw = mw_map.get(ex_id)
        if mw is None:
            raise ValueError(f"No MW for {ex_id} in --mw_table; required for g/day conversion.")
        mmol_per_day = (lb_model_units * 1000.0) / mw
        return mmol_per_day / max(gDW, 1e-9) / max(hours_per_day, 1e-9)
    raise ValueError(f"Unsupported units: {units}")

def apply_diet_bounds_layer1(model, *, code: str, diet_bounds: Dict, diet_units: str,
                             mw_map: Dict[str, float] = None, gDW: float = 1.0, hours_per_day: float = 24.0,
                             logger=print):
    code = canonical_code(code)
    if not diet_bounds or code not in diet_bounds:
        logger(f"[INFO] No diet bounds for {code}. Skipping Layer 1.")
        return model, 0
    changed = 0
    for rxn_id, bounds in diet_bounds[code].items():
        if isinstance(bounds, (list, tuple)) and len(bounds) == 2:
            lb, ub = bounds[0], bounds[1]
        elif isinstance(bounds, dict) and "lb" in bounds and "ub" in bounds:
            lb, ub = bounds["lb"], bounds["ub"]
        else:
            raise ValueError(f"Invalid bounds format for {rxn_id}: {bounds}")
        try:
            rxn = model.reactions.get_by_id(rxn_id)
        except KeyError:
            logger(f"[WARN] Diet lists {rxn_id} but it is not in model. Skipping.")
            continue
        if not rxn.id.startswith("EX_"):
            continue
        new_lb = convert_bound_value(lb, units=diet_units, ex_id=rxn_id, mw_map=(mw_map or {}), gDW=gDW, hours_per_day=hours_per_day)
        new_ub = convert_bound_value(ub, units=diet_units, ex_id=rxn_id, mw_map=(mw_map or {}), gDW=gDW, hours_per_day=hours_per_day) if diet_units != "model" else float(ub)
        old_lb, old_ub = rxn.lower_bound, rxn.upper_bound
        rxn.lower_bound, rxn.upper_bound = float(new_lb), float(new_ub)
        if (rxn.lower_bound != old_lb) or (rxn.upper_bound != old_ub):
            changed += 1
    logger(f"[INFO] Layer 1 (Diet) applied to {changed} EX_ reactions for code={code}.")
    return model, changed


def validate_model(model, logger=print, tolerance=1e-12):
    """
    Validate reaction bounds after Layer 1-3 constraint application.

    Scientifically important: do NOT silently repair materially inconsistent
    bounds, because that changes the model without a declared rule. Tiny
    floating-point inconsistencies within `tolerance` are collapsed to the
    midpoint; larger inconsistencies raise an error with reaction IDs.
    """
    bad = []
    tiny = []
    for rxn in model.reactions:
        if rxn.lower_bound > rxn.upper_bound:
            gap = float(rxn.lower_bound - rxn.upper_bound)
            if gap <= tolerance:
                mid = 0.5 * (rxn.lower_bound + rxn.upper_bound)
                rxn.lower_bound = mid
                rxn.upper_bound = mid
                tiny.append(rxn.id)
            else:
                bad.append((rxn.id, float(rxn.lower_bound), float(rxn.upper_bound)))
    if tiny:
        logger(f"[INFO] Collapsed {len(tiny)} tiny numerical bound inconsistencies (<= {tolerance:g}).")
    if bad:
        preview = "; ".join(f"{rid}[lb={lb:g},ub={ub:g}]" for rid, lb, ub in bad[:10])
        raise ValueError(
            f"{len(bad)} reactions have lower_bound > upper_bound after constraint application. "
            f"Refusing to alter them silently. First examples: {preview}"
        )
    return model

def _transform_for_pca(X: np.ndarray, mode: str = "none"):
    mode = (mode or "none").lower()
    Y = X.copy().astype(float)
    def sign_log1p(a):
        return np.sign(a) * np.log1p(np.abs(a))
    if mode in ("log", "log1p", "log-scale", "log_scale"):
        return sign_log1p(Y)
    if mode in ("zscore", "standardize", "standardise", "standardize_features"):
        mu = np.nanmean(Y, axis=0); sd = np.nanstd(Y, axis=0); sd[sd == 0] = 1.0
        return (Y - mu) / sd
    if mode in ("log_zscore", "log1p_zscore", "log_then_zscore"):
        Y = sign_log1p(Y)
        mu = np.nanmean(Y, axis=0); sd = np.nanstd(Y, axis=0); sd[sd == 0] = 1.0
        return (Y - mu) / sd
    return Y

def analyze_flux_distributions(model, solutions_dict, baseline_key, results_dir="results", pca_scale="none", write_replicates_long=False, logger=print, replicate_labels=None):
    os.makedirs(results_dir, exist_ok=True)
    out_dir = os.path.join(results_dir, "flux_analysis")
    os.makedirs(out_dir, exist_ok=True)
    cond_names = list(solutions_dict.keys())
    flux_data = {cond: np.array([sol.fluxes.values for sol in sols]) for cond, sols in solutions_dict.items()}
    rep_name_map = {}
    for cond in cond_names:
        n = flux_data[cond].shape[0]
        if replicate_labels and cond in replicate_labels and len(replicate_labels[cond]) == n:
            rep_name_map[cond] = list(replicate_labels[cond])
        else:
            rep_name_map[cond] = [f"{cond}_rep{i+1}" for i in range(n)]
    cols = ["ReactionID", "ReactionName", "Subsystem"]
    for cond in cond_names:
        cols += [f"{cond}_MeanFlux", f"{cond}_StdFlux", f"{cond}_N"]
        for rn in rep_name_map[cond]:
            cols.append(f"{rn}_Flux")
    for c in [x for x in cond_names if x != baseline_key]:
        cols += [f"Diff({c}-{baseline_key})", f"Ratio({c}/{baseline_key})"]
    df = pd.DataFrame(columns=cols)
    for rxn in model.reactions:
        i = model.reactions.index(rxn)
        row = {"ReactionID": rxn.id, "ReactionName": rxn.name, "Subsystem": rxn.subsystem or ""}
        for cond in cond_names:
            arr = flux_data[cond][:, i]
            row[f"{cond}_MeanFlux"] = float(np.mean(arr)) if arr.size else float('nan')
            row[f"{cond}_StdFlux"]  = float(np.std(arr)) if arr.size else float('nan')
            row[f"{cond}_N"]        = int(arr.shape[0])
            for r, v in enumerate(arr, start=0):
                row[f"{rep_name_map[cond][r]}_Flux"] = float(v)
        ref_vals = flux_data[baseline_key][:, i]
        ref_mean = float(np.mean(ref_vals)) if ref_vals.size else float('nan')
        ref_mean_abs = float(np.mean(np.abs(ref_vals)) + 1e-12) if ref_vals.size else float('nan')
        for c in [x for x in cond_names if x != baseline_key]:
            c_vals = flux_data[c][:, i]
            c_mean = float(np.mean(c_vals)) if c_vals.size else float('nan')
            row[f"Diff({c}-{baseline_key})"] = c_mean - ref_mean if (not np.isnan(c_mean) and not np.isnan(ref_mean)) else float('nan')
            row[f"Ratio({c}/{baseline_key})"] = (float(np.mean(np.abs(c_vals))) / ref_mean_abs) if (c_vals.size and not np.isnan(ref_mean_abs) and ref_mean_abs != 0.0) else float('nan')
        df.loc[len(df)] = row
    out_csv = os.path.join(out_dir, "reaction_flux_comparison_extended.csv")
    df.to_csv(out_csv, index=False); logger(f"[INFO] Extended flux comparison saved to: {out_csv}")
    if write_replicates_long:
        long_records = []
        for cond, arr in flux_data.items():
            for r in range(arr.shape[0]):
                rep_label = rep_name_map[cond][r]
                for i, rxn in enumerate(model.reactions):
                    long_records.append({
                        "ReactionID": rxn.id, "ReactionName": rxn.name, "Subsystem": rxn.subsystem or "",
                        "Group": cond, "Replicate": r + 1, "ReplicateName": rep_label, "Flux": float(arr[r, i]),
                    })
        long_df = pd.DataFrame(long_records)
        long_csv = os.path.join(out_dir, "flux_replicates_long.csv")
        long_df.to_csv(long_csv, index=False); logger(f"[INFO] Replicate-level long table saved to: {long_csv}")
    all_vectors, all_labels = [], []
    for cond, arr in flux_data.items():
        for r in range(arr.shape[0]):
            all_vectors.append(arr[r, :])
            label = rep_name_map[cond][r] if arr.shape[0] > 0 else cond
            all_labels.append(label)
    all_vectors = np.array(all_vectors)
    if all_vectors.shape[0] >= 2:
        Xp_in = _transform_for_pca(all_vectors, mode=pca_scale)
        pca = PCA(n_components=2)
        Xp = pca.fit_transform(Xp_in)
        plt.figure(figsize=(7, 6))
        plt.scatter(Xp[:, 0], Xp[:, 1], s=50, alpha=0.8)
        for idx, label in enumerate(all_labels):
            plt.text(Xp[idx, 0], Xp[idx, 1], label, fontsize=8)
        plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        plt.title("PCA of Flux Distributions")
        plt.grid(True)
        pca_fig = os.path.join(out_dir, "flux_pca.png")
        plt.savefig(pca_fig); plt.close()
        logger(f"[INFO] PCA plot saved to: {pca_fig}")

# =====================================
# Stats with BH-FDR
# =====================================
def bh_fdr(pvals: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(pvals), dtype=float)
    m = len(p)
    order = np.argsort(p)
    ranks = np.arange(1, m+1, dtype=float)
    q = np.empty_like(p)
    q[order] = p[order] * (m / ranks)
    # monotone
    for i in range(m-2, -1, -1):
        q[order[i]] = min(q[order[i]], q[order[i+1]])
    q = np.clip(q, 0.0, 1.0)
    return q


def _cohen_d(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return float("nan")
    nx, ny = len(x), len(y)
    dof = nx + ny - 2
    pooled_var = ((nx - 1) * np.var(x, ddof=1) + (ny - 1) * np.var(y, ddof=1)) / max(dof, 1)
    diff = float(np.mean(x) - np.mean(y))
    if pooled_var <= 0:
        if abs(diff) < 1e-15:
            return 0.0
        return float(np.sign(diff) * np.inf)
    return diff / np.sqrt(pooled_var)


def compare_flux_statistically(solutions_dict, model, results_dir="results", test_type="t-test", logger=print):
    """
    Pairwise per-reaction tests.

    This table is suitable as primary inference for a single-study input
    (e.g. one RQ2 strain from GSE182668). For multi-dataset RQ1 input it is
    explicitly labelled UNADJUSTED; the revised comprehensive_flux_analysis.py
    performs the primary dataset-aware inference on the original feasible fluxes.

    BH-FDR is applied separately within each contrast (not across all contrasts).
    """
    out_dir = os.path.join(results_dir, "stats_comparison")
    os.makedirs(out_dir, exist_ok=True)
    cond_names = list(solutions_dict.keys())
    flux_data = {cond: np.array([sol.fluxes.values for sol in sols]) for cond, sols in solutions_dict.items()}

    rep_labels = globals().get('_GLOBAL_REP_NAME_MAP_FOR_STATS', None) or {}
    ds_re = re.compile(r"(?:GSE|GSM)\d+", re.I)
    datasets = set()
    for cond, names in rep_labels.items():
        for nm in names:
            m = ds_re.search(str(nm))
            if m:
                datasets.add(m.group(0).upper())
    multi_dataset = len({d for d in datasets if d.startswith("GSE")}) > 1

    records = []
    for c_idx, cA in enumerate(cond_names):
        for cB in cond_names[c_idx+1:]:
            contrast_rows = []
            for i, rxn in enumerate(model.reactions):
                arrA = flux_data[cA][:, i]
                arrB = flux_data[cB][:, i]
                arrA = arrA[np.isfinite(arrA)]
                arrB = arrB[np.isfinite(arrB)]
                if len(arrA) < 2 or len(arrB) < 2:
                    pval = float('nan'); stat = float('nan')
                else:
                    with np.errstate(all="ignore"):
                        if test_type.lower() == "wilcoxon":
                            stat, pval = ranksums(arrA, arrB)
                        elif test_type.lower() == "mann-whitney":
                            stat, pval = mannwhitneyu(arrA, arrB, alternative='two-sided')
                        else:
                            stat, pval = ttest_ind(arrA, arrB, equal_var=False, nan_policy="omit")
                mA = float(np.mean(arrA)) if arrA.size else float('nan')
                mB = float(np.mean(arrB)) if arrB.size else float('nan')
                contrast_rows.append({
                    "ReactionID": rxn.id,
                    "ReactionName": rxn.name,
                    "Subsystem": rxn.subsystem or "",
                    "GroupA": cA,
                    "GroupB": cB,
                    "MeanA": mA,
                    "MeanB": mB,
                    "Diff(A-B)": mA - mB if np.isfinite(mA) and np.isfinite(mB) else float('nan'),
                    "Cohen_d": _cohen_d(arrA, arrB),
                    "TestStat": float(stat) if np.isfinite(stat) else stat,
                    "PValue": float(pval) if np.isfinite(pval) else pval,
                    "N_A": int(len(arrA)),
                    "N_B": int(len(arrB)),
                    "InferenceScope": "UNADJUSTED_MULTI_DATASET" if multi_dataset else "SINGLE_DATASET"
                })
            cdf = pd.DataFrame(contrast_rows)
            valid = cdf["PValue"].notna() & np.isfinite(cdf["PValue"].astype(float))
            cdf["FDR_BH"] = np.nan
            if valid.any():
                cdf.loc[valid, "FDR_BH"] = bh_fdr(cdf.loc[valid, "PValue"].astype(float).values)
            records.append(cdf)

    df = pd.concat(records, ignore_index=True) if records else pd.DataFrame()
    out_name = "flux_pairwise_stats_unadjusted.csv" if multi_dataset else "flux_pairwise_stats.csv"
    out_csv = os.path.join(out_dir, out_name)
    df.to_csv(out_csv, index=False)
    if multi_dataset:
        logger("[WARN] Multi-dataset input detected: pooled pairwise tests are saved only as UNADJUSTED diagnostics. "
               "Use comprehensive_flux_analysis.py for dataset-aware RQ1 inference.")
    logger(f"[INFO] Pairwise stats saved to: {out_csv}")

    all_vectors, all_names = [], []
    for cond in cond_names:
        arr = flux_data[cond]
        names_list = rep_labels.get(cond, [f"{cond}_rep{r+1}" for r in range(arr.shape[0])])
        for r in range(arr.shape[0]):
            all_vectors.append(arr[r, :])
            all_names.append(names_list[r])
    all_vectors = np.asarray(all_vectors, dtype=float)
    if all_vectors.shape[0] > 1:
        dmat = squareform(pdist(all_vectors, metric='euclidean'))
        pd.DataFrame(dmat, index=all_names, columns=all_names).to_csv(
            os.path.join(out_dir, "flux_global_distance_matrix.csv")
        )
    return df

def compute_rank_product(flux_csv: str, baseline: str, targets: List[str], out_csv: str, logger=print):
    df = pd.read_csv(flux_csv)
    df = df.set_index("ReactionID")
    diff_cols = [c for c in df.columns if c.startswith("Diff(")]
    valid_diff_cols = [c for c in diff_cols if any(t in c for t in targets)]
    if not valid_diff_cols:
        logger(f"[WARN] No diff columns found for targets={targets} in flux_csv.")
        return pd.DataFrame()
    sub_df = df[valid_diff_cols].copy()
    # Rank each diff col
    ranks = {}
    for col in sub_df.columns:
        # s = sub_df[col].dropna() 
        s = sub_df[col].dropna().abs()        # If your goal was to find the most responsive reactions in either direction (up or down)
        ranked = s.rank(method="average", ascending=False) #  you might need a second rank product run for downregulated reactions using ascending=True
        ranks[col] = ranked
    rank_df = pd.DataFrame(ranks)
    # Geometric mean of ranks
    rp = rank_df.apply(lambda row: np.exp(np.log(row + 1e-9).mean()), axis=1)
    rp_sorted = rp.sort_values(ascending=True)
    rp_df = pd.DataFrame({
        "ReactionID": rp_sorted.index,
        "RankProduct": rp_sorted.values,
        "RankProductRank": range(1, len(rp_sorted)+1)
    })
    for col in valid_diff_cols:
        rp_df[col] = df.loc[rp_sorted.index, col].values
    rp_df.to_csv(out_csv, index=False)
    logger(f"[INFO] Rank-product saved to: {out_csv}")
    return rp_df

# =====================================
# Extended visuals
# =====================================
def extended_visuals(model, solutions_dict, baseline_key, results_dir="results", logger=print):
    out_dir = os.path.join(results_dir, "extended_visuals")
    os.makedirs(out_dir, exist_ok=True)
    cond_names = list(solutions_dict.keys())
    flux_data = {cond: np.array([sol.fluxes.values for sol in sols]) for cond, sols in solutions_dict.items()}
    # Heatmap
    mean_fluxes = []
    for cond in cond_names:
        arr = flux_data[cond]
        mean_fluxes.append(np.mean(arr, axis=0))
    heatmap_data = np.array(mean_fluxes).T
    fig, ax = plt.subplots(figsize=(8, 10))
    cax = ax.imshow(heatmap_data, aspect='auto', cmap='viridis', interpolation='nearest')
    ax.set_xlabel("Condition")
    ax.set_ylabel("Reaction Index")
    ax.set_title("Mean Flux Heatmap")
    ax.set_xticks(range(len(cond_names)))
    ax.set_xticklabels(cond_names)
    fig.colorbar(cax, ax=ax, label="Flux")
    out_fig = os.path.join(out_dir, "mean_flux_heatmap.png")
    plt.savefig(out_fig, dpi=150, bbox_inches='tight'); plt.close()
    logger(f"[INFO] Mean flux heatmap saved to: {out_fig}")

# =====================================
# Cytoscape edges
# =====================================
def build_cytoscape_edges_for_comparison(model, diff_series, ratio_series, cond, baseline, abs_diff_threshold=0.0, out_path=None, logger=print):
    edges = []
    for rxn in model.reactions:
        rid = rxn.id
        if rid not in diff_series.index:
            continue
        diff_val = diff_series.loc[rid]
        if np.isnan(diff_val) or abs(diff_val) < abs_diff_threshold:
            continue
        ratio_val = ratio_series.loc[rid] if (rid in ratio_series.index) else float('nan')
        # Build edges from reaction metabolites
        for met, coeff in rxn.metabolites.items():
            mid = met.id
            edges.append({
                "ReactionID": rid,
                "ReactionName": rxn.name or "",
                "Subsystem": rxn.subsystem or "",
                "MetaboliteID": mid,
                "MetaboliteName": met.name or "",
                "Coefficient": float(coeff),
                f"Diff({cond}-{baseline})": float(diff_val),
                f"Ratio({cond}/{baseline})": float(ratio_val)
            })
    df = pd.DataFrame(edges)
    if out_path:
        df.to_csv(out_path, index=False)
        logger(f"[INFO] Cytoscape edges ({len(df)} rows) saved to: {out_path}")
    return df

# =====================================
# Report files
# =====================================
def write_report_files(report: dict, results_dir: str):
    report_json = os.path.join(results_dir, "analysis_report.json")
    with open(report_json, "w") as f:
        json.dump(report, f, indent=2)
    print(f"[INFO] Analysis report written to: {report_json}")

# =====================================
# Main pipeline
# =====================================

# --------------------------------------------------------------------------- #
# Solver step: FBA vs pFBA
#
# WHY THIS EXISTS
# Plain FBA returns SOME optimal solution, not a unique one. When the optimum is
# degenerate (an optimal face rather than a single vertex) the reported per-reaction
# fluxes depend on which vertex the solver happens to land on, even though the
# objective value is identical every time. FVA on this model shows ~68% of reactions
# have a non-zero feasible range at the optimum, so this is not hypothetical.
#
# pFBA removes that ambiguity: it fixes the biological objective at its optimum
# (or a fraction of it) and then minimises total absolute flux, which selects one
# specific, reproducible solution.
#
# IMPORTANT — objective_value semantics:
# cobra's pfba() returns a Solution whose .objective_value is the MINIMISED TOTAL
# FLUX, not the biomass/objective flux. Every downstream consumer here (console
# output, analysis_report.json, group means) expects the biological objective. So
# after a pFBA solve we capture the biological objective separately and write it
# back into .objective_value, keeping the total flux available as .total_flux.
# Without this, "objective_values" in the report would silently become a completely
# different quantity and would not be comparable to any previous run.
# --------------------------------------------------------------------------- #

def configure_solver(model, solver="gurobi", deterministic=True, logger=print):
    """Select solver and, when available, apply deterministic Gurobi LP settings."""
    try:
        model.solver = solver
    except Exception as exc:
        logger(f"[WARN] Could not assign solver '{solver}' directly: {exc}")
    if deterministic and "gurobi" in str(model.solver.interface.__name__).lower():
        params = {
            "Threads": 1,
            "Seed": 1,
            "NumericFocus": 3,
            "Method": 1,
        }
        problem = getattr(model.solver, "problem", None)
        gparams = getattr(problem, "Params", None)
        if gparams is not None:
            for key, value in params.items():
                try:
                    setattr(gparams, key, value)
                except Exception as exc:
                    logger(f"[WARN] Could not set Gurobi parameter {key}={value}: {exc}")
            logger("[INFO] Deterministic Gurobi settings applied: Threads=1, Seed=1, NumericFocus=3, Method=1.")
    return model


def solve_flux(mdl, solve_mode="pfba", pfba_fraction=1.0):
    """
    Solve one constrained model.

    pFBA first maximizes the declared biological objective and then minimizes
    total absolute flux while retaining the requested fraction of that optimum.
    pFBA reduces solution degeneracy but does not guarantee mathematical
    uniqueness; FVA is used separately to quantify remaining alternative optima.
    """
    if not (0 < float(pfba_fraction) <= 1.0):
        raise ValueError("pfba_fraction must satisfy 0 < fraction <= 1.")

    if solve_mode == "fba":
        sol = mdl.optimize()
        if getattr(sol, "status", None) != "optimal":
            raise RuntimeError(f"FBA did not reach optimality (status={getattr(sol, 'status', None)}).")
        return sol
    if solve_mode != "pfba":
        raise ValueError(f"unknown solve_mode {solve_mode!r}; expected 'fba' or 'pfba'")

    biological_objective = mdl.slim_optimize(error_value=np.nan)
    if not np.isfinite(biological_objective):
        raise RuntimeError("Biological objective optimization failed before pFBA.")
    sol = _cobra_pfba(mdl, fraction_of_optimum=float(pfba_fraction))
    if getattr(sol, "status", None) != "optimal":
        raise RuntimeError(f"pFBA did not reach optimality (status={getattr(sol, 'status', None)}).")
    # COBRApy's pFBA Solution objective may refer to the secondary total-flux
    # objective. Preserve both quantities explicitly.
    secondary_total_flux = float(sol.objective_value)
    try:
        sol.total_flux = secondary_total_flux
        sol.biological_objective_value = float(biological_objective)
        sol.objective_value = float(biological_objective)
    except Exception:
        pass
    return sol


def _load_fva_targets(model, targets=None, targets_file=None):
    ids = []
    if targets:
        ids.extend([x.strip() for x in re.split(r"[,;\s]+", str(targets)) if x.strip()])
    if targets_file:
        tdf = pd.read_csv(targets_file)
        if "ReactionID" in tdf.columns:
            ids.extend(tdf["ReactionID"].dropna().astype(str).tolist())
        elif len(tdf.columns):
            ids.extend(tdf.iloc[:, 0].dropna().astype(str).tolist())
    if not ids:
        return None
    available = {r.id for r in model.reactions}
    unique = []
    seen = set()
    for rid in ids:
        if rid in available and rid not in seen:
            unique.append(rid); seen.add(rid)
    missing = sorted(set(ids) - available)
    if missing:
        print(f"[WARN] {len(missing)} requested FVA reactions are absent from the model; skipped.")
    return unique


def _run_fva_and_write(mdl, name, results_dir, fraction, reaction_ids=None, processes=1):
    out_dir = os.path.join(results_dir, "fva")
    os.makedirs(out_dir, exist_ok=True)
    kwargs = {"fraction_of_optimum": float(fraction), "processes": int(max(1, processes))}
    if reaction_ids:
        kwargs["reaction_list"] = reaction_ids
    res = flux_variability_analysis(mdl, **kwargs)
    res = res.copy()
    res.index.name = "ReactionID"
    res["range"] = res["maximum"] - res["minimum"]
    out = os.path.join(out_dir, f"FVA_{name}.csv")
    res.reset_index().to_csv(out, index=False)
    print(f"[INFO] FVA saved: {out} ({len(res)} reactions; fraction={fraction})")
    return out


def run_pipeline(
    rnaseq_file,
    met_data_file=None,
    model_file="iMM1415.json",
    results_dir="results",
    aggregate=False,
    infer_groups_from_filename=True,
    baseline_code=None,
    explicit_groups=None,
    mapping_file=None,
    column_regex=None,
    eflux_quantile=0.95,
    eflux_floor=0.1,
    eflux_cap=1000.0,
    diet_bounds_json=None,
    diet_bounds_units="model",
    mw_table=None,
    gDW=1.0,
    hours_per_day=24.0,
    pca_scale="none",
    write_replicates_long=False,
    test_type="t-test",
    objective_id=None,
    objective_regex=None,
    objective_sense="max",
    transporter_strategy="e_to_non_e",
    transporter_regex=None,
    transporter_subsystem_regex=None,
    transporter_compartments="e",
    edge_abs_diff_threshold=0.0,
    solve_mode="pfba",
    pfba_fraction=1.0,
    solver="gurobi",
    deterministic_solver=True,
    fva_mode="none",
    fva_fraction=1.0,
    fva_targets=None,
    fva_targets_file=None,
    fva_processes=1,
):
    if fva_mode not in {"none", "representative", "per_sample"}:
        raise ValueError("fva_mode must be one of: none, representative, per_sample")
    if not (0 < float(fva_fraction) <= 1.0):
        raise ValueError("fva_fraction must satisfy 0 < fraction <= 1.")

    symbol_to_entrez = None
    if mapping_file:
        if not os.path.exists(mapping_file):
            raise FileNotFoundError(f"mapping_file not found: {mapping_file}")
        print(f"[INFO] Loading gene symbol-to-entrez mapping from: {mapping_file}")
        symbol_to_entrez = load_symbol_to_entrez_mapping(mapping_file)
        print(f"[INFO] Loaded {len(symbol_to_entrez)} symbol->entrez mappings")

    if infer_groups_from_filename:
        tags = parse_groups_from_filename_multi(rnaseq_file)
    elif explicit_groups:
        tags = [canonical_code(x) for x in re.split(r"[,\s]+", explicit_groups) if x]
    else:
        raise ValueError("Must infer groups from filename or supply explicit_groups.")
    baseline = canonical_code(baseline_code) if baseline_code else choose_default_baseline(tags)
    if baseline not in tags:
        raise ValueError(f"Baseline '{baseline}' is not among detected groups {tags}.")
    print(f"[INFO] Detected groups: {tags} (baseline={baseline})")

    os.makedirs(results_dir, exist_ok=True)
    report = {
        "groups": {"all_groups": tags, "baseline": baseline, "sample_counts": {}},
        "rna_seq": {"file": rnaseq_file, "total_genes": 0, "genes_with_expr": 0},
        "parameters": {
            "aggregate": aggregate,
            "eflux_quantile": eflux_quantile,
            "eflux_floor": eflux_floor,
            "eflux_cap": eflux_cap,
            "diet_bounds_json": diet_bounds_json,
            "diet_bounds_units": diet_bounds_units,
            "test_type": test_type,
            "objective_id": objective_id,
            "objective_regex": objective_regex,
            "objective_sense": objective_sense,
            "transporter_strategy": transporter_strategy,
            "solve_mode": solve_mode,
            "pfba_fraction": pfba_fraction,
            "solver": solver,
            "deterministic_solver": deterministic_solver,
            "fva_mode": fva_mode,
            "fva_fraction": fva_fraction,
            "fva_targets_file": fva_targets_file,
        },
        "per_group": {}, "global": {}, "outputs": {}
    }

    print("[STEP] Loading model...")
    model = load_json_model(model_file)
    model = configure_solver(model, solver=solver, deterministic=deterministic_solver)
    model, obj_id = set_objective_reaction(
        model, objective_id=objective_id, objective_regex=objective_regex,
        sense=objective_sense
    )
    report["model"] = {
        "num_reactions": len(model.reactions),
        "num_metabolites": len(model.metabolites),
        "num_genes": len(model.genes),
        "objective_id": obj_id,
        "solver": str(model.solver.interface.__name__),
    }

    rna = pd.read_csv(rnaseq_file)
    if "Gene_Symbol" not in rna.columns:
        raise ValueError("RNA-seq CSV must contain a 'Gene_Symbol' column.")
    gene_symbols = rna["Gene_Symbol"].astype(str).str.lower().tolist()
    gene_ids_for_matching = gene_symbols
    if symbol_to_entrez:
        gene_ids_for_matching = [
            symbol_to_entrez.get(sym.lower().strip(), sym) for sym in gene_symbols
        ]

    # Identify sample columns from their diet labels, NOT from the proportion of
    # non-missing gene values.  The merged RQ1 matrix is an outer merge of six GEO
    # studies and therefore legitimately contains dataset-specific NA values.
    # A completeness filter (e.g. >=50% non-missing genes) silently drops entire
    # cohorts/diets and must not be used here.
    annotation_cols = {"Gene_Symbol", "Gene_ID"}
    candidate_cols = [c for c in rna.columns if c not in annotation_cols]

    def _column_matches_tag(col, tag):
        aliases = GROUP_ALIASES.get(tag, {tag}) | {tag}
        alias_alt = "(?:" + "|".join(sorted(re.escape(a) for a in aliases)) + ")"
        if column_regex:
            pat = column_regex.replace("{TAG}", re.escape(tag)).replace("{ALIASES}", alias_alt)
        else:
            pat = rf"(?<![A-Za-z0-9]){alias_alt}(?![A-Za-z0-9])"
        return re.search(pat, str(col), flags=re.IGNORECASE) is not None

    tag_to_cols = {tag: [c for c in candidate_cols if _column_matches_tag(c, tag)] for tag in tags}

    # A sample column must map to exactly one requested diet group.
    assigned = {}
    for tag, cols in tag_to_cols.items():
        for c in cols:
            if c in assigned and assigned[c] != tag:
                raise ValueError(
                    f"Sample column '{c}' matched multiple diet groups: {assigned[c]} and {tag}. "
                    "Use --column_regex to disambiguate."
                )
            assigned[c] = tag

    sample_cols = [c for c in candidate_cols if c in assigned]
    if not sample_cols:
        raise ValueError("No sample-expression columns were detected from the requested diet labels.")

    # Convert expression columns to numeric while PRESERVING NA as missing/unmeasured.
    # Do not fill NA with zero: the subsequent E-Flux floor would turn a zero into
    # a positive capacity and would incorrectly treat missing measurements as low
    # measured expression.  NA genes are instead left unconstrained for that sample.
    expr_df = rna[sample_cols].apply(pd.to_numeric, errors="coerce")
    numeric_counts = expr_df.notna().sum(axis=0)
    empty_cols = numeric_counts[numeric_counts == 0].index.tolist()
    if empty_cols:
        raise ValueError(f"Sample-expression columns contain no numeric values: {empty_cols[:5]}")
    expr = expr_df.to_numpy(dtype=float)

    coverage = numeric_counts / max(len(rna), 1)
    print(
        "[INFO] Per-sample numeric gene coverage: "
        f"min={coverage.min():.1%}, median={coverage.median():.1%}, max={coverage.max():.1%}. "
        "Missing values are preserved as unmeasured (not converted to zero)."
    )

    report["rna_seq"]["total_genes"] = int(len(gene_symbols))
    report["rna_seq"]["genes_with_expr"] = int((np.nansum(expr, axis=1) > 0).sum())
    report["rna_seq"]["sample_numeric_coverage"] = {
        str(c): float(coverage.loc[c]) for c in sample_cols
    }

    tag_to_indices = {tag: [sample_cols.index(c) for c in tag_to_cols[tag]] for tag in tags}
    for tag, idxs in tag_to_indices.items():
        if not idxs:
            raise ValueError(f"No sample columns found for group '{tag}'.")
        report["groups"]["sample_counts"][tag] = len(idxs)
        print(f"[INFO] {tag}: {len(idxs)} sample columns")

    unmatched = [c for c in candidate_cols if c not in assigned]
    if unmatched:
        print(f"[INFO] Ignoring {len(unmatched)} non-sample/unmatched columns: {unmatched[:5]}")

    def _nanmean_expression(matrix, indices):
        """Mean expression across replicates using available observations only.

        Genes missing in every replicate remain NaN (unmeasured) rather than being
        converted to zero/low expression. This is used only for aggregate mode and
        representative FVA models; primary RQ1/RQ2 inference remains replicate-level.
        """
        sub = np.asarray(matrix[:, indices], dtype=float)
        counts = np.isfinite(sub).sum(axis=1)
        sums = np.nansum(sub, axis=1)
        out = np.full(sub.shape[0], np.nan, dtype=float)
        np.divide(sums, counts, out=out, where=counts > 0)
        return out

    diet_bounds = None
    if diet_bounds_json:
        with open(diet_bounds_json, "r", encoding="utf-8") as handle:
            raw = json.load(handle)
        diet_bounds = {
            canonical_code(k): {rid: [float(v[0]), float(v[1])] for rid, v in d.items()}
            for k, d in raw.items()
        }
    if diet_bounds_units not in {"model", "g_per_day", "mmol_per_day"}:
        raise ValueError("diet_bounds_units must be model, g_per_day, or mmol_per_day")
    mw_map = load_mw_table(mw_table) if (mw_table and diet_bounds_units == "g_per_day") else {}

    comps = tuple(c.strip() for c in str(transporter_compartments).split(",") if c.strip())
    ex_set, trans_set, int_set = classify_reactions(
        model,
        transporter_strategy=transporter_strategy,
        transporter_regex=transporter_regex,
        transporter_subsystem_regex=transporter_subsystem_regex,
        compartments_for_transport=comps,
    )
    print(f"[INFO] Reaction classes: EX={len(ex_set)} | Transporters={len(trans_set)} | Internal={len(int_set)}")

    def build_constrained_model(vec, tag):
        mdl = model.copy()
        mdl = configure_solver(mdl, solver=solver, deterministic=deterministic_solver, logger=lambda *_: None)
        mdl, L1_count = apply_diet_bounds_layer1(
            mdl, code=tag, diet_bounds=diet_bounds, diet_units=diet_bounds_units,
            mw_map=mw_map, gDW=gDW, hours_per_day=hours_per_day
        )
        mdl, L2_changed, _, _, genes_set2, rxn_set_L2 = apply_expression_constraints_scoped(
            mdl, gene_ids_for_matching, vec, trans_set,
            eflux_quantile=eflux_quantile, eflux_floor=eflux_floor,
            eflux_cap=eflux_cap, label="L2", symbol_to_entrez=symbol_to_entrez
        )
        mdl, L3_changed, _, _, genes_set3, rxn_set_L3 = apply_expression_constraints_scoped(
            mdl, gene_ids_for_matching, vec, int_set,
            eflux_quantile=eflux_quantile, eflux_floor=eflux_floor,
            eflux_cap=eflux_cap, label="L3", symbol_to_entrez=symbol_to_entrez
        )
        mdl = validate_model(mdl)
        info = {
            "L1": int(L1_count), "L2": int(L2_changed), "L3": int(L3_changed),
            "genes": set(genes_set2) | set(genes_set3),
            "rxn_L2": set(rxn_set_L2), "rxn_L3": set(rxn_set_L3)
        }
        return mdl, info

    fva_target_ids = _load_fva_targets(model, fva_targets, fva_targets_file)
    solutions_dict = {}
    replicate_labels = {}

    for tag in tags:
        idxs = tag_to_indices[tag]
        reps = [_nanmean_expression(expr, idxs)] if aggregate else [expr[:, si] for si in idxs]
        labels = [tag] if aggregate else [str(sample_cols[si]) for si in idxs]
        sols, obj_vals = [], []
        rxn_union_L2, rxn_union_L3, gene_union = set(), set(), set()

        for rep_i, (vec, rep_label) in enumerate(zip(reps, labels), start=1):
            mdl, info = build_constrained_model(vec, tag)
            sol = solve_flux(mdl, solve_mode=solve_mode, pfba_fraction=pfba_fraction)
            sols.append(sol)
            obj_vals.append(float(sol.objective_value))
            rxn_union_L2 |= info["rxn_L2"]; rxn_union_L3 |= info["rxn_L3"]; gene_union |= info["genes"]
            print(
                f"[{tag}] {rep_label} ({solve_mode}) objective={sol.objective_value:.6g} | "
                f"L1={info['L1']} | L2={info['L2']} | L3={info['L3']}"
            )
            if fva_mode == "per_sample":
                safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", rep_label)
                _run_fva_and_write(
                    mdl, safe, results_dir, fva_fraction,
                    reaction_ids=fva_target_ids, processes=fva_processes
                )

        solutions_dict[tag] = sols
        replicate_labels[tag] = labels
        report["per_group"][tag] = {
            "samples": int(len(idxs)),
            "modeled_flux_vectors": int(len(sols)),
            "objective_values": obj_vals,
            "objective_mean": float(np.mean(obj_vals)) if obj_vals else None,
            "objective_std": float(np.std(obj_vals, ddof=1)) if len(obj_vals) > 1 else 0.0,
            "rxns_constrained_unique_L2": len(rxn_union_L2),
            "rxns_constrained_unique_L3": len(rxn_union_L3),
            "genes_mapped_unique": len(gene_union),
        }

        if fva_mode == "representative":
            rep_vec = _nanmean_expression(expr, idxs)
            rep_mdl, _ = build_constrained_model(rep_vec, tag)
            # Confirm feasibility/optimum before FVA.
            _ = solve_flux(rep_mdl, solve_mode=solve_mode, pfba_fraction=pfba_fraction)
            _run_fva_and_write(
                rep_mdl, f"{tag}_representative", results_dir, fva_fraction,
                reaction_ids=fva_target_ids, processes=fva_processes
            )

    globals()['_GLOBAL_REP_NAME_MAP_FOR_STATS'] = replicate_labels

    analyze_flux_distributions(
        model, solutions_dict, baseline_key=baseline, results_dir=results_dir,
        pca_scale=pca_scale, write_replicates_long=write_replicates_long,
        replicate_labels=replicate_labels
    )
    flux_table_csv = os.path.join(results_dir, "flux_analysis", "reaction_flux_comparison_extended.csv")
    report["outputs"]["flux_comparison_csv"] = flux_table_csv
    report["outputs"]["pca_png"] = os.path.join(results_dir, "flux_analysis", "flux_pca.png")

    extended_visuals(model, solutions_dict, baseline_key=baseline, results_dir=results_dir)
    stats_df = compare_flux_statistically(
        solutions_dict, model, results_dir=results_dir, test_type=test_type
    )
    report["outputs"]["stats_dir"] = os.path.join(results_dir, "stats_comparison")

    # Descriptive Cytoscape edges from raw mean differences. Primary significant
    # RQ1 networks are generated later by comprehensive_flux_analysis.py.
    try:
        ft = pd.read_csv(flux_table_csv).set_index("ReactionID")
        edge_dir = os.path.join(results_dir, "cytoscape_edges")
        os.makedirs(edge_dir, exist_ok=True)
        for cond in [t for t in tags if t != baseline]:
            diff_col = f"Diff({cond}-{baseline})"
            ratio_col = f"Ratio({cond}/{baseline})"
            if diff_col not in ft.columns:
                continue
            ratio_series = ft[ratio_col] if ratio_col in ft.columns else pd.Series(index=ft.index, data=np.nan)
            build_cytoscape_edges_for_comparison(
                model, ft[diff_col], ratio_series, cond, baseline,
                abs_diff_threshold=float(edge_abs_diff_threshold),
                out_path=os.path.join(edge_dir, f"edges_{cond}_vs_{baseline}.csv")
            )
        report["outputs"]["cytoscape_edge_dir"] = edge_dir
    except Exception as exc:
        print(f"[WARN] Descriptive Cytoscape edge build failed: {exc}")

    write_report_files(report, results_dir)
    print("[INFO] Layered replicate-level modeling complete.")
    return report


def main():
    p = argparse.ArgumentParser(
        description="Three-layer E-Flux modeling with replicate-preserving pFBA and optional FVA."
    )
    p.add_argument("rnaseq_file")
    p.add_argument("met_data_file", nargs="?", default=None)
    p.add_argument("--model_file", default="iMM1415.json")
    p.add_argument("--results_dir", default="results")
    p.add_argument("--aggregate", action="store_true",
                   help="Average expression within groups before modeling. Do NOT use for primary RQ1/RQ2 inference.")
    p.add_argument("--pca_scale", default="none",
                   choices=["none","zscore","log","log_zscore","log1p","log1p_zscore"])
    p.add_argument("--write_replicates_long", action="store_true")
    p.add_argument("--infer_groups_from_filename", action="store_true", default=True)
    p.add_argument("--baseline_code", default=None)
    p.add_argument("--explicit_groups", default=None)
    p.add_argument("--mapping_file", default=None)
    p.add_argument("--column_regex", default=None)
    p.add_argument("--eflux_quantile", type=float, default=0.95)
    p.add_argument("--eflux_floor", type=float, default=0.1)
    p.add_argument("--eflux_cap", type=float, default=1000.0)
    p.add_argument("--diet_bounds_json", default=None)
    p.add_argument("--diet_bounds_units", default="model",
                   choices=["model", "g_per_day", "mmol_per_day"])
    p.add_argument("--mw_table", default=None)
    p.add_argument("--gDW", type=float, default=1.0)
    p.add_argument("--hours_per_day", type=float, default=24.0)
    p.add_argument("--test_type", default="t-test",
                   choices=["t-test", "wilcoxon", "mann-whitney"])
    p.add_argument("--objective_id", default=None)
    p.add_argument("--objective_regex", default=None)
    p.add_argument("--objective_sense", default="max", choices=["max","min"])
    p.add_argument("--transporter_strategy", default="e_to_non_e",
                   choices=["e_to_non_e","regex","either"])
    p.add_argument("--transporter_regex", default=None)
    p.add_argument("--transporter_subsystem_regex", default=None)
    p.add_argument("--transporter_compartments", default="e")
    p.add_argument("--edge_abs_diff_threshold", type=float, default=0.0)
    p.add_argument("--solve_mode", choices=["fba", "pfba"], default="pfba",
                   help="Primary revision default is pFBA. Use fba only for sensitivity comparison.")
    p.add_argument("--pfba_fraction", type=float, default=1.0)
    p.add_argument("--solver", default="gurobi")
    p.add_argument("--non_deterministic_solver", action="store_true",
                   help="Disable deterministic Gurobi settings.")
    p.add_argument("--fva_mode", choices=["none","representative","per_sample"], default="none")
    p.add_argument("--fva_fraction", type=float, default=1.0)
    p.add_argument("--fva_targets", default=None,
                   help="Optional comma/space-separated reaction IDs. Empty means all reactions.")
    p.add_argument("--fva_targets_file", default=None,
                   help="CSV whose ReactionID (or first) column lists reactions for FVA.")
    p.add_argument("--fva_processes", type=int, default=1)
    # Backward-compatible flag: --no_fva forces fva_mode=none.
    p.add_argument("--no_fva", action="store_true")
    args = p.parse_args()

    if not (0 < args.eflux_quantile <= 1):
        p.error("--eflux_quantile must be in (0, 1].")
    if args.eflux_floor < 0 or args.eflux_cap <= 0 or args.eflux_floor > args.eflux_cap:
        p.error("Require 0 <= eflux_floor <= eflux_cap.")
    if args.no_fva:
        args.fva_mode = "none"

    try:
        run_pipeline(
            rnaseq_file=args.rnaseq_file,
            met_data_file=args.met_data_file,
            model_file=args.model_file,
            results_dir=args.results_dir,
            aggregate=args.aggregate,
            infer_groups_from_filename=args.infer_groups_from_filename,
            baseline_code=args.baseline_code,
            explicit_groups=args.explicit_groups,
            mapping_file=args.mapping_file,
            column_regex=args.column_regex,
            eflux_quantile=args.eflux_quantile,
            eflux_floor=args.eflux_floor,
            eflux_cap=args.eflux_cap,
            diet_bounds_json=args.diet_bounds_json,
            diet_bounds_units=args.diet_bounds_units,
            mw_table=args.mw_table,
            gDW=args.gDW,
            hours_per_day=args.hours_per_day,
            pca_scale=args.pca_scale,
            write_replicates_long=args.write_replicates_long,
            test_type=args.test_type,
            objective_id=args.objective_id,
            objective_regex=args.objective_regex,
            objective_sense=args.objective_sense,
            transporter_strategy=args.transporter_strategy,
            transporter_regex=args.transporter_regex,
            transporter_subsystem_regex=args.transporter_subsystem_regex,
            transporter_compartments=args.transporter_compartments,
            edge_abs_diff_threshold=args.edge_abs_diff_threshold,
            solve_mode=args.solve_mode,
            pfba_fraction=args.pfba_fraction,
            solver=args.solver,
            deterministic_solver=not args.non_deterministic_solver,
            fva_mode=args.fva_mode,
            fva_fraction=args.fva_fraction,
            fva_targets=args.fva_targets,
            fva_targets_file=args.fva_targets_file,
            fva_processes=args.fva_processes,
        )
    except Exception as exc:
        print("ERROR:", exc)
        traceback.print_exc()
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
RQ4: Microbiome Community-Level Metabolic Modeling using MICOM
================================================================

Purpose:
--------
Builds microbial community models using MICOM framework to predict
community-level metabolite production (SCFAs, bile acids, amino acids)
that enter hepatic portal circulation.

Scientific Rationale:
--------------------
Individual species modeling (from map_agorav3.py) captures species-specific
responses, but real gut microbiomes exhibit metabolic crosstalk where:
1. Species compete for/share metabolites
2. Cross-feeding networks emerge (one species' waste = another's nutrient)
3. Community-level emergent properties arise

MICOM cooperative tradeoff modeling identifies community solutions that maintain a specified fraction of maximum community growth while distributing growth across taxa.

Integration with Research Questions:
-----------------------------------
- Links to RQ1: Compare microbiome metabolite production across diets
- Links to RQ2: Model strain-specific microbiome compositions
- Links to RQ3: Identify hepatic cell types most responsive to microbial metabolites

Author: Research Collaboration (Roland + AI-assisted review)
Date: 2026-07-08

BIOLOGICAL CORRECTIONS IN THIS VERSION:
---------------------------------------
1. Removed biologically unsafe DCA proxy mapping to straight-chain fatty acids.
2. Removed IMP as an imidazole-propionate proxy; IMP is usually inosine
   monophosphate in metabolic model namespaces unless model annotations prove
   otherwise.
3. Treats metatranscriptome-derived totals as activity-weighted abundance
   proxies, not measured biomass abundances.
4. Supports explicit species-to-AGORA mapping files and labels fuzzy matches as
   provisional.
5. Renormalizes retained community abundance after unmapped species are removed
   and writes retained-abundance audit tables.
6. Loads medium per condition instead of silently using the first nested medium.
7. Interprets MICOM cooperative_tradeoff fraction correctly as the minimum
   fraction of maximum community growth retained.
8. Uses pFBA by default for more parsimonious, auditable exchange fluxes.
9. Writes portal-metabolite audit tables documenting included, absent, and
   intentionally excluded metabolites.
10. Optionally validates portal hepatic exchange reactions against a hepatic
    model before writing the integration JSON.
"""

import os
import sys
import re
import json
import argparse
import warnings
import traceback
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# COBRApy and MICOM
try:
    import cobra
    from cobra.io import load_json_model, load_matlab_model
    COBRA_AVAILABLE = True
except ImportError:
    cobra = None
    load_json_model = None
    load_matlab_model = None
    COBRA_AVAILABLE = False
    print("[WARNING] COBRApy not installed. Install with: pip install cobra")

# MICOM imports (with graceful fallback)
try:
    from micom import Community, load_pickle
    from micom.workflows import grow
    from micom.workflows.core import workflow
    MICOM_AVAILABLE = True
except ImportError:
    Community = object  # keeps type annotations import-safe when MICOM is absent
    print("[WARNING] MICOM not installed. Install with: pip install micom")
    print("[WARNING] Community modeling requires MICOM; the script will exit before modeling if MICOM is unavailable.")
    MICOM_AVAILABLE = False


################################################################################
# CONFIGURATION AND CONSTANTS
################################################################################

# =============================================================================
# PORTAL METABOLITES: Evidence-Based Gut-Liver Axis Metabolites
# =============================================================================
# This dictionary defines metabolites that:
# 1. Are produced by gut microbiota (documented in literature)
# 2. Enter hepatic portal circulation (physiological concentrations)
# 3. Have known metabolic effects on liver (mechanistic evidence)
# 4. Are available in iMM1415 mouse metabolic model (practical constraint)
#
# Selection rule:
# - Include only metabolite IDs with plausible gut-liver biology and model IDs
#   that are not known namespace collisions.
# - Do not use fatty-acid identifiers as bile-acid proxies.
# - If --hepatic_model is supplied, metabolites are included in the portal JSON
#   only when the corresponding hepatic exchange reaction exists.
#
# Literature support in metadata is a guide for interpretation, not evidence that
# the model predicted a given exchange. Interpret only metabolites with active,
# exact-matched exchange fluxes and preferably measured metabolite validation.
# =============================================================================

PORTAL_METABOLITES = {
    # =========================================================================
    # TIER 1: SHORT-CHAIN FATTY ACIDS (SCFAs) - HIGH PRIORITY
    # =========================================================================
    # Major microbial fermentation products. In portal circulation, acetate and
    # propionate are typically more directly liver-exposed than butyrate, which
    # is often consumed by colonocytes; all three are retained as mechanistically
    # relevant candidates when the model predicts exchange.
    'ac': {
        'name': 'acetate',
        'hepatic_rxn': 'EX_ac_e',
        'importance': 'high',
        'category': 'SCFA',
        'portal_conc_range': '0.1-0.5 mM',
        'mechanisms': 'Lipogenesis substrate, acetyl-CoA supply, HDAC-linked signaling',
        'evidence': 'Canfora 2015; Nogal 2021; Pant 2023',
        'interpretation_note': 'Only report as active if exact exchange flux is detected; never infer from nac/lactate substring matches.'
    },
    'ppa': {
        'name': 'propionate',
        'hepatic_rxn': 'EX_ppa_e',
        'importance': 'high',
        'category': 'SCFA',
        'portal_conc_range': '0.01-0.1 mM',
        'mechanisms': 'Gluconeogenic/anaplerotic substrate; cholesterol and GPCR-linked signaling',
        'evidence': 'Canfora 2015; Tirosh 2019; Nogal 2021'
    },
    'but': {
        'name': 'butyrate',
        'hepatic_rxn': 'EX_but_e',
        'importance': 'high',
        'category': 'SCFA',
        'portal_conc_range': '0.001-0.01 mM',
        'mechanisms': 'HDAC inhibition; mitochondrial beta-oxidation; anti-inflammatory signaling',
        'evidence': 'Canfora 2015; Pant 2023',
        'interpretation_note': 'Portal liver exposure is usually lower than acetate/propionate because colonocytes consume substantial butyrate.'
    },

    # =========================================================================
    # BILE ACIDS - INCLUDE ONLY TRUE BILE ACID IDS, NO FATTY-ACID PROXIES
    # =========================================================================
    # Do not substitute hdca/ocdca/ptdca for deoxycholate. Those identifiers are
    # fatty acids in common constraint-based-model namespaces and are not valid
    # C24 steroid bile-acid proxies.
    'dca': {
        'name': 'deoxycholate / deoxycholic acid',
        'hepatic_rxn': 'EX_dca_e',
        'importance': 'medium',
        'category': 'bile_acid_secondary',
        'portal_conc_range': 'low µM; diet and microbiome dependent',
        'mechanisms': 'Secondary bile acid; FXR/TGR5 and hepatobiliary signaling',
        'evidence': 'Wahlstrom 2016; Chiang and Ferrell 2020; Ridlon 2014',
        'interpretation_note': 'Retained only when the exact exchange is present and, if --hepatic_model is supplied, the hepatic exchange exists.'
    },
    'tdchola': {
        'name': 'taurodeoxycholate',
        'hepatic_rxn': 'EX_tdchola_e',
        'importance': 'medium',
        'category': 'bile_acid_conjugated',
        'portal_conc_range': '1-10 µM',
        'mechanisms': 'Taurine-conjugated bile acid; hepatobiliary and receptor-linked signaling',
        'evidence': 'Wahlstrom 2016; Chiang and Ferrell 2020'
    },

    # =========================================================================
    # NITROGEN METABOLISM
    # =========================================================================
    'nh4': {
        'name': 'ammonia / ammonium',
        'hepatic_rxn': 'EX_nh4_e',
        'importance': 'high',
        'category': 'nitrogen',
        'portal_conc_range': '0.05-0.2 mM',
        'mechanisms': 'Urea-cycle substrate; nitrogen load; glutamine/glutamate metabolism',
        'evidence': 'Haeussinger 1990; Shawcross 2023; Zhu 2023',
        'interpretation_note': 'Negative community exchange means microbial uptake from the medium/lumen, not microbiome production for hepatic delivery.'
    },

    # =========================================================================
    # ORGANIC ACIDS AND TCA-RELATED METABOLITES
    # =========================================================================
    'succ': {
        'name': 'succinate',
        'hepatic_rxn': 'EX_succ_e',
        'importance': 'high',
        'category': 'organic_acid',
        'portal_conc_range': '1-10 µM',
        'mechanisms': 'TCA anaplerosis; SUCNR1/GPR91 signaling; inflammatory/metabolic signaling',
        'evidence': 'Fernandez-Veledo 2019; Mills 2016; Wei 2023'
    },
    'lac__L': {
        'name': 'L-lactate',
        'hepatic_rxn': 'EX_lac__L_e',
        'importance': 'medium',
        'category': 'organic_acid',
        'portal_conc_range': '0.5-2 mM',
        'mechanisms': 'Gluconeogenesis substrate; redox-linked carbon shuttle',
        'evidence': 'Brooks 2018; Gladden 2004'
    },
    'lac__D': {
        'name': 'D-lactate',
        'hepatic_rxn': 'EX_lac__D_e',
        'importance': 'low',
        'category': 'organic_acid',
        'portal_conc_range': '0.01-0.1 mM',
        'mechanisms': 'Dysbiosis marker; slower mammalian clearance than L-lactate',
        'evidence': 'Ewaschuk 2005; Chu 2020'
    },
    'pyr': {
        'name': 'pyruvate',
        'hepatic_rxn': 'EX_pyr_e',
        'importance': 'low',
        'category': 'organic_acid',
        'portal_conc_range': '0.05-0.15 mM',
        'mechanisms': 'Gluconeogenesis and acetyl-CoA precursor; community cross-feeding intermediate',
        'evidence': 'Jeoung 2015; Magnusdottir 2017'
    },

    # =========================================================================
    # ETHANOL AND VITAMINS/COFACTORS
    # =========================================================================
    'etoh': {
        'name': 'ethanol',
        'hepatic_rxn': 'EX_etoh_e',
        'importance': 'medium',
        'category': 'alcohol',
        'portal_conc_range': '<0.1 mM; can be higher in dysbiosis',
        'mechanisms': 'ADH/ALDH metabolism; acetaldehyde, redox stress, lipogenesis',
        'evidence': 'Meijnikman 2022; Yuan 2019; Bishehsari 2017'
    },
    'fol': {
        'name': 'folate',
        'hepatic_rxn': 'EX_fol_e',
        'importance': 'low',
        'category': 'vitamin',
        'portal_conc_range': '0.01-0.1 µM',
        'mechanisms': 'One-carbon metabolism; methionine cycle',
        'evidence': 'Rossi 2011; Tarracchini 2024'
    },
    'ribflv': {
        'name': 'riboflavin / vitamin B2',
        'hepatic_rxn': 'EX_ribflv_e',
        'importance': 'low',
        'category': 'vitamin',
        'portal_conc_range': '0.01-0.1 µM',
        'mechanisms': 'FAD/FMN precursor; redox metabolism',
        'evidence': 'Thakur 2015; Tarracchini 2024'
    },
    'thm': {
        'name': 'thiamine / vitamin B1',
        'hepatic_rxn': 'EX_thm_e',
        'importance': 'low',
        'category': 'vitamin',
        'portal_conc_range': '<0.1 µM',
        'mechanisms': 'TPP precursor; pyruvate dehydrogenase and TCA-linked metabolism',
        'evidence': 'Begley 2007; Magnusdottir 2015'
    },
    'nac': {
        'name': 'niacin / nicotinate / vitamin B3',
        'hepatic_rxn': 'EX_nac_e',
        'importance': 'low',
        'category': 'vitamin',
        'portal_conc_range': '<0.1 µM',
        'mechanisms': 'NAD+/NADP+ precursor; redox homeostasis',
        'evidence': 'Magnusdottir 2015; Grozio 2019'
    },
}

# Intentionally excluded identifiers that were biologically unsafe in the prior
# version. They may still appear in unfiltered exploratory exchange tables, but
# they are not used as curated portal constraints.
EXCLUDED_PORTAL_METABOLITES = {
    'imp': {
        'reason': 'Model ID imp usually denotes inosine monophosphate, not imidazole propionate. Do not use as a histidine-derived insulin-resistance metabolite without explicit model annotation validation.',
        'previous_mislabel': 'imidazole propionate',
    },
    'hdca': {
        'reason': 'Commonly denotes hexadecanoate/palmitate, a straight-chain fatty acid, not a deoxycholate bile-acid proxy.',
        'previous_mislabel': 'DCA analog',
    },
    'ocdca': {
        'reason': 'Commonly denotes octadecanoate/stearate, a straight-chain fatty acid, not a deoxycholate bile-acid proxy.',
        'previous_mislabel': 'DCA analog',
    },
    'ptdca': {
        'reason': 'Straight-chain fatty acid identifier; not a validated secondary bile-acid proxy.',
        'previous_mislabel': 'DCA analog',
    },
}

# Metabolite importance tiers for reporting only. These are not statistical
# significance tiers.
PRIORITY_TIERS = {
    'HIGH': ['ac', 'ppa', 'but', 'nh4', 'succ'],
    'MEDIUM': ['lac__L', 'etoh', 'tdchola', 'dca'],
    'LOW': ['lac__D', 'pyr', 'fol', 'ribflv', 'thm', 'nac'],
}


################################################################################
# DATA LOADING AND PREPROCESSING
################################################################################

def load_metatranscriptome_data(
    filepath: str,
    expression_cols: List[str] = ['ND_SCD', 'DD_HFD']
) -> pd.DataFrame:
    """
    Load metatranscriptomic data from CSV file.
    
    Scientific Context:
    ------------------
    Metatranscriptomics measures community transcriptional activity. In this
    script it is used as an activity-weighted proxy for community composition,
    not as a direct measurement of biomass abundance. True expression-constrained
    AGORA modeling would require mapping microbial genes to each species model's
    GPR rules and constraining reaction bounds; this script does not claim that.
    
    Parameters
    ----------
    filepath : str
        Path to metatranscriptome CSV file
    expression_cols : list of str
        Column names containing expression values for each condition
        
    Returns
    -------
    pd.DataFrame
        Metatranscriptome data with species and gene expression
        
    Engineering Benefit:
    -------------------
    Standardized data loading ensures reproducibility across different
    metatranscriptome datasets. Explicit column specification allows
    easy adaptation to different experimental designs.
    """
    print(f"[INFO] Loading metatranscriptome data from: {filepath}")
    
    # Load data
    df = pd.read_csv(filepath)
    print(f"[INFO] Loaded {len(df)} gene expression records")
    print(f"[INFO] Columns: {df.columns.tolist()}")
    
    # Validate expression columns exist
    missing_cols = [col for col in expression_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Expression columns not found: {missing_cols}")
    
    # Infer species from Genome_name or similar column
    if 'Genome_name' in df.columns:
        species_col = 'Genome_name'
    elif 'Genome' in df.columns:
        species_col = 'Genome'
    elif 'species' in df.columns:
        species_col = 'species'
    else:
        raise ValueError("No species column found (expected: Genome_name, Genome, or species)")
    
    # Standardize species column name
    df = df.rename(columns={species_col: 'species'})
    
    # Count unique species
    n_species = df['species'].nunique()
    print(f"[INFO] Found {n_species} unique species")
    
    # Show top species by total expression
    species_totals = df.groupby('species')[expression_cols].sum().sum(axis=1).sort_values(ascending=False)
    print(f"[INFO] Top 10 species by total expression:")
    print(species_totals.head(10))
    
    return df


def infer_species_abundances(
    metatranscriptome_df: pd.DataFrame,
    expression_cols: List[str],
    method: str = 'total_expression'
) -> pd.DataFrame:
    """
    Infer relative species abundances from metatranscriptomic data.
    
    Scientific Rationale:
    --------------------
    Total metatranscriptomic signal per species is an activity-weighted
    composition proxy. It should not be described as measured biomass abundance.
    This proxy can reflect both cell abundance and per-cell transcriptional
    activity, so downstream claims should be phrased as model predictions from
    metatranscriptome-derived community activity.

    If true metagenomic/16S abundance is available, use that for community
    composition and reserve metatranscriptomics for pathway/reaction activity.
    
    Parameters
    ----------
    metatranscriptome_df : pd.DataFrame
        Metatranscriptome data with species and expression columns
    expression_cols : list of str
        Column names to use for abundance estimation
    method : str
        Method for abundance inference:
        - 'total_expression': Sum all gene expression per species
        - 'housekeeping': Use constitutively expressed genes (if available)
        
    Returns
    -------
    pd.DataFrame
        Activity-weighted relative composition proxies for each condition
        (sum to 1.0 before AGORA mapping). Columns: species,
        {condition}_abundance
        
    Engineering Benefit:
    -------------------
    Multiple abundance estimation methods allow sensitivity analysis.
    Normalized abundances (sum=1.0) ensure consistent MICOM input regardless
    of sequencing depth differences between samples.
    """
    print(f"[INFO] Inferring metatranscriptome-derived activity-weighted composition using method: {method}")
    
    if method == 'total_expression':
        # Sum all gene expression per species per condition
        abundance_df = metatranscriptome_df.groupby('species')[expression_cols].sum()
        
        # Normalize to relative activity proportions (sum = 1.0 per condition)
        for col in expression_cols:
            total = abundance_df[col].sum()
            if total <= 0 or not np.isfinite(total):
                raise ValueError(f"Cannot infer abundance proxy for {col}: total expression is {total}")
            abundance_df[f'{col}_abundance'] = abundance_df[col] / total
        
        # Keep only abundance columns
        abundance_df = abundance_df[[f'{col}_abundance' for col in expression_cols]].reset_index()
        
    elif method == 'housekeeping':
        # TODO: Implement housekeeping gene method
        # Requires pre-defined list of universally expressed genes
        raise NotImplementedError("Housekeeping gene method not yet implemented")
    
    else:
        raise ValueError(f"Unknown abundance method: {method}")
    
    print(f"[INFO] Inferred activity-weighted composition for {len(abundance_df)} species")
    
    # Show top 5 most abundant species per condition
    for col in expression_cols:
        print(f"\n{col}_abundance:")
        print(abundance_df.nlargest(5, f'{col}_abundance')[['species', f'{col}_abundance']])
    
    return abundance_df


def _normalize_name_for_matching(name: str) -> str:
    """Normalize taxon/model names for deterministic matching."""
    return re.sub(r'[^a-z0-9]+', '_', str(name).lower()).strip('_')


def map_species_to_agora_models(
    species_list: List[str],
    agora_dir: str,
    mapping_file: Optional[str] = None,
    allow_fuzzy: bool = False,
    strict_mapping_file: bool = False,
    min_prefix_words: int = 2,
    results_dir: Optional[str] = None,
) -> Tuple[Dict[str, str], pd.DataFrame]:
    """
    Map species names from the metatranscriptome to AGORA model files.

    Scientific correction:
    ----------------------
    Strain-level AGORA models can differ biologically. Therefore, exact or
    user-curated mappings are preferred. Fuzzy prefix matching is disabled by
    default because automatically selecting the first matching strain can create
    silent biological errors. If --allow_fuzzy_model_matching is used, all fuzzy
    matches are labeled as provisional in the mapping report.

    Optional mapping_file format:
    -----------------------------
    CSV/TSV with at least one species column and one model column. Accepted
    species columns: species, Genome, Genome_name, taxon. Accepted model columns:
    model_file, file, agora_model, model_path.

    Returns
    -------
    species_to_model : dict
        Mapping {species_name: absolute_model_filepath}
    mapping_report : pd.DataFrame
        One row per input species with match_status and confidence fields.
    """
    print(f"[INFO] Mapping species to AGORA models in: {agora_dir}")
    agora_path = Path(agora_dir)
    model_files = list(agora_path.glob('*.mat')) + list(agora_path.glob('*.json'))
    if not model_files:
        raise FileNotFoundError(f"No .mat or .json AGORA models found in {agora_dir}")

    print(f"[INFO] Found {len(model_files)} AGORA model files")
    model_by_stem = {f.stem: str(f.resolve()) for f in model_files}
    model_by_norm = {_normalize_name_for_matching(f.stem): str(f.resolve()) for f in model_files}

    explicit_map = {}
    if mapping_file:
        print(f"[INFO] Loading explicit species-model mapping: {mapping_file}")
        sep = '\t' if str(mapping_file).lower().endswith(('.tsv', '.txt')) else ','
        map_df = pd.read_csv(mapping_file, sep=sep)
        species_cols = [c for c in ['species', 'Genome', 'Genome_name', 'taxon'] if c in map_df.columns]
        model_cols = [c for c in ['model_file', 'file', 'agora_model', 'model_path'] if c in map_df.columns]
        if not species_cols or not model_cols:
            raise ValueError(
                "Mapping file must contain a species column "
                "(species/Genome/Genome_name/taxon) and a model column "
                "(model_file/file/agora_model/model_path)."
            )
        sp_col, model_col = species_cols[0], model_cols[0]
        for _, row in map_df.iterrows():
            sp = str(row[sp_col])
            sp_norm = _normalize_name_for_matching(sp)
            model_value = str(row[model_col])
            candidate = Path(model_value)
            if not candidate.is_absolute():
                # Accept either a basename/stem or a path relative to agora_dir.
                by_stem = model_by_stem.get(candidate.stem)
                by_name = model_by_norm.get(_normalize_name_for_matching(candidate.stem))
                candidate = Path(by_stem or by_name or (agora_path / model_value))
            if candidate.exists():
                if sp_norm in explicit_map and explicit_map[sp_norm] != str(candidate.resolve()):
                    raise ValueError(
                        f"Duplicate normalized explicit mapping key {sp_norm!r} points to different models."
                    )
                explicit_map[sp_norm] = str(candidate.resolve())
            else:
                print(f"[WARNING] Explicit mapping for {sp} points to missing model: {model_value}")

    species_to_model = {}
    report_rows = []

    for species in species_list:
        species_str = str(species)
        norm_species = _normalize_name_for_matching(species_str)
        status = 'unmatched'
        confidence = 'none'
        model_path = None
        match_basis = ''

        if norm_species in explicit_map:
            model_path = explicit_map[norm_species]
            status = 'matched'
            confidence = 'curated_explicit_normalized'
            match_basis = 'mapping_file_normalized'
        elif mapping_file and strict_mapping_file:
            status = 'unmatched'
            confidence = 'excluded_by_strict_mapping_file'
            match_basis = 'strict_mapping_file'
        elif norm_species in model_by_norm:
            model_path = model_by_norm[norm_species]
            status = 'matched'
            confidence = 'exact_normalized'
            match_basis = Path(model_path).stem
        elif allow_fuzzy:
            words = [w for w in norm_species.split('_') if w]
            for n_words in range(len(words), max(min_prefix_words, 1) - 1, -1):
                prefix = '_'.join(words[:n_words])
                matches = sorted(k for k in model_by_norm if k.startswith(prefix))
                if len(matches) == 1:
                    model_path = model_by_norm[matches[0]]
                    status = 'matched_provisional'
                    confidence = f'fuzzy_unique_prefix_{n_words}_words'
                    match_basis = matches[0]
                    break
                elif len(matches) > 1:
                    status = 'ambiguous_fuzzy'
                    confidence = f'{len(matches)}_matches_for_prefix_{prefix}'
                    match_basis = ';'.join(matches[:5])
                    break

        if model_path:
            species_to_model[species_str] = model_path

        report_rows.append({
            'species': species_str,
            'model_file': os.path.basename(model_path) if model_path else '',
            'model_path': model_path or '',
            'match_status': status,
            'match_confidence': confidence,
            'match_basis': match_basis,
        })

    mapping_report = pd.DataFrame(report_rows)
    n_matched = mapping_report['match_status'].isin(['matched', 'matched_provisional']).sum()
    n_provisional = (mapping_report['match_status'] == 'matched_provisional').sum()
    n_ambiguous = (mapping_report['match_status'] == 'ambiguous_fuzzy').sum()
    n_unmatched = (mapping_report['match_status'] == 'unmatched').sum()

    print(f"[INFO] Successfully mapped {n_matched}/{len(species_list)} species to AGORA models")
    if n_provisional:
        print(f"[WARNING] {n_provisional} mappings are fuzzy/provisional. Review mapping report before publication.")
    if n_ambiguous:
        print(f"[WARNING] {n_ambiguous} species had ambiguous fuzzy matches and were not mapped.")
    if n_unmatched:
        print(f"[WARNING] {n_unmatched} species could not be matched.")

    if results_dir:
        path = os.path.join(results_dir, 'species_to_model_mapping_audit.csv')
        mapping_report.to_csv(path, index=False)
        print(f"[SAVED] Species-model mapping audit: {path}")

    return species_to_model, mapping_report


# -----------------------------------------------------------------------------
# Medium loading/selection helpers
# -----------------------------------------------------------------------------

def load_medium_config(medium_json: Optional[str]) -> Dict:
    """Load a medium JSON or return a conservative default if none is supplied."""
    if medium_json is None:
        print("[INFO] No medium JSON supplied; using default generic gut medium.")
        return {
            'EX_glc__D_e': 10.0,
            'EX_fru_e': 5.0,
            'EX_gal_e': 2.0,
            'EX_mal__L_e': 5.0,
            'EX_pyr_e': 1.0,
            'EX_ac_e': 0.1,
            'EX_o2_e': 20.0,
            'EX_pi_e': 10.0,
            'EX_so4_e': 10.0,
            'EX_nh4_e': 10.0,
        }

    with open(medium_json, 'r') as f:
        medium_dict = json.load(f)
    if not isinstance(medium_dict, dict) or len(medium_dict) == 0:
        raise ValueError(f"Medium JSON must contain a non-empty dictionary, got {type(medium_dict)}")
    return medium_dict


def _is_bounds_or_number(value) -> bool:
    """Return True for a single medium entry value."""
    return isinstance(value, (int, float)) or (
        isinstance(value, list) and len(value) >= 1 and all(isinstance(v, (int, float)) for v in value[:2])
    )


def _condition_medium_keys(condition: str) -> List[str]:
    """Candidate diet/media keys for a condition label such as ND_SCD or DD_HFD."""
    condition = str(condition)
    upper = condition.upper()
    keys = [condition]
    aliases = {
        'ND_SCD': ['ND_SCD', 'SCD', 'ND', 'NORMAL', 'CONTROL'],
        'SCD': ['SCD', 'ND_SCD', 'CONTROL', 'NORMAL'],
        'DD_HFD': ['DD_HFD', 'HFD', 'DD', 'HIGH_FAT'],
        'HFD': ['HFD', 'DD_HFD', 'HIGH_FAT'],
        'KD': ['KD', 'KETOGENIC'],
        'WD': ['WD', 'WESTERN'],
    }
    for pattern, candidates in aliases.items():
        if pattern in upper:
            keys.extend(candidates)
    # Also try final token after underscore.
    if '_' in condition:
        keys.append(condition.split('_')[-1])
    # Preserve order but deduplicate case-insensitively later.
    deduped = []
    seen = set()
    for k in keys:
        if k and k.upper() not in seen:
            deduped.append(k)
            seen.add(k.upper())
    return deduped


def select_medium_for_condition(medium_config: Dict, condition: str) -> Tuple[Dict, str]:
    """
    Select the correct medium block for a condition.

    Critical correction: if the medium JSON is nested by diet/condition, this
    function matches the current condition instead of silently taking the first
    block. A flat dictionary is treated as shared across all conditions.
    """
    first_value = next(iter(medium_config.values()))
    if _is_bounds_or_number(first_value):
        return medium_config, 'flat_shared_medium'

    if not isinstance(first_value, dict):
        raise ValueError(f"Unexpected medium JSON value type: {type(first_value)}")

    key_lookup = {str(k).upper(): k for k in medium_config.keys()}
    for candidate in _condition_medium_keys(condition):
        if candidate.upper() in key_lookup:
            selected_key = key_lookup[candidate.upper()]
            selected = medium_config[selected_key]
            print(f"[INFO] Using medium block '{selected_key}' for condition '{condition}'")
            return selected, str(selected_key)

    available = ', '.join(map(str, medium_config.keys()))
    raise ValueError(
        f"No medium block matched condition '{condition}'. Available medium keys: {available}. "
        "Provide matching keys or use a flat shared medium if this is intentional."
    )


def medium_bounds_to_uptake(medium_raw: Dict[str, object]) -> Dict[str, float]:
    """
    Convert medium entries to MICOM uptake capacities.

    Accepts either:
    - {'EX_glc__D_e': 10.0}
    - {'EX_glc__D_e': [-10.0, 1000.0]}

    Returns positive uptake rates expected by MICOM's medium property.
    """
    medium = {}
    for rxn_id, value in medium_raw.items():
        if isinstance(value, list):
            if len(value) < 1 or not isinstance(value[0], (int, float)):
                raise ValueError(f"Bounds list for {rxn_id} must start with numeric lower bound: {value}")
            lower_bound = value[0]
            uptake_rate = abs(float(lower_bound))
        elif isinstance(value, (int, float)):
            uptake_rate = abs(float(value))
        else:
            raise ValueError(f"Unexpected medium value for {rxn_id}: {type(value)}")

        # Skip blocked medium entries; keeping zeros in MICOM medium is harmless
        # but noisy and may hide accidental missing diet constraints.
        if uptake_rate > 0:
            medium[str(rxn_id)] = uptake_rate
    return medium


def get_condition_medium(
    medium_config: Dict,
    condition: str,
    target_compartment: str = 'e',
) -> Tuple[Dict[str, float], str]:
    """Select, convert, and compartment-standardize medium for one condition."""
    medium_raw, source_key = select_medium_for_condition(medium_config, condition)
    medium = medium_bounds_to_uptake(medium_raw)
    medium = standardize_medium_compartments(medium, target_compartment=target_compartment)
    print(f"[INFO] Medium for {condition}: {len(medium)} exchange reactions from {source_key}")
    print(f"[INFO] Sample medium reactions: {list(medium.keys())[:5]}")
    return medium, source_key


def standardize_medium_compartments(
    medium: Dict[str, float],
    target_compartment: str = 'e'
) -> Dict[str, float]:
    """
    Standardize exchange reaction compartment suffixes for MICOM.
    
    Scientific Rationale:
    --------------------
    AGORA models use different compartment conventions:
    - Some use '_m' (medium/external)
    - Others use '_e' (extracellular)
    - MICOM community models expect consistent naming
    
    This function ensures all exchange reactions use the same compartment suffix.
    
    Parameters
    ----------
    medium : dict
        Medium composition {reaction_id: flux_value}
    target_compartment : str
        Target compartment suffix ('e' or 'm')
        
    Returns
    -------
    dict
        Standardized medium dictionary
        
    Engineering Benefit:
    -------------------
    Prevents subtle bugs from compartment naming inconsistencies.
    Makes medium definitions portable across different AGORA versions.
    """
    standardized = {}
    
    for rxn_id, flux_val in medium.items():
        # Replace any compartment suffix with target
        # E.g., EX_glc__D_m → EX_glc__D_e
        if rxn_id.endswith('_m'):
            new_id = rxn_id[:-2] + f'_{target_compartment}'
        elif rxn_id.endswith('_e'):
            new_id = rxn_id[:-2] + f'_{target_compartment}'
        else:
            # No compartment suffix, add one
            new_id = f'{rxn_id}_{target_compartment}'
        
        standardized[new_id] = flux_val
    
    return standardized


def validate_medium(medium: Dict[str, float]) -> bool:
    """
    Validate medium composition before using with MICOM.
    
    Checks:
    ------
    1. Medium is a dictionary (not an integer or other type)
    2. All values are numeric
    3. Essential nutrients present (glucose, oxygen, etc.)
    
    Parameters
    ----------
    medium : dict
        Medium composition to validate
        
    Returns
    -------
    bool
        True if medium is valid
        
    Raises
    ------
    ValueError
        If medium fails validation
    """
    # Check type
    if not isinstance(medium, dict):
        raise ValueError(f"Medium must be dict, got {type(medium)}")
    
    if len(medium) == 0:
        raise ValueError("Medium is empty")
    
    # Check all values are numeric
    for rxn_id, flux_val in medium.items():
        if not isinstance(flux_val, (int, float)):
            raise ValueError(f"Medium flux must be numeric, got {type(flux_val)} for {rxn_id}")
    
    # Check for essential nutrients
    essential = ['glc', 'o2']  # Glucose and oxygen
    found_essential = {nutrient: False for nutrient in essential}
    
    for rxn_id in medium.keys():
        for nutrient in essential:
            if nutrient in rxn_id.lower():
                found_essential[nutrient] = True
    
    missing = [n for n, found in found_essential.items() if not found]
    if missing:
        print(f"[WARNING] Essential nutrients may be missing: {missing}")
        print(f"[WARNING] Available reactions: {list(medium.keys())[:10]}")
    
    print(f"[INFO] Medium validation passed: {len(medium)} exchange reactions")
    return True


def build_micom_community(
    abundance_df: pd.DataFrame,
    species_to_model: Dict[str, str],
    condition: str,
    medium: Dict[str, float],
    solver: str = 'gurobi',
    results_dir: Optional[str] = None,
    min_retained_abundance: float = 0.80,
) -> Optional[Community]:
    """
    Build a MICOM community model for one condition.

    Biological correction:
    ----------------------
    Abundances are renormalized after unmapped species are removed. The retained
    fraction is saved and printed because low retained abundance means the
    community no longer represents the measured/transcriptome-derived community.
    """
    abundance_col = f'{condition}_abundance'
    if abundance_col not in abundance_df.columns:
        print(f"[ERROR] Abundance column not found: {abundance_col}")
        return None

    species_with_abundance = abundance_df[abundance_df[abundance_col] > 0].copy()
    species_with_abundance['has_agora_model'] = species_with_abundance['species'].isin(species_to_model)

    retained_before_norm = float(
        species_with_abundance.loc[species_with_abundance['has_agora_model'], abundance_col].sum()
    )
    dropped_abundance = float(
        species_with_abundance.loc[~species_with_abundance['has_agora_model'], abundance_col].sum()
    )

    if retained_before_norm <= 0:
        print("[ERROR] No positive-abundance species with AGORA models available for community building!")
        return None

    if retained_before_norm < min_retained_abundance:
        print(
            f"[WARNING] Only {retained_before_norm:.1%} of activity-weighted community "
            f"retained after AGORA mapping for {condition}. Interpret results cautiously."
        )
    else:
        print(f"[INFO] Retained {retained_before_norm:.1%} of activity-weighted community after AGORA mapping")

    taxonomy = species_with_abundance[species_with_abundance['has_agora_model']].copy()
    taxonomy['original_abundance_proxy'] = taxonomy[abundance_col]
    taxonomy['abundance'] = taxonomy[abundance_col] / retained_before_norm
    taxonomy['id'] = taxonomy['species'].astype(str)
    taxonomy['file'] = taxonomy['species'].map(species_to_model)
    taxonomy = taxonomy[['id', 'species', 'abundance', 'original_abundance_proxy', 'file']]

    print(f"\n[INFO] Building MICOM community for condition: {condition}")
    print(f"[INFO] Community contains {len(taxonomy)} modeled species")
    print(f"[INFO] Renormalized abundance sum: {taxonomy['abundance'].sum():.6f}")

    if results_dir:
        os.makedirs(results_dir, exist_ok=True)
        taxonomy_path = os.path.join(results_dir, f'taxonomy_{condition}.csv')
        taxonomy.to_csv(taxonomy_path, index=False)
        coverage_path = os.path.join(results_dir, f'mapping_coverage_{condition}.json')
        with open(coverage_path, 'w') as f:
            json.dump({
                'condition': condition,
                'retained_activity_weighted_abundance_before_renormalization': retained_before_norm,
                'dropped_activity_weighted_abundance': dropped_abundance,
                'n_species_total_positive': int(len(species_with_abundance)),
                'n_species_modeled': int(len(taxonomy)),
                'n_species_unmapped_positive': int((~species_with_abundance['has_agora_model']).sum()),
                'abundance_column': abundance_col,
                'abundance_interpretation': 'metatranscriptome-derived activity-weighted proxy, renormalized after AGORA mapping',
            }, f, indent=2)
        print(f"[SAVED] MICOM taxonomy table: {taxonomy_path}")
        print(f"[SAVED] Mapping coverage summary: {coverage_path}")

    try:
        validate_medium(medium)
    except ValueError as e:
        print(f"[ERROR] Medium validation failed: {e}")
        return None

    try:
        print("[INFO] Loading AGORA models and building community...")
        community = Community(
            taxonomy=taxonomy[['id', 'abundance', 'file']],
            model_db=None,
            solver=solver,
            progress=True
        )

        print("[INFO] Setting community medium...")
        community_exchanges = [r.id for r in community.exchanges]
        community_exchange_set = set(community_exchanges)
        print(f"[INFO] Community has {len(community_exchanges)} exchange reactions")

        matched_medium = {}
        unmatched = []
        for rxn_id, flux_val in medium.items():
            candidates = [rxn_id]
            if rxn_id.startswith('EX_'):
                base_id = re.sub(r'(_[a-z]\d?|\([a-z]\d?\)|\[[a-z]\d?\])$', '', rxn_id)
                candidates.extend([base_id + s for s in ['_e', '_m', '[e]', '[m]', '(e)', '(m)']])
            for test_id in candidates:
                if test_id in community_exchange_set:
                    matched_medium[test_id] = float(flux_val)
                    break
            else:
                unmatched.append(rxn_id)

        if unmatched:
            print(f"[WARNING] {len(unmatched)}/{len(medium)} medium reactions could not be matched to community")
            if results_dir:
                pd.DataFrame({'unmatched_medium_reaction': unmatched}).to_csv(
                    os.path.join(results_dir, f'unmatched_medium_{condition}.csv'), index=False
                )

        if not matched_medium:
            print("[ERROR] No medium reactions matched the community exchanges. Refusing to use hidden default medium.")
            return None

        community.medium = matched_medium
        print(f"[INFO] Matched and set {len(matched_medium)}/{len(medium)} medium reactions")
        print(f"[INFO] Community objective: {community.objective.expression}")
        return community

    except Exception as e:
        print(f"[ERROR] Failed to build community: {e}")
        traceback.print_exc()
        return None


def _parse_exchange_metabolite_id(rxn_id: str) -> str:
    """Parse metabolite ID from exchange reaction ID without substring matching."""
    met_id = str(rxn_id)
    if met_id.startswith('EX_'):
        met_id = met_id[3:]
    met_id = re.sub(r'(_[a-z]\d?|\([a-z]\d?\)|\[[a-z]\d?\])$', '', met_id)
    return met_id


def _aggregate_exchange_fluxes(solution, community: Community, threshold: float = 1e-9) -> pd.DataFrame:
    """
    Aggregate community exchange fluxes robustly across MICOM/Cobrapy variants.

    Handles:
    - solution.exchanges DataFrame if available
    - solution.fluxes as taxa x reactions DataFrame
    - solution.fluxes as reaction Series/dict
    """
    exchange_rxn_ids = [r.id for r in community.exchanges]
    exchange_set = set(exchange_rxn_ids)
    exchange_data = []
    all_fluxes = []

    if hasattr(solution, 'exchanges') and isinstance(solution.exchanges, pd.DataFrame):
        ex = solution.exchanges.copy()
        # MICOM versions differ in column names. Normalize to reaction/flux.
        if 'reaction' not in ex.columns:
            if ex.index.name and 'reaction' in ex.index.name.lower():
                ex = ex.reset_index()
            elif ex.index.dtype == object:
                ex = ex.reset_index().rename(columns={'index': 'reaction'})
        flux_col = 'flux' if 'flux' in ex.columns else None
        if flux_col is None:
            numeric_cols = ex.select_dtypes(include=[np.number]).columns.tolist()
            if numeric_cols:
                flux_col = numeric_cols[0]
        if 'reaction' in ex.columns and flux_col:
            for _, row in ex.iterrows():
                rxn_id = str(row['reaction'])
                flux_val = float(row[flux_col])
                all_fluxes.append((rxn_id, abs(flux_val)))
                if abs(flux_val) > threshold:
                    exchange_data.append({
                        'reaction': rxn_id,
                        'flux': flux_val,
                        'metabolite': _parse_exchange_metabolite_id(rxn_id),
                        'direction': 'export' if flux_val > 0 else 'import'
                    })

    if not exchange_data and hasattr(solution, 'fluxes'):
        fluxes = solution.fluxes
        if isinstance(fluxes, pd.DataFrame):
            # Expected MICOM format: rows taxa, columns reactions.
            if len(set(fluxes.columns) & exchange_set) > 0:
                matched = [rxn_id for rxn_id in exchange_rxn_ids if rxn_id in fluxes.columns]
                print(f"[INFO] Found {len(matched)}/{len(exchange_rxn_ids)} exchange reactions in solution columns")
                for rxn_id in matched:
                    flux_val = float(pd.to_numeric(fluxes[rxn_id], errors='coerce').fillna(0).sum())
                    all_fluxes.append((rxn_id, abs(flux_val)))
                    if abs(flux_val) > threshold:
                        exchange_data.append({
                            'reaction': rxn_id,
                            'flux': flux_val,
                            'metabolite': _parse_exchange_metabolite_id(rxn_id),
                            'direction': 'export' if flux_val > 0 else 'import'
                        })
            # Alternate orientation: rows reactions, columns taxa or flux.
            elif len(set(map(str, fluxes.index)) & exchange_set) > 0:
                matched = [rxn_id for rxn_id in exchange_rxn_ids if rxn_id in fluxes.index]
                print(f"[INFO] Found {len(matched)}/{len(exchange_rxn_ids)} exchange reactions in solution index")
                for rxn_id in matched:
                    vals = fluxes.loc[rxn_id]
                    flux_val = float(pd.to_numeric(vals, errors='coerce').fillna(0).sum())
                    all_fluxes.append((rxn_id, abs(flux_val)))
                    if abs(flux_val) > threshold:
                        exchange_data.append({
                            'reaction': rxn_id,
                            'flux': flux_val,
                            'metabolite': _parse_exchange_metabolite_id(rxn_id),
                            'direction': 'export' if flux_val > 0 else 'import'
                        })
        else:
            # Series/dict-like flux vector.
            for rxn_id in exchange_rxn_ids:
                try:
                    flux_val = float(fluxes.get(rxn_id, 0.0))
                except AttributeError:
                    flux_val = 0.0
                all_fluxes.append((rxn_id, abs(flux_val)))
                if abs(flux_val) > threshold:
                    exchange_data.append({
                        'reaction': rxn_id,
                        'flux': flux_val,
                        'metabolite': _parse_exchange_metabolite_id(rxn_id),
                        'direction': 'export' if flux_val > 0 else 'import'
                    })

    if not exchange_data:
        all_fluxes.sort(key=lambda x: x[1], reverse=True)
        print(f"[DIAGNOSTIC] No exchange fluxes above threshold {threshold}")
        print("[DIAGNOSTIC] Top 10 exchange flux magnitudes:")
        for rxn_id, flux_mag in all_fluxes[:10]:
            print(f"  {rxn_id}: {flux_mag:.2e}")
        return pd.DataFrame(columns=['reaction', 'flux', 'metabolite', 'direction'])

    exchange_fluxes = pd.DataFrame(exchange_data)
    print(f"[SUCCESS] Found {len(exchange_fluxes)} active community exchanges")
    print("[INFO] Top 5 active exchanges by |flux|:")
    for _, ex in exchange_fluxes.reindex(exchange_fluxes['flux'].abs().sort_values(ascending=False).index).head(5).iterrows():
        print(f"  {ex['reaction']:30s}: {ex['flux']:10.6f} ({ex['direction']})")
    return exchange_fluxes


def simulate_community_growth(
    community: Community,
    tradeoff: float = 0.5,
    strategy: str = 'linear',
    save_dir: Optional[str] = None,
    condition: Optional[str] = None,
    pfba: bool = True,
    exchange_threshold: float = 1e-9,
) -> Optional[pd.DataFrame]:
    """
    Simulate community growth and extract exchange fluxes.

    MICOM interpretation correction:
    --------------------------------
    ``cooperative_tradeoff(fraction=x)`` does not mean x is a linear weight
    between individual and community objectives. It means the solution must
    retain at least x of the maximum community growth rate while finding a
    cooperative tradeoff across taxa. Use sensitivity analysis across fraction
    values before making strong biological claims.
    """
    if not (0 < tradeoff <= 1):
        raise ValueError("MICOM cooperative_tradeoff fraction must be in (0, 1].")

    print(f"\n[INFO] Simulating community growth...")
    print(f"[INFO] MICOM cooperative_tradeoff fraction: {tradeoff}")
    print(f"[INFO] pFBA for exchange fluxes: {pfba}")

    try:
        solution = community.cooperative_tradeoff(
            fraction=tradeoff,
            fluxes=True,
            pfba=pfba
        )
        print("[INFO] Using MICOM cooperative_tradeoff optimization")
    except Exception as e:
        if "only supports linear" in str(e).lower() or "quadratic" in str(e).lower():
            print("[WARNING] Cooperative tradeoff failed, likely due to solver limitations.")
            print(f"[WARNING] Original error: {e}")
            print("[INFO] Falling back to simple community.optimize(); interpret exchange fluxes cautiously.")
            try:
                solution = community.optimize()
            except Exception as e2:
                print(f"[ERROR] Linear optimization also failed: {e2}")
                return None
        elif "Unable to retrieve attribute" in str(e):
            print("[ERROR] Optimization failed, likely infeasible model or invalid medium.")
            print("[INFO] Trying relaxed tolerances for cooperative_tradeoff...")
            try:
                solution = community.cooperative_tradeoff(
                    fraction=tradeoff,
                    fluxes=True,
                    pfba=pfba,
                    atol=1e-4,
                    rtol=1e-4
                )
            except Exception as e2:
                print(f"[ERROR] Relaxed optimization also failed: {e2}")
                return None
        else:
            print(f"[ERROR] Community simulation failed: {e}")
            traceback.print_exc()
            return None

    if solution is None:
        print("[WARNING] Community simulation returned None")
        return None

    # Extract a scalar growth summary robustly.
    total_growth = np.nan
    if hasattr(solution, 'growth_rate'):
        try:
            total_growth = float(solution.growth_rate.sum()) if hasattr(solution.growth_rate, 'sum') else float(solution.growth_rate)
        except Exception:
            pass
    if not np.isfinite(total_growth) and hasattr(solution, 'objective_value'):
        try:
            total_growth = float(solution.objective_value)
        except Exception:
            pass

    if np.isfinite(total_growth):
        print(f"[INFO] Community growth summary: {total_growth:.6g}")
        if total_growth < 1e-6:
            print("[WARNING] Community growth is near zero. Exchange predictions may not be biologically meaningful.")
    else:
        print("[WARNING] Could not retrieve community growth summary from solution object.")

    if save_dir is not None and condition is not None and hasattr(solution, 'fluxes'):
        try:
            _sol_path = os.path.join(save_dir, f'community_solution_fluxes_{condition}.csv')
            if hasattr(solution.fluxes, 'to_csv'):
                solution.fluxes.to_csv(_sol_path)
                print(f"[SAVED] Per-taxon/raw solution fluxes: {_sol_path}")
        except Exception as _e:
            print(f"[WARNING] Could not save raw solution fluxes: {_e}")

    exchange_fluxes = _aggregate_exchange_fluxes(solution, community, threshold=exchange_threshold)
    print(f"[INFO] Number of active exchanges: {len(exchange_fluxes)}")
    return exchange_fluxes


def _match_metabolite(exchange_fluxes: pd.DataFrame, met_id: str) -> pd.DataFrame:
    """
    Exact-match a portal metabolite to its exchange reaction / metabolite ID.

    Replaces the previous unanchored ``str.contains`` matching, under which the
    short ID ``ac`` (acetate) matched every longer ID sharing that substring --
    e.g. ``nac`` (niacin), ``lac__L`` / ``lac__D`` (lactate) -- silently aliasing
    niacin's flux onto acetate. Metabolite IDs form a namespace in which short
    IDs are substrings of longer ones, so exact matching is the only
    collision-safe contract.

    For a given ``met_id`` this matches:
      * reaction IDs of the form ``EX_<met_id>`` optionally followed by a single
        compartment tag (``_e`` / ``_m`` / ``_c`` / ... , optionally one digit)
        or a bracketed compartment (``(e)`` / ``[e]`` / ...); anchored via
        ``fullmatch`` so ``EX_ac`` can never absorb ``EX_nac``.
      * metabolite IDs equal to ``met_id`` or ``met_id`` plus a compartment tag,
        via exact ``isin`` (never substring).
    """
    esc = re.escape(met_id)
    rxn_pat = rf'EX_{esc}(?:_[a-z]\d?|\([a-z]\d?\)|\[[a-z]\d?\])?'
    rxn_ok = exchange_fluxes['reaction'].str.fullmatch(rxn_pat, na=False)
    met_ok = exchange_fluxes['metabolite'].isin({
        met_id,
        f'{met_id}_e', f'{met_id}_m', f'{met_id}_c', f'{met_id}_p',
        f'{met_id}[e]', f'{met_id}(e)', f'{met_id}[m]', f'{met_id}(m)',
    })
    return exchange_fluxes[rxn_ok | met_ok]


def _load_cobra_model_any(path: str):
    """Load a COBRA model from JSON or MATLAB .mat."""
    if not COBRA_AVAILABLE:
        raise ImportError("COBRApy is required to load/validate hepatic models. Install with: pip install cobra")
    suffix = Path(path).suffix.lower()
    if suffix == '.json':
        return load_json_model(path)
    if suffix == '.mat':
        return load_matlab_model(path)
    raise ValueError(f"Unsupported model format for {path}. Expected .json or .mat")


def get_hepatic_exchange_set(hepatic_model_path: Optional[str]) -> Optional[set]:
    """Return hepatic reaction IDs if a hepatic model is supplied."""
    if hepatic_model_path is None:
        return None
    print(f"[INFO] Validating portal hepatic reactions against hepatic model: {hepatic_model_path}")
    model = _load_cobra_model_any(hepatic_model_path)
    rxns = {rxn.id for rxn in model.reactions}
    print(f"[INFO] Hepatic model loaded with {len(rxns)} reactions")
    return rxns


def build_portal_metabolite_audit(
    exchange_fluxes: pd.DataFrame,
    portal_fluxes: Dict[str, float],
    hepatic_exchange_set: Optional[set] = None,
) -> pd.DataFrame:
    """Build an audit table for curated, absent, excluded, and invalid metabolites."""
    rows = []

    for met_id, meta in PORTAL_METABOLITES.items():
        flux_rows = _match_metabolite(exchange_fluxes, met_id) if exchange_fluxes is not None and len(exchange_fluxes) else pd.DataFrame()
        total_flux = float(flux_rows['flux'].sum()) if len(flux_rows) else 0.0
        hepatic_rxn = meta.get('hepatic_rxn', '')
        hepatic_valid = None if hepatic_exchange_set is None else hepatic_rxn in hepatic_exchange_set
        included = met_id in portal_fluxes
        status = 'included' if included else 'absent_or_below_threshold'
        if hepatic_valid is False:
            status = 'excluded_missing_hepatic_exchange'
        rows.append({
            'metabolite_id': met_id,
            'name': meta.get('name', met_id),
            'category': meta.get('category', ''),
            'importance': meta.get('importance', ''),
            'matched_exchange_rows': int(len(flux_rows)),
            'summed_exchange_flux': total_flux,
            'direction_if_active': 'export' if total_flux > 0 else ('import' if total_flux < 0 else 'none'),
            'hepatic_rxn': hepatic_rxn,
            'hepatic_rxn_validated': hepatic_valid,
            'included_in_portal_json': included,
            'status': status,
            'interpretation_note': meta.get('interpretation_note', ''),
        })

    for met_id, meta in EXCLUDED_PORTAL_METABOLITES.items():
        flux_rows = _match_metabolite(exchange_fluxes, met_id) if exchange_fluxes is not None and len(exchange_fluxes) else pd.DataFrame()
        total_flux = float(flux_rows['flux'].sum()) if len(flux_rows) else 0.0
        rows.append({
            'metabolite_id': met_id,
            'name': meta.get('previous_mislabel', met_id),
            'category': 'excluded_previous_candidate',
            'importance': 'excluded',
            'matched_exchange_rows': int(len(flux_rows)),
            'summed_exchange_flux': total_flux,
            'direction_if_active': 'export' if total_flux > 0 else ('import' if total_flux < 0 else 'none'),
            'hepatic_rxn': '',
            'hepatic_rxn_validated': False,
            'included_in_portal_json': False,
            'status': 'intentionally_excluded_even_if_detected',
            'interpretation_note': meta.get('reason', ''),
        })

    return pd.DataFrame(rows)


def extract_portal_metabolites(
    exchange_fluxes: pd.DataFrame,
    threshold: float = 1e-6,
    filter_metabolites: bool = True,
    hepatic_exchange_set: Optional[set] = None,
    audit_path: Optional[str] = None,
) -> Dict[str, float]:
    """
    Extract biologically curated portal metabolites from community exchanges.

    Corrections relative to the previous version:
    ---------------------------------------------
    - Uses exact metabolite matching only; no substring matching.
    - Does not map deoxycholate to hdca/ocdca/ptdca fatty acids.
    - Does not treat IMP as imidazole propionate.
    - If a hepatic model is supplied, excludes metabolites whose hepatic
      exchange reactions are absent.

    Flux sign convention:
    ---------------------
    Positive community exchange flux is treated as microbial export/secretion
    to the shared medium; negative flux is microbial import/consumption from the
    medium/lumen. Downstream hepatic integration should normally constrain the
    liver only with microbial production/export events.
    """
    portal_fluxes = {}

    if exchange_fluxes is None or len(exchange_fluxes) == 0:
        print("[WARNING] No exchange fluxes provided - returning empty portal metabolite set")
        if audit_path:
            build_portal_metabolite_audit(pd.DataFrame(), portal_fluxes, hepatic_exchange_set).to_csv(audit_path, index=False)
            print(f"[SAVED] Portal metabolite audit: {audit_path}")
        return portal_fluxes

    required = {'reaction', 'flux', 'metabolite'}
    if not required.issubset(exchange_fluxes.columns):
        print(f"[WARNING] Exchange flux table missing required columns {required}. Has: {exchange_fluxes.columns.tolist()}")
        return portal_fluxes

    if filter_metabolites:
        print(f"[INFO] Extracting curated portal metabolites: {len(PORTAL_METABOLITES)} candidates")
        for met_id, meta in PORTAL_METABOLITES.items():
            flux_rows = _match_metabolite(exchange_fluxes, met_id)
            if len(flux_rows) == 0:
                continue
            total_flux = float(flux_rows['flux'].sum())
            if abs(total_flux) <= threshold:
                continue

            hepatic_rxn = meta.get('hepatic_rxn')
            if hepatic_exchange_set is not None and hepatic_rxn not in hepatic_exchange_set:
                print(
                    f"[WARNING] Excluding {met_id} ({meta.get('name', met_id)}): "
                    f"hepatic reaction {hepatic_rxn} not found in hepatic model"
                )
                continue

            portal_fluxes[met_id] = total_flux

        # Notify if excluded unsafe IDs are active.
        for met_id, meta in EXCLUDED_PORTAL_METABOLITES.items():
            flux_rows = _match_metabolite(exchange_fluxes, met_id)
            if len(flux_rows):
                total_flux = float(flux_rows['flux'].sum())
                if abs(total_flux) > threshold:
                    print(
                        f"[BIOLOGY WARNING] Detected active exchange for excluded ID '{met_id}' "
                        f"({total_flux:+.4g}) but did not include it: {meta['reason']}"
                    )
    else:
        print("[INFO] Extracting ALL metabolites above threshold (unfiltered exploratory mode)")
        for _, row in exchange_fluxes.iterrows():
            if abs(row['flux']) > threshold:
                met_id = str(row['metabolite'])
                if met_id in EXCLUDED_PORTAL_METABOLITES:
                    print(f"[INFO] Unfiltered mode includes {met_id}; it remains biologically excluded from curated interpretation.")
                portal_fluxes[met_id] = float(row['flux'])
        print(f"[INFO] Found {len(portal_fluxes)} metabolites above threshold")

    if filter_metabolites:
        print("\n[SUMMARY] Curated portal metabolites extracted by priority:")
        for tier, met_list in PRIORITY_TIERS.items():
            tier_metabolites = [m for m in met_list if m in portal_fluxes]
            if tier_metabolites:
                print(f"  {tier}: {len(tier_metabolites)} metabolites")
                for met_id in tier_metabolites:
                    flux = portal_fluxes[met_id]
                    direction = "microbial export/production" if flux > 0 else "microbial import/consumption"
                    print(f"    {met_id:12s} ({PORTAL_METABOLITES[met_id]['name']:35s}): {flux:+10.4f} ({direction})")

    audit_df = build_portal_metabolite_audit(exchange_fluxes, portal_fluxes, hepatic_exchange_set)
    if audit_path:
        audit_df.to_csv(audit_path, index=False)
        print(f"[SAVED] Portal metabolite audit: {audit_path}")

    return portal_fluxes


################################################################################
# VISUALIZATION AND REPORTING
################################################################################

def compare_diet_metabolite_production(
    portal_metabolites_dict: Dict[str, Dict[str, float]],
    results_dir: str
):
    """
    Generate comparative visualizations of portal metabolite production.
    
    Figures:
    -------
    1. Heatmap: Metabolites × Diets
    2. Bar plots: Per-metabolite comparisons
    3. Network diagram: Metabolite relationships
    
    Parameters
    ----------
    portal_metabolites_dict : dict
        {condition: {metabolite_id: flux}}
    results_dir : str
        Output directory for figures
    """
    if len(portal_metabolites_dict) == 0:
        print("[WARNING] No portal metabolite data to plot")
        return
    
    # Create results dataframe
    rows = []
    for condition, metabolites in portal_metabolites_dict.items():
        for met_id, flux in metabolites.items():
            # Handle both filtered and unfiltered metabolites
            if met_id in PORTAL_METABOLITES:
                met_name = PORTAL_METABOLITES[met_id]['name']
                importance = PORTAL_METABOLITES[met_id]['importance']
            else:
                met_name = met_id  # Use ID as name for unfiltered metabolites
                importance = 'unknown'
            
            rows.append({
                'Condition': condition,
                'Metabolite': met_name,
                'Metabolite_ID': met_id,
                'Flux': flux,
                'Importance': importance
            })
    
    df = pd.DataFrame(rows)
    
    if len(df) == 0:
        print("[WARNING] No metabolite fluxes to visualize")
        return
    
    # Save detailed table
    table_path = os.path.join(results_dir, 'portal_metabolite_production.csv')
    df.to_csv(table_path, index=False)
    print(f"[SAVED] Portal metabolite table: {table_path}")
    
    # Figure 1: Heatmap
    plt.figure(figsize=(10, 8))
    pivot_df = df.pivot(index='Metabolite', columns='Condition', values='Flux')
    sns.heatmap(
        pivot_df,
        annot=True,
        fmt='.3f',
        cmap='RdYlGn',
        center=0,
        cbar_kws={'label': 'Community exchange flux (positive = export)'}
    )
    plt.title('Curated Portal Metabolite Exchange Across Conditions')
    plt.tight_layout()
    heatmap_path = os.path.join(results_dir, 'portal_metabolites_heatmap.png')
    plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[SAVED] Heatmap: {heatmap_path}")
    
    # Figure 2: Bar plot comparison (high importance metabolites only)
    high_importance = df[df['Importance'] == 'high']
    if len(high_importance) > 0:
        plt.figure(figsize=(12, 6))
        sns.barplot(
            data=high_importance,
            x='Metabolite',
            y='Flux',
            hue='Condition'
        )
        plt.ylabel('Community exchange flux (positive = export)')
        plt.title('Key Curated Portal Metabolite Exchange')
        plt.xticks(rotation=45, ha='right')
        plt.legend(title='Diet')
        plt.tight_layout()
        barplot_path = os.path.join(results_dir, 'portal_metabolites_barplot.png')
        plt.savefig(barplot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[SAVED] Bar plot: {barplot_path}")


################################################################################
# COMMAND LINE INTERFACE
################################################################################

def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description='RQ4: Biologically corrected microbiome community metabolic modeling with MICOM',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
---------
# Recommended: explicit species mapping + hepatic reaction validation
python rq4_microbiome_community_modeling_v2026_biocorrected.py \
    --metatranscriptome Meta_GSE104913.csv \
    --agora_dir /path/to/AGORA2 \
    --species_model_map curated_species_to_agora.csv \
    --hepatic_model iMM1415.json \
    --expression_cols ND_SCD DD_HFD \
    --medium_json expanded_diet_bounds_flat.json \
    --results_dir results_community_biocorrected

# Sensitivity example
python rq4_microbiome_community_modeling_v2026_biocorrected.py \
    --metatranscriptome Meta_GSE104913.csv \
    --agora_dir /path/to/AGORA2 \
    --expression_cols ND_SCD DD_HFD \
    --tradeoff 0.7 \
    --no-pfba \
    --results_dir results_tradeoff_0p7_no_pfba
        """
    )

    parser.add_argument('--metatranscriptome', required=True, help='Path to metatranscriptome CSV file')
    parser.add_argument('--agora_dir', required=True, help='Directory containing AGORA model files (.mat or .json)')
    parser.add_argument('--expression_cols', nargs='+', default=['ND_SCD', 'DD_HFD'], help='Column names for expression/activity values')
    parser.add_argument('--medium_json', default=None, help='Optional flat or condition/diet-specific medium JSON file')
    parser.add_argument('--hepatic_model', default=None, help='Optional hepatic COBRA model (.json/.mat) used to validate hepatic exchange reaction IDs before writing portal JSON')
    parser.add_argument('--species_model_map', default=None, help='Optional curated CSV/TSV mapping from input species to AGORA model files')
    parser.add_argument('--allow_fuzzy_model_matching', action='store_true', help='Allow provisional fuzzy species-to-model matching when no explicit/exact match is available')
    parser.add_argument('--strict_species_map', action='store_true',
                        help='Treat --species_model_map as a strict whitelist; recommended for final publication reruns.')
    parser.add_argument('--expected_mapped_species', type=int, default=None,
                        help='Hard QC: require exactly this many species to map before community modeling.')
    parser.add_argument('--min_retained_abundance', type=float, default=0.80, help='Warn if retained mapped community activity before renormalization is below this fraction')
    parser.add_argument('--abundance_method', choices=['total_expression', 'housekeeping'], default='total_expression', help='Method for inferring activity-weighted community composition')
    parser.add_argument('--tradeoff', type=float, default=0.5, help='MICOM cooperative_tradeoff fraction: minimum fraction of maximum community growth retained; use sensitivity analysis')
    parser.add_argument('--pfba', dest='pfba', action='store_true', default=True, help='Use pFBA in cooperative_tradeoff for parsimonious exchange fluxes (default)')
    parser.add_argument('--no-pfba', dest='pfba', action='store_false', help='Disable pFBA in cooperative_tradeoff for sensitivity analysis')
    parser.add_argument('--solver', choices=['glpk', 'gurobi', 'cplex'], default='gurobi', help='LP/QP solver backend')
    parser.add_argument('--portal_threshold', type=float, default=1e-6, help='Minimum absolute curated exchange flux to include in portal JSON')
    parser.add_argument('--exchange_threshold', type=float, default=1e-9, help='Minimum absolute exchange flux to keep in active exchange table')
    parser.add_argument('--filter_metabolites', action='store_true', default=True, help='Filter portal metabolites to curated biologically meaningful candidates (default)')
    parser.add_argument('--no-filter_metabolites', dest='filter_metabolites', action='store_false', help='Extract all exchange metabolites above threshold; exploratory only')
    parser.add_argument('--results_dir', default='results_microbiome_community_biocorrected', help='Output directory for results')

    args = parser.parse_args()
    if not (0 < args.tradeoff <= 1):
        parser.error('--tradeoff must be in (0, 1]. In MICOM it is a fraction of maximum community growth retained.')
    if not (0 <= args.min_retained_abundance <= 1):
        parser.error('--min_retained_abundance must be between 0 and 1')
    return args


################################################################################
# MAIN ANALYSIS WORKFLOW
################################################################################

def main():
    """
    Main workflow for biologically corrected microbiome community modeling.
    """
    args = parse_arguments()
    os.makedirs(args.results_dir, exist_ok=True)

    print("="*80)
    print("RQ4: MICROBIOME COMMUNITY METABOLIC MODELING (BIOLOGICALLY CORRECTED)")
    print("="*80)
    print("[METHOD NOTE] Species composition is a metatranscriptome-derived activity proxy, not measured biomass abundance.")
    print("[METHOD NOTE] Curated portal extraction uses exact metabolite matching and excludes unsafe DCA/IMP proxy assumptions.")

    hepatic_exchange_set = get_hepatic_exchange_set(args.hepatic_model)

    # Step 1: Load metatranscriptome/activity table.
    metatranscriptome_df = load_metatranscriptome_data(
        args.metatranscriptome,
        args.expression_cols
    )

    # Step 2: Infer activity-weighted community composition.
    abundance_df = infer_species_abundances(
        metatranscriptome_df,
        args.expression_cols,
        args.abundance_method
    )
    abundance_path = os.path.join(args.results_dir, 'species_activity_weighted_composition.csv')
    abundance_df.to_csv(abundance_path, index=False)
    print(f"[SAVED] Activity-weighted species composition: {abundance_path}")

    # Step 3: Map species to AGORA models.
    species_list = abundance_df['species'].astype(str).tolist()
    species_to_model, mapping_report = map_species_to_agora_models(
        species_list,
        args.agora_dir,
        mapping_file=args.species_model_map,
        allow_fuzzy=args.allow_fuzzy_model_matching,
        strict_mapping_file=args.strict_species_map,
        results_dir=args.results_dir,
    )

    mapping_path = os.path.join(args.results_dir, 'species_to_model_mapping.csv')
    mapping_report.to_csv(mapping_path, index=False)
    print(f"[SAVED] Species-model mapping: {mapping_path}")

    n_mapped = int(mapping_report['match_status'].isin(['matched','matched_provisional']).sum())
    if args.expected_mapped_species is not None and n_mapped != int(args.expected_mapped_species):
        raise RuntimeError(
            f"Species mapping QC failed: mapped {n_mapped}, expected {args.expected_mapped_species}. "
            "Community modeling aborted before MICOM."
        )
    if args.strict_species_map and n_mapped == 0:
        raise RuntimeError("Strict species map produced zero mapped species; aborting before MICOM.")

    # Step 4: Load medium configuration once; select per condition inside loop.
    medium_config = load_medium_config(args.medium_json)

    if not MICOM_AVAILABLE:
        print("\n[ERROR] MICOM is required for community modeling. Install with: pip install micom")
        sys.exit(1)

    portal_metabolites_dict = {}
    medium_sources = {}

    for condition in args.expression_cols:
        print(f"\n{'='*80}")
        print(f"Processing condition: {condition}")
        print(f"{'='*80}")

        try:
            condition_medium, medium_source = get_condition_medium(
                medium_config,
                condition,
                target_compartment='e'
            )
            medium_sources[condition] = medium_source
        except Exception as e:
            print(f"[ERROR] Could not prepare medium for {condition}: {e}")
            continue

        community = build_micom_community(
            abundance_df,
            species_to_model,
            condition,
            condition_medium,
            args.solver,
            results_dir=args.results_dir,
            min_retained_abundance=args.min_retained_abundance,
        )
        if community is None:
            print(f"[WARNING] Skipping {condition} - community building failed")
            continue

        exchange_fluxes = simulate_community_growth(
            community,
            tradeoff=args.tradeoff,
            save_dir=args.results_dir,
            condition=condition,
            pfba=args.pfba,
            exchange_threshold=args.exchange_threshold,
        )
        if exchange_fluxes is None:
            print(f"[WARNING] Skipping {condition} - simulation failed")
            continue

        raw_path = os.path.join(args.results_dir, f'community_exchange_fluxes_{condition}.csv')
        exchange_fluxes.to_csv(raw_path, index=False)
        print(f"[SAVED] Active community exchange fluxes ({len(exchange_fluxes)} rows): {raw_path}")

        audit_path = os.path.join(args.results_dir, f'portal_metabolite_audit_{condition}.csv')
        portal_fluxes = extract_portal_metabolites(
            exchange_fluxes,
            threshold=args.portal_threshold,
            filter_metabolites=args.filter_metabolites,
            hepatic_exchange_set=hepatic_exchange_set,
            audit_path=audit_path,
        )
        portal_metabolites_dict[condition] = portal_fluxes

    # Step 6: Save portal metabolite fluxes for hepatic integration.
    portal_output = {}
    for condition, metabolites in portal_metabolites_dict.items():
        portal_output[condition] = {}
        for met_id, flux in metabolites.items():
            if met_id in PORTAL_METABOLITES:
                meta = PORTAL_METABOLITES[met_id]
                portal_output[condition][met_id] = {
                    'flux': float(flux),
                    'name': meta['name'],
                    'hepatic_rxn': meta['hepatic_rxn'],
                    'importance': meta['importance'],
                    'category': meta['category'],
                    'flux_sign_convention': 'positive = microbial export/production to shared medium; negative = microbial import/consumption',
                    'medium_source': medium_sources.get(condition, ''),
                    'hepatic_rxn_validated': None if hepatic_exchange_set is None else meta['hepatic_rxn'] in hepatic_exchange_set,
                }
            else:
                # Exploratory unfiltered mode. Keep metadata conservative.
                portal_output[condition][met_id] = {
                    'flux': float(flux),
                    'name': met_id,
                    'hepatic_rxn': f'EX_{met_id}_e',
                    'importance': 'unknown',
                    'category': 'unfiltered_exploratory',
                    'flux_sign_convention': 'positive = microbial export/production to shared medium; negative = microbial import/consumption',
                    'medium_source': medium_sources.get(condition, ''),
                    'hepatic_rxn_validated': None,
                }

    portal_json_path = os.path.join(args.results_dir, 'portal_metabolites_for_hepatic_model.json')
    with open(portal_json_path, 'w') as f:
        json.dump(portal_output, f, indent=2)
    print(f"\n[SAVED] Portal metabolites for hepatic integration: {portal_json_path}")

    # Save run metadata.
    metadata_path = os.path.join(args.results_dir, 'run_metadata_biocorrected.json')
    with open(metadata_path, 'w') as f:
        json.dump({
            'script': os.path.basename(__file__),
            'metatranscriptome': args.metatranscriptome,
            'agora_dir': args.agora_dir,
            'species_model_map': args.species_model_map,
            'allow_fuzzy_model_matching': args.allow_fuzzy_model_matching,
            'hepatic_model': args.hepatic_model,
            'expression_cols': args.expression_cols,
            'medium_json': args.medium_json,
            'medium_sources_by_condition': medium_sources,
            'abundance_interpretation': 'metatranscriptome-derived activity-weighted composition proxy, not measured biomass abundance',
            'tradeoff_fraction': args.tradeoff,
            'pfba': args.pfba,
            'portal_threshold': args.portal_threshold,
            'exchange_threshold': args.exchange_threshold,
            'excluded_previous_candidates': EXCLUDED_PORTAL_METABOLITES,
        }, f, indent=2)
    print(f"[SAVED] Run metadata: {metadata_path}")

    # Print summary statistics.
    print(f"\n{'='*80}")
    print("CURATED PORTAL METABOLITE EXTRACTION SUMMARY")
    print(f"{'='*80}")
    for condition, metabolites in portal_output.items():
        print(f"\n{condition}:")
        print(f"  Total metabolites in portal JSON: {len(metabolites)}")
        by_importance = {}
        for met_data in metabolites.values():
            importance = met_data['importance']
            by_importance[importance] = by_importance.get(importance, 0) + 1
        print(f"  By importance: {dict(by_importance)}")
        production = sum(1 for m in metabolites.values() if m['flux'] > 0)
        consumption = sum(1 for m in metabolites.values() if m['flux'] < 0)
        print(f"  Microbial export/production: {production}, microbial import/consumption: {consumption}")
        top_metabolites = sorted(metabolites.items(), key=lambda x: abs(x[1]['flux']), reverse=True)[:5]
        print("  Top 5 by |exchange flux|:")
        for met_id, met_data in top_metabolites:
            direction = "export" if met_data['flux'] > 0 else "import"
            print(f"    {met_id:12s} ({met_data['name'][:35]:35s}): {met_data['flux']:+10.4f} ({direction})")

    if len(portal_metabolites_dict) > 0:
        compare_diet_metabolite_production(portal_metabolites_dict, args.results_dir)

    print("\n" + "="*80)
    print("RQ4 MICROBIOME COMMUNITY MODELING COMPLETE")
    print("="*80)
    print(f"Results saved to: {args.results_dir}")
    print("\nRecommended manuscript wording:")
    print("- 'metatranscriptome-derived activity-weighted community composition proxy'")
    print("- 'model-predicted community exchange fluxes'")
    print("- 'curated portal candidate metabolites with exact ID matching'")
    print("Avoid claiming measured abundance, imidazole propionate, or DCA proxy effects unless independently validated.")


if __name__ == '__main__':
    main()

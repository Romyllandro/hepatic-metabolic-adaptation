#!/usr/bin/env python3
"""
Cross-RQ Integration Analysis Script
=====================================
Generates all cross-RQ integration tables for the manuscript.

This script produces:
1. Cross-RQ Integration Tables (RQ1→RQ2, RQ1→RQ3, RQ2↔RQ3, RQ1→RQ4, RQ3→RQ4)
2. Reaction Tracing Tables (15 RQ1-HFD reactions through pipeline)
3. Microbiome/Diet/Antagonistic Reaction Tables
4. Summary Statistics Table
5. Microbiome Community Tables

Author: Roland Chen
Date: 2026
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

DATA_DIR = Path('Final_run_3142026/RQ_outputs_files')
# OUTPUT_DIR = Path('Final_run_3142026/ALL_RQs_Tables_CrossRQ')
OUTPUT_DIR = Path('Final_run_3142026/ALL_RQs_Tables_CrossRQ_Enhanced_4222026')
OUTPUT_DIR.mkdir(exist_ok=True)

STRAINS = ['129S1SvImJ', 'AJ', 'C57BL6J', 'CASTEiJ', 'DBA2J', 
           'NODShiLtJ', 'NZOHlLtJ', 'PWKPhJ', 'WSBEiJ']

# ============================================================================
# DATA LOADING
# ============================================================================

def load_all_data():
    """Load all required data files."""
    print("Loading data files...")
    data = {}
    
    """
    # RQ1 data
    data['rq1_flux'] = pd.read_csv(DATA_DIR / 'RQ1_multidataset_reaction_flux_comparison_extended.csv')
    data['rq1_rank'] = pd.read_csv(DATA_DIR / 'RQ1_multidataset_rank_product.csv')
    data['rq1_edges_hfd'] = pd.read_csv(DATA_DIR / 'RQ1_multidataset_edges_HFD_vs_SCD.csv')
    data['rq1_edges_kd'] = pd.read_csv(DATA_DIR / 'RQ1_multidataset_edges_KD_vs_SCD.csv')
    data['rq1_edges_wd'] = pd.read_csv(DATA_DIR / 'RQ1_multidataset_edges_WD_vs_SCD.csv')
    data['rq1_stats'] = pd.read_csv(DATA_DIR / 'RQ1_multidataset_flux_pairwise_stats.csv')
    """
    
    # RQ1 data
    data['rq1_flux'] = pd.read_csv(DATA_DIR / 'RQ1_reaction_flux_comparison_extended.csv')
    data['rq1_rank'] = pd.read_csv(DATA_DIR / 'RQ1_rank_product.csv')
    data['rq1_edges_hfd'] = pd.read_csv(DATA_DIR / 'RQ1_edges_HFD_vs_SCD.csv')
    data['rq1_edges_kd'] = pd.read_csv(DATA_DIR / 'RQ1_edges_KD_vs_SCD.csv')
    data['rq1_edges_wd'] = pd.read_csv(DATA_DIR / 'RQ1_edges_WD_vs_SCD.csv')
    data['rq1_stats'] = pd.read_csv(DATA_DIR / 'RQ1_flux_pairwise_stats.csv')
    
    
    # RQ2 data - strain-specific
    data['rq2_edges'] = {}
    data['rq2_rank'] = {}
    data['rq2_flux'] = {}
    for strain in STRAINS:
        data['rq2_edges'][strain] = pd.read_csv(DATA_DIR / f'RQ2_{strain}_edges_HFD_vs_SCD.csv')
        data['rq2_rank'][strain] = pd.read_csv(DATA_DIR / f'RQ2_{strain}_rank_product.csv')
        data['rq2_flux'][strain] = pd.read_csv(DATA_DIR / f'RQ2_{strain}_reaction_flux_comparison_extended.csv')
    
    # RQ3 data
    data['rq3_contrib'] = pd.read_csv(DATA_DIR / 'RQ3_phase2_contribution_analysis.csv')
    data['rq3_pathway'] = pd.read_csv(DATA_DIR / 'RQ3_phase4_pathway_attribution.csv')
    data['rq3_overlap'] = pd.read_csv(DATA_DIR / 'RQ3_phase1_bulk_cellular_overlap.csv')
    data['rq3_conservation'] = pd.read_csv(DATA_DIR / 'RQ3_conservation_analysis.csv')
    data['rq3_driver'] = pd.read_csv(DATA_DIR / 'RQ3_driver_consistency.csv')
    data['rq3_attribution'] = pd.read_csv(DATA_DIR / 'RQ3_phase3_reaction_attribution.csv')
    data['rq3_all_contrib'] = pd.read_csv(DATA_DIR / 'RQ3_all_strain_contributions.csv')
    data['rq3_flux'] = pd.read_csv(DATA_DIR / 'RQ3_flux_comparison.csv')
    
    # RQ4 data
    data['rq4_attrib'] = pd.read_csv(DATA_DIR / 'RQ4_flux_attribution_analysis.csv')
    data['rq4_synergy'] = pd.read_csv(DATA_DIR / 'RQ4_pathway_synergy_analysis.csv')
    data['rq4_species'] = pd.read_csv(DATA_DIR / 'RQ4_species_abundances.csv')
    data['rq4_portal'] = pd.read_csv(DATA_DIR / 'RQ4_portal_metabolite_production.csv')
    data['rq4_pathway_enrich'] = pd.read_csv(DATA_DIR / 'RQ4_pathway_enrichment_results.csv')
    data['rq4_annotations'] = pd.read_csv(DATA_DIR / 'RQ4_reaction_annotations.csv')
    
    print(f"Loaded {len(data)} data categories")
    return data


# ============================================================================
# CROSS-RQ INTEGRATION TABLES
# ============================================================================

def generate_rq1_to_rq2_table(data):
    """
    RQ1 → RQ2: Genetic Conservation of Diet Signatures
    """
    print("\n" + "="*70)
    print("RQ1 → RQ2: Genetic Conservation of Diet Signatures")
    print("="*70)
    
    rq1_hfd_rxns = set(data['rq1_edges_hfd']['ReactionID'].unique())
    
    # Calculate conservation for each RQ1-HFD reaction
    results = []
    for rxn in rq1_hfd_rxns:
        strain_count = 0
        strains_present = []
        for strain in STRAINS:
            strain_rxns = set(data['rq2_edges'][strain]['ReactionID'].unique())
            if rxn in strain_rxns:
                strain_count += 1
                strains_present.append(strain)
        
        # Get reaction name from RQ1 flux data
        rxn_info = data['rq1_flux'][data['rq1_flux']['ReactionID'] == rxn]
        rxn_name = rxn_info['ReactionName'].iloc[0] if len(rxn_info) > 0 else rxn
        subsystem = rxn_info['Subsystem'].iloc[0] if len(rxn_info) > 0 else 'Unknown'
        
        results.append({
            'ReactionID': rxn,
            'ReactionName': rxn_name[:40],
            'Subsystem': subsystem[:30],
            'Strains_Conserved': strain_count,
            'Conservation_Level': 'Universal' if strain_count == 9 else 
                                  'High (≥5)' if strain_count >= 5 else
                                  'Moderate (3-4)' if strain_count >= 3 else 'Low (<3)'
        })
    
    df = pd.DataFrame(results).sort_values('Strains_Conserved', ascending=False)
    
    # Print table
    print(f"\n{'Metric':<40} {'Value':<20} {'Interpretation':<30}")
    print("-"*90)
    print(f"{'RQ1-HFD reactions':<40} {len(rq1_hfd_rxns):<20} {'Reference set':<30}")
    print(f"{'Conserved in ≥5/9 strains':<40} {(df['Strains_Conserved'] >= 5).sum():<20} {'Genetically robust':<30}")
    print(f"{'Universal (9/9 strains)':<40} {(df['Strains_Conserved'] == 9).sum():<20} {'Fundamental response':<30}")
    
    print(f"\n\nDetailed Reaction Conservation:")
    print("-"*120)
    print(f"{'ReactionID':<15} {'ReactionName':<40} {'Subsystem':<30} {'Strains':<8} {'Level':<15}")
    print("-"*120)
    for _, row in df.iterrows():
        print(f"{row['ReactionID']:<15} {row['ReactionName']:<40} {row['Subsystem']:<30} {row['Strains_Conserved']:<8} {row['Conservation_Level']:<15}")
    
    return df


def generate_rq1_to_rq3_table(data):
    """
    RQ1 → RQ3: Cell-Type Drivers of Diet Signatures
    """
    print("\n" + "="*70)
    print("RQ1 → RQ3: Cell-Type Drivers of Diet Signatures")
    print("="*70)
    
    contrib = data['rq3_contrib'].copy()
    contrib = contrib.sort_values('contribution_percent', ascending=False)
    
    print(f"\n{'Cell Type':<30} {'Contribution (%)':<20} {'Response Rate (%)':<20}")
    print("-"*70)
    for _, row in contrib.iterrows():
        print(f"{row['cell_type']:<30} {row['contribution_percent']:<20.1f} {row['response_rate']*100:<20.1f}")
    
    # Summary statistics
    lec_contrib = contrib[contrib['cell_type'] == 'LECs']['contribution_percent'].values[0]
    hep_contrib = contrib[contrib['cell_type'] == 'Hepatocytes']['contribution_percent'].values[0]
    
    print(f"\n\nSummary:")
    print("-"*70)
    print(f"{'LECs contribution:':<40} {lec_contrib:.1f}%")
    print(f"{'Hepatocytes contribution:':<40} {hep_contrib:.1f}%")
    print(f"{'LEC/Hepatocyte ratio:':<40} {lec_contrib/hep_contrib:.1f}x")
    
    return contrib


def generate_rq2_to_rq3_table(data):
    """
    RQ2 ↔ RQ3: Conservation of Cell-Type Hierarchy Across Strains
    """
    print("\n" + "="*70)
    print("RQ2 ↔ RQ3: Conservation of Cell-Type Hierarchy Across Strains")
    print("="*70)
    
    conservation = data['rq3_conservation'].copy()
    conservation = conservation.sort_values('mean_contribution', ascending=False)
    
    print(f"\n{'Cell Type':<25} {'Mean Contrib':<15} {'Std':<10} {'CV (%)':<10} {'Min Strain':<15} {'Max Strain':<15}")
    print("-"*100)
    for _, row in conservation.iterrows():
        print(f"{row['cell_type']:<25} {row['mean_contribution']:<15.4f} {row['std_contribution']:<10.4f} {row['cv']*100:<10.1f} {row['min_strain']:<15} {row['max_strain']:<15}")
    
    # Key insight
    lec_cv = conservation[conservation['cell_type'] == 'LECs']['cv'].values[0]
    print(f"\n\nKey Insight:")
    print("-"*70)
    print(f"LEC contribution CV across strains: {lec_cv*100:.1f}% (very low = robust)")
    print("LEC dominance is conserved across genetically diverse strains")
    
    return conservation


def generate_rq1_to_rq4_table(data):
    """
    RQ1 → RQ4: Diet vs Microbiome Attribution of Diet Signatures
    """
    print("\n" + "="*70)
    print("RQ1 → RQ4: Diet vs Microbiome Attribution of Diet Signatures")
    print("="*70)
    
    rq1_hfd_rxns = set(data['rq1_edges_hfd']['ReactionID'].unique())
    rq4_attrib = data['rq4_attrib'].copy()
    
    # Filter to RQ1-HFD reactions
    rq4_for_rq1 = rq4_attrib[rq4_attrib['reaction_id'].isin(rq1_hfd_rxns)].copy()
    
    # Add reaction names
    rxn_names = data['rq1_flux'][['ReactionID', 'ReactionName', 'Subsystem']].drop_duplicates()
    rq4_for_rq1 = rq4_for_rq1.merge(rxn_names, left_on='reaction_id', right_on='ReactionID', how='left')
    
    print(f"\n{'ReactionID':<15} {'ReactionName':<35} {'Diet Var%':<12} {'Micro Var%':<12} {'Driver':<25}")
    print("-"*110)
    for _, row in rq4_for_rq1.iterrows():
        rxn_name = str(row.get('ReactionName', row['reaction_id']))[:35]
        print(f"{row['reaction_id']:<15} {rxn_name:<35} {row['variance_explained_diet']:<12.1f} {row['variance_explained_microbiome']:<12.1f} {row['dominant_driver']:<25}")
    
    # Summary
    diet_driven = rq4_for_rq1['dominant_driver'].str.contains('Diet', na=False).sum()
    micro_driven = rq4_for_rq1['dominant_driver'].str.contains('Microbiome', na=False).sum()
    stable = rq4_for_rq1['dominant_driver'].str.contains('Stable', na=False).sum()
    
    print(f"\n\nSummary of {len(rq4_for_rq1)} RQ1-HFD reactions:")
    print("-"*70)
    print(f"{'Diet-driven:':<40} {diet_driven}")
    print(f"{'Microbiome-driven:':<40} {micro_driven}")
    print(f"{'Stable:':<40} {stable}")
    
    return rq4_for_rq1


def generate_rq3_to_rq4_table(data):
    """
    RQ3 → RQ4: Cell Types Predict Microbiome Synergy
    """
    print("\n" + "="*70)
    print("RQ3 → RQ4: Cell Types Predict Microbiome Synergy")
    print("="*70)
    
    synergy = data['rq4_synergy'].copy()
    synergy = synergy.sort_values('synergistic_pct', ascending=False)
    
    print(f"\n{'Pathway':<40} {'N Rxns':<10} {'Synerg%':<12} {'Antag%':<12} {'Class':<20}")
    print("-"*100)
    for _, row in synergy.iterrows():
        print(f"{row['subsystem'][:40]:<40} {row['n_reactions']:<10} {row['synergistic_pct']:<12.1f} {row['antagonistic_pct']:<12.1f} {row['pathway_class']:<20}")
    
    # Focus on transport pathways (LEC-dominated)
    transport = synergy[synergy['subsystem'].str.contains('Transport', na=False)]
    if len(transport) > 0:
        mean_synergy = transport['synergistic_pct'].mean()
        print(f"\n\nKey Insight:")
        print("-"*70)
        print(f"Transport pathways (LEC-dominated): {mean_synergy:.1f}% synergistic")
        print("LEC-dominated pathways show coordinated diet-microbiome responses")
    
    return synergy


# ============================================================================
# REACTION TRACING TABLES
# ============================================================================

def generate_rq1_hfd_traced_table(data):
    """
    RQ1-HFD Reactions Traced Through Full Pipeline (15 reactions)
    """
    print("\n" + "="*70)
    print("RQ1-HFD Reactions Traced Through Full Pipeline (15 reactions)")
    print("="*70)
    
    rq1_hfd = data['rq1_edges_hfd'].copy()
    rq1_hfd_rxns = rq1_hfd['ReactionID'].unique()
    
    # Get flux data
    rq1_flux = data['rq1_flux']
    rq4_attrib = data['rq4_attrib']
    rq3_attr = data['rq3_attribution']
    
    results = []
    for rxn in rq1_hfd_rxns:
        # RQ1 info
        rq1_info = rq1_flux[rq1_flux['ReactionID'] == rxn]
        rxn_name = rq1_info['ReactionName'].iloc[0] if len(rq1_info) > 0 else rxn
        subsystem = rq1_info['Subsystem'].iloc[0] if len(rq1_info) > 0 else 'Unknown'
        diff_hfd = rq1_info['Diff(HFD-SCD)'].iloc[0] if len(rq1_info) > 0 else 0
        
        # RQ2: Count strains conserved
        strain_count = sum(1 for s in STRAINS if rxn in set(data['rq2_edges'][s]['ReactionID'].unique()))
        
        # RQ3: Primary driver cell type
        rq3_info = rq3_attr[rq3_attr['reaction_id'] == rxn]
        primary_driver = rq3_info['primary_driver'].iloc[0] if len(rq3_info) > 0 else 'Unknown'
        
        # RQ4: Diet vs microbiome
        rq4_info = rq4_attrib[rq4_attrib['reaction_id'] == rxn]
        if len(rq4_info) > 0:
            dominant = rq4_info['dominant_driver'].iloc[0]
            diet_var = rq4_info['variance_explained_diet'].iloc[0]
            micro_var = rq4_info['variance_explained_microbiome'].iloc[0]
        else:
            dominant = 'Unknown'
            diet_var = 0
            micro_var = 0
        
        results.append({
            'ReactionID': rxn,
            'ReactionName': rxn_name[:30],
            'Subsystem': subsystem[:25],
            'Diff(HFD-SCD)': diff_hfd,
            'Strains(RQ2)': f"{strain_count}/9",
            'CellDriver(RQ3)': primary_driver[:15] if primary_driver else 'Unknown',
            'DietVar%(RQ4)': diet_var,
            'MicroVar%(RQ4)': micro_var,
            'Dominant(RQ4)': dominant[:20]
        })
    
    df = pd.DataFrame(results)
    
    print(f"\n{'RxnID':<12} {'Name':<30} {'Subsystem':<25} {'Δ(HFD-SCD)':<12} {'Strains':<8} {'CellDriver':<15} {'DietVar%':<10} {'MicroVar%':<10} {'Dominant':<20}")
    print("-"*160)
    for _, row in df.iterrows():
        print(f"{row['ReactionID']:<12} {row['ReactionName']:<30} {row['Subsystem']:<25} {row['Diff(HFD-SCD)']:<12.3f} {row['Strains(RQ2)']:<8} {row['CellDriver(RQ3)']:<15} {row['DietVar%(RQ4)']:<10.1f} {row['MicroVar%(RQ4)']:<10.1f} {row['Dominant(RQ4)']:<20}")
    
    return df


def generate_microbiome_dominated_table(data):
    """
    Microbiome-Dominated Reactions (72 total)
    Top 20 by microbiome variance explained
    """
    print("\n" + "="*70)
    print("Microbiome-Dominated Reactions (72 total)")
    print("Top 20 by microbiome variance explained:")
    print("="*70)
    
    rq4_attrib = data['rq4_attrib'].copy()
    
    # Filter to microbiome-dominated
    micro_dom = rq4_attrib[rq4_attrib['dominant_driver'].str.contains('Microbiome', na=False)].copy()
    micro_dom = micro_dom.sort_values('variance_explained_microbiome', ascending=False).head(20)
    
    # Add reaction names
    annotations = data['rq4_annotations'][['reaction_id', 'reaction_name', 'subsystem']].drop_duplicates()
    micro_dom = micro_dom.merge(annotations, on='reaction_id', how='left')
    
    print(f"\n{'Rank':<5} {'ReactionID':<15} {'ReactionName':<35} {'Subsystem':<25} {'MicroVar%':<12} {'DietVar%':<12} {'Driver':<25}")
    print("-"*140)
    for i, (_, row) in enumerate(micro_dom.iterrows(), 1):
        rxn_name = str(row.get('reaction_name', row['reaction_id']))[:35]
        subsys = str(row.get('subsystem', 'Unknown'))[:25]
        print(f"{i:<5} {row['reaction_id']:<15} {rxn_name:<35} {subsys:<25} {row['variance_explained_microbiome']:<12.1f} {row['variance_explained_diet']:<12.1f} {row['dominant_driver']:<25}")
    
    print(f"\nTotal microbiome-dominated reactions: {len(rq4_attrib[rq4_attrib['dominant_driver'].str.contains('Microbiome', na=False)])}")
    
    return micro_dom


def generate_diet_dominated_table(data):
    """
    Diet-Dominated Reactions (308 total)
    Top reactions (100% diet variance)
    """
    print("\n" + "="*70)
    print("Diet-Dominated Reactions (308 total)")
    print("Top reactions (100% diet variance):")
    print("="*70)
    
    rq4_attrib = data['rq4_attrib'].copy()
    
    # Filter to diet-dominated
    diet_dom = rq4_attrib[rq4_attrib['dominant_driver'].str.contains('Diet', na=False)].copy()
    
    # Sort by diet variance and get top 20
    diet_dom = diet_dom.sort_values('variance_explained_diet', ascending=False).head(20)
    
    # Add reaction names
    annotations = data['rq4_annotations'][['reaction_id', 'reaction_name', 'subsystem']].drop_duplicates()
    diet_dom = diet_dom.merge(annotations, on='reaction_id', how='left')
    
    print(f"\n{'Rank':<5} {'ReactionID':<15} {'ReactionName':<35} {'Subsystem':<25} {'DietVar%':<12} {'MicroVar%':<12} {'Driver':<25}")
    print("-"*140)
    for i, (_, row) in enumerate(diet_dom.iterrows(), 1):
        rxn_name = str(row.get('reaction_name', row['reaction_id']))[:35]
        subsys = str(row.get('subsystem', 'Unknown'))[:25]
        print(f"{i:<5} {row['reaction_id']:<15} {rxn_name:<35} {subsys:<25} {row['variance_explained_diet']:<12.1f} {row['variance_explained_microbiome']:<12.1f} {row['dominant_driver']:<25}")
    
    total_diet = len(rq4_attrib[rq4_attrib['dominant_driver'].str.contains('Diet', na=False)])
    print(f"\nTotal diet-dominated reactions: {total_diet}")
    
    return diet_dom


def generate_antagonistic_table(data):
    """
    Antagonistic Reactions (75 total)
    Top by effect magnitude
    """
    print("\n" + "="*70)
    print("Antagonistic Reactions (75 total)")
    print("Top by effect magnitude:")
    print("="*70)
    
    rq4_attrib = data['rq4_attrib'].copy()
    
    # Filter to antagonistic
    antagonistic = rq4_attrib[rq4_attrib['opposing_effects'] == True].copy()
    
    # Calculate total effect magnitude
    antagonistic['effect_magnitude'] = antagonistic['abs_diet_contribution'].abs() + antagonistic['abs_microbiome_contribution'].abs()
    antagonistic = antagonistic.sort_values('effect_magnitude', ascending=False).head(20)
    
    # Add reaction names
    annotations = data['rq4_annotations'][['reaction_id', 'reaction_name', 'subsystem']].drop_duplicates()
    antagonistic = antagonistic.merge(annotations, on='reaction_id', how='left')
    
    print(f"\n{'Rank':<5} {'ReactionID':<15} {'ReactionName':<30} {'Δ_Diet':<12} {'Δ_Micro':<12} {'Interaction':<25} {'Driver':<25}")
    print("-"*140)
    for i, (_, row) in enumerate(antagonistic.iterrows(), 1):
        rxn_name = str(row.get('reaction_name', row['reaction_id']))[:30]
        print(f"{i:<5} {row['reaction_id']:<15} {rxn_name:<30} {row['delta_diet']:<12.3f} {row['delta_microbiome']:<12.3f} {row['interaction']:<25} {row['dominant_driver']:<25}")
    
    print(f"\nTotal antagonistic reactions: {len(rq4_attrib[rq4_attrib['opposing_effects'] == True])}")
    
    return antagonistic


# ============================================================================
# SUMMARY STATISTICS TABLE
# ============================================================================

def generate_summary_table(data):
    """
    Generate comprehensive summary statistics table
    """
    print("\n" + "="*70)
    print("SUMMARY STATISTICS TABLE")
    print("="*70)
    
    # Collect all statistics
    stats = {}
    
    # RQ1
    rq1_hfd_rxns = set(data['rq1_edges_hfd']['ReactionID'].unique())
    rq1_kd_rxns = set(data['rq1_edges_kd']['ReactionID'].unique())
    rq1_wd_rxns = set(data['rq1_edges_wd']['ReactionID'].unique())
    
    stats['RQ1_Total_Reactions'] = len(data['rq1_flux'])
    stats['RQ1_HFD_Significant'] = len(rq1_hfd_rxns)
    stats['RQ1_KD_Significant'] = len(rq1_kd_rxns)
    stats['RQ1_WD_Significant'] = len(rq1_wd_rxns)
    stats['RQ1_Cohorts'] = 6
    stats['RQ1_Samples'] = 99
    
    # RQ2
    stats['RQ2_Strains'] = 9
    universal = sum(1 for rxn in rq1_hfd_rxns 
                   if sum(1 for s in STRAINS if rxn in set(data['rq2_edges'][s]['ReactionID'].unique())) == 9)
    conserved_5plus = sum(1 for rxn in rq1_hfd_rxns 
                         if sum(1 for s in STRAINS if rxn in set(data['rq2_edges'][s]['ReactionID'].unique())) >= 5)
    stats['RQ2_Universal_Rxns'] = universal
    stats['RQ2_Conserved_5plus'] = conserved_5plus
    
    # RQ3
    contrib = data['rq3_contrib']
    stats['RQ3_Cell_Types'] = len(contrib)
    stats['RQ3_LEC_Contribution'] = contrib[contrib['cell_type'] == 'LECs']['contribution_percent'].values[0]
    stats['RQ3_Hepatocyte_Contribution'] = contrib[contrib['cell_type'] == 'Hepatocytes']['contribution_percent'].values[0]
    stats['RQ3_LEC_CV'] = data['rq3_conservation'][data['rq3_conservation']['cell_type'] == 'LECs']['cv'].values[0] * 100
    
    # RQ4
    rq4_attrib = data['rq4_attrib']
    stats['RQ4_Diet_Dominated'] = rq4_attrib['dominant_driver'].str.contains('Diet', na=False).sum()
    stats['RQ4_Microbiome_Dominated'] = rq4_attrib['dominant_driver'].str.contains('Microbiome', na=False).sum()
    stats['RQ4_Antagonistic'] = rq4_attrib['opposing_effects'].sum()
    stats['RQ4_Stable'] = rq4_attrib['dominant_driver'].str.contains('Stable', na=False).sum()
    stats['RQ4_Species'] = len(data['rq4_species'])
    stats['RQ4_Mean_Diet_Var'] = rq4_attrib['variance_explained_diet'].mean()
    stats['RQ4_Mean_Micro_Var'] = rq4_attrib['variance_explained_microbiome'].mean()
    
    # Synergy
    transport = data['rq4_synergy'][data['rq4_synergy']['subsystem'].str.contains('Transport', na=False)]
    stats['RQ4_Transport_Synergy'] = transport['synergistic_pct'].mean() if len(transport) > 0 else 0
    
    # Print summary table
    print(f"\n{'Cross-RQ Connection':<45} {'Key Statistic':<30} {'Verdict':<25}")
    print("-"*100)
    print(f"{'RQ1→RQ2':<45} {f'{conserved_5plus} ({conserved_5plus/len(rq1_hfd_rxns)*100:.1f}%) conserved in ≥5 strains':<30} {'Partially upheld':<25}")
    
    rq1_rq3_str = (
        f"LECs={stats['RQ3_LEC_Contribution']:.1f}%, "
        f"Hep={stats['RQ3_Hepatocyte_Contribution']:.1f}%"
    )
    
    rq2_rq3_str = f"LEC CV = {stats['RQ3_LEC_CV']:.1f}% across strains"
    rq3_rq4_str = f"{stats['RQ4_Transport_Synergy']:.1f}% synergy in transport"


    print(f"{'RQ1→RQ3':<45} {rq1_rq3_str:<30} {'Upheld with surprise':<25}")
    print(f"{'RQ2→RQ3':<45} {rq2_rq3_str:<30} {'Strongly upheld':<25}")
    print(f"{'RQ1→RQ4':<45} {f'7 diet-driven, 4 microbiome-driven':<30} {'Nuanced':<25}")
    print(f"{'RQ3→RQ4':<45} {rq3_rq4_str:<30} {'Upheld':<25}")

    
    
    
    
    #print(f"{'RQ1→RQ3':<45} {f'LECs={stats['RQ3_LEC_Contribution']:.1f}%, Hep={stats['RQ3_Hepatocyte_Contribution']:.1f}%':<30} {'Upheld with surprise':<25}")
    #print(f"{'RQ2→RQ3':<45} {f'LEC CV = {stats['RQ3_LEC_CV']:.1f}% across strains':<30} {'Strongly upheld':<25}")
    #print(f"{'RQ1→RQ4':<45} {f'7 diet-driven, 4 microbiome-driven':<30} {'Nuanced':<25}")
    #print(f"{'RQ3→RQ4':<45} {f'{stats['RQ4_Transport_Synergy']:.1f}% synergy in transport':<30} {'Upheld':<25}")
    
    print(f"\n\nDetailed Statistics:")
    print("-"*70)
    for key, value in stats.items():
        if isinstance(value, float):
            print(f"  {key:<40} {value:.2f}")
        else:
            print(f"  {key:<40} {value}")
    
    return stats


# ============================================================================
# MICROBIOME COMMUNITY TABLES
# ============================================================================

def generate_microbiome_tables(data):
    """
    Generate microbiome community composition tables
    """
    print("\n" + "="*70)
    print("MICROBIOME COMMUNITY TABLES")
    print("="*70)
    
    species = data['rq4_species'].copy()
    portal = data['rq4_portal'].copy()
    
    # Table 1: Dominant Species Under SCD
    print("\n" + "-"*50)
    print("Dominant Species Under Standard Chow Diet (SCD)")
    print("-"*50)
    scd_species = species.nlargest(10, 'ND_SCD_abundance')
    print(f"\n{'Rank':<5} {'Species':<40} {'Abundance':<15}")
    print("-"*60)
    for i, (_, row) in enumerate(scd_species.iterrows(), 1):
        print(f"{i:<5} {row['species']:<40} {row['ND_SCD_abundance']:<15.4f}")
    
    # Table 2: Dominant Species Under HFD
    print("\n" + "-"*50)
    print("Dominant Species Under High-Fat Diet (HFD)")
    print("-"*50)
    hfd_species = species.nlargest(10, 'DD_HFD_abundance')
    print(f"\n{'Rank':<5} {'Species':<40} {'Abundance':<15}")
    print("-"*60)
    for i, (_, row) in enumerate(hfd_species.iterrows(), 1):
        print(f"{i:<5} {row['species']:<40} {row['DD_HFD_abundance']:<15.4f}")
    
    # Table 3: Portal Metabolites
    print("\n" + "-"*50)
    print("Portal Metabolites (What Reaches the Liver)")
    print("-"*50)
    print(f"\n{'Condition':<12} {'Metabolite':<25} {'Metabolite_ID':<15} {'Flux':<15} {'Importance':<15}")
    print("-"*80)
    for _, row in portal.iterrows():
        print(f"{row['Condition']:<12} {row['Metabolite']:<25} {row['Metabolite_ID']:<15} {row['Flux']:<15.4f} {row['Importance']:<15}")
    
    # Table 4: Functional Guild Shifts
    print("\n" + "-"*50)
    print("Functional Guild Shifts (SCD → HFD)")
    print("-"*50)
    
    # Calculate shifts
    species['abundance_change'] = species['DD_HFD_abundance'] - species['ND_SCD_abundance']
    species['fold_change'] = species['DD_HFD_abundance'] / (species['ND_SCD_abundance'] + 1e-10)
    
    # Top increasers
    print("\nTop 5 Increasers (SCD → HFD):")
    increasers = species.nlargest(5, 'abundance_change')
    print(f"{'Species':<40} {'SCD':<12} {'HFD':<12} {'Change':<12} {'Fold':<10}")
    print("-"*90)
    for _, row in increasers.iterrows():
        print(f"{row['species']:<40} {row['ND_SCD_abundance']:<12.4f} {row['DD_HFD_abundance']:<12.4f} {row['abundance_change']:<12.4f} {row['fold_change']:<10.2f}")
    
    # Top decreasers
    print("\nTop 5 Decreasers (SCD → HFD):")
    decreasers = species.nsmallest(5, 'abundance_change')
    print(f"{'Species':<40} {'SCD':<12} {'HFD':<12} {'Change':<12} {'Fold':<10}")
    print("-"*90)
    for _, row in decreasers.iterrows():
        print(f"{row['species']:<40} {row['ND_SCD_abundance']:<12.4f} {row['DD_HFD_abundance']:<12.4f} {row['abundance_change']:<12.4f} {row['fold_change']:<10.2f}")
    
    return species, portal


# ============================================================================
# CSV EXPORT
# ============================================================================

def export_all_tables(data):
    """Export all tables to CSV files"""
    print("\n" + "="*70)
    print("EXPORTING TABLES TO CSV")
    print("="*70)
    
    # Generate all tables
    rq1_rq2 = generate_rq1_to_rq2_table(data)
    rq1_rq3 = generate_rq1_to_rq3_table(data)
    rq2_rq3 = generate_rq2_to_rq3_table(data)
    rq1_rq4 = generate_rq1_to_rq4_table(data)
    rq3_rq4 = generate_rq3_to_rq4_table(data)
    
    traced = generate_rq1_hfd_traced_table(data)
    micro_dom = generate_microbiome_dominated_table(data)
    diet_dom = generate_diet_dominated_table(data)
    antagonistic = generate_antagonistic_table(data)
    
    summary = generate_summary_table(data)
    species, portal = generate_microbiome_tables(data)
    
    # Export to CSV
    rq1_rq2.to_csv(OUTPUT_DIR / 'Table_RQ1_to_RQ2_Conservation.csv', index=False)
    rq1_rq3.to_csv(OUTPUT_DIR / 'Table_RQ1_to_RQ3_CellType.csv', index=False)
    rq2_rq3.to_csv(OUTPUT_DIR / 'Table_RQ2_to_RQ3_Conservation.csv', index=False)
    rq1_rq4.to_csv(OUTPUT_DIR / 'Table_RQ1_to_RQ4_Attribution.csv', index=False)
    rq3_rq4.to_csv(OUTPUT_DIR / 'Table_RQ3_to_RQ4_Synergy.csv', index=False)
    
    traced.to_csv(OUTPUT_DIR / 'Table_RQ1_HFD_Full_Trace.csv', index=False)
    micro_dom.to_csv(OUTPUT_DIR / 'Table_Microbiome_Dominated_Top20.csv', index=False)
    diet_dom.to_csv(OUTPUT_DIR / 'Table_Diet_Dominated_Top20.csv', index=False)
    antagonistic.to_csv(OUTPUT_DIR / 'Table_Antagonistic_Top20.csv', index=False)
    
    pd.DataFrame([summary]).to_csv(OUTPUT_DIR / 'Table_Summary_Statistics.csv', index=False)
    species.to_csv(OUTPUT_DIR / 'Table_Microbiome_Species_Shifts.csv', index=False)
    portal.to_csv(OUTPUT_DIR / 'Table_Portal_Metabolites.csv', index=False)
    
    print(f"\nExported all tables to {OUTPUT_DIR}")
    print("Files created:")
    for f in sorted(OUTPUT_DIR.glob('Table_*.csv')):
        print(f"  - {f.name}")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function"""
    print("="*70)
    print("CROSS-RQ INTEGRATION ANALYSIS")
    print("="*70)
    
    # Load data
    data = load_all_data()
    
    # Export all tables
    export_all_tables(data)
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)


if __name__ == "__main__":
    main()

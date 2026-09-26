"""
generate_Fig4_5_RQ3.py
======================
Reproduces Figures 4A–4D and 5A–5D for Section 2.3 and 2.4
(Cellular Heterogeneity Redistributes Metabolic Responses & Abundance-weighted Validation).

Required input files:
    RQ3_Section23_Supplementary_Tables.xlsx
    Section24_Validation_Supplementary_Tables.xlsx (or SJ2 CSV)

Dependencies:
    pip install pandas numpy matplotlib openpyxl scipy

Usage:
    python generate_Fig4_5_RQ3.py
    python generate_Fig4_5_RQ3.py --supp_xlsx /path/to/rq3.xlsx --val_file /path/to/sec24_SJ2.csv
    python generate_Fig4_5_RQ3.py --panel  # Generates composite panels (Fig 4 and Fig 5)
"""

import argparse
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy.cluster import hierarchy

warnings.filterwarnings('ignore')

# ============================================================
# CONFIG
# ============================================================
DEFAULTS = {
    'supp_xlsx': 'RQ3_Section23_Supplementary_Tables.xlsx',
    'val_file':  'Section24_Validation_Supplementary_Tables.xlsx',
    'out_dir':   '.',
}

STRAIN_COLS = [
    '129S1', 'A/J', 'C57BL/6J', 'CAST/EiJ', 'DBA/2J',
    'NOD/ShiLtJ', 'NZO/HlLtJ', 'PWK/PhJ', 'WSB/EiJ',
]

CAT_COLORS = {
    'Sinusoidal/Endothelial': '#1565C0',
    'Stellate/Mesenchymal':   '#7B2D8B',
    'Kupffer/Macrophage':     '#BF360C',
    'Monocyte/DC':            '#E65100',
    'Granulocyte/Neutrophil': '#2E7D32',
    'Cholangiocyte':          '#00838F',
    'Lymphocyte':             '#AD1457',
    'Parenchymal':            '#F57F17',
    'Other':                  '#78909C',
}

DRIVER_COLORS = {
    'qHSCs':                 '#7B2D8B',
    'LECs':                  '#1565C0',
    'Neutrophiles':          '#2E7D32',
    'Transitioning Mo':      '#E65100',
    'cAMP qHSCs':            '#9C27B0',
    'Cholangiocytes':        '#00838F',
    'Cycling':               '#78909C',
    'NONE':                  '#BDBDBD',
    'Timd4+ resKC':          '#BF360C',
    'Cd207-, Trem2+ Mo-KC':  '#FF7043',
    'Hepatocytes':           '#F57F17',
}

plt.rcParams.update({
    'font.family':       'sans-serif',
    'font.sans-serif':   ['DejaVu Sans'],
    'font.size':         9,
    'axes.titlesize':    10,
    'axes.labelsize':    9,
    'xtick.labelsize':   8,
    'ytick.labelsize':   8,
    'legend.fontsize':   8,
    'figure.dpi':        150,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'axes.linewidth':    0.8,
})

# ============================================================
# Data loading helpers
# ============================================================
def load_sa(xl: pd.ExcelFile) -> pd.DataFrame:
    sa = xl.parse('SA_Contribution_Summary', header=1).dropna(subset=['Cell Type'])
    # Column was renamed 'N Sig. Reactions' -> 'N Threshold-\nResponsive Reactions' in the
    # 2026-09-25 supplementary-tables refresh (deterministic threshold criteria, not inferential
    # significance testing) -- support both so this script works against either workbook.
    n_col = 'N Threshold-\nResponsive Reactions' if 'N Threshold-\nResponsive Reactions' in sa.columns else 'N Sig. Reactions'
    for col in ['Contribution (%)', 'Abundance (%)', 'Response Rate',
                'Per-Cell\nContribution Index', n_col]:
        sa[col] = pd.to_numeric(sa[col], errors='coerce')
    return sa.dropna(subset=['Contribution (%)'])

def load_sb(xl: pd.ExcelFile) -> pd.DataFrame:
    sb = xl.parse('SB_CrossStrain_Matrix', header=1).dropna(subset=['Cell Type'])
    for col in STRAIN_COLS + ['Mean (%)', 'Std (%)']:
        sb[col] = pd.to_numeric(sb[col], errors='coerce')
    return sb

def load_sf(xl: pd.ExcelFile) -> pd.DataFrame:
    return xl.parse('SF_Pathway_Attribution', header=1).dropna(subset=['Pathway'])

def _cat_color(cat: str) -> str:
    return CAT_COLORS.get(str(cat), '#888888')

def _row_bands(ax, n: int) -> None:
    for i in range(0, n, 2):
        ax.axhspan(i - 0.45, i + 0.45, alpha=0.03, color='gray', zorder=0)

# ============================================================
# Core Plotting Functions (Fig 4)
# ============================================================

def plot_fig4a(ax: plt.Axes, sa: pd.DataFrame, title: str) -> None:
    df = sa.sort_values('Contribution (%)', ascending=False).reset_index(drop=True)
    colors = [_cat_color(c) for c in df['Category']]
    y = np.arange(len(df))

    ax.barh(y, df['Contribution (%)'], height=0.65, color=colors, alpha=0.88, edgecolor='white', linewidth=0.5)

    for i, abund in enumerate(df['Abundance (%)']):
        ax.scatter(abund, i, marker='|', s=25, color='#333333', zorder=4, linewidths=1)

    for i, row in df.iterrows():
        v = row['Contribution (%)']
        abund = row['Abundance (%)']
        text_x = v + 0.3
        if (v - 0.5) <= abund <= (v + 3.5):
            text_x = abund + 0.8
        ax.text(text_x, i, f'{v:.1f}%', va='center', ha='left', fontsize=7.5, color='#222222')

    ax.set_yticks(y)
    ax.set_yticklabels(df['Cell Type'], fontsize=8)
    ax.set_xlabel('Abundance-weighted contribution to modeled flux changes (%)')
    ax.set_title(title, fontsize=10.5, pad=9, fontweight='bold')
    ax.set_xlim(0, df['Contribution (%)'].max() * 1.18)
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, alpha=0.25, linestyle='--', linewidth=0.5)
    _row_bands(ax, len(df))

    marker_leg = Line2D([0], [0], marker='|', color='#333333', lw=0, markersize=5, markeredgewidth=1, label='Cell abundance (%)')
    cat_patches = [mpatches.Patch(color=CAT_COLORS[c], label=c, alpha=0.88) for c in CAT_COLORS if c in df['Category'].values]
    ax.legend(handles=cat_patches + [marker_leg], frameon=False, fontsize=7.2, loc='center right', ncol=1)

def plot_fig4b(ax: plt.Axes, sa: pd.DataFrame, title: str) -> None:
    df = sa.dropna(subset=['Abundance (%)'])
    for _, row in df.iterrows():
        c = _cat_color(row['Category'])
        ax.scatter(row['Abundance (%)'], row['Contribution (%)'], c=c, s=60, alpha=0.85, edgecolors='white', linewidths=0.5, zorder=3)

    maxval = max(df['Abundance (%)'].max(), df['Contribution (%)'].max()) * 1.08
    ax.plot([0, maxval], [0, maxval], '--', color='#aaaaaa', lw=1.2, zorder=1, label='Proportional contribution (y = x)')
    ax.fill_between([0, maxval], [0, maxval], [maxval, maxval], alpha=0.04, color='#1565C0')
    ax.fill_between([0, maxval], [0, 0], [0, maxval], alpha=0.04, color='#E65100')

    MAIN_LABELS = {'LECs': (0.5, 2.0)}
    for _, row in df.iterrows():
        if row['Cell Type'] in MAIN_LABELS:
            dx, dy = MAIN_LABELS[row['Cell Type']]
            ax.annotate(row['Cell Type'], xy=(row['Abundance (%)'], row['Contribution (%)']), 
                        xytext=(row['Abundance (%)']+dx, row['Contribution (%)']+dy), 
                        fontsize=7.2, color='#222222', arrowprops=dict(arrowstyle='-', color='#aaaaaa', lw=0.6))

    axins = ax.inset_axes([0.48, 0.08, 0.48, 0.48])
    for _, row in df.iterrows():
        c = _cat_color(row['Category'])
        axins.scatter(row['Abundance (%)'], row['Contribution (%)'], c=c, s=45, alpha=0.85, edgecolors='white', linewidths=0.5, zorder=3)
        
    axins.plot([0, 12], [0, 12], '--', color='#aaaaaa', lw=1.0, zorder=1)
    axins.fill_between([0, 12], [0, 12], [12, 12], alpha=0.04, color='#1565C0')
    axins.fill_between([0, 12], [0, 0], [12, 12], alpha=0.04, color='#E65100')
    axins.set_xlim(-0.2, 10.5)
    axins.set_ylim(-0.2, 9.0)
    axins.tick_params(labelsize=7)
    axins.grid(True, alpha=0.25, linestyle='--', linewidth=0.5)
    
    INSET_LABELS = {'qHSCs': (0.2, 0.5), 'aHSCs': (-1.5, 0.4), 'Hepatocytes': (0.3, -0.8), 'Neutrophiles': (0.2, 0.4), 'Cholangiocytes': (0.2, -0.8), 'cAMP qHSCs': (0.2, -0.4)}
    for _, row in df.iterrows():
        if row['Cell Type'] in INSET_LABELS:
            dx, dy = INSET_LABELS[row['Cell Type']]
            axins.annotate(row['Cell Type'], xy=(row['Abundance (%)'], row['Contribution (%)']), 
                           xytext=(row['Abundance (%)']+dx, row['Contribution (%)']+dy), 
                           fontsize=7.0, color='#222222', arrowprops=dict(arrowstyle='-', color='#aaaaaa', lw=0.6))

    ax.indicate_inset_zoom(axins, edgecolor="gray", alpha=0.6)
    ax.set_xlabel('Cell abundance (%)')
    ax.set_ylabel('Abundance-weighted contribution (%)')
    ax.set_title(title, fontsize=10.5, pad=9, fontweight='bold')
    ax.set_xlim(-0.5, maxval)
    ax.set_ylim(-1.0, max(df['Contribution (%)'].max() * 1.12, maxval))
    
    cat_patches = [mpatches.Patch(color=CAT_COLORS[c], label=c, alpha=0.88) for c in CAT_COLORS if c in df['Category'].values]
    diag_handle = Line2D([0], [0], linestyle='--', color='#aaaaaa', lw=1.2, label='Proportional contribution (y = x)')
    ax.legend(handles=cat_patches + [diag_handle], frameon=False, fontsize=7.2, loc='upper left', bbox_to_anchor=(0.0, 1.0))

def plot_fig4c(ax: plt.Axes, sa: pd.DataFrame, title: str) -> None:
    df = sa.dropna(subset=['Response Rate']).sort_values('Response Rate', ascending=True).reset_index(drop=True)
    colors = [_cat_color(c) for c in df['Category']]
    y = np.arange(len(df))

    ax.barh(y, df['Response Rate'] * 100, height=0.65, color=colors, alpha=0.88, edgecolor='white', linewidth=0.5)

    hep = df[df['Cell Type'] == 'Hepatocytes']
    if len(hep):
        hep_rate = hep.iloc[0]['Response Rate'] * 100
        ax.axvline(hep_rate, color='#F57F17', lw=1.2, linestyle='--', alpha=0.7,
                   label=f'Hepatocyte response rate ({hep_rate:.1f}%)')

    for i, v in enumerate(df['Response Rate'] * 100):
        ax.text(v + 0.05, i, f'{v:.1f}%', va='center', ha='left', fontsize=7.3, color='#222222')

    ax.set_yticks(y)
    ax.set_yticklabels(df['Cell Type'], fontsize=8)
    ax.set_xlabel('Response rate (% of reactions classified as threshold-responsive)')
    ax.set_title(title, fontsize=10.5, pad=9, fontweight='bold')
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, alpha=0.25, linestyle='--', linewidth=0.5)
    _row_bands(ax, len(df))
    
    cat_patches = [mpatches.Patch(color=CAT_COLORS[c], label=c, alpha=0.88) for c in CAT_COLORS if c in df['Category'].values]
    hep_handle = Line2D([0], [0], linestyle='--', color='#F57F17', lw=1.2, alpha=0.7,
                        label=f'Hepatocyte response rate ({hep_rate:.1f}%)')
    ax.legend(handles=cat_patches + [hep_handle], frameon=False, fontsize=7.2, loc='lower right', ncol=1)

def plot_fig4d(ax: plt.Axes, sa: pd.DataFrame, title: str) -> None:
    df = sa.dropna(subset=['Per-Cell\nContribution Index']).sort_values('Per-Cell\nContribution Index', ascending=True).reset_index(drop=True)
    colors = [_cat_color(c) for c in df['Category']]
    y = np.arange(len(df))

    ax.barh(y, df['Per-Cell\nContribution Index'], height=0.65, color=colors, alpha=0.88, edgecolor='white', linewidth=0.5)
    ax.axvline(1.0, color='#555555', lw=1.2, linestyle='--', label='Proportional (index = 1.0)')

    for i, v in enumerate(df['Per-Cell\nContribution Index']):
        ax.text(v + 0.02, i, f'{v:.2f}', va='center', ha='left', fontsize=7.3, color='#222222')

    ax.set_yticks(y)
    ax.set_yticklabels(df['Cell Type'], fontsize=8)
    ax.set_xlabel('Per-cell contribution index (contribution% ÷ abundance%)')
    ax.set_title(title, fontsize=10.5, pad=9, fontweight='bold')
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, alpha=0.25, linestyle='--', linewidth=0.5)
    _row_bands(ax, len(df))
    ax.legend(frameon=False, fontsize=8, loc='lower right')

# ============================================================
# Core Plotting Functions (Fig 5)
# ============================================================

def plot_fig5a(ax_left: plt.Axes, ax_right: plt.Axes, sa: pd.DataFrame, title: str) -> None:
    cat_contrib = sa.groupby('Category')['Contribution (%)'].sum().sort_values(ascending=False)
    cat_abund   = sa.groupby('Category')['Abundance (%)'].sum().reindex(cat_contrib.index)

    for ax, values, subtitle in [(ax_left, cat_contrib, 'Abundance-weighted\ncontribution (%)'), (ax_right, cat_abund, 'Cell abundance (%)')]:
        wedge_colors = [CAT_COLORS.get(c, '#888888') for c in values.index]
        _, _, autotexts = ax.pie(values.values, labels=None, colors=wedge_colors, autopct=lambda p: f'{p:.1f}%' if p > 2 else '', startangle=90, pctdistance=0.75, wedgeprops=dict(width=0.55, edgecolor='white', linewidth=1.5))
        for at in autotexts: at.set_fontsize(7.5)
        ax.add_patch(plt.Circle((0, 0), 0.4, fc='white'))
        ax.set_title(subtitle, fontsize=10, pad=6)

    cat_patches = [mpatches.Patch(color=CAT_COLORS.get(c, '#888888'), label=c, alpha=0.88) for c in cat_contrib.index]
    ax_left.legend(handles=cat_patches, frameon=False, fontsize=8, loc='upper center', bbox_to_anchor=(1.0, -0.05), ncol=3)
    ax_left.text(1.0, 1.15, title, ha='center', va='bottom', transform=ax_left.transAxes, fontsize=11, fontweight='bold')

def plot_fig5b(ax: plt.Axes, sb: pd.DataFrame, title: str) -> None:
    """Plots the heatmap with 2D hierarchical clustering."""
    sb_plot = sb.dropna(subset=['Mean (%)']).sort_values('Mean (%)', ascending=False).head(15).reset_index(drop=True)
    matrix = sb_plot[STRAIN_COLS].values.astype(float)
    cat_map = dict(zip(sb_plot['Cell Type'], sb_plot['Category']))

    # Perform 2D Hierarchical Clustering
    row_linkage = hierarchy.linkage(matrix, method='ward')
    col_linkage = hierarchy.linkage(matrix.T, method='ward')

    row_order = hierarchy.leaves_list(row_linkage)
    col_order = hierarchy.leaves_list(col_linkage)

    # Reorder the matrix and labels
    matrix_clustered = matrix[row_order, :][:, col_order]
    ordered_cells = [sb_plot['Cell Type'].iloc[i] for i in row_order]
    ordered_strains = [STRAIN_COLS[j] for j in col_order]
    ordered_means = [sb_plot['Mean (%)'].iloc[i] for i in row_order]
    # Population SD (ddof=0) computed directly from the same per-strain values shown in the heatmap,
    # to match the manuscript's CV convention (verified against RQ3_conservation_analysis.csv) rather
    # than the sheet's own 'Std (%)' column, which is a sample SD (ddof=1) and doesn't reconcile.
    raw_stds = matrix.std(axis=1, ddof=0)
    ordered_stds = [raw_stds[i] for i in row_order]

    im = ax.imshow(matrix_clustered, cmap='Blues', aspect='auto', vmin=0, vmax=float(np.nanmax(matrix_clustered)))
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('Contribution (%)', fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    for i in range(matrix_clustered.shape[0]):
        for j in range(matrix_clustered.shape[1]):
            v = matrix_clustered[i, j]
            tc = 'white' if v > np.nanmax(matrix_clustered) * 0.6 else '#222222'
            ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=7.5, color=tc, fontweight='bold')

    ax.set_xticks(range(len(ordered_strains)))
    ax.set_xticklabels(ordered_strains, rotation=30, ha='right', fontsize=8.5)
    ax.set_yticks(range(len(ordered_cells)))
    ax.set_yticklabels(ordered_cells, fontsize=8.5)

    for tick in ax.get_yticklabels():
        cat = cat_map.get(tick.get_text(), '')
        tick.set_color(CAT_COLORS.get(cat, '#333333'))

    for i in range(len(ordered_cells)):
        ax.text(len(ordered_strains) + 0.1, i, f"{ordered_means[i]:.1f}±{ordered_stds[i]:.1f}", ha='left', va='center', fontsize=7, color='#555555')
    
    ax.text(len(ordered_strains) + 0.1, -0.75, 'Mean±SD', ha='left', va='center', fontsize=7, color='#555555', style='italic')
    ax.set_xlim(-0.5, len(ordered_strains) + 1.6)
    ax.set_title(title, fontsize=10.5, pad=10, fontweight='bold')

def plot_fig5c(ax: plt.Axes, sf: pd.DataFrame, title: str) -> None:
    sf_plot = sf[sf['N Reactions'] >= 3].sort_values('N Reactions', ascending=True).reset_index(drop=True)
    y = np.arange(len(sf_plot))
    STACK_SPECS = [('N Multicellular', 'Multi-cellular', '#1565C0'), ('N Cooperative', 'Cooperative', '#7B2D8B'), ('N Unique\n(1 cell type)', 'Cell-type unique', '#2E7D32'), ('N Non-Cellular', 'Non-cellular', '#BDBDBD')]

    left = np.zeros(len(sf_plot))
    for col, label, color in STACK_SPECS:
        vals = sf_plot[col].fillna(0).values.astype(float)
        ax.barh(y, vals, height=0.65, left=left, color=color, alpha=0.88, edgecolor='white', linewidth=0.5, label=label)
        left += vals

    for i, (_, row) in enumerate(sf_plot.iterrows()):
        driver = str(row['Primary Driver'])
        pct    = row['Primary Driver (%)']
        dc     = DRIVER_COLORS.get(driver, '#555555')
        ax.text(row['N Reactions'] + 0.2, i, f"{driver} ({pct:.0f}%)", va='center', ha='left', fontsize=7, color=dc)

    ax.set_yticks(y)
    labels_y = [p[:35] + '…' if len(p) > 37 else p for p in sf_plot['Pathway']]
    ax.set_yticklabels(labels_y, fontsize=8)
    ax.set_xlabel('Number of bulk-reference reactions with cell-type attribution')
    ax.set_title(title, fontsize=10.5, pad=9, fontweight='bold')
    ax.set_xlim(0, sf_plot['N Reactions'].max() * 1.9)
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, alpha=0.25, linestyle='--', linewidth=0.5)
    for i in range(0, len(sf_plot), 2):
        ax.axhspan(i - 0.45, i + 0.45, alpha=0.03, color='gray', zorder=0)
    ax.legend(frameon=False, fontsize=8, loc='lower right')

def plot_fig5d(ax: plt.Axes, sj2: pd.DataFrame, title: str) -> None:
    """Plots Baseline vs Delta Flux Concordance based on Section 2.4 data."""
    if sj2.empty:
        ax.text(0.5, 0.5, "Validation data missing.\nPlease pass --val_file", ha='center', va='center')
        ax.set_title(title, fontweight='bold')
        return
        
    cols_lower = sj2.columns.astype(str).str.lower()

    # Robust column matching. The 2026-09-25 canonical supplementary workbook's SJ2 sheet uses
    # explicit correlation-column names, e.g. 'r (Composite Chow\nvs Bulk SCD)' for baseline and
    # 'r (Composite Δ\nvs Bulk Δ)' for the differential/delta correlation, rather than a bare
    # 'baseline'/'delta' substring -- match on the 'r (' correlation-column prefix plus a
    # baseline-vs-delta keyword so this works against either workbook layout.
    is_corr_col = lambda cl: cl.strip().startswith('r (') or cl.strip().startswith('r(') or 'corr' in cl
    base_col = [c for c, cl in zip(sj2.columns, cols_lower)
                if is_corr_col(cl) and ('baseline' in cl or 'chow' in cl or 'scd' in cl)]
    delta_col = [c for c, cl in zip(sj2.columns, cols_lower)
                 if is_corr_col(cl) and ('delta' in cl or 'diff' in cl or 'Δ' in cl or 'δ' in cl)]
    sub_col = [c for c, cl in zip(sj2.columns, cols_lower) if 'subsystem' in cl or 'pathway' in cl]
    # Prefer an explicit "active reactions" count (the n actually used in each subsystem's
    # correlation) over a bare total-reactions-in-subsystem column, when both are present.
    n_active_col = [c for c, cl in zip(sj2.columns, cols_lower) if 'active' in cl and 'reaction' in cl]
    n_col = n_active_col or [c for c, cl in zip(sj2.columns, cols_lower) if 'n' in cl and ('reaction' in cl or 'count' in cl)]
    
    if not (base_col and delta_col and sub_col):
        ax.text(0.5, 0.5, f"Could not map columns.\nFound: {list(sj2.columns)}", ha='center', va='center', fontsize=8)
        ax.set_title(title, fontweight='bold')
        return

    bc, dc, sc = base_col[0], delta_col[0], sub_col[0]
    nc = n_col[0] if n_col else None
    
    # Draw quadrant guides
    ax.axhline(0, color='#999999', lw=1.0, zorder=1)
    ax.axvline(0, color='#999999', lw=1.0, zorder=1)
    
    # Highlight zones (fixed data-coordinate rectangles, independent of axis limits)
    ax.add_patch(mpatches.Rectangle((0.5, 0.525), 0.55, 0.525, alpha=0.05, color='#1565C0', zorder=0, linewidth=0))  # Strong both
    ax.add_patch(mpatches.Rectangle((0.5, -1.05), 0.55, 0.525, alpha=0.05, color='#E65100', zorder=0, linewidth=0))  # Divergent

    # Zone labels anchored at the LEFT edge of each band: all 5 hardcoded pathway callouts below
    # sit at x=0.96-1.00 (the right edge of these bands), so a left-anchored label doesn't compete
    # with them for space.
    ax.text(0.55, 1.0, 'Strong concordance\nin both states', fontsize=7.5, color='#1565C0', alpha=0.9,
            ha='left', va='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1.5), zorder=2)
    ax.text(0.55, -1.0, 'Baseline concordance only\n(Divergent ΔFlux)', fontsize=7.5, color='#E65100', alpha=0.9,
            ha='left', va='bottom', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1.5), zorder=2)

    # Plot points
    for _, row in sj2.iterrows():
        try:
            x, y = float(row[bc]), float(row[dc])
            sub = str(row[sc])
            s = float(row[nc]) * 1.5 if nc and not pd.isna(row[nc]) else 35
            
            # Color logic
            color = '#1565C0' if (x > 0.5 and y > 0.5) else ('#E65100' if (x > 0.5 and y < -0.5) else '#9E9E9E')
                
            ax.scatter(x, y, s=s, color=color, alpha=0.8, edgecolor='white', linewidth=0.5, zorder=3)
            
            # Annotate specific mentioned pathways. These 5 all sit within x=0.96-1.00 (the
            # extreme right edge of the plot), two of them (cholesterol/carnitine) essentially on
            # top of each other -- so labels use an ABSOLUTE external position in the right-hand
            # margin, stacked top-to-bottom in the same order as the points, rather than a small
            # offset relative to each point (which caused label-on-label and label-on-zone-text
            # collisions).
            targets = {
                'cholesterol': (1.28, 1.15),
                'carnitine': (1.28, 0.80),
                'urea': (1.28, -0.65),
                'sphingolipid': (1.28, -0.95),
                'lysine': (1.28, -1.25),
            }

            for t, label_pos in targets.items():
                if t in sub.lower():
                    ax.annotate(sub, (x, y), xytext=label_pos, ha='left', va='center',
                                fontsize=8, color='#222222', zorder=5, fontweight='bold',
                                annotation_clip=False,
                                bbox=dict(facecolor='white', alpha=0.85, edgecolor='none', pad=1),
                                arrowprops=dict(arrowstyle='-', color='#888888', lw=0.6, zorder=4,
                                                 shrinkA=2, shrinkB=2))
        except:
            continue

    ax.set_xlabel('Baseline Concordance (Pearson r)', fontsize=9)
    ax.set_ylabel('Diet-Induced ΔFlux Concordance (Pearson r)', fontsize=9)
    ax.set_title(title, fontsize=10.5, pad=10, fontweight='bold')
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_axisbelow(True)
    ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.5)
    ax.text(0.02, 0.02, 'Point size ∝ number of active reactions', transform=ax.transAxes,
            fontsize=7, color='#555555', ha='left', va='bottom', style='italic')


# ============================================================
# Panel Generation Functions
# ============================================================

def make_panel_fig4(sa: pd.DataFrame, out_path: str) -> None:
    fig = plt.figure(figsize=(16, 14))
    gs = gridspec.GridSpec(2, 2, hspace=0.3, wspace=0.25)

    plot_fig4a(fig.add_subplot(gs[0, 0]), sa, 'a | Cell-type contributions to modeled WD-versus-chow responses')
    plot_fig4b(fig.add_subplot(gs[0, 1]), sa, 'b | Contribution vs. abundance across populations')
    plot_fig4c(fig.add_subplot(gs[1, 0]), sa, 'c | Per-cell-type response rates')
    plot_fig4d(fig.add_subplot(gs[1, 1]), sa, 'd | Per-cell contribution index')

    plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  Saved Composite Panel: {out_path}')

def make_panel_fig5(sa: pd.DataFrame, sb: pd.DataFrame, sf: pd.DataFrame, sj2: pd.DataFrame, out_path: str) -> None:
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1.2], wspace=0.2, hspace=0.35)

    gs_5a = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0, 0], wspace=0.1)
    plot_fig5a(fig.add_subplot(gs_5a[0, 0]), fig.add_subplot(gs_5a[0, 1]), sa, 'a | Category-level modeled contribution versus estimated abundance')
    
    plot_fig5b(fig.add_subplot(gs[0, 1]), sb, 'b | Heuristic cross-strain projection of cell-type contributions (Clustered)')
    plot_fig5c(fig.add_subplot(gs[1, 0]), sf, 'c | Pathway-level cell-type attribution of bulk-reference reactions')
    plot_fig5d(fig.add_subplot(gs[1, 1]), sj2, 'd | Independent bulk comparison: baseline and differential concordance')

    plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  Saved Composite Panel: {out_path}')

# ============================================================
# Entry point
# ============================================================
def main() -> None:
    parser = argparse.ArgumentParser(description='Generate Figures 4A–4D and 5A–5D for RQ3.')
    parser.add_argument('--supp_xlsx', default=DEFAULTS['supp_xlsx'], help='Path to RQ3_Section23_Supplementary_Tables.xlsx')
    parser.add_argument('--val_file', default=DEFAULTS['val_file'], help='Path to Section24 CSV or XLSX')
    parser.add_argument('--out_dir', default=DEFAULTS['out_dir'], help='Directory for output PNG files')
    parser.add_argument('--panel', action='store_true', help='Generate composite panels')
    args = parser.parse_args()

    if not os.path.isfile(args.supp_xlsx):
        raise FileNotFoundError(f"Missing base file: {args.supp_xlsx}")

    os.makedirs(args.out_dir, exist_ok=True)
    def o(fname: str) -> str: return os.path.join(args.out_dir, fname)

    # Load 2.3 data
    xl = pd.ExcelFile(args.supp_xlsx)
    sa = load_sa(xl)
    sb = load_sb(xl)
    sf = load_sf(xl)

    # Load Section 2.4 Validation Data dynamically (CSV or XLSX)
    sj2 = pd.DataFrame()
    val_file = args.val_file
    
    if os.path.isfile(val_file):
        try:
            if val_file.endswith('.csv'):
                # Read CSV directly
                sj2 = pd.read_csv(val_file)
                # If headers are weird (e.g. metadata on row 1), shift down
                if not any(c.lower() in ['subsystem', 'pathway'] for c in sj2.columns.astype(str)):
                    sj2 = pd.read_csv(val_file, header=1)
            else:
                # Read from Excel
                val_xl = pd.ExcelFile(val_file)
                sj2_sheets = [s for s in val_xl.sheet_names if 'SJ2' in s or 'Subsystem' in s]
                if sj2_sheets:
                    sj2 = val_xl.parse(sj2_sheets[0])
                    if not any(c.lower() in ['subsystem', 'pathway'] for c in sj2.columns.astype(str)):
                        sj2 = val_xl.parse(sj2_sheets[0], header=1)
        except Exception as e:
            print(f"Error reading validation file {val_file}: {e}")
    else:
        print(f"Validation file not found: {val_file}")

    if args.panel:
        print('\n=== Generating Composite Panels (Fig 4 and Fig 5) ===')
        make_panel_fig4(sa, o('Fig4_Composite_Panel.png'))
        make_panel_fig5(sa, sb, sf, sj2, o('Fig5_Composite_Panel.png'))
    else:
        print('\n=== Generating Individual Figures ===')
        make_panel_fig5(sa, sb, sf, sj2, o('Fig5_Composite_Panel.png'))

    print(f'\nFinished writing to: {os.path.abspath(args.out_dir)}')

if __name__ == '__main__':
    main()
#!/usr/bin/env python3
"""
Figure 6 (hybrid redesign) — Gut microbiome community modeling and hepatic
flux attribution, RQ4.

Panels keep the CURRENT, corrected 27-species / 47-reaction scientific base
(matching the in-manuscript media/image6.png and the current Section 3.6
text). Presentation borrows from the old (superseded) Figure 6 only where
explicitly requested:

  (a) current 27-species log2FC composition-shift bars, now red/green
      directionally colored, + "27 of 57 retained" on-panel annotation.
      (31 decreased / 19 increased / 7 unchanged goes in the CAPTION, not
      the panel — per instruction.)
  (b) current acetate-vs-tradeoff line plot, unchanged data, + "positive
      exchange = microbial export" note + primary-tradeoff annotation.
  (c) REPLACED: old donut/pie style, but with CURRENT counts
      (260 / 115 / 27 / 2 / 3,322).
  (d) current subsystem/pathway panel, but as a stacked horizontal bar by
      dominant-driver category (data supports it — see verification notes),
      with Cholesterol Metabolism (6 rxns, 100% microbiome-dominant)
      highlighted.

Explicitly NOT restored: old portal-metabolite values/panel, old
synergy/antagonism classification (superseded).

Panel (a) caveat (disclosed to Roland): the exact script that rendered the
currently-embedded panel (a) could not be located on disk despite an
extensive search, so its log2FC pseudocount convention is not reproduced
pixel-for-pixel. The same 10 species and the same sign/direction for every
one of them were verified against the primary 27-species composition data;
only the exact bar magnitudes use a disclosed pseudocount (0.001) rather
than whatever the original convention was.
"""

import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

plt.rcParams.update({
    "font.family":        "Arial",
    "font.size":          9,
    "axes.labelsize":     10,
    "axes.titlesize":     10,
    "xtick.labelsize":    8,
    "ytick.labelsize":    8,
    "legend.fontsize":    8,
    "figure.dpi":         300,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "savefig.pad_inches": 0.1,
    "axes.linewidth":     0.8,
})

DATA_DIR = Path(os.environ.get("RQ4_RESULTS", "results/RQ4_final_full_v3"))  # [release] configurable
ANNOT_FILE = Path(os.environ.get("RQ4_ANNOTATIONS", "inputs/bulk/RQ4_reaction_annotations.csv"))  # [release] configurable
OUT_DIR  = Path(os.environ.get("FIG_OUT", "results/figures/figure6"))  # [release] configurable
OUT_DIR.mkdir(parents=True, exist_ok=True)

COL_UP   = "#D6604D"   # red   = increased under DD/HFD
COL_DOWN = "#2CA02C"   # green = decreased under DD/HFD
COL_ACETATE_ND = "#E69F00"
COL_ACETATE_DD = "#2166AC"


def _clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


PANEL_A_SPECIES = [
    "Lactobacillus_acidophilus_NCFM",
    "Enterobacter_aerogenes_KCTC_2190",
    "Streptococcus_pneumoniae_R6",
    "Staphylococcus_aureus_subsp_aureus_NCTC_8325",
    "Bacteroides_fragilis_YCH46",
    "Lactococcus_lactis_subsp_lactis_Il1403",
    "Enterococcus_faecalis_V583",
    "Streptococcus_agalactiae_2603V/R",
    "Enterococcus_faecium_DO",
    "Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586",
]

DISPLAY_NAME = {
    "Lactobacillus_acidophilus_NCFM": "Lactobacillus acidophilus NCFM",
    "Enterobacter_aerogenes_KCTC_2190": "Enterobacter aerogenes KCTC 2190",
    "Streptococcus_pneumoniae_R6": "Streptococcus pneumoniae R6",
    "Staphylococcus_aureus_subsp_aureus_NCTC_8325": "Staphylococcus aureus subsp aureus NCTC 8325",
    "Bacteroides_fragilis_YCH46": "Bacteroides fragilis YCH46",
    "Lactococcus_lactis_subsp_lactis_Il1403": "Lactococcus lactis subsp lactis Il1403",
    "Enterococcus_faecalis_V583": "Enterococcus faecalis V583",
    "Streptococcus_agalactiae_2603V/R": "Streptococcus agalactiae 2603V/R",
    "Enterococcus_faecium_DO": "Enterococcus faecium DO",
    "Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586": "Fusobacterium nucleatum subsp nucleatum ATCC 25586",
}

PSEUDOCOUNT_A = 0.001


def panel_a(ax):
    comp = pd.read_csv(DATA_DIR / "community_primary_tradeoff_0p5" /
                        "species_activity_weighted_composition.csv")
    comp["norm"] = comp["species"].apply(_norm)
    lut = {row.norm: row for row in comp.itertuples()}

    rows = []
    for sp in PANEL_A_SPECIES:
        r = lut[_norm(sp)]
        nd, dd = r.ND_SCD_abundance, r.DD_HFD_abundance
        log2fc = np.log2((dd + PSEUDOCOUNT_A) / (nd + PSEUDOCOUNT_A))
        rows.append((sp, nd, dd, log2fc))
    df = pd.DataFrame(rows, columns=["species", "ND", "DD", "log2fc"]).sort_values("log2fc")

    y = np.arange(len(df))
    colors = [COL_UP if v > 0 else COL_DOWN for v in df["log2fc"]]

    ax.barh(y, df["log2fc"], color=colors, edgecolor="black", linewidth=0.4, height=0.7)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([DISPLAY_NAME[s] for s in df["species"]], fontsize=7.2, style="italic")
    ax.set_xlabel("log2 fold change (DD/HFD vs ND/SCD)")
    ax.set_title("a | Largest modeled-community composition shifts\namong retained species",
                  fontweight="bold", loc="left", fontsize=10)

    handles = [mpatches.Patch(color=COL_UP, label="Increased under DD/HFD"),
               mpatches.Patch(color=COL_DOWN, label="Decreased under DD/HFD")]
    ax.legend(handles=handles, fontsize=7, loc="lower right")

    ax.text(0.02, 0.03, "27 of 57 species retained\nin primary community",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5,
            style="italic", color="#333333",
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="#999999", linewidth=0.5, pad=2))
    _clean(ax)


ACETATE = {
    0.3: {"ND_SCD": 0.5791, "DD_HFD": 0.1196},
    0.5: {"ND_SCD": 0.9652, "DD_HFD": 0.1994},
    0.7: {"ND_SCD": 1.6863, "DD_HFD": 0.2802},
}


def panel_b(ax):
    tradeoffs = sorted(ACETATE)
    nd_vals = [ACETATE[t]["ND_SCD"] for t in tradeoffs]
    dd_vals = [ACETATE[t]["DD_HFD"] for t in tradeoffs]

    ax.plot(tradeoffs, dd_vals, "o-", color=COL_ACETATE_DD, label="DD_HFD")
    ax.plot(tradeoffs, nd_vals, "o-", color=COL_ACETATE_ND, label="ND_SCD")

    ax.set_xlabel("MICOM cooperative tradeoff")
    ax.set_ylabel("Acetate community exchange flux")
    ax.set_title("b | Acetate export remains lower under DD/HFD",
                  fontweight="bold", loc="left", fontsize=10)

    ax.text(0.02, 0.96, "Positive exchange = microbial export",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=7.5, style="italic", color="#333333")

    nd_ref, dd_ref = ACETATE[0.5]["ND_SCD"], ACETATE[0.5]["DD_HFD"]
    pct = (dd_ref - nd_ref) / nd_ref * 100
    ax.annotate(f"{nd_ref:.3f} → {dd_ref:.3f}; {pct:.1f}%",
                xy=(0.5, dd_ref), xytext=(0.55, 0.75),
                fontsize=7.5, color="#333333",
                arrowprops=dict(arrowstyle="-", color="#888888", lw=0.6),
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.5))

    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(0.02, 0.88))
    _clean(ax)


PANEL_C_GROUPS = [
    ("Diet-dominant",                    260,  "#E69F00"),
    ("Diet-dominant\n(opposed by microbiome)", 115, "#F5C842"),
    ("Microbiome-dominant",               27,  "#56B4E9"),
    ("Microbiome-dominant\n(opposed by diet)",  2,  "#89CFF0"),
    ("Stable",                          3322,  "#AAAAAA"),
]


def panel_c(ax):
    labels = [g[0] for g in PANEL_C_GROUPS]
    sizes  = [g[1] for g in PANEL_C_GROUPS]
    colors = [g[2] for g in PANEL_C_GROUPS]
    total = sum(sizes)
    explode = [0.04, 0.04, 0.06, 0.06, 0.0]

    wedges, texts, autotexts = ax.pie(
        sizes, labels=None, colors=colors,
        autopct=lambda p: f"{p:.1f}%" if p > 0.3 else "",
        startangle=90, explode=explode,
        wedgeprops={"edgecolor": "white", "linewidth": 1.0},
        pctdistance=0.72,
    )
    for at in autotexts:
        at.set_fontsize(7.5)
        at.set_color("white")
        at.set_fontweight("bold")
    for i, at in zip(range(len(sizes)), autotexts):
        if sizes[i] / total * 100 < 2:
            at.set_color("black")
            at.set_position((at.get_position()[0] * 1.55, at.get_position()[1] * 1.55))

    centre = plt.Circle((0, 0), 0.45, fc="white")
    ax.add_artist(centre)
    ax.text(0, 0.08, f"{total:,}", ha="center", va="center", fontsize=11, fontweight="bold")
    ax.text(0, -0.14, "reactions", ha="center", va="center", fontsize=8, color="#555555")

    ax.legend(wedges, labels, loc="upper center", bbox_to_anchor=(0.5, -0.02),
              fontsize=6.5, framealpha=0.8, ncol=2, columnspacing=1.0, handletextpad=0.5)
    ax.set_title("c | Primary diet–microbiome attribution",
                 fontweight="bold", loc="left", fontsize=10)


DRIVER_COLORS = {
    "Diet": "#E69F00",
    "Diet (opposed by microbiome)": "#F5C842",
    "Microbiome": "#56B4E9",
    "Microbiome (opposed by diet)": "#89CFF0",
    "Stable": "#AAAAAA",
}
DRIVER_ORDER = ["Diet", "Diet (opposed by microbiome)", "Microbiome",
                "Microbiome (opposed by diet)", "Stable"]


def panel_d(ax):
    attr = pd.read_csv(DATA_DIR / "host_scenarios" / "primary" / "attribution" /
                        "flux_attribution_analysis.csv")
    annot = pd.read_csv(ANNOT_FILE)
    m = attr.merge(annot[["reaction_id", "subsystem"]], on="reaction_id", how="left")
    thresh = m[m["delta_microbiome"].abs() >= 0.01].copy()
    assert len(thresh) == 47, f"expected 47 threshold reactions, got {len(thresh)}"
    # One threshold reaction (PIt2m_2) has no subsystem annotation in
    # RQ4_reaction_annotations.csv. Keep it visible as its own row rather
    # than silently dropping it from the pathway breakdown, so the 15
    # annotated-subsystem bars plus this row sum to the full 47.
    thresh["subsystem"] = thresh["subsystem"].fillna("Unannotated")

    ct = pd.crosstab(thresh["subsystem"], thresh["dominant_driver"])
    for d in DRIVER_ORDER:
        if d not in ct.columns:
            ct[d] = 0
    ct = ct[DRIVER_ORDER]
    ct["Total"] = ct.sum(axis=1)
    ct = ct.sort_values("Total", ascending=False)

    y = np.arange(len(ct))
    left = np.zeros(len(ct))
    for d in DRIVER_ORDER:
        vals = ct[d].values
        ax.barh(y, vals, left=left, color=DRIVER_COLORS[d], label=d,
                edgecolor="white", linewidth=0.5, height=0.72)
        left = left + vals

    for yi, total in zip(y, ct["Total"]):
        ax.text(total + 0.15, yi, str(int(total)), va="center", fontsize=7.5)

    ax.set_yticks(y)
    labels = list(ct.index)
    ax.set_yticklabels(labels, fontsize=7.2)
    for lbl in ax.get_yticklabels():
        if lbl.get_text() == "Cholesterol Metabolism":
            lbl.set_fontweight("bold")
            lbl.set_color("#0B5394")

    if "Cholesterol Metabolism" in labels:
        idx = labels.index("Cholesterol Metabolism")
        ax.axhspan(idx - 0.42, idx + 0.42, color="#0B5394", alpha=0.06, zorder=0)

    ax.set_xlabel("HFD reactions with |Δv$_{microbiome}$| ≥ 0.01")
    ax.set_title("d | Pathway distribution of the 47 threshold-attributable reactions",
                 fontweight="bold", loc="left", fontsize=10)

    handles = [mpatches.Patch(color=DRIVER_COLORS[d], label=d) for d in DRIVER_ORDER]
    ax.legend(handles=handles, fontsize=6.3, loc="upper center", bbox_to_anchor=(0.5, -0.16),
              ncol=3, framealpha=0.85, columnspacing=0.8, handletextpad=0.4)
    ax.set_xlim(0, 9.5)
    _clean(ax)


def make_figure6():
    fig = plt.figure(figsize=(13, 10.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.35], hspace=0.60, wspace=0.32,
                           left=0.07, right=0.97, top=0.88, bottom=0.09)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    fig.suptitle("Figure 6 | Gut microbiome community modeling and hepatic flux attribution",
                 fontsize=13, fontweight="bold", y=0.995)

    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_c)
    panel_d(ax_d)

    fig.savefig(OUT_DIR / "Figure6_hybrid.png")
    fig.savefig(OUT_DIR / "Figure6_hybrid.pdf")
    plt.close(fig)
    print(f"Saved to {OUT_DIR.resolve()}/Figure6_hybrid.png")


if __name__ == "__main__":
    make_figure6()

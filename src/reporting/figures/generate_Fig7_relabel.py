#!/usr/bin/env python3
"""
Figure 7 (relabeled) — Reaction-level cross-scale integration.

Same underlying data/counts as the currently-approved Figure 7 (Roland
confirmed no numeric issues); only the RQ1-4 labels are replaced with
descriptive names, per his request, since RQ1-4 notation is being removed
from the manuscript.

Row-level set membership for the 15 bulk-HFD reference reactions was
derived directly from the currently-approved figure (each anchor's
dot pattern across the RQ2/RQ3/RQ4 columns), then cross-checked against
every summary statistic printed on both panels of that figure:
  - panel (a) subtitle: 4/15 genetic core, 6/15 single-cell, 7/15 microbiome
  - panel (b) bars: 6 / 3 / 2 / 2 / 2 (sums to 15)
All checks passed exactly, so this membership table is treated as a
faithful transcription of the approved figure's data, not a re-derivation
from scratch.

Note on panel (b) category labels: Roland's requested relabeling
("Bulk + single-cell" and "Bulk + genetic + microbiome") does not match
the actual combinations present in the data. The two middle categories
are actually "Bulk + genetic + single-cell" (RQ2 ∩ RQ3, no RQ4) and
"Bulk + single-cell + microbiome" (RQ3 ∩ RQ4, no RQ2) — there is no
"genetic + microbiome only" (RQ2 ∩ RQ4) category in the data at all.
Using the corrected compound labels here; flagged to Roland.
"""

from pathlib import Path

import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

plt.rcParams.update({
    "font.family":        "Arial",
    "font.size":          9,
    "axes.labelsize":     9,
    "axes.titlesize":     10,
    "xtick.labelsize":    8,
    "ytick.labelsize":    8,
    "legend.fontsize":    8,
    "figure.dpi":         300,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "savefig.pad_inches": 0.15,
    "axes.linewidth":     0.8,
})

OUT_DIR = Path(os.environ.get("FIG_OUT", "results/figures/figure7"))  # [release] configurable
OUT_DIR.mkdir(parents=True, exist_ok=True)

COL_BULK   = "black"
COL_GENETIC = "#4C72B0"   # blue  (was RQ2 core)
COL_SC      = "#DD8452"   # orange (was RQ3 attrib.)
COL_MICRO   = "#55A868"   # green (was RQ4 attrib.)
COL_GRAY    = "#B0B0B0"
COL_RED     = "#C44E52"

# 15 bulk-HFD reference reactions, transcribed from the currently-approved
# figure's dot pattern. Each maps to the subset of {genetic, singlecell,
# microbiome} it also belongs to (bulk membership is implicit/universal).
ANCHORS = [
    ("EX_co2_e", {"genetic", "singlecell", "microbiome"}),
    ("EX_nh4_e", {"singlecell", "microbiome"}),
    ("EX_o2_e",  {"genetic", "singlecell"}),
    ("EX_o2s_e", set()),
    ("CO2t",     {"genetic", "singlecell", "microbiome"}),
    ("CO2tm",    {"singlecell", "microbiome"}),
    ("H2O2tm",   set()),
    ("AKGDm",    {"microbiome"}),
    ("AKGMALtm", {"microbiome"}),
    ("NAt",      set()),
    ("O2St",     set()),
    ("O2Stm",    set()),
    ("O2t",      {"genetic", "singlecell"}),
    ("O2tm",     {"microbiome"}),
    ("STCOAtx",  set()),
]

COLUMNS = [
    ("bulk",       "Bulk HFD\nreference",            COL_BULK),
    ("genetic",    "Genetic\nuniversal core",         COL_GENETIC),
    ("singlecell", "Single-cell\nWD–chow response",   COL_SC),
    ("microbiome", "Microbiome-\nassociated response", COL_MICRO),
]


def panel_a(ax):
    n = len(ANCHORS)
    y = np.arange(n)[::-1]  # top-to-bottom in listed order

    for xi, (key, _, color) in enumerate(COLUMNS):
        for yi, (rxn, members) in zip(y, ANCHORS):
            present = (key == "bulk") or (key in members)
            if present:
                ax.scatter(xi, yi, s=55, color=color, zorder=3,
                           edgecolor="black" if key == "bulk" else "none", linewidth=0.4)

    ax.set_yticks(y)
    ax.set_yticklabels([rxn for rxn, _ in ANCHORS], fontsize=8)
    ax.set_xticks(range(len(COLUMNS)))
    ax.set_xticklabels([lbl for _, lbl, _ in COLUMNS], fontsize=8)
    ax.set_xlim(-0.6, len(COLUMNS) - 0.4)
    ax.set_ylim(-0.8, n - 0.2)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(left=False)
    ax.grid(axis="y", color="#EEEEEE", linewidth=0.6, zorder=0)

    n_genetic = sum(1 for _, m in ANCHORS if "genetic" in m)
    n_sc = sum(1 for _, m in ANCHORS if "singlecell" in m)
    n_micro = sum(1 for _, m in ANCHORS if "microbiome" in m)
    ax.text(0.0, 1.22, "a | Reaction-level tracing of 15 bulk HFD reference reactions",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=10, color="black",
            fontweight="bold")
    ax.text(0.0, 1.16,
            f"{n_genetic}/15 overlap the genetic core; {n_sc}/15 the single-cell response;\n"
            f"{n_micro}/15 the microbiome-associated set",
            transform=ax.transAxes, ha="left", va="top", fontsize=8.5, color="#333333")


def panel_b(ax):
    categories = [
        ("Bulk only", 6, COL_GRAY),
        ("Bulk + microbiome", 3, COL_GENETIC),
        ("Bulk + genetic\n+ single-cell", 2, COL_GENETIC),
        ("Bulk + single-cell\n+ microbiome", 2, COL_GENETIC),
        ("Bulk + genetic + single-cell\n+ microbiome", 2, COL_RED),
    ]
    x = np.arange(len(categories))
    heights = [c[1] for c in categories]
    colors = [c[2] for c in categories]

    ax.bar(x, heights, color=colors, edgecolor="black", linewidth=0.5, width=0.6, zorder=3)
    for xi, h in zip(x, heights):
        ax.text(xi, h + 0.15, str(h), ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels([c[0] for c in categories], fontsize=7, rotation=20, ha="right",
                        rotation_mode="anchor")
    ax.set_ylabel("Number of the 15 reference reactions")
    ax.set_ylim(0, 7.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#EEEEEE", linewidth=0.6, zorder=0)

    ax.text(0.0, 1.22, "b | Interpretive boundary",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=10, color="black",
            fontweight="bold")
    ax.text(0.0, 1.16, "Only 2/15 reference reactions occur in all three additional contexts",
            transform=ax.transAxes, ha="left", va="top", fontsize=8.5, color="#333333")
    ax.text(0.98, 0.97, "Contextual overlap,\nnot condition-matched validation",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            style="italic", color="#444444")


def make_figure7():
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(15.5, 6.4),
                                      gridspec_kw={"width_ratios": [1.4, 1], "wspace": 0.38,
                                                   "bottom": 0.20, "top": 0.72})
    fig.suptitle("Figure 7 | Reaction-level cross-scale integration.",
                 fontsize=12, y=0.985, fontweight="bold")

    panel_a(ax_a)
    panel_b(ax_b)

    fig.savefig(OUT_DIR / "Figure7_relabel.png")
    fig.savefig(OUT_DIR / "Figure7_relabel.pdf")
    plt.close(fig)
    print(f"Saved to {OUT_DIR.resolve()}/Figure7_relabel.png")


if __name__ == "__main__":
    make_figure7()

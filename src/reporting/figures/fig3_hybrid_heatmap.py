#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 3 hybrid (2026-09-25d): panel (b) unchanged from the current manuscript;
panel (c) rebuilt as an 11-reaction x 9-strain heatmap (2026-09-25c pass);
panel (a) now also rebuilt (this pass) as stacked horizontal bars breaking
each strain's magnitude-responsive count into three conservation tiers --
universal (9/9), shared (2-8/9), strain-unique (1/9) -- using CURRENT
values/terminology throughout (not the older, now-incorrect per-strain counts
or "significantly altered" language from the earlier figure version this
layout is borrowed from).

Why stacked bars for panel (a): a plain per-strain total (the previous
version) shows breadth of response but nothing about composition. Results 3.3
now specifically discusses strain-unique counts and the universal core, so
showing each strain's bar broken into universal/shared/unique ties the figure
directly to the text it illustrates.

Why a heatmap for panel (c): it shows, in one view, all 11 core reactions x
all 9 strains x the actual signed magnitude, so the reader can see directly
that O2t/EX_o2_e flip sign specifically in C57BL/6J, rather than only reading
it in the caption. The mean+/-SD bar chart it replaces collapses across
strains and hides exactly that fact.

Data: <data_dir>/<STRAIN>.csv, one file per strain, columns include
ReactionID, ReactionName, Subsystem, Diff(HFD-SCD). Same source files used by
fig3_regenerate_directionfix.py (2026-09-25b pass) -- panel (b) is regenerated
with identical code/data so it is byte-for-byte the same content as currently
embedded in the manuscript, with bold titles added per Roland's standing
"bold every figure caption" rule.

Usage: python3 fig3_hybrid_heatmap.py [data_dir] [out_png]
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "data"
OUT_PNG = sys.argv[2] if len(sys.argv) > 2 else "Figure3_hybrid.png"

STRAINS = {
    "129S1/SvImJ": f"{DATA_DIR}/RQ2_129S1SvImJ_edges_HFD_vs_SCD.csv",
    "A/J": f"{DATA_DIR}/RQ2_AJ_edges_HFD_vs_SCD.csv",
    "C57BL/6J": f"{DATA_DIR}/RQ2_C57BL6J_edges_HFD_vs_SCD.csv",
    "CAST/EiJ": f"{DATA_DIR}/RQ2_CASTEiJ_edges_HFD_vs_SCD.csv",
    "DBA/2J": f"{DATA_DIR}/RQ2_DBA2J_edges_HFD_vs_SCD.csv",
    "NOD/ShiLtJ": f"{DATA_DIR}/RQ2_NODShiLtJ_edges_HFD_vs_SCD.csv",
    "NZO/HlLtJ": f"{DATA_DIR}/RQ2_NZOHlLtJ_edges_HFD_vs_SCD.csv",
    "PWK/PhJ": f"{DATA_DIR}/RQ2_PWKPhJ_edges_HFD_vs_SCD.csv",
    "WSB/EiJ": f"{DATA_DIR}/RQ2_WSBEiJ_edges_HFD_vs_SCD.csv",
}
THRESH = 0.20
NZO_STRAIN = "NZO/HlLtJ"

SUBSYS_COLOR = {
    "Tyrosine metabolism": "#7D3C98",
    "Extracellular exchange": "#229954",
    "Transport, Extracellular": "#2874A6",
}
NAME_OVERRIDE = {
    # Cleaned display names for readability; underlying data/ID unaffected.
    "34HPPOR": "4-Hydroxyphenylpyruvate oxidoreductase",
    "HGNTOR": "Homogentisate 1,2-dioxygenase",
    "O2t": "O2 transport via diffusion",
}

# ---------------------------------------------------------------- load data
per_strain_sets, per_strain_counts, diff_series = {}, {}, {}
meta = {}  # ReactionID -> (ReactionName, Subsystem)
for strain, fname in STRAINS.items():
    df = pd.read_csv(fname).dropna(subset=["Diff(HFD-SCD)"]).drop_duplicates(subset=["ReactionID"])
    df["absdiff"] = df["Diff(HFD-SCD)"].abs()
    hit = df[df["absdiff"] >= THRESH]
    per_strain_sets[strain] = set(hit["ReactionID"])
    per_strain_counts[strain] = len(per_strain_sets[strain])
    diff_series[strain] = df.set_index("ReactionID")["Diff(HFD-SCD)"]
    for _, row in df.iterrows():
        meta.setdefault(row["ReactionID"], (row["ReactionName"], row["Subsystem"]))

all_ids = set().union(*per_strain_sets.values())
union_n = len(all_ids)
core = sorted(set.intersection(*per_strain_sets.values()))
core_n = len(core)
assert union_n == 243 and core_n == 11, (
    f"Recomputed union={union_n}, core={core_n} -- does not match manuscript "
    f"(243/11); STOP, do not proceed on mismatched data."
)

strain_labels = list(STRAINS.keys())
matrix = np.array([[diff_series[s].get(rid, np.nan) for s in strain_labels] for rid in core])

# direction consistency per core reaction (sign-consistent across all 9 strains?)
consistent = {}
for rid, row in zip(core, matrix):
    signs = np.sign(row)
    consistent[rid] = bool(np.all(signs == signs[0]))
mixed = [rid for rid in core if not consistent[rid]]
assert set(mixed) == {"O2t", "EX_o2_e"}, (
    f"Expected exactly O2t/EX_o2_e to show mixed directionality; got {mixed}. "
    f"STOP -- this contradicts the corrected Results 3.3 text, do not proceed."
)

from collections import Counter
share_count = {}  # ReactionID -> number of strains (of 9) sharing it
for rid in all_ids:
    share_count[rid] = sum(1 for ids in per_strain_sets.values() if rid in ids)
spectrum = Counter(share_count.values())
spectrum_x = sorted(spectrum)
spectrum_y = [spectrum[k] for k in spectrum_x]

# --- per-strain conservation-tier breakdown, for panel (a) ---
def tier(n):
    if n == 9:
        return "universal"
    if n == 1:
        return "unique"
    return "shared"

tier_counts = {s: {"universal": 0, "shared": 0, "unique": 0} for s in strain_labels}
for s in strain_labels:
    for rid in per_strain_sets[s]:
        tier_counts[s][tier(share_count[rid])] += 1
    tot = sum(tier_counts[s].values())
    assert tot == per_strain_counts[s], (
        f"{s}: tier breakdown sums to {tot}, expected {per_strain_counts[s]}")
    assert tier_counts[s]["universal"] == core_n, (
        f"{s}: universal-tier count {tier_counts[s]['universal']} != core_n {core_n} "
        f"-- every strain must contain all {core_n} universal-core reactions by definition")

# ------------------------------------------------------------------- layout
fig = plt.figure(figsize=(19, 15), dpi=200)
outer = fig.add_gridspec(2, 1, height_ratios=[1, 1.35], hspace=0.32)
top = outer[0].subgridspec(1, 2, wspace=0.28)
ax_a = fig.add_subplot(top[0, 0])
ax_b = fig.add_subplot(top[0, 1])
ax_c = fig.add_subplot(outer[1])

# --- Panel (a): per-strain magnitude-responsive reactions, by conservation tier ---
TIER_COLOR = {"universal": "#C0392B", "shared": "#4472A8", "unique": "#229954"}
TIER_LABEL = {
    "universal": "Universal magnitude-responsive (9/9)",
    "shared": "Shared magnitude-responsive (2–8/9)",
    "unique": "Strain-unique (1/9)",
}
y_pos = np.arange(len(strain_labels))
universal_vals = [tier_counts[s]["universal"] for s in strain_labels]
shared_vals = [tier_counts[s]["shared"] for s in strain_labels]
unique_vals = [tier_counts[s]["unique"] for s in strain_labels]
left_shared = universal_vals
left_unique = [u + sh for u, sh in zip(universal_vals, shared_vals)]

bars_universal = ax_a.barh(y_pos, universal_vals, color=TIER_COLOR["universal"],
                            edgecolor="black", linewidth=0.6, label=TIER_LABEL["universal"])
bars_shared = ax_a.barh(y_pos, shared_vals, left=left_shared, color=TIER_COLOR["shared"],
                         edgecolor="black", linewidth=0.6, label=TIER_LABEL["shared"])
bars_unique = ax_a.barh(y_pos, unique_vals, left=left_unique, color=TIER_COLOR["unique"],
                         edgecolor="black", linewidth=0.6, label=TIER_LABEL["unique"])

# hatch overlay on the NZO/HlLtJ row only (within-strain testing non-estimable),
# applied on top of all three tier segments so the tier breakdown stays visible
nzo_idx = strain_labels.index(NZO_STRAIN)
for group in (bars_universal, bars_shared, bars_unique):
    group[nzo_idx].set_hatch("//")

totals = [per_strain_counts[s] for s in strain_labels]
for y, tot in zip(y_pos, totals):
    ax_a.text(tot + 1.5, y, str(tot), va="center", ha="left", fontsize=9)
mean_v = np.mean(totals)
ax_a.axvline(mean_v, color="gray", linestyle="--", linewidth=1)
ax_a.text(mean_v, len(strain_labels) - 0.35, f"mean = {mean_v:.1f}", color="gray",
          fontsize=9, ha="center", va="top")

ax_a.set_yticks(y_pos)
ax_a.set_yticklabels(strain_labels, fontsize=9)
ax_a.invert_yaxis()
ax_a.set_xlabel("Magnitude-responsive reactions (|Δv| ≥ 0.20)")
ax_a.set_xlim(0, max(totals) + 10)
ax_a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=1, fontsize=8, frameon=False)
ax_a.set_title(
    "(a) Per-strain magnitude-responsive reactions, by conservation tier\n"
    "hatched = NZO/HlLtJ (within-strain testing non-estimable, HFD n=1)",
    fontsize=10, fontweight="bold")

# --- Panel (b): conservation spectrum (unchanged content) ---
bar_colors = ["#c0392b" if k == 9 else "#4472A8" for k in spectrum_x]
ax_b.bar([str(k) for k in spectrum_x], spectrum_y, color=bar_colors, edgecolor="black", linewidth=0.6)
for x, v in zip(spectrum_x, spectrum_y):
    ax_b.text(str(x), v + 1.2, str(v), ha="center", va="bottom", fontsize=10)
ax_b.set_xlabel("Number of strains sharing the reaction (of 9)")
ax_b.set_ylabel("Number of reactions")
ax_b.set_title(
    f"(b) Conservation spectrum of the {union_n}-reaction union\n"
    f"{core_n} reactions in the universal core (red)",
    fontsize=10, fontweight="bold")

# --- Panel (c): NEW heatmap, 11 core reactions x 9 strains ---
vmax = np.nanmax(np.abs(matrix))
vmax = np.ceil(vmax * 10) / 10  # round up to nearest 0.1, e.g. 3.391 -> 3.4
im = ax_c.imshow(matrix, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")

for i, row in enumerate(matrix):
    for j, val in enumerate(row):
        txt_color = "white" if abs(val) > vmax * 0.6 else "black"
        ax_c.text(j, i, f"{val:+.2f}", ha="center", va="center", fontsize=8.5, color=txt_color)

ax_c.set_xticks(range(len(strain_labels)))
ax_c.set_xticklabels(strain_labels, rotation=45, ha="right", fontsize=9)
ax_c.set_yticks(range(len(core)))
ylabels = []
for rid in core:
    rname = NAME_OVERRIDE.get(rid, meta[rid][0])
    ylabels.append(f"{rid}  ({rname})")
ax_c.set_yticklabels(ylabels, fontsize=9)
for tick, rid in zip(ax_c.get_yticklabels(), core):
    tick.set_color(SUBSYS_COLOR[meta[rid][1]])
ax_c.set_xlim(-0.5, len(strain_labels) - 0.5)

# direction marker column, just right of the heatmap
dir_x = len(strain_labels) - 0.5 + 0.6
for i, rid in enumerate(core):
    if consistent[rid]:
        ax_c.plot(dir_x, i, marker="o", markersize=9, markerfacecolor="black",
                  markeredgecolor="black", clip_on=False)
    else:
        ax_c.plot(dir_x, i, marker="o", markersize=9, fillstyle="right",
                  markerfacecolor="#E67E22", markerfacecoloralt="white",
                  markeredgecolor="black", clip_on=False)
ax_c.text(dir_x, -1.0, "Dir.", ha="center", va="bottom", fontsize=9, style="italic", clip_on=False)

cbar = fig.colorbar(im, ax=ax_c, fraction=0.025, pad=0.1)
cbar.set_label("Δ Flux: HFD − SCD (mmol/gDW/h)", fontsize=9)

ax_c.set_title(
    f"(c) Universal magnitude-responsive reactions across nine strains (n = {core_n})\n"
    "Nine reactions show concordant directionality; O2t and EX_o2_e reverse direction in C57BL/6J.",
    fontsize=10, fontweight="bold")

# legend: subsystem colors + direction markers
subsys_handles = [Patch(facecolor=c, edgecolor="none", label=s) for s, c in SUBSYS_COLOR.items()]
dir_handles = [
    Line2D([0], [0], marker="o", linestyle="None", markersize=8, markerfacecolor="black",
           markeredgecolor="black", label="Consistent direction (9/9)"),
    Line2D([0], [0], marker="o", linestyle="None", markersize=8, fillstyle="right",
           markerfacecolor="#E67E22", markerfacecoloralt="white", markeredgecolor="black",
           label="Mixed direction"),
]
ax_c.legend(handles=subsys_handles + dir_handles, loc="upper center",
            bbox_to_anchor=(0.42, -0.14), ncol=5, frameon=False, fontsize=8.5,
            title="Row colours = Subsystem                              Direction markers",
            title_fontsize=8.5)

fig.suptitle(
    "Figure 3 | Genetic-background analysis using standard FBA with aggregate (mean-collapsed) expression per group.",
    fontsize=13, y=0.995, fontweight="bold")

plt.tight_layout(rect=[0, 0.02, 1, 0.97])
plt.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
print("Saved", OUT_PNG)
print("union:", union_n, "core:", core_n, "core reactions:", core)
print("mixed-direction reactions:", mixed)
print("per-strain counts:", per_strain_counts)
print("mean per-strain:", round(mean_v, 1))

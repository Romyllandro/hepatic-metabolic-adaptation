#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2, hybrid composition (2026-09-25): all data/analysis from the current
standard-FBA regeneration (fig2_regenerate.py, built for the 2026-09-21
RQ1-primary revert); panel layout for the pathway-response panel borrowed from
the pre-revert original figure, per Roland's review.

Panels:
  (a) PCA of batch-adjusted flux profiles (no confidence ellipses -- points and
      the restricted-PERMANOVA statement are sufficient; ellipses would imply
      more about group uncertainty than was actually tested, especially at
      n=12 for KD/WD).
  (b) Euclidean distance between diet centroids (PC1/PC2).
  (c) FDR-significant reaction counts per pairwise contrast; KD-vs-WD marked
      NE (no dataset contains both groups -- confounded with batch). Fixes the
      duplicate "(b)" labeling bug in fig2_regenerate.py's output (the
      centroid-distance and significant-reaction panels were both titled "(b)").
  (d) Dot plot of effect sizes (Cohen's d) for the 10 reactions significant in
      all three of HFD-vs-SCD, KD-vs-SCD, WD-vs-SCD (the cross-diet convergent
      signature). Unchanged from fig2_regenerate.py's panel (c), renumbered.
  (e) Directional pathway responses across dietary interventions: three
      side-by-side small multiples (HFD / KD / WD), each vs. SCD, sharing one
      subsystem list (the top 12 by combined significant-reaction count) and
      ONE shared x-axis scale across all three panels so magnitudes are
      directly, honestly comparable diet-to-diet -- this is the one place the
      pre-revert original figure's layout read better than a single densely
      overlaid bar chart, but its data (and titles/colors) are replaced here
      with the current, correct standard-FBA counts and this figure's own
      Okabe-Ito diet-color mapping (for consistency with panels a and d).
      NOTE (found while adding count labels for this revision): this panel
      counts each diet's FDR-significant reactions by Subsystem x Direction
      from RQ1_reaction_stats_{d}_vs_SCD.csv, asserted to sum to exactly
      95 / 154 / 134 (HFD/KD/WD vs SCD, matching panel c). It deliberately
      does NOT use the N_positive/N_negative columns in
      RQ1_subsystem_analysis_{d}_vs_SCD.csv the way fig2_regenerate.py's own
      single-panel (d) does -- those columns count ALL reactions in a
      subsystem by effect sign regardless of significance (962 total for HFD
      alone, vs. 95 truly significant), which silently inflates every bar in
      that panel. See the data-hygiene note in this repo's history for the
      precedent of a stale/mis-scoped source file causing this kind of
      quiet error.

Source data: identical to fig2_regenerate.py (see that script's docstring /
the bundled README.md). No new source files are required.

Run: python3 fig2_hybrid_v2.py --data-dir data_fig2 --out Figure2_hybrid.png
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Okabe-Ito colorblind-safe palette (identical to fig2_regenerate.py, kept
# consistent across every panel in this hybrid figure, including the new
# panel (e), rather than reusing the old figure's arbitrary title colors)
DIET_COLOR = {
    "SCD": "#000000",   # black -- reference/control
    "HFD": "#E69F00",   # orange
    "KD":  "#56B4E9",   # sky blue
    "WD":  "#CC79A7",   # reddish purple
}
UP_COLOR = "#009E73"    # bluish green
DOWN_COLOR = "#D55E00"  # vermillion
NE_COLOR = "#BBBBBB"    # neutral gray for non-estimable


def load(data_dir, name):
    return pd.read_csv(os.path.join(data_dir, name))


def panel_a(ax, data_dir):
    scores = load(data_dir, "RQ1_pca_scores.csv")
    var = load(data_dir, "RQ1_pca_variance_explained.csv").set_index("PC")
    pc1_pct = var.loc["PC1", "Variance"] * 100
    pc2_pct = var.loc["PC2", "Variance"] * 100

    for diet in ["SCD", "HFD", "KD", "WD"]:
        sub = scores[scores["Group"] == diet]
        ax.scatter(sub["PC1"], sub["PC2"], s=26, alpha=0.85,
                   color=DIET_COLOR[diet], edgecolor="white", linewidth=0.4,
                   label=f"{diet} (n={len(sub)})")
    ax.set_xlabel(f"PC1 ({pc1_pct:.1f}%)")
    ax.set_ylabel(f"PC2 ({pc2_pct:.1f}%)")
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.axhline(0, color="#DDDDDD", lw=0.6, zorder=0)
    ax.axvline(0, color="#DDDDDD", lw=0.6, zorder=0)
    ax.set_title("(a) PCA of batch-adjusted flux profiles\nPERMANOVA p = 0.001 (restricted by dataset)",
                 fontsize=9, loc="left", fontweight="bold")


def panel_b(ax, data_dir):
    """Centroid distance only (was the left half of fig2_regenerate.py's panel_b)."""
    cent = load(data_dir, "RQ1_centroid_distances_PC12.csv")
    order = ["HFD_vs_SCD", "KD_vs_SCD", "WD_vs_SCD", "HFD_vs_KD", "HFD_vs_WD", "KD_vs_WD"]
    cent = cent.set_index("Contrast").reindex(order)
    labels = [c.replace("_vs_", " vs ") for c in order]

    ax.barh(labels, cent["Distance"], color="#4C72B0", edgecolor="black", linewidth=0.5)
    ax.invert_yaxis()
    ax.set_xlabel("Euclidean distance\n(diet centroids, PC1/PC2)")
    ax.set_title("(b) Centroid distance", fontsize=9, loc="left", fontweight="bold")
    return order, labels


def panel_c(ax, data_dir, order, labels):
    """Significant-reaction counts only (was the right half of fig2_regenerate.py's
    panel_b, mislabeled "(b)" a second time there -- fixed to "(c)" here)."""
    sig_counts = {}
    for c in order:
        if c == "KD_vs_WD":
            sig_counts[c] = None  # non-estimable: no dataset contains both KD and WD
            continue
        df = load(data_dir, f"RQ1_reaction_stats_{c}.csv")
        sig_counts[c] = int((df["Significant"] == True).sum())

    bar_vals = [sig_counts[c] if sig_counts[c] is not None else 0 for c in order]
    colors = ["#4C72B0" if sig_counts[c] is not None else NE_COLOR for c in order]
    ax.barh(labels, bar_vals, color=colors, edgecolor="black", linewidth=0.5)
    ax.invert_yaxis()
    ax.set_xlabel("FDR-significant reactions\n(q < 0.05)")
    ax.set_title("(c) Significant reactions", fontsize=9, loc="left", fontweight="bold")
    for i, c in enumerate(order):
        if sig_counts[c] is None:
            ax.text(2, i, "NE", va="center", ha="left", fontsize=8, color="#555555",
                    fontstyle="italic")
        else:
            ax.text(sig_counts[c] + 2, i, str(sig_counts[c]), va="center", ha="left", fontsize=8)
    legend_elems = [Patch(facecolor=NE_COLOR, edgecolor="black",
                          label="NE: no dataset contains\nboth groups (confounded\nwith batch)")]
    ax.legend(handles=legend_elems, frameon=False, fontsize=7, loc="lower right")
    return sig_counts


def panel_d(ax, data_dir):
    """Cross-diet convergent signature (was fig2_regenerate.py's panel_c; renumbered)."""
    hfd = load(data_dir, "RQ1_reaction_stats_HFD_vs_SCD.csv").set_index("ReactionID")
    kd = load(data_dir, "RQ1_reaction_stats_KD_vs_SCD.csv").set_index("ReactionID")
    wd = load(data_dir, "RQ1_reaction_stats_WD_vs_SCD.csv").set_index("ReactionID")

    sig_hfd = set(hfd[hfd["Significant"] == True].index)
    sig_kd = set(kd[kd["Significant"] == True].index)
    sig_wd = set(wd[wd["Significant"] == True].index)
    core = sorted(sig_hfd & sig_kd & sig_wd)
    assert len(core) == 10, f"expected 10 cross-diet convergent reactions, got {len(core)}: {core}"

    diets = ["HFD", "KD", "WD"]
    tables = {"HFD": hfd, "KD": kd, "WD": wd}
    y_positions = np.arange(len(core))

    for j, diet in enumerate(diets):
        vals = [tables[diet].loc[r, "Cohen_d"] for r in core]
        ax.scatter(vals, y_positions + (j - 1) * 0.22, s=55,
                   color=DIET_COLOR[diet], edgecolor="black", linewidth=0.4,
                   label=diet, zorder=3)

    ax.axvline(0, color="#999999", lw=0.8, zorder=0)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(core, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Effect size (Cohen's d) vs. SCD")
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.set_title("(d) Cross-diet convergent signature (n=10)\nGLNALANaEx is the sole direction-switching reaction",
                 fontsize=9, loc="left", fontweight="bold")
    if "GLNALANaEx" in core:
        i = core.index("GLNALANaEx")
        ax.annotate("direction\nswitches", xy=(0, i), xytext=(1.0, i - 1.3),
                    fontsize=7, ha="center", color="#555555",
                    arrowprops=dict(arrowstyle="-", color="#999999", lw=0.6))


def panel_e(axes, data_dir):
    """Directional pathway responses across dietary interventions: three
    side-by-side small multiples (HFD / KD / WD vs. SCD), current standard-FBA
    counts, laid out like the pre-revert original figure's panel (d) -- but
    with a shared subsystem list AND a shared x-axis scale across all three
    sub-panels, so e.g. 'KD showed stronger mitochondrial-transport remodeling
    than HFD' is a claim the reader can verify directly by bar length, not
    just by eye across differently-scaled axes.

    IMPORTANT data-source note: this counts each diet's FDR-significant
    reactions (Significant == True in RQ1_reaction_stats_{d}_vs_SCD.csv),
    split by Direction and grouped by Subsystem -- NOT the N_positive/
    N_negative columns in RQ1_subsystem_analysis_{d}_vs_SCD.csv, which count
    ALL reactions in a subsystem by the sign of their effect regardless of
    significance (verified: those two columns sum to 962 for HFD alone, far
    above the 95 truly significant HFD-vs-SCD reactions from panel (c)).
    fig2_regenerate.py's original single-panel (d) uses the N_positive/
    N_negative columns and is therefore over-counting non-significant
    reactions in that panel -- flagged separately, not fixed here since that
    script is a distinct prior deliverable."""
    diets = ["HFD", "KD", "WD"]
    frames = {}
    for d in diets:
        df = load(data_dir, f"RQ1_reaction_stats_{d}_vs_SCD.csv")
        sig = df[df["Significant"] == True]
        counts = sig.groupby(["Subsystem", "Direction"]).size().unstack(fill_value=0)
        for col in ("Up", "Down"):
            if col not in counts.columns:
                counts[col] = 0
        frames[d] = counts[["Up", "Down"]]
        # cross-check against the panel (c) / RQ1_reaction_stats totals directly
        # (95 / 154 / 134 for HFD / KD / WD vs SCD) rather than trusting the
        # groupby: catches exactly the kind of over-count this replaced.
        expected = {"HFD": 95, "KD": 154, "WD": 134}[d]
        got = int(frames[d].values.sum())
        assert got == expected, f"{d}: subsystem-summed significant count {got} != {expected} from panel (c)"

    # Same "top 12 by combined significant count" subsystem selection as the
    # single-panel version, so this remains the identical analysis -- only the
    # layout changes.
    combined = pd.concat([frames[d] for d in diets], axis=1, keys=diets).fillna(0)
    combined["total"] = combined.sum(axis=1)
    top = combined.sort_values("total", ascending=False).head(12).index
    combined = combined.loc[top]
    n = len(top)
    y = np.arange(n)

    # Shared x-axis range across all three sub-panels (deliberate departure
    # from the original figure, whose three facets happened to share a range
    # only because the pre-revert magnitudes were all small; with the current,
    # larger counts a shared range must be set explicitly or the panels would
    # silently mislead on relative magnitude).
    xmax = max((combined[(d, "Up")].max(), combined[(d, "Down")].max()) for d in diets)
    xmax = max(max(pair) for pair in [(combined[(d, "Up")].max(), combined[(d, "Down")].max()) for d in diets])
    xlim = xmax * 1.30

    # small, constant offset (in data units) so the count label sits just
    # outside the bar tip regardless of that bar's own length
    label_offset = xlim * 0.015

    for j, (ax, d) in enumerate(zip(axes, diets)):
        up = combined[(d, "Up")].values
        down = combined[(d, "Down")].values
        ax.barh(y, up, height=0.62, color=UP_COLOR, edgecolor="black", linewidth=0.3)
        ax.barh(y, -down, height=0.62, color=DOWN_COLOR, edgecolor="black", linewidth=0.3)
        for i in range(n):
            if up[i] > 0:
                ax.text(up[i] + label_offset, i, str(int(up[i])), va="center", ha="left",
                        fontsize=7.5, color=UP_COLOR, fontweight="bold")
            if down[i] > 0:
                ax.text(-down[i] - label_offset, i, str(int(down[i])), va="center", ha="right",
                        fontsize=7.5, color=DOWN_COLOR, fontweight="bold")
        ax.axvline(0, color="black", lw=0.8)
        ax.set_xlim(-xlim, xlim)
        ax.invert_yaxis()
        ax.set_yticks(y)
        if j == 0:
            ax.set_yticklabels(top, fontsize=7.5)
        else:
            ax.set_yticklabels([])
        ax.set_title(d, fontsize=10, color=DIET_COLOR[d], fontweight="bold")
        ax.set_xlabel("N reactions" if j == 1 else "", fontsize=8)

    legend_elems = [Patch(facecolor=UP_COLOR, edgecolor="black", label="Up"),
                    Patch(facecolor=DOWN_COLOR, edgecolor="black", label="Down")]
    axes[-1].legend(handles=legend_elems, frameon=False, fontsize=8, loc="lower right")
    # Row-level title, analogous to the single-panel version's title, spanning
    # the three sub-panels via the middle axis (title text is left-anchored
    # relative to axes[0] using the figure's transform so it reads as one
    # heading for the whole row rather than three separate panel titles).
    axes[0].annotate("(e) Directional pathway responses across dietary interventions\n"
                      "top 12 subsystems by combined significant-reaction count; HFD / KD / WD each vs. SCD, shared scale",
                      xy=(0, 1.12), xycoords=("axes fraction"), fontsize=9, fontweight="bold",
                      ha="left", va="bottom", annotation_clip=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out", default="Figure2_hybrid.png")
    args = ap.parse_args()

    fig = plt.figure(figsize=(20.0, 13.5), dpi=300)
    outer = fig.add_gridspec(2, 1, height_ratios=[1, 1.15], hspace=0.42)
    top_gs = outer[0].subgridspec(1, 3, wspace=0.55)
    # extra-wide gap right after panel (d) so panel (e)'s subsystem labels
    # (rendered outside its own left spine) never overlap panel (d)'s content
    bottom_gs = outer[1].subgridspec(1, 5, width_ratios=[1.35, 0.55, 1, 1, 1], wspace=0.35)

    ax_a = fig.add_subplot(top_gs[0, 0])
    ax_b = fig.add_subplot(top_gs[0, 1])
    ax_c = fig.add_subplot(top_gs[0, 2])
    ax_d = fig.add_subplot(bottom_gs[0, 0])
    # bottom_gs[0, 1] is left empty on purpose (spacer for panel-e labels)
    ax_e1 = fig.add_subplot(bottom_gs[0, 2])
    ax_e2 = fig.add_subplot(bottom_gs[0, 3])
    ax_e3 = fig.add_subplot(bottom_gs[0, 4])

    panel_a(ax_a, args.data_dir)
    order, labels = panel_b(ax_b, args.data_dir)
    panel_c(ax_c, args.data_dir, order, labels)
    panel_d(ax_d, args.data_dir)
    panel_e([ax_e1, ax_e2, ax_e3], args.data_dir)

    fig.suptitle("Figure 2 | Global structure and dietary specificity of hepatic metabolic flux reprogramming "
                 "(standard FBA, primary)", fontsize=11, y=0.995, fontweight="bold")
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    print(f"[SAVED] {args.out}")


if __name__ == "__main__":
    main()

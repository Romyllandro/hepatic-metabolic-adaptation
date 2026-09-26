#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rq2_conservation_analysis.py

Recompute the RQ2 cross-strain conservation analysis using the new replicate-level
run, combining THREE independent filters that the original analysis did not have:

  1. MAGNITUDE   |Delta flux| > threshold            (the original criterion, alone)
  2. STATISTICS  replicate-based t-test, FDR-corrected (was impossible before:
                 the aggregate run had n=1 per group and every p-value was NaN)
  3. RESOLVABILITY |Delta flux| > pooled FVA range    (is the difference larger
                 than the width of the solution space it was drawn from?)

WHY ALL THREE
-------------
They fail in different ways and none is sufficient alone:

  * Magnitude alone cannot distinguish a real shift from the solver picking a
    different vertex of a degenerate optimum. On this model most reactions have a
    non-zero FVA range, so this is the common case.
  * Statistics alone says the three HFD animals differ from the three SCD animals,
    but with n=3 the smallest attainable p is small only for a t-test, and a
    difference can be statistically consistent yet still smaller than solver
    ambiguity.
  * Resolvability alone ignores biological variance entirely -- it bounds solver
    ambiguity, nothing else.

A reaction that passes all three is one where the animals differ consistently AND
the difference is bigger than the model's own indeterminacy. That is the set the
manuscript can defend at reaction level.

CONSERVATION TIERS ARE RECOMPUTED UNDER EACH DEFINITION
------------------------------------------------------
so you can see directly how much of the published "11 of 242 conserved across all
nine strains" survives each additional requirement. Expect the counts to fall.
That is the point of the exercise, not a failure of it.

EXPECTED LAYOUT (from run_layered_executor.py)
----------------------------------------------
    <rq2_root>/results_<STRAIN>_GSE182668/
        flux_analysis/reaction_flux_comparison_extended.csv
        stats_comparison/flux_pairwise_stats.csv
        fva/FVA_<GROUP>_rep<N>.txt

USAGE
-----
    python rq2_conservation_analysis.py \
        --rq2_root 972026_Step_2_RQ2 \
        --annotations RQ4_reaction_annotations.csv \
        --out_prefix RQ2_revised

    # match the original edge threshold exactly
    python rq2_conservation_analysis.py --rq2_root ... --abs_diff_threshold 0.2
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd


# ----------------------------------------------------------------------------- #
# FVA loading (mirrors merge_fva_into_flux.py so the two agree)
# ----------------------------------------------------------------------------- #
def read_fva_file(path):
    try:
        df = pd.read_csv(path, sep="\t", comment="#")
        if not {"Rxn", "Min", "Max"}.issubset(df.columns):
            raise ValueError
    except Exception:
        df = pd.read_csv(path, sep=r"\s+", comment="#", engine="python")
    df.columns = [str(c).strip() for c in df.columns]
    for c in ("Min", "Max"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["Min", "Max"]).drop_duplicates("Rxn").set_index("Rxn")
    df["Range"] = df["Max"] - df["Min"]
    return df[["Min", "Max", "Range"]]


def group_range(fva_dir, group, rxn_index):
    """Mean FVA range across that group's replicates."""
    pat = re.compile(rf"FVA_{re.escape(group)}_rep(\d+)\.(txt|csv|tsv)$", re.IGNORECASE)
    files = [p for p in sorted(glob.glob(os.path.join(fva_dir, "*")))
             if pat.search(os.path.basename(p))]
    if not files:
        return None, 0
    rngs = [read_fva_file(p).reindex(rxn_index)["Range"] for p in files]
    return pd.concat(rngs, axis=1).mean(axis=1), len(files)


def sample_cols(flux, group):
    return [c for c in flux.columns
            if c.startswith(group + "_") and c.endswith("_Flux")
            and "MeanFlux" not in c and "StdFlux" not in c and "SEMFlux" not in c]


# ----------------------------------------------------------------------------- #
def analyse_strain(strain, path, args):
    fpath = os.path.join(path, "flux_analysis", "reaction_flux_comparison_extended.csv")
    spath = os.path.join(path, "stats_comparison", "flux_pairwise_stats.csv")
    fva_dir = os.path.join(path, "fva")
    if not os.path.exists(fpath):
        print(f"[WARN] {strain}: no flux table at {fpath}; skipped")
        return None

    flux = pd.read_csv(fpath).set_index("ReactionID")
    idx = flux.index
    out = pd.DataFrame(index=idx)
    out["strain"] = strain

    tc, bc = f"{args.test_group}_MeanFlux", f"{args.baseline}_MeanFlux"
    if tc not in flux.columns or bc not in flux.columns:
        print(f"[WARN] {strain}: missing {tc}/{bc}; skipped")
        return None
    out["delta"] = pd.to_numeric(flux[tc], errors="coerce") - pd.to_numeric(flux[bc], errors="coerce")

    # 1. magnitude
    out["pass_magnitude"] = out["delta"].abs() > args.abs_diff_threshold

    # 2. statistics
    out["pvalue"] = np.nan
    out["fdr"] = np.nan
    if os.path.exists(spath):
        st = pd.read_csv(spath).drop_duplicates("ReactionID").set_index("ReactionID")
        for src, dst in (("PValue", "pvalue"), ("FDR_BH", "fdr")):
            if src in st.columns:
                out[dst] = pd.to_numeric(st[src], errors="coerce").reindex(idx)
    else:
        print(f"[WARN] {strain}: no pairwise stats at {spath}")
    out["pass_stats"] = out["fdr"] < args.fdr_alpha          # NaN -> False
    out["testable"] = out["pvalue"].notna()

    # 3. resolvability
    rt, nt = group_range(fva_dir, args.test_group, idx)
    rb, nb = group_range(fva_dir, args.baseline, idx)
    if rt is None or rb is None:
        print(f"[WARN] {strain}: FVA missing (test={nt}, base={nb} files); resolvability skipped")
        out["fva_pooled_range"] = np.nan
        out["pass_resolvable"] = False
        out["has_fva"] = False
    else:
        pooled = (rt + rb) / 2.0
        out["fva_pooled_range"] = pooled
        determined = pooled.abs() < args.zero_tol
        out["pass_resolvable"] = np.where(determined,
                                          out["delta"].abs() > args.zero_tol,
                                          out["delta"].abs() > pooled)
        out["has_fva"] = True

        # containment check: every point estimate inside its own interval
        chk = viol = 0
        for grp in (args.test_group, args.baseline):
            cols = sample_cols(flux, grp)
            pat = re.compile(rf"FVA_{re.escape(grp)}_rep(\d+)\.", re.IGNORECASE)
            fs = sorted([p for p in glob.glob(os.path.join(fva_dir, "*"))
                         if pat.search(os.path.basename(p))],
                        key=lambda p: int(pat.search(os.path.basename(p)).group(1)))
            for i, p in enumerate(fs):
                if i >= len(cols):
                    break
                fv = read_fva_file(p).reindex(idx)
                pt = pd.to_numeric(flux[cols[i]], errors="coerce")
                ok = pt.between(fv["Min"] - 1e-6, fv["Max"] + 1e-6)
                v = fv["Min"].notna() & pt.notna()
                chk += int(v.sum()); viol += int((~ok & v).sum())
        rate = 1.0 - (viol / chk if chk else 1.0)
        if rate < args.min_containment:
            print(f"[WARN] {strain}: FVA containment only {rate:.1%} -- replicate-to-sample "
                  f"mapping may be wrong for this strain. Treat its resolvability with caution.")
        else:
            print(f"[INFO] {strain}: FVA containment {rate:.1%} ({nt}+{nb} files)")

    out["pass_all"] = out["pass_magnitude"] & out["pass_stats"] & out["pass_resolvable"]
    return out.reset_index()


def tiers(sets_by_strain, n_strains, label):
    """Count reactions by how many strains they appear in."""
    from collections import Counter
    c = Counter()
    for s in sets_by_strain.values():
        for r in s:
            c[r] += 1
    dist = Counter(c.values())
    union = len(c)
    core = dist.get(n_strains, 0)
    partial = sum(dist.get(k, 0) for k in (4, 5, 6))
    high = sum(dist.get(k, 0) for k in (7, 8, 9))
    print(f"\n--- {label} ---")
    print(f"  union across strains : {union}")
    print(f"  core (all {n_strains})        : {core}")
    print(f"  high (7-9 strains)   : {high}")
    print(f"  partial (4-6)        : {partial}")
    print("  n_strains : count  " + "  ".join(f"{k}:{dist.get(k,0)}" for k in range(n_strains, 0, -1)))
    return {"criterion": label, "union": union, "core": core,
            "high_7_9": high, "partial_4_6": partial}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rq2_root", required=True)
    ap.add_argument("--annotations", default=None)
    ap.add_argument("--out_prefix", default="RQ2_revised")
    ap.add_argument("--test_group", default="HFD")
    ap.add_argument("--baseline", default="SCD")
    ap.add_argument("--abs_diff_threshold", type=float, default=0.2,
                    help="matches --edge_abs_diff_threshold from the production run")
    ap.add_argument("--fdr_alpha", type=float, default=0.05)
    ap.add_argument("--zero_tol", type=float, default=1e-9)
    ap.add_argument("--min_containment", type=float, default=0.95)
    args = ap.parse_args()

    dirs = sorted(glob.glob(os.path.join(args.rq2_root, "results_*")))
    if not dirs:
        sys.exit(f"[FATAL] no results_* directories under {args.rq2_root}")

    frames = []
    for d in dirs:
        m = re.search(r"results_(.+?)_GSE", os.path.basename(d))
        strain = m.group(1) if m else os.path.basename(d)
        r = analyse_strain(strain, d, args)
        if r is not None:
            frames.append(r)
    if not frames:
        sys.exit("[FATAL] no strain produced usable output.")

    allr = pd.concat(frames, ignore_index=True)
    allr.to_csv(f"{args.out_prefix}_per_strain_reactions.csv", index=False)
    print(f"\n[WRITE] {args.out_prefix}_per_strain_reactions.csv")

    strains = sorted(allr.strain.unique())
    n = len(strains)
    print(f"\n[INFO] {n} strains: {', '.join(strains)}")

    print("\n" + "=" * 78)
    print("PER-STRAIN COUNTS")
    print("=" * 78)
    per = allr.groupby("strain").agg(
        testable=("testable", "sum"),
        magnitude=("pass_magnitude", "sum"),
        fdr_sig=("pass_stats", "sum"),
        resolvable=("pass_resolvable", "sum"),
        all_three=("pass_all", "sum"),
    )
    print(per.to_string())

    print("\n" + "=" * 78)
    print("CONSERVATION TIERS UNDER EACH CRITERION")
    print("=" * 78)
    rows = []
    for col, lab in [("pass_magnitude", "A. magnitude only (original criterion)"),
                     ("pass_stats", "B. FDR-significant only"),
                     ("pass_resolvable", "C. resolvable vs FVA only"),
                     ("pass_all", "D. all three (defensible set)")]:
        sets = {s: set(allr.loc[(allr.strain == s) & (allr[col] == True), "ReactionID"])
                for s in strains}
        rows.append(tiers(sets, n, lab))
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}_tier_comparison.csv", index=False)
    print(f"\n[WRITE] {args.out_prefix}_tier_comparison.csv")

    # the defensible core, named
    core_sets = {s: set(allr.loc[(allr.strain == s) & (allr.pass_all == True), "ReactionID"])
                 for s in strains}
    from collections import Counter
    cc = Counter(r for s in core_sets.values() for r in s)
    if cc:
        top = pd.DataFrame(sorted(cc.items(), key=lambda kv: -kv[1]),
                           columns=["ReactionID", "n_strains"])
        if args.annotations and os.path.exists(args.annotations):
            ann = pd.read_csv(args.annotations)
            idc = "reaction_id" if "reaction_id" in ann.columns else ann.columns[0]
            keep = [c for c in ("reaction_name", "subsystem") if c in ann.columns]
            top = top.merge(ann.set_index(idc)[keep], left_on="ReactionID",
                            right_index=True, how="left")
        top.to_csv(f"{args.out_prefix}_defensible_reactions.csv", index=False)
        print(f"[WRITE] {args.out_prefix}_defensible_reactions.csv")
        print("\nReactions passing all three filters, most widely shared first:")
        print(top.head(25).to_string(index=False))
    else:
        print("\nNo reaction passed all three filters in any strain.")

    print("\nNOTE: criterion A reproduces the published analysis. B, C and D add the "
          "replicate statistics and the FVA width. Report D as the reaction-level "
          "result and A only as the historical comparison.")


if __name__ == "__main__":
    sys.exit(main())

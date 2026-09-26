#!/usr/bin/env python3
"""
rq1_signature_summary.py
========================
Summarise the primary RQ1 (bulk multi-cohort, standard FBA) differential results that
are reported in the manuscript (Results 3.2, Fig. 2, Supplementary Tables S1-S2).

Input : the six reaction_stats_<A>_vs_<B>.csv files written by
        comprehensive_flux_analysis.py (column `Significant` = BH-FDR < 0.05 on Welch tests).
Output: per-contrast counts, the union of the three SCD-referenced contrasts, the union
        across all estimable contrasts, and the cross-diet signature (intersection of
        HFD/KD/WD vs SCD) with per-diet direction.

Expected manuscript values: HFD 95, KD 154, WD 134, HFD-KD 32, HFD-WD 227
(KD-WD is computed by the script but is NOT estimable because no dataset contains both
diets; it is excluded from the 450 union), SCD-referenced union 262, all-contrast union 450,
signature 10 (9 directionally concordant; GLNALANaEx switches).
NOTE: the manuscript's "450 across all estimable pairwise comparisons" equals the union over
all SIX contrasts including KD_vs_WD; the union over the five estimable contrasts is 360.
Both values are reported so the text can be reconciled.

Usage
-----
python rq1_signature_summary.py --stats_dir <dir with reaction_stats_*.csv> --output <dir>
    [--prefix "reaction_stats_"]   # use "RQ1_reaction_stats_" for the frozen reference files
"""
import argparse, json
from pathlib import Path
import pandas as pd

SCD_REF = ["HFD_vs_SCD", "KD_vs_SCD", "WD_vs_SCD"]
ESTIMABLE = SCD_REF + ["HFD_vs_KD", "HFD_vs_WD"]   # KD_vs_WD: no dataset contains both diets


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stats_dir", required=True)
    ap.add_argument("--prefix", default="reaction_stats_")
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    d = Path(a.stats_dir); out = Path(a.output); out.mkdir(parents=True, exist_ok=True)

    sig, tabs, counts = {}, {}, {}
    for c in ESTIMABLE + ["KD_vs_WD"]:
        f = d / f"{a.prefix}{c}.csv"
        if not f.exists():
            continue
        t = pd.read_csv(f); tabs[c] = t
        s = set(t.loc[t["Significant"].astype(str).str.lower() == "true", "ReactionID"])
        sig[c] = s; counts[c] = len(s)

    union_scd = set().union(*[sig[c] for c in SCD_REF])
    union_all = set().union(*[sig[c] for c in ESTIMABLE if c in sig])
    core = set.intersection(*[sig[c] for c in SCD_REF])

    rows = []
    for rid in sorted(core):
        row = {"ReactionID": rid}
        for c in SCD_REF:
            r = tabs[c].set_index("ReactionID").loc[rid]
            row["ReactionName"] = r.get("ReactionName", ""); row["Subsystem"] = r.get("Subsystem", "")
            row[f"{c}_Cohen_d"] = r.get("Cohen_d"); row[f"{c}_q"] = r.get("q_value"); row[f"{c}_dir"] = r.get("Direction")
        dirs = {row[f"{c}_dir"] for c in SCD_REF}
        row["direction_concordant"] = len(dirs) == 1
        rows.append(row)
    sigdf = pd.DataFrame(rows)
    sigdf.to_csv(out / "RQ1_cross_diet_signature.csv", index=False)
    pd.DataFrame(sorted(union_scd), columns=["ReactionID"]).to_csv(out / "RQ1_union_SCD_referenced.csv", index=False)

    summary = {"significant_per_contrast": counts,
               "KD_vs_WD_note": "not estimable (no dataset contains both diets); excluded from union",
               "union_SCD_referenced": len(union_scd),
               "union_all_estimable_contrasts": len(union_all),
               "union_all_six_contrasts_including_KD_vs_WD": len(union_all | sig.get("KD_vs_WD", set())),
               "cross_diet_signature": len(core),
               "signature_directionally_concordant": int(sigdf["direction_concordant"].sum()) if len(sigdf) else 0,
               "signature_reactions": sorted(core)}
    (out / "RQ1_signature_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

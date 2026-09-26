#!/usr/bin/env python3
"""
rq2_conservation_tiers.py
=========================
Primary RQ2 magnitude-based conservation analysis (Results 3.3, Fig. 3,
Supplementary Tables S3-S5, SP5) from the nine strain-level standard-FBA runs
(aggregate design: expression averaged per strain x diet before FBA).

For every strain s and reaction j:  dv_sj = v_sj(HFD) - v_sj(SCD)  (column "Diff(HFD-SCD)")
A reaction is magnitude-responsive in strain s when |dv_sj| >= theta (primary theta = 0.20).
K_j(theta) = number of strains in which reaction j is responsive (Eq. 7).

Outputs (for the primary theta):
  RQ2_union_membership.csv   reaction x strain responsive matrix + K_j + tier
  RQ2_per_strain_summary.csv per-strain responsive / universal / shared / unique counts
  RQ2_universal_core.csv     signed dv for the universal (9/9) reactions + direction check
  RQ2_theta_sweep.csv        union and universal-core size for each theta
  RQ2_tiers_summary.json     headline numbers

Expected manuscript values at theta = 0.20: union 243, universal 11 (9 sign-consistent;
O2t and EX_o2_e reverse in C57BL6J), strain-unique 80, 7-8 strains 14, 4-6 strains 60,
per-strain 65-96 (mean 84.4); theta sweep 0.10/0.15/0.20/0.25/0.30 -> 14/13/11/9/4.

Usage
-----
# from a pipeline run (results_<STRAIN>_GSE182668/flux_analysis/reaction_flux_comparison_extended.csv)
python rq2_conservation_tiers.py --rq2_root results/RQ2_aggregate_FBA --output results/RQ2_tiers
# from the frozen reference files (RQ2_<STRAIN>_reaction_flux_comparison_extended.csv)
python rq2_conservation_tiers.py --flat_dir reference_outputs/RQ2 --output results/RQ2_tiers_check
"""
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

STRAINS = ["129S1SvImJ", "AJ", "C57BL6J", "CASTEiJ", "DBA2J", "NODShiLtJ", "NZOHlLtJ", "PWKPhJ", "WSBEiJ"]
DIFF = "Diff(HFD-SCD)"


def load(a):
    d = {}
    for s in STRAINS:
        if a.flat_dir:
            f = Path(a.flat_dir) / f"RQ2_{s}_reaction_flux_comparison_extended.csv"
        else:
            f = Path(a.rq2_root) / f"results_{s}_GSE182668" / "flux_analysis" / "reaction_flux_comparison_extended.csv"
        t = pd.read_csv(f, usecols=lambda c: c in ("ReactionID", "ReactionName", "Subsystem", DIFF))
        d[s] = t.set_index("ReactionID")
    return d


def tiers(diff, theta):
    resp = diff.abs() >= theta   # exact comparison, identical to the edge filter in the flux script (e.g. EX_cys__L_e at 0.19999999999999996 is excluded)
    k = resp.sum(axis=1)
    return resp, k


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--rq2_root"); g.add_argument("--flat_dir")
    ap.add_argument("--theta", type=float, default=0.20)
    ap.add_argument("--sweep", default="0.10,0.15,0.20,0.25,0.30")
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)

    d = load(a)
    ann = pd.concat([d[s][["ReactionName", "Subsystem"]] for s in STRAINS]).groupby(level=0).first()
    diff = pd.DataFrame({s: d[s][DIFF] for s in STRAINS}).fillna(0.0)

    resp, k = tiers(diff, a.theta)
    n = len(STRAINS)
    tier = pd.Series("not_responsive", index=diff.index)
    tier[k == 1] = "strain_unique"; tier[(k >= 2) & (k <= 3)] = "shared_2_3"
    tier[(k >= 4) & (k <= 6)] = "shared_4_6"; tier[(k >= 7) & (k <= 8)] = "shared_7_8"; tier[k == n] = "universal"
    union = k[k >= 1].index
    mem = resp.loc[union].astype(int).join(ann).assign(K=k.loc[union], tier=tier.loc[union])
    mem.sort_values(["K", "Subsystem"], ascending=[False, True]).to_csv(out / "RQ2_union_membership.csv")

    per = []
    for s in STRAINS:
        r = resp[s]
        per.append({"strain": s, "responsive": int(r.sum()),
                    "universal": int((r & (k == n)).sum()),
                    "shared_2_8": int((r & (k >= 2) & (k < n)).sum()),
                    "strain_unique": int((r & (k == 1)).sum())})
    per = pd.DataFrame(per); per.to_csv(out / "RQ2_per_strain_summary.csv", index=False)

    core = diff.loc[k[k == n].index]
    sign = np.sign(core)
    core_tab = core.join(ann)
    core_tab["sign_consistent"] = (sign.nunique(axis=1) == 1)
    core_tab["reversing_strains"] = [",".join(s for s in STRAINS if sign.loc[i, s] != sign.loc[i].mode().iloc[0])
                                     for i in core.index]
    core_tab.to_csv(out / "RQ2_universal_core.csv")

    sw = []
    for th in [float(x) for x in a.sweep.split(",")]:
        _, kk = tiers(diff, th)
        sw.append({"theta": th, "union": int((kk >= 1).sum()), "universal": int((kk == n).sum()),
                   "strain_unique": int((kk == 1).sum())})
    sw = pd.DataFrame(sw); sw.to_csv(out / "RQ2_theta_sweep.csv", index=False)

    summary = {"theta": a.theta, "union": int(len(union)), "universal": int((k == n).sum()),
               "universal_sign_consistent": int(core_tab["sign_consistent"].sum()),
               "universal_reactions": sorted(core.index.tolist()),
               "strain_unique": int((k == 1).sum()), "shared_7_8": int(((k >= 7) & (k <= 8)).sum()),
               "shared_4_6": int(((k >= 4) & (k <= 6)).sum()),
               "per_strain_min": int(per.responsive.min()), "per_strain_max": int(per.responsive.max()),
               "per_strain_mean": round(float(per.responsive.mean()), 1),
               "theta_sweep_universal": dict(zip(sw.theta.astype(str), sw.universal.astype(int)))}
    (out / "RQ2_tiers_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2)); print(per.to_string(index=False))


if __name__ == "__main__":
    main()

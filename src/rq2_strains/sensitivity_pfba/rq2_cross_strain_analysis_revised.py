#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cross-strain RQ2 analysis for replicate-preserving pFBA outputs.

Separates:
1) magnitude-defined conservation: |mean HFD - mean SCD| >= theta
2) within-strain statistical support: BH-FDR + effect-size criterion
3) FVA robustness of representative HFD/SCD models, when available

This avoids calling a reaction "significant in all strains" merely because it
crosses a fixed magnitude threshold.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def bh(p):
    from statsmodels.stats.multitest import multipletests
    p=np.asarray(p,float); q=np.full_like(p,np.nan)
    ok=np.isfinite(p)
    if ok.any(): q[ok]=multipletests(p[ok],method="fdr_bh")[1]
    return q

def load_config(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))

def read_strain(job,root):
    rdir=root/Path(job["results_dir"])
    f=rdir/"stats_comparison"/"flux_pairwise_stats.csv"
    if not f.exists(): raise FileNotFoundError(f)
    df=pd.read_csv(f)
    # Retain the HFD-SCD row regardless of internal GroupA/B ordering.
    mask=((df.GroupA=="HFD")&(df.GroupB=="SCD"))|((df.GroupA=="SCD")&(df.GroupB=="HFD"))
    df=df.loc[mask].copy()
    if df.empty: raise ValueError(f"No HFD/SCD contrast in {f}")
    sign=np.where((df.GroupA=="HFD")&(df.GroupB=="SCD"),1.0,-1.0)
    df["HFD_minus_SCD"]=sign*pd.to_numeric(df["Diff(A-B)"],errors="coerce")
    df["Cohen_d_HFD_vs_SCD"]=sign*pd.to_numeric(df["Cohen_d"],errors="coerce")
    df["q_value"]=pd.to_numeric(df["FDR_BH"],errors="coerce")
    df["Strain"]=job["name"]
    return df,rdir

def fva_diff(rdir):
    h=rdir/"fva"/"FVA_HFD_representative.csv"
    s=rdir/"fva"/"FVA_SCD_representative.csv"
    if not h.exists() or not s.exists(): return None
    H=pd.read_csv(h).set_index("ReactionID")
    S=pd.read_csv(s).set_index("ReactionID")
    idx=H.index.intersection(S.index)
    out=pd.DataFrame(index=idx)
    out["FVA_DiffMin"]=H.loc[idx,"minimum"].to_numpy()-S.loc[idx,"maximum"].to_numpy()
    out["FVA_DiffMax"]=H.loc[idx,"maximum"].to_numpy()-S.loc[idx,"minimum"].to_numpy()
    out["FVA_RobustDirection"]=np.where(out["FVA_DiffMin"]>0,"Up",
                                 np.where(out["FVA_DiffMax"]<0,"Down","Ambiguous"))
    return out.reset_index()

def main():
    ap=argparse.ArgumentParser(description="Combine nine RQ2 strain runs.")
    ap.add_argument("--config",required=True)
    ap.add_argument("--output",default="RQ2_cross_strain_revision")
    ap.add_argument("--theta",type=float,default=0.20)
    ap.add_argument("--thresholds",default="0.05,0.10,0.15,0.20,0.25,0.30,0.50")
    ap.add_argument("--fdr",type=float,default=0.05)
    ap.add_argument("--effect_size",type=float,default=0.5)
    a=ap.parse_args()
    cfg=load_config(a.config)
    root=Path(cfg.get("runner",{}).get("project_root",".")).resolve()
    out=Path(a.output); out.mkdir(parents=True,exist_ok=True)

    frames=[]; fva_frames=[]
    for job in cfg["jobs"]:
        df,rdir=read_strain(job,root); frames.append(df)
        fv=fva_diff(rdir)
        if fv is not None:
            fv["Strain"]=job["name"]; fva_frames.append(fv)
    long=pd.concat(frames,ignore_index=True)
    long.to_csv(out/"RQ2_all_strain_reaction_stats_long.csv",index=False)
    strains=[j["name"] for j in cfg["jobs"]]
    nstr=len(strains)

    # Reaction metadata and matrices.
    meta=long.drop_duplicates("ReactionID").set_index("ReactionID")[["ReactionName","Subsystem"]]
    diff=long.pivot(index="ReactionID",columns="Strain",values="HFD_minus_SCD").reindex(columns=strains)
    q=long.pivot(index="ReactionID",columns="Strain",values="q_value").reindex(columns=strains)
    d=long.pivot(index="ReactionID",columns="Strain",values="Cohen_d_HFD_vs_SCD").reindex(columns=strains)

    responsive=diff.abs()>=a.theta
    stat_supported=(q<=a.fdr)&(d.abs()>=a.effect_size)

    summary=meta.copy()
    summary["N_magnitude_responsive"]=responsive.sum(axis=1)
    summary["N_statistically_supported"]=stat_supported.sum(axis=1)
    summary["MeanDelta"]=diff.mean(axis=1)
    summary["MedianAbsDelta"]=diff.abs().median(axis=1)
    summary["DirectionAgreement"]=np.maximum((diff>0).mean(axis=1),(diff<0).mean(axis=1))
    summary["MagnitudeUniversal"]=summary["N_magnitude_responsive"]==nstr
    summary["StatisticallyUniversal"]=summary["N_statistically_supported"]==nstr
    summary["ConservationTier"]=pd.cut(
        summary["N_magnitude_responsive"],bins=[-1,1,3,6,8,9],
        labels=["low/strain-specific","partial","moderate","high","universal"]
    )
    summary.reset_index().to_csv(out/"RQ2_cross_strain_reaction_summary.csv",index=False)

    # Per-strain counts and unique responses.
    per=[]
    for s in strains:
        mag=responsive[s]
        sig=stat_supported[s]
        unique=mag & (responsive.sum(axis=1)==1)
        per.append({"Strain":s,"N_magnitude_responsive":int(mag.sum()),
                    "N_statistically_supported":int(sig.sum()),
                    "N_strain_unique_magnitude":int(unique.sum())})
    pd.DataFrame(per).to_csv(out/"RQ2_per_strain_summary.csv",index=False)

    # Threshold sensitivity.
    th=[float(x) for x in a.thresholds.split(",") if x.strip()]
    sens=[]
    for t in th:
        R=diff.abs()>=t
        counts=R.sum(axis=1)
        sens.append({"theta":t,"N_union":int((counts>=1).sum()),
                     "N_universal":int((counts==nstr).sum()),
                     "N_high_7to9":int((counts>=7).sum()),
                     "N_unique":int((counts==1).sum())})
    sens=pd.DataFrame(sens)
    sens.to_csv(out/"RQ2_threshold_sensitivity.csv",index=False)
    fig,ax=plt.subplots(figsize=(7,5))
    ax.plot(sens.theta,sens.N_universal,marker="o",label="Universal (9/9)")
    ax.plot(sens.theta,sens.N_union,marker="o",label="Union (>=1)")
    ax.set_xlabel("Magnitude threshold θ")
    ax.set_ylabel("Number of reactions")
    ax.legend(); fig.tight_layout(); fig.savefig(out/"RQ2_threshold_sensitivity.png",dpi=300); plt.close(fig)

    # Rank product across strains using normalized ranks of |delta|.
    ranks=diff.abs().rank(axis=0,ascending=False,method="average")
    norm=ranks.divide(ranks.count(axis=0),axis=1)
    rp=np.exp(np.log(norm.clip(lower=1e-12)).mean(axis=1))
    rpdf=pd.DataFrame({"ReactionID":rp.index,"RankProduct":rp.values}).sort_values("RankProduct")
    rpdf["RankProductRank"]=np.arange(1,len(rpdf)+1)
    rpdf["DirectionAgreement"]=summary.loc[rpdf.ReactionID,"DirectionAgreement"].to_numpy()
    rpdf["MedianAbsDelta"]=summary.loc[rpdf.ReactionID,"MedianAbsDelta"].to_numpy()
    rpdf.to_csv(out/"RQ2_rank_product_across_strains.csv",index=False)

    # Universal magnitude core.
    core=summary[summary["MagnitudeUniversal"]].copy().reset_index()
    core=core.sort_values(["DirectionAgreement","MedianAbsDelta"],ascending=False)
    core.to_csv(out/"RQ2_universal_magnitude_core.csv",index=False)

    # FVA robustness, if representative FVA was run.
    if fva_frames:
        fvlong=pd.concat(fva_frames,ignore_index=True)
        fvlong.to_csv(out/"RQ2_FVA_all_strains_long.csv",index=False)
        fvs=fvlong.groupby("ReactionID").agg(
            N_FVA_strains=("Strain","nunique"),
            N_FVA_robust=("FVA_RobustDirection",lambda x:int((x!="Ambiguous").sum())),
            N_FVA_ambiguous=("FVA_RobustDirection",lambda x:int((x=="Ambiguous").sum()))
        ).reset_index()
        fvs.to_csv(out/"RQ2_FVA_robustness_summary.csv",index=False)

    report={
        "theta":a.theta,"fdr":a.fdr,"effect_size":a.effect_size,
        "n_strains":nstr,
        "n_union_magnitude":int((summary.N_magnitude_responsive>=1).sum()),
        "n_universal_magnitude":int(summary.MagnitudeUniversal.sum()),
        "n_union_statistically_supported":int((summary.N_statistically_supported>=1).sum()),
        "n_universal_statistically_supported":int(summary.StatisticallyUniversal.sum()),
        "interpretation":"Conservation is defined by modeled effect magnitude; replicate-level FDR support is reported separately."
    }
    (out/"RQ2_summary.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))
    return 0
if __name__=="__main__":
    raise SystemExit(main())

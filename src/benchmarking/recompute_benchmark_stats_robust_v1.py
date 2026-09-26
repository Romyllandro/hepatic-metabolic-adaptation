#!/usr/bin/env python3
"""
recompute_benchmark_stats_robust_v1.py
--------------------------------------
Post-hoc statistical cleanup for benchmark_liver_methods.py outputs.

Purpose
-------
The optimization runs are valid, but many reactions are exactly/near-constant
within one or both diet groups. scipy.stats.ttest_ind emits precision-loss
warnings for such reactions. This script recomputes Welch statistics explicitly
and labels degenerate-variance cases instead of treating them as ordinary tests.

No metabolic optimization is rerun.

Outputs
-------
- reaction_stats_robust_<method>_<objective>.csv
- degeneracy_summary.csv
- significant_set_overlap_robust.csv

Interpretation
--------------
- both_constant_equal: p=1, Cohen_d=0
- both_constant_different: p/q left NA; effect is structural/deterministic and
  should not be treated as replicate-level inferential evidence
- ordinary / one_group_constant: explicit Welch t test is computed
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t as student_t


def bh(p):
    p=np.asarray(p,float)
    q=np.full_like(p,np.nan)
    ok=np.isfinite(p)
    vals=p[ok]
    if not len(vals):
        return q
    order=np.argsort(vals)
    ranked=vals[order]
    adj=ranked*len(ranked)/np.arange(1,len(ranked)+1)
    adj=np.minimum.accumulate(adj[::-1])[::-1]
    out=np.empty_like(vals)
    out[order]=np.minimum(adj,1.0)
    q[np.where(ok)[0]]=out
    return q


def cohen_d(a,b,var_tol=1e-14):
    a=np.asarray(a,float); b=np.asarray(b,float)
    a=a[np.isfinite(a)]; b=b[np.isfinite(b)]
    if len(a)<2 or len(b)<2:
        return np.nan
    va=np.var(a,ddof=1); vb=np.var(b,ddof=1)
    pooled=((len(a)-1)*va+(len(b)-1)*vb)/(len(a)+len(b)-2)
    md=float(np.mean(a)-np.mean(b))
    if pooled <= var_tol:
        return 0.0 if abs(md) <= np.sqrt(var_tol) else np.nan
    return md/np.sqrt(pooled)


def welch_explicit(a,b,var_tol=1e-14,mean_tol=1e-12):
    a=np.asarray(a,float); b=np.asarray(b,float)
    a=a[np.isfinite(a)]; b=b[np.isfinite(b)]
    if len(a)<2 or len(b)<2:
        return np.nan,np.nan,"insufficient_n"

    ma=float(np.mean(a)); mb=float(np.mean(b))
    va=float(np.var(a,ddof=1)); vb=float(np.var(b,ddof=1))
    ca=va <= var_tol
    cb=vb <= var_tol

    if ca and cb:
        if abs(ma-mb) <= mean_tol:
            return 0.0,1.0,"both_constant_equal"
        return np.nan,np.nan,"both_constant_different"

    se2=va/len(a)+vb/len(b)
    if se2 <= var_tol:
        return np.nan,np.nan,"degenerate_se"

    stat=(ma-mb)/np.sqrt(se2)
    num=se2**2
    den=0.0
    if va>var_tol:
        den += (va/len(a))**2/(len(a)-1)
    if vb>var_tol:
        den += (vb/len(b))**2/(len(b)-1)
    if den <= 0:
        return np.nan,np.nan,"degenerate_df"
    df=num/den
    p=2*student_t.sf(abs(stat),df)

    case="one_group_constant" if (ca or cb) else "ordinary"
    return float(stat),float(p),case


def stats_from_flux(df, test, baseline, var_tol, mean_tol):
    tcols=[c for c in df.columns if str(c).startswith(test+"_")]
    bcols=[c for c in df.columns if str(c).startswith(baseline+"_")]
    if not tcols or not bcols:
        raise ValueError(
            f"Could not identify columns for {test}/{baseline}. "
            f"Example columns: {list(df.columns[:8])}"
        )

    rows=[]
    for rid,row in df.iterrows():
        a=pd.to_numeric(row[tcols],errors="coerce").to_numpy(float)
        b=pd.to_numeric(row[bcols],errors="coerce").to_numpy(float)
        aa=a[np.isfinite(a)]; bb=b[np.isfinite(b)]
        ma=float(np.mean(aa)) if len(aa) else np.nan
        mb=float(np.mean(bb)) if len(bb) else np.nan
        stat,p,case=welch_explicit(aa,bb,var_tol=var_tol,mean_tol=mean_tol)
        rows.append({
            "ReactionID":rid,
            f"{test}_mean":ma,
            f"{baseline}_mean":mb,
            "MeanDiff":ma-mb if np.isfinite(ma) and np.isfinite(mb) else np.nan,
            "Welch_t":stat,
            "Cohen_d":cohen_d(aa,bb,var_tol=var_tol),
            "p_value":p,
            "variance_case":case,
            f"{test}_sd":float(np.std(aa,ddof=1)) if len(aa)>=2 else np.nan,
            f"{baseline}_sd":float(np.std(bb,ddof=1)) if len(bb)>=2 else np.nan,
        })
    out=pd.DataFrame(rows)
    out["q_value"]=bh(out["p_value"].to_numpy(float))
    out["Significant"]=(out["q_value"]<0.05)&(out["Cohen_d"].abs()>=0.5)
    out["Direction"]=np.where(out.MeanDiff>0,"Up",np.where(out.MeanDiff<0,"Down","Zero"))
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--benchmark_dir",required=True)
    ap.add_argument("--test",default="HFD")
    ap.add_argument("--baseline",default="SCD")
    ap.add_argument("--methods",default="three_layer,eflux,gimme,pfba")
    ap.add_argument("--objectives",default="BIOMASS_mm_1_no_glygln,ATPM")
    ap.add_argument("--var_tol",type=float,default=1e-14)
    ap.add_argument("--mean_tol",type=float,default=1e-12)
    args=ap.parse_args()

    root=Path(args.benchmark_dir)
    methods=[x.strip() for x in args.methods.split(",") if x.strip()]
    objectives=[x.strip() for x in args.objectives.split(",") if x.strip()]

    robust={}
    diag=[]
    for obj in objectives:
        for method in methods:
            f=root/f"fluxes_{method}_{obj}.csv"
            if not f.exists():
                print(f"[WARN] missing {f}; skipping")
                continue
            flux=pd.read_csv(f,index_col=0)
            st=stats_from_flux(flux,args.test,args.baseline,args.var_tol,args.mean_tol)
            out=root/f"reaction_stats_robust_{method}_{obj}.csv"
            st.to_csv(out,index=False)
            robust[(method,obj)]=st
            vc=st["variance_case"].value_counts()
            for case,n in vc.items():
                diag.append({"method":method,"objective":obj,"variance_case":case,"n_reactions":int(n)})
            print(f"[OK] {out}")

    pd.DataFrame(diag).to_csv(root/"degeneracy_summary.csv",index=False)

    overlaps=[]
    for obj in objectives:
        key=("three_layer",obj)
        if key not in robust:
            continue
        ref=robust[key].set_index("ReactionID")
        ref_sig=set(ref.index[ref["Significant"].fillna(False)])
        for method in methods:
            if method=="three_layer" or (method,obj) not in robust:
                continue
            q=robust[(method,obj)].set_index("ReactionID")
            sig=set(q.index[q["Significant"].fillna(False)])
            overlaps.append({
                "objective":obj,
                "method":method,
                "reference":"three_layer",
                "reference_n_significant":len(ref_sig),
                "method_n_significant":len(sig),
                "intersection":len(ref_sig & sig),
                "union":len(ref_sig | sig),
                "jaccard":len(ref_sig & sig)/len(ref_sig | sig) if (ref_sig|sig) else np.nan,
            })
    pd.DataFrame(overlaps).to_csv(root/"significant_set_overlap_robust.csv",index=False)
    print(f"[DONE] robust statistics written under {root}")


if __name__=="__main__":
    main()

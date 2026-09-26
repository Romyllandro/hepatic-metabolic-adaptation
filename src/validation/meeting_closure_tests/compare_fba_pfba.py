#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

def response(df, test="HFD", baseline="SCD"):
    diff=f"Diff({test}-{baseline})"
    if diff in df.columns:
        return pd.to_numeric(df.set_index("ReactionID")[diff], errors="coerce")
    t=f"{test}_MeanFlux"; b=f"{baseline}_MeanFlux"
    if t in df.columns and b in df.columns:
        x=df.set_index("ReactionID")
        return pd.to_numeric(x[t],errors="coerce")-pd.to_numeric(x[b],errors="coerce")
    raise ValueError(f"Cannot find {diff} or {t}/{b} columns.")

def top_set(s,n):
    return set(s.abs().sort_values(ascending=False).head(min(n,len(s))).index)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--fba",required=True)
    ap.add_argument("--pfba",required=True)
    ap.add_argument("--test",default="HFD")
    ap.add_argument("--baseline",default="SCD")
    ap.add_argument("--top_n",type=int,default=100)
    ap.add_argument("--tol",type=float,default=1e-9)
    ap.add_argument("--output",default="reviewer_results/meeting_closure/fba_vs_pfba")
    a=ap.parse_args()
    out=Path(a.output); out.mkdir(parents=True,exist_ok=True)
    A=pd.read_csv(a.fba); B=pd.read_csv(a.pfba)
    ra=response(A,a.test,a.baseline); rb=response(B,a.test,a.baseline)
    idx=ra.index.intersection(rb.index)
    ra=ra.loc[idx]; rb=rb.loc[idx]
    ok=ra.notna()&rb.notna()
    ra=ra[ok]; rb=rb[ok]
    union=(ra.abs()>a.tol)|(rb.abs()>a.tol)
    both=(ra.abs()>a.tol)&(rb.abs()>a.tol)
    def rho(x,y):
        if len(x)<3:return np.nan
        return float(spearmanr(x,y).correlation)
    topa=top_set(ra,a.top_n); topb=top_set(rb,a.top_n)
    result={
        "n_reactions_compared":int(len(ra)),
        "n_nontrivial_union":int(union.sum()),
        "n_nontrivial_both":int(both.sum()),
        "response_spearman_all":rho(ra,rb),
        "response_spearman_nontrivial_union":rho(ra[union],rb[union]),
        "response_spearman_both_active":rho(ra[both],rb[both]),
        "absolute_magnitude_spearman_both_active":rho(ra[both].abs(),rb[both].abs()),
        "directional_agreement_union":float((np.sign(ra[union])==np.sign(rb[union])).mean()) if union.any() else np.nan,
        "directional_agreement_both_active":float((np.sign(ra[both])==np.sign(rb[both])).mean()) if both.any() else np.nan,
        "top_n":a.top_n,
        "top_n_jaccard":len(topa&topb)/len(topa|topb) if topa|topb else np.nan,
        "fba_total_L1_response":float(ra.abs().sum()),
        "pfba_total_L1_response":float(rb.abs().sum()),
        "pfba_to_fba_L1_ratio":float(rb.abs().sum()/ra.abs().sum()) if ra.abs().sum()>0 else np.nan,
        "median_absolute_response_difference":float(np.median(np.abs(rb-ra))),
        "n_direction_flips_both_active":int(((np.sign(ra[both])!=np.sign(rb[both]))).sum()),
    }
    pd.DataFrame({"ReactionID":idx, "FBA_response":ra.reindex(idx).values,
                  "pFBA_response":rb.reindex(idx).values,
                  "absolute_difference":(rb-ra).abs().reindex(idx).values}).to_csv(
                      out/"reaction_level_FBA_vs_pFBA.csv",index=False)
    (out/"FBA_vs_pFBA_summary.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
if __name__=="__main__":
    main()

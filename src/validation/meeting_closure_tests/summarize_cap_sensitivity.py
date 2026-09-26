#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from pathlib import Path
import pandas as pd, numpy as np
from scipy.stats import spearmanr

def response(df,test,baseline):
    t=[c for c in df.columns if c.startswith(test+"_")]
    b=[c for c in df.columns if c.startswith(baseline+"_")]
    return df[t].mean(axis=1)-df[b].mean(axis=1)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--dir",required=True)
    ap.add_argument("--primary_cap",type=float,default=1000)
    ap.add_argument("--test",default="HFD")
    ap.add_argument("--baseline",default="SCD")
    a=ap.parse_args()
    p=Path(a.dir)
    files=list(p.glob("flux_q*_floor*_cap*_*.csv"))
    rows=[]
    by={}
    pat=re.compile(r"flux_q([^_]+)_floor([^_]+)_cap([^_]+)_(.+)\.csv$")
    for f in files:
        m=pat.match(f.name)
        if not m: continue
        q=float(m.group(1)); fl=float(m.group(2)); cap=float(m.group(3)); obj=m.group(4)
        if abs(q-.95)>1e-12 or abs(fl-.1)>1e-12: continue
        d=pd.read_csv(f,index_col=0)
        r=response(d,a.test,a.baseline)
        by[(cap,obj)]=r
    for (cap,obj),r in by.items():
        ref=by.get((a.primary_cap,obj))
        if ref is None: continue
        idx=r.index.intersection(ref.index)
        x=ref.loc[idx]; y=r.loc[idx]
        union=(x.abs()>1e-9)|(y.abs()>1e-9)
        rho=float(spearmanr(x[union],y[union]).correlation) if union.sum()>=3 else np.nan
        direction=float((np.sign(x[union])==np.sign(y[union])).mean()) if union.any() else np.nan
        ta=set(x.abs().nlargest(100).index); tb=set(y.abs().nlargest(100).index)
        rows.append({"objective":obj,"cap":cap,"reference_cap":a.primary_cap,
                     "response_spearman":rho,"directional_agreement":direction,
                     "top100_jaccard":len(ta&tb)/len(ta|tb),
                     "total_L1_response_ratio_vs_primary":float(y.abs().sum()/x.abs().sum()) if x.abs().sum()>0 else np.nan})
    out=pd.DataFrame(rows).sort_values(["objective","cap"])
    out.to_csv(p/"cap_sensitivity_summary.csv",index=False)
    print(out.to_string(index=False))
if __name__=="__main__":
    main()

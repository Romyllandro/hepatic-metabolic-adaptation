#!/usr/bin/env python3
"""
run_rq3_sensitivity_suite.py
============================
Targeted RQ3 sensitivity analysis requested by Reviewer 2.

This does NOT invent biological replication. RQ3 remains one deterministic
pseudo-bulk model per cell type × condition. The purpose is to test whether the
reported cell-type hierarchy and pathway attribution are stable to:
1. normalization percentile (P95 vs P99 by default),
2. biological objective (biomass-associated demand vs ATPM),
3. post-hoc response threshold (|Δv| = 0.05, 0.10, 0.20).

The underlying RQ3 modeling script is run unchanged for the optimization step;
the suite then re-summarizes its deterministic reaction-level outputs.

Recommended primary reference:
P99 + per_condition normalization + biomass objective.
"""
from __future__ import annotations
import argparse, itertools, json, os, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

def run(cmd,log,dry=False):
    print(" ".join(map(str,cmd)))
    if dry:return 0
    with open(log,"w",encoding="utf-8") as f:
        p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,text=True)
    if p.returncode:
        tail=Path(log).read_text(errors="replace").splitlines()[-80:]
        print("\n".join(tail))
        raise RuntimeError(f"Command failed: {cmd}")
    return p.returncode

def threshold_mask(df,fold_thr,abs_thr,rel_thr):
    # Mirrors RQ3 perform_statistical_tests: responsive if >=2 of 3 criteria.
    a=(pd.to_numeric(df["abs_fold_change"],errors="coerce")>fold_thr)
    b=(pd.to_numeric(df["abs_flux_change"],errors="coerce")>abs_thr)
    c=(pd.to_numeric(df["relative_magnitude"],errors="coerce").abs()>rel_thr)
    return (a.astype(int)+b.astype(int)+c.astype(int))>=2

def abundance_from_aggregation(agg):
    q=agg.groupby("cell_type",as_index=False)["n_cells"].sum()
    total=q.n_cells.sum()
    return dict(zip(q.cell_type,q.n_cells/total))

def contribution_table(stats,agg,comparison,fold_thr=1.5,abs_thr=0.1,rel_thr=0.5):
    s=stats[stats.comparison==comparison].copy()
    s["responsive"]=threshold_mask(s,fold_thr,abs_thr,rel_thr)
    abund=abundance_from_aggregation(agg)
    rows=[]
    for ct,sub in s.groupby("cell_type"):
        sig=sub[sub.responsive]
        n=len(sig); total=len(sub)
        rate=n/total if total else 0
        mean_change=float(pd.to_numeric(sig.abs_flux_change,errors="coerce").mean()) if n else 0.0
        abundance=float(abund.get(ct,0.0))
        score=abundance*rate*mean_change
        rows.append({"cell_type":ct,"abundance":abundance,"n_responsive":n,
                     "n_total":total,"response_rate":rate,
                     "mean_abs_flux_change":mean_change,"contribution_score":score})
    out=pd.DataFrame(rows).sort_values("contribution_score",ascending=False)
    denom=out.contribution_score.sum()
    out["contribution_percent"]=100*out.contribution_score/denom if denom else 0
    out["rank"]=np.arange(1,len(out)+1)
    return out,s

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--rq3_script",required=True)
    ap.add_argument("--sc_data",required=True)
    ap.add_argument("--sc_metadata",required=True)
    ap.add_argument("--diet_bounds_file",required=True)
    ap.add_argument("--condition_mapping",required=True,
                    help='JSON string, e.g. {"Chow":"SCD","WesternDiet":"WD"}')
    ap.add_argument("--model",default="iMM1415.json")
    ap.add_argument("--quantiles",default="0.95,0.99")
    ap.add_argument("--objectives",default="BIOMASS_mm_1_no_glygln,ATPM")
    ap.add_argument("--normalization_strategy",default="per_condition")
    ap.add_argument("--eflux_floor",type=float,default=0.001)
    ap.add_argument("--eflux_cap",type=float,default=10000.0)
    ap.add_argument("--baseline",default="Chow")
    ap.add_argument("--test_condition",default="WesternDiet")
    ap.add_argument("--fold_threshold",type=float,default=1.5)
    ap.add_argument("--relative_threshold",type=float,default=0.5)
    ap.add_argument("--abs_thresholds",default="0.05,0.10,0.20")
    ap.add_argument("--min_cells",type=int,default=20)
    ap.add_argument("--output",default="RQ3_sensitivity")
    ap.add_argument("--dry-run",action="store_true")
    args=ap.parse_args()

    out=Path(args.output); out.mkdir(parents=True,exist_ok=True)
    qs=[float(x) for x in args.quantiles.split(",")]
    objs=[x.strip() for x in args.objectives.split(",")]
    abs_thrs=[float(x) for x in args.abs_thresholds.split(",")]
    settings=[]
    for q,obj in itertools.product(qs,objs):
        label=f"q{int(round(q*100))}_{obj}"
        rdir=out/label
        cmd=[sys.executable,args.rq3_script,
             "--sc_data",args.sc_data,
             "--sc_metadata",args.sc_metadata,
             "--diet_bounds_file",args.diet_bounds_file,
             "--condition_mapping",args.condition_mapping,
             "--model",args.model,
             "--objective",obj,
             "--eflux_quantile",str(q),
             "--eflux_floor",str(args.eflux_floor),
             "--eflux_cap",str(args.eflux_cap),
             "--normalization_strategy",args.normalization_strategy,
             "--baseline",args.baseline,
             "--test_conditions",args.test_condition,
             "--fold_change_threshold",str(args.fold_threshold),
             "--abs_change_threshold","0.1",
             "--min_cells",str(args.min_cells),
             "--results_dir",str(rdir)]
        run(cmd,out/f"{label}.log",args.dry_run)
        settings.append({"label":label,"quantile":q,"objective":obj,"dir":str(rdir)})
    if args.dry_run:
        pd.DataFrame(settings).to_csv(out/"settings.csv",index=False)
        return

    comparison=f"{args.test_condition}_vs_{args.baseline}"
    contribution_all=[]; threshold_rows=[]; pathway_rows=[]
    setting_tables={}
    for st in settings:
        rdir=Path(st["dir"])
        stats=pd.read_csv(rdir/"statistics/statistical_tests.csv")
        agg=pd.read_csv(rdir/"aggregation_summary.csv")
        # Sensitivity across absolute thresholds with all other criteria fixed.
        for athr in abs_thrs:
            q=stats[stats.comparison==comparison].copy()
            q["responsive"]=threshold_mask(q,args.fold_threshold,athr,args.relative_threshold)
            threshold_rows.append({
                **{k:st[k] for k in ("label","quantile","objective")},
                "abs_threshold":athr,
                "n_responsive_rows":int(q.responsive.sum()),
                "n_unique_reactions":int(q.loc[q.responsive,"reaction_id"].nunique()),
                "n_cell_types_with_response":int(q.loc[q.responsive,"cell_type"].nunique())
            })
        contrib,marked=contribution_table(
            stats,agg,comparison,args.fold_threshold,0.1,args.relative_threshold
        )
        contrib["setting"]=st["label"]; contrib["quantile"]=st["quantile"]; contrib["objective"]=st["objective"]
        contribution_all.append(contrib)
        setting_tables[st["label"]]=contrib.set_index("cell_type")
        # pathway counts at the primary |Δv| threshold
        sig=marked[marked.responsive]
        p=(sig.groupby("subsystem")
           .agg(n_responsive_rows=("reaction_id","size"),
                n_unique_reactions=("reaction_id","nunique"),
                n_cell_types=("cell_type","nunique"))
           .reset_index())
        p["setting"]=st["label"]; pathway_rows.append(p)

    contrib_all=pd.concat(contribution_all,ignore_index=True)
    contrib_all.to_csv(out/"rq3_contribution_sensitivity.csv",index=False)
    pd.DataFrame(threshold_rows).to_csv(out/"rq3_threshold_sensitivity.csv",index=False)
    pd.concat(pathway_rows,ignore_index=True).to_csv(out/"rq3_pathway_sensitivity.csv",index=False)

    # Compare all settings to the manuscript reference P99+biomass when available.
    ref_label=f"q99_{objs[0]}"
    if ref_label not in setting_tables:
        ref_label=settings[0]["label"]
    ref=setting_tables[ref_label]
    rows=[]
    key_cts=["LECs","qHSCs","aHSCs","Hepatocytes"]
    for st in settings:
        lab=st["label"]; q=setting_tables[lab]
        idx=sorted(set(ref.index)|set(q.index))
        a=ref.reindex(idx).contribution_percent.fillna(0)
        b=q.reindex(idx).contribution_percent.fillna(0)
        rho=float(spearmanr(a,b).correlation)
        top5_ref=set(ref.sort_values("rank").head(5).index)
        top5=set(q.sort_values("rank").head(5).index)
        row={
            "reference_setting":ref_label,"setting":lab,
            "quantile":st["quantile"],"objective":st["objective"],
            "contribution_spearman":rho,
            "top5_celltype_jaccard":len(top5_ref&top5)/len(top5_ref|top5)
        }
        for ct in key_cts:
            row[f"{ct}_rank"]=float(q.loc[ct,"rank"]) if ct in q.index else np.nan
            row[f"{ct}_contribution_percent"]=float(q.loc[ct,"contribution_percent"]) if ct in q.index else np.nan
        rows.append(row)
    pd.DataFrame(rows).to_csv(out/"rq3_setting_concordance.csv",index=False)
    pd.DataFrame(settings).to_csv(out/"settings.csv",index=False)
    print(f"[DONE] {out}")
    print(pd.DataFrame(rows).to_string(index=False))

if __name__=="__main__":
    main()

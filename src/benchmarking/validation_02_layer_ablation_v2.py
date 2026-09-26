#!/usr/bin/env python3
"""
validation_02_layer_ablation_v2.py
==================================
Current-production layer ablation for the revised pFBA pipeline.

Variants
--------
full       L1 diet + L2 transport + L3 intracellular
vanilla    conventional single-scope E-Flux on all non-exchange reactions, no diet
no_diet    L2 + L3
no_inter   L1 + L2
no_trans   L1 + L3

Unlike the historical script, all variants use pFBA at fraction_of_optimum=1.0.
For a clean methodological comparison, use one within-study cohort (default GSE101657).
"""
from __future__ import annotations
import argparse, importlib.util, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def condition_means(rna,dataset,conditions):
    sample_cols=list(rna.columns[2:])
    out={}; counts={}
    for c in conditions:
        hits=[x for x in sample_cols if x.startswith(f"{c}_{dataset}_")] if dataset else [x for x in sample_cols if x.startswith(c+"_")]
        if not hits and dataset:
            hits=[x for x in sample_cols if x.startswith(c+"_") and dataset in x]
        if hits:
            X=rna[hits].apply(pd.to_numeric,errors="coerce").to_numpy(dtype=float)
            out[c]=np.nanmean(X,axis=1); counts[c]=len(hits)
    return out,counts

def solve_variant(M,base,variant,code,gene_names,vec,diet,sym2ent,trans,intset,args):
    mdl=base.copy()
    if variant in ("full","no_inter","no_trans"):
        mdl,_=M.apply_diet_bounds_layer1(mdl,code=code,diet_bounds=diet,diet_units="model")
    if variant in ("full","no_inter","no_diet"):
        mdl,*_=M.apply_expression_constraints_scoped(
            mdl,gene_names,vec,trans,eflux_quantile=args.eflux_quantile,
            eflux_floor=args.eflux_floor,eflux_cap=args.eflux_cap,label="L2",
            symbol_to_entrez=None)
    if variant in ("full","no_trans","no_diet"):
        mdl,*_=M.apply_expression_constraints_scoped(
            mdl,gene_names,vec,intset,eflux_quantile=args.eflux_quantile,
            eflux_floor=args.eflux_floor,eflux_cap=args.eflux_cap,label="L3",
            symbol_to_entrez=None)
    if variant=="vanilla":
        mdl,*_=M.apply_expression_constraints_scoped(
            mdl,gene_names,vec,trans|intset,eflux_quantile=args.eflux_quantile,
            eflux_floor=args.eflux_floor,eflux_cap=args.eflux_cap,label="vanilla",
            symbol_to_entrez=None)
    mdl=M.validate_model(mdl)
    sol=M.solve_flux(mdl,solve_mode="pfba",pfba_fraction=1.0)
    if sol.status!="optimal": raise RuntimeError(f"{variant}/{code}: {sol.status}")
    return sol.fluxes

def jac(a,b,k):
    A=set(np.argsort(-np.abs(a))[:k]); B=set(np.argsort(-np.abs(b))[:k])
    return len(A&B)/len(A|B) if A|B else np.nan

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--flux_script",required=True)
    ap.add_argument("--rnaseq",required=True)
    ap.add_argument("--dataset",default="GSE101657")
    ap.add_argument("--conditions",default="SCD,HFD,KD")
    ap.add_argument("--baseline",default="SCD")
    ap.add_argument("--model",required=True)
    ap.add_argument("--diet_bounds",required=True)
    ap.add_argument("--mapping",required=True)
    ap.add_argument("--objective",default="BIOMASS_mm_1_no_glygln")
    ap.add_argument("--eflux_quantile",type=float,default=0.95)
    ap.add_argument("--eflux_floor",type=float,default=0.1)
    ap.add_argument("--eflux_cap",type=float,default=1000.0)
    ap.add_argument("--topk",type=int,default=100)
    ap.add_argument("--out",default="layer_ablation_v2.csv")
    args=ap.parse_args()

    from cobra.io import load_json_model
    M=load_module(args.flux_script,"production")
    rna=pd.read_csv(args.rnaseq)
    if "Gene_Symbol" not in rna.columns:
        raise ValueError(f"Gene_Symbol column not found: {list(rna.columns[:8])}")
    gene_symbols=rna["Gene_Symbol"].astype(str).str.lower().tolist()
    conds=[x.strip() for x in args.conditions.split(",")]
    cexpr,counts=condition_means(rna,args.dataset,conds)
    print("[INFO] cohort condition sample counts:",counts)
    if args.baseline not in cexpr:raise ValueError("Baseline missing.")
    sym2ent=M.load_symbol_to_entrez_mapping(args.mapping)
    genes=[sym2ent.get(str(s).lower().strip(),str(s)) for s in gene_symbols]
    with open(args.diet_bounds) as f:raw=json.load(f)
    diet={M.canonical_code(k):{rid:[float(v[0]),float(v[1])] for rid,v in d.items()} for k,d in raw.items()}
    base=load_json_model(args.model)
    base,_=M.set_objective_reaction(base,objective_id=args.objective,sense="max")
    _,trans,intset=M.classify_reactions(base,transporter_strategy="either")
    model_gene_ids={str(g.id) for g in base.genes}
    n_model_hits=len(set(genes)&model_gene_ids)
    print(f"[PRECHECK] production-aligned gene IDs present in iMM1415: {n_model_hits}")
    if n_model_hits < 100:
        raise RuntimeError("Too few genes map to iMM1415; aborting layer ablation.")

    variants=["full","vanilla","no_diet","no_inter","no_trans"]
    flux={v:{} for v in variants}
    for c,vec in cexpr.items():
        for v in variants:
            flux[v][c]=solve_variant(M,base,v,c,genes,vec,diet,sym2ent,trans,intset,args)
        print(f"[INFO] solved {c}")

    rows=[]
    for c in cexpr:
        if c==args.baseline:continue
        fd=(flux["full"][c]-flux["full"][args.baseline]).values
        for v in variants[1:]:
            vd=(flux[v][c]-flux[v][args.baseline]).reindex(flux["full"][c].index).values
            # compare condition solutions and response vectors separately
            fv=flux["full"][c].values
            vv=flux[v][c].reindex(flux["full"][c].index).values
            rows.append({
                "contrast":f"{c}_vs_{args.baseline}","variant":v,
                "condition_flux_spearman":float(spearmanr(fv,vv).correlation),
                "condition_flux_directional_agreement":float(np.mean(np.sign(fv)==np.sign(vv))),
                "response_spearman":float(spearmanr(fd,vd).correlation),
                "response_directional_agreement":float(np.mean(np.sign(fd)==np.sign(vd))),
                f"response_top{args.topk}_jaccard":jac(fd,vd,args.topk)
            })
    pd.DataFrame(rows).to_csv(args.out,index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"[OK] {args.out}")

if __name__=="__main__":
    main()

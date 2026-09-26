#!/usr/bin/env python3
"""
benchmark_liver_methods.py
==========================
Reviewer-facing comparative benchmark on one within-study mouse-liver cohort.

This complements the external E. coli 13C benchmark. Here all methods receive the
same liver expression samples, same iMM1415 model, same objective, and the same
biological replicates. The three-layer method is compared with:
- conventional single-scope E-Flux + pFBA
- GIMME + pFBA representative solution
- optional unconstrained pFBA baseline

Primary recommended cohort: GSE101657, HFD vs SCD.
This cohort is used as a proof-of-concept benchmark, not as a replacement for the
multi-cohort RQ1 inference.

Outputs
-------
fluxes_<method>_<objective>.csv
reaction_stats_<method>_<objective>.csv
method_concordance.csv
significant_set_overlap.csv
positive_control_results.csv
objective_sensitivity_three_layer.csv
run_manifest.json
"""
from __future__ import annotations
import argparse, importlib.util, json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, spearmanr

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
sys.path.insert(0,str(HERE.parent/"common"))
from benchmark_methods_core import gimme_solve, safe_spearman, directional_agreement, topk_jaccard
from deterministic_solver_v2 import deterministic_optimize, configure_solver_deterministic

def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def bh(p):
    p=np.asarray(p,float)
    q=np.full_like(p,np.nan)
    ok=np.isfinite(p)
    vals=p[ok]
    if not len(vals):return q
    order=np.argsort(vals); ranked=vals[order]
    adj=ranked*len(ranked)/np.arange(1,len(ranked)+1)
    adj=np.minimum.accumulate(adj[::-1])[::-1]
    out=np.empty_like(vals); out[order]=np.minimum(adj,1.0)
    q[np.where(ok)[0]]=out
    return q

def cohen_d(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    a=a[np.isfinite(a)]; b=b[np.isfinite(b)]
    if len(a)<2 or len(b)<2:return np.nan
    va=np.var(a,ddof=1); vb=np.var(b,ddof=1)
    pooled=((len(a)-1)*va+(len(b)-1)*vb)/(len(a)+len(b)-2)
    return (np.mean(a)-np.mean(b))/np.sqrt(pooled) if pooled>0 else 0.0

def stats_from_flux(df, groups, test, baseline):
    rows=[]
    for rid,row in df.iterrows():
        a=row[[c for c in df.columns if groups[c]==test]].astype(float).values
        b=row[[c for c in df.columns if groups[c]==baseline]].astype(float).values
        aa=a[np.isfinite(a)]; bb=b[np.isfinite(b)]
        p=ttest_ind(aa,bb,equal_var=False,nan_policy="omit").pvalue if len(aa)>=2 and len(bb)>=2 else np.nan
        rows.append({
            "ReactionID":rid,
            f"{test}_mean":np.nanmean(aa) if len(aa) else np.nan,
            f"{baseline}_mean":np.nanmean(bb) if len(bb) else np.nan,
            "MeanDiff":(np.nanmean(aa)-np.nanmean(bb)) if len(aa) and len(bb) else np.nan,
            "Cohen_d":cohen_d(aa,bb),"p_value":p
        })
    out=pd.DataFrame(rows)
    out["q_value"]=bh(out.p_value.values)
    out["Significant"]=(out.q_value<0.05)&(out.Cohen_d.abs()>=0.5)
    out["Direction"]=np.where(out.MeanDiff>0,"Up",np.where(out.MeanDiff<0,"Down","Zero"))
    return out

def load_expression(args,M):
    rna=pd.read_csv(args.rnaseq)
    if "Gene_Symbol" not in rna.columns:
        raise ValueError(f"Gene_Symbol column not found. First columns: {list(rna.columns[:8])}")
    gene_symbols=rna["Gene_Symbol"].astype(str).str.lower().tolist()
    sample_cols=[c for c in rna.columns if c not in {"Gene_Symbol","Gene_ID"}]
    groups={}
    selected=[]
    for c in args.conditions.split(","):
        c=c.strip()
        prefix=f"{c}_{args.dataset}_"
        hits=[x for x in sample_cols if x.startswith(prefix)]
        if not hits:
            # fallback for slightly different naming
            hits=[x for x in sample_cols if x.startswith(c+"_") and args.dataset in x]
        for x in hits: groups[x]=c
        selected.extend(hits)
    if len(selected)==0:
        raise ValueError(f"No sample columns matched dataset={args.dataset}.")
    mat=rna[selected].apply(pd.to_numeric,errors="coerce")
    sym2ent=M.load_symbol_to_entrez_mapping(args.mapping) if args.mapping else None
    gene_ids=[sym2ent.get(str(s).lower().strip(),str(s)) for s in gene_symbols] if sym2ent else gene_symbols
    return gene_symbols,gene_ids,mat,groups,sym2ent

def expr_dict(gene_symbols,vec,sym2ent):
    out={}
    for g,v in zip(gene_symbols,vec):
        if v is None or not np.isfinite(v):continue
        key=sym2ent.get(str(g).lower().strip(),str(g)) if sym2ent else str(g)
        out[str(key)]=float(v)
    return out

def solve_one(M,base,method,condition,gene_ids,vec,sym2ent,diet_bounds,
              trans_set,int_set,args):
    mdl=base.copy()
    configure_solver_deterministic(mdl,solver=args.solver,verbose=False)
    # All methods receive the same condition-specific environmental (diet) bounds.
    # This isolates the transcriptome-integration method rather than confounding it
    # with presence/absence of the known dietary environment.
    mdl,_=M.apply_diet_bounds_layer1(mdl,code=condition,diet_bounds=diet_bounds,diet_units="model")

    if method=="three_layer":
        mdl,l2_changed,l2_matched,*_=M.apply_expression_constraints_scoped(
            mdl,gene_ids,vec,trans_set,eflux_quantile=args.eflux_quantile,
            eflux_floor=args.eflux_floor,eflux_cap=args.eflux_cap,label="L2",
            symbol_to_entrez=None)
        mdl,l3_changed,l3_matched,*_=M.apply_expression_constraints_scoped(
            mdl,gene_ids,vec,int_set,eflux_quantile=args.eflux_quantile,
            eflux_floor=args.eflux_floor,eflux_cap=args.eflux_cap,label="L3",
            symbol_to_entrez=None)
        if (l2_matched+l3_matched)==0 or (l2_changed+l3_changed)==0:
            raise RuntimeError(
                f"Expression mapping failed for three_layer: L2 matched={l2_matched}, "
                f"L3 matched={l3_matched}, L2 changed={l2_changed}, L3 changed={l3_changed}"
            )
        mdl=M.validate_model(mdl)
        return deterministic_optimize(mdl,use_pfba=True,fraction_of_optimum=1.0)

    if method=="eflux":
        scope=trans_set|int_set
        mdl,changed,matched,*_=M.apply_expression_constraints_scoped(
            mdl,gene_ids,vec,scope,eflux_quantile=args.eflux_quantile,
            eflux_floor=args.eflux_floor,eflux_cap=args.eflux_cap,label="vanilla",
            symbol_to_entrez=None)
        if matched==0 or changed==0:
            raise RuntimeError(f"Expression mapping failed for E-Flux: matched={matched}, changed={changed}")
        mdl=M.validate_model(mdl)
        return deterministic_optimize(mdl,use_pfba=True,fraction_of_optimum=1.0)

    if method=="gimme":
        # gene_ids are already production-aligned Entrez IDs.
        ex={str(g):float(v) for g,v in zip(gene_ids,vec) if np.isfinite(v)}
        sol,_=gimme_solve(mdl,ex,args.current_objective,
                          threshold_quantile=args.gimme_threshold_quantile,
                          objective_fraction=args.gimme_objective_fraction,
                          solver=args.solver)
        return sol

    if method=="pfba":
        # Diet-only pFBA baseline.
        return deterministic_optimize(mdl,use_pfba=True,fraction_of_optimum=1.0)
    raise ValueError(method)

def eval_controls(flux,groups,controls,test,baseline):
    rows=[]
    for ctl in controls.get("controls",[]):
        comp=f"{test}_vs_{baseline}"
        if comp not in ctl.get("comparisons",[]):continue
        ids=[r for r in ctl["reaction_ids"] if r in flux.index]
        if not ids:continue
        t=[c for c in flux.columns if groups[c]==test]
        b=[c for c in flux.columns if groups[c]==baseline]
        # mean absolute pathway activity, then group mean
        tv=float(np.nanmean(np.abs(flux.loc[ids,t].values)))
        bv=float(np.nanmean(np.abs(flux.loc[ids,b].values)))
        observed="increase" if tv>bv else ("decrease" if tv<bv else "no_change")
        expected=ctl["expected_direction"]
        rows.append({"control":ctl["name"],"comparison":comp,
                     "n_reactions":len(ids),"test_mean_abs":tv,"baseline_mean_abs":bv,
                     "fold_test_vs_baseline":tv/bv if bv else np.nan,
                     "expected":expected,"observed":observed,"passed":observed==expected})
    return rows

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--flux_script",required=True,
                    help="Revised production map_fix...manuscript_rerun.py")
    ap.add_argument("--rnaseq",required=True)
    ap.add_argument("--dataset",default="GSE101657")
    ap.add_argument("--conditions",default="SCD,HFD")
    ap.add_argument("--baseline",default="SCD")
    ap.add_argument("--test",default="HFD")
    ap.add_argument("--model",required=True)
    ap.add_argument("--diet_bounds",required=True)
    ap.add_argument("--mapping",required=True)
    ap.add_argument("--positive_controls")
    ap.add_argument("--methods",default="three_layer,eflux,gimme,pfba")
    ap.add_argument("--objectives",default="BIOMASS_mm_1_no_glygln,ATPM")
    ap.add_argument("--eflux_quantile",type=float,default=0.95)
    ap.add_argument("--eflux_floor",type=float,default=0.1)
    ap.add_argument("--eflux_cap",type=float,default=1000.0)
    ap.add_argument("--gimme_threshold_quantile",type=float,default=0.25)
    ap.add_argument("--gimme_objective_fraction",type=float,default=0.90)
    ap.add_argument("--topk",type=int,default=100)
    ap.add_argument("--solver",default="gurobi")
    ap.add_argument("--output",default="benchmark_liver_methods")
    args=ap.parse_args()

    from cobra.io import load_json_model
    out=Path(args.output); out.mkdir(parents=True,exist_ok=True)
    M=load_module(args.flux_script,"production")
    gene_symbols,gene_ids,expr,groups,sym2ent=load_expression(args,M)
    if args.baseline not in groups.values() or args.test not in groups.values():
        raise ValueError("Requested baseline/test not present in selected cohort.")
    print("[INFO] sample counts:",pd.Series(groups).value_counts().to_dict())

    with open(args.diet_bounds) as f: raw=json.load(f)
    diet_bounds={M.canonical_code(k):{rid:[float(v[0]),float(v[1])] for rid,v in d.items()}
                 for k,d in raw.items()}
    base0=load_json_model(args.model)
    _,trans_set,int_set=M.classify_reactions(base0,transporter_strategy="either")
    model_gene_ids={str(g.id) for g in base0.genes}
    n_model_hits=len(set(str(g) for g in gene_ids) & model_gene_ids)
    print(f"[PRECHECK] production-aligned gene IDs present in iMM1415: {n_model_hits}")
    if n_model_hits < 100:
        raise RuntimeError("Too few expression genes map to iMM1415; aborting benchmark.")
    methods=[x.strip() for x in args.methods.split(",") if x.strip()]
    objectives=[x.strip() for x in args.objectives.split(",") if x.strip()]
    controls=json.loads(Path(args.positive_controls).read_text()) if args.positive_controls else {"controls":[]}

    all_flux={}; all_stats={}; pc_rows=[]
    t0=time.time()
    for objective in objectives:
        base,_=M.set_objective_reaction(base0.copy(),objective_id=objective,sense="max")
        args.current_objective=objective
        for method in methods:
            # Objective sensitivity is most interpretable for three-layer; comparators
            # are still solved under both objectives when requested.
            cols={}
            for j,sample in enumerate(expr.columns,1):
                cond=groups[sample]
                vec=expr[sample].to_numpy(dtype=float)
                try:
                    sol=solve_one(M,base,method,cond,gene_ids,vec,sym2ent,diet_bounds,
                                  trans_set,int_set,args)
                    if sol.status!="optimal":
                        raise RuntimeError(f"status={sol.status}")
                    cols[sample]=sol.fluxes
                except Exception as exc:
                    raise RuntimeError(f"{method}/{objective}/{sample}: {exc}") from exc
                print(f"[{objective}] {method}: {j}/{len(expr.columns)} {sample}")
            flux=pd.DataFrame(cols)
            all_flux[(method,objective)]=flux
            flux.to_csv(out/f"fluxes_{method}_{objective}.csv")
            st=stats_from_flux(flux,groups,args.test,args.baseline)
            all_stats[(method,objective)]=st
            st.to_csv(out/f"reaction_stats_{method}_{objective}.csv",index=False)
            for row in eval_controls(flux,groups,controls,args.test,args.baseline):
                row.update({"method":method,"objective":objective}); pc_rows.append(row)

    # Within-objective comparator concordance to three-layer
    conc=[]; overlaps=[]
    for objective in objectives:
        if ("three_layer",objective) not in all_stats:continue
        ref=all_stats[("three_layer",objective)].set_index("ReactionID")
        ref_sig=set(ref.index[ref.Significant])
        for method in methods:
            if method=="three_layer" or (method,objective) not in all_stats:continue
            q=all_stats[(method,objective)].set_index("ReactionID").reindex(ref.index)
            conc.append({
                "objective":objective,"method":method,
                "reference":"three_layer",
                "spearman_mean_diff":safe_spearman(ref.MeanDiff,q.MeanDiff),
                "directional_agreement_mean_diff":directional_agreement(ref.MeanDiff,q.MeanDiff),
                f"top{args.topk}_jaccard":topk_jaccard(ref.MeanDiff.values,q.MeanDiff.values,args.topk)
            })
            sig=set(q.index[q.Significant.fillna(False)])
            overlaps.append({
                "objective":objective,"method":method,
                "reference":"three_layer",
                "n_sig_reference":len(ref_sig),"n_sig_method":len(sig),
                "n_intersection":len(ref_sig&sig),
                "jaccard_significant":len(ref_sig&sig)/len(ref_sig|sig) if ref_sig|sig else np.nan
            })
    pd.DataFrame(conc).to_csv(out/"method_concordance.csv",index=False)
    pd.DataFrame(overlaps).to_csv(out/"significant_set_overlap.csv",index=False)
    pd.DataFrame(pc_rows).to_csv(out/"positive_control_results.csv",index=False)

    # Objective sensitivity of the primary three-layer method
    osens=[]
    if len(objectives)>=2:
        ref_obj=objectives[0]
        ref=all_stats[("three_layer",ref_obj)].set_index("ReactionID")
        for objective in objectives[1:]:
            q=all_stats[("three_layer",objective)].set_index("ReactionID").reindex(ref.index)
            osens.append({
                "reference_objective":ref_obj,"comparison_objective":objective,
                "spearman_mean_diff":safe_spearman(ref.MeanDiff,q.MeanDiff),
                "directional_agreement_mean_diff":directional_agreement(ref.MeanDiff,q.MeanDiff),
                f"top{args.topk}_jaccard":topk_jaccard(ref.MeanDiff.values,q.MeanDiff.values,args.topk),
                "n_sig_reference":int(ref.Significant.sum()),
                "n_sig_comparison":int(q.Significant.fillna(False).sum()),
                "significant_jaccard":len(set(ref.index[ref.Significant])&set(q.index[q.Significant.fillna(False)]))/
                    max(1,len(set(ref.index[ref.Significant])|set(q.index[q.Significant.fillna(False)])))
            })
    pd.DataFrame(osens).to_csv(out/"objective_sensitivity_three_layer.csv",index=False)
    manifest={
        "dataset":args.dataset,"sample_counts":pd.Series(groups).value_counts().to_dict(),
        "methods":methods,"objectives":objectives,
        "eflux":{"quantile":args.eflux_quantile,"floor":args.eflux_floor,"cap":args.eflux_cap},
        "gimme":{"threshold_quantile":args.gimme_threshold_quantile,
                 "objective_fraction":args.gimme_objective_fraction},
        "elapsed_seconds":time.time()-t0
    }
    (out/"run_manifest.json").write_text(json.dumps(manifest,indent=2))
    print(f"[DONE] {out}")

if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""
sweep_threelayer_pipeline_v2.py
===============================
Reviewer-facing normalization / scaling / objective sensitivity for the production
three-layer pipeline.

Key corrections relative to the earlier sweep
----------------------------------------------
- Uses the revised production pFBA solver, not the historical FBA pipeline.
- Rejects the old bug where the string "fba" was passed into a boolean use_pfba
  argument and therefore always executed pFBA.
- Uses np.nanmean for the outer-merged multi-cohort expression matrix.
- Can restrict the sweep to one within-study cohort (recommended: GSE101657) so
  normalization sensitivity is not confounded by study-specific missingness.
- Does not tune/select the best setting. The primary q=0.95, floor=0.1, cap=1000,
  biomass configuration is marked in the output and all sensitivity settings are
  reported.

Suggested revision run
----------------------
python sweep_threelayer_pipeline_v2.py ^
  --flux_script map_fixv5_multigroupsv8_layered_manuscript_rerun.py ^
  --model iMM1415.json ^
  --rnaseq GSEMERGED_SCD_HFD_KD_WD_gene_expression.csv ^
  --dataset GSE101657 ^
  --mapping mouse_entrez_to_symbol.csv ^
  --diet_bounds expanded_diet_bounds_flat.json ^
  --positive_controls positive_controls.json ^
  --quantiles 0.90,0.95,0.99 ^
  --floors 0.05,0.1 ^
  --caps 1000 ^
  --objectives BIOMASS_mm_1_no_glygln,ATPM ^
  --results_dir threelayer_sensitivity_GSE101657
"""
from __future__ import annotations
import argparse, importlib.util, itertools, json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def condition_means(rna,dataset,conditions):
    sample_cols=list(rna.columns[2:])
    out={}
    counts={}
    for c in conditions:
        hits=[x for x in sample_cols if x.startswith(f"{c}_{dataset}_")] if dataset else [x for x in sample_cols if x.startswith(c+"_")]
        if not hits and dataset:
            hits=[x for x in sample_cols if x.startswith(c+"_") and dataset in x]
        if hits:
            vals=rna[hits].apply(pd.to_numeric,errors="coerce").to_numpy(dtype=float)
            out[c]=np.nanmean(vals,axis=1)
            counts[c]=len(hits)
    return out,counts

def evaluate_controls(fluxes,config,baseline):
    rows=[]
    for ctl in config.get("controls",[]):
        ids=[r for r in ctl.get("reaction_ids",[]) if r in next(iter(fluxes.values())).index]
        if not ids:continue
        for comp in ctl.get("comparisons",[]):
            if "_vs_" not in comp:continue
            t,b=comp.split("_vs_",1)
            if t not in fluxes or b not in fluxes:continue
            tv=float(np.nanmean(np.abs(fluxes[t].loc[ids].values)))
            bv=float(np.nanmean(np.abs(fluxes[b].loc[ids].values)))
            obs="increase" if tv>bv else ("decrease" if tv<bv else "no_change")
            rows.append({"control":ctl["name"],"comparison":comp,"n_reactions":len(ids),
                         "fold":tv/bv if bv else np.nan,"expected":ctl["expected_direction"],
                         "observed":obs,"passed":obs==ctl["expected_direction"]})
    return pd.DataFrame(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--flux_script",required=True)
    ap.add_argument("--model",required=True)
    ap.add_argument("--rnaseq",required=True)
    ap.add_argument("--dataset",default="GSE101657")
    ap.add_argument("--conditions",default="SCD,HFD,KD")
    ap.add_argument("--baseline",default="SCD")
    ap.add_argument("--mapping",required=True)
    ap.add_argument("--diet_bounds",required=True)
    ap.add_argument("--positive_controls",required=True)
    ap.add_argument("--quantiles",default="0.90,0.95,0.99")
    ap.add_argument("--floors",default="0.05,0.1")
    ap.add_argument("--caps",default="1000")
    ap.add_argument("--objectives",default="BIOMASS_mm_1_no_glygln,ATPM")
    ap.add_argument("--solver",default="gurobi")
    ap.add_argument("--results_dir",default="threelayer_sensitivity")
    args=ap.parse_args()

    from cobra.io import load_json_model
    M=load_module(args.flux_script,"production")
    out=Path(args.results_dir); out.mkdir(parents=True,exist_ok=True)
    rna=pd.read_csv(args.rnaseq)
    if "Gene_Symbol" not in rna.columns:
        raise ValueError(f"Gene_Symbol column not found: {list(rna.columns[:8])}")
    gene_symbols=rna["Gene_Symbol"].astype(str).str.lower().tolist()
    conds=[x.strip() for x in args.conditions.split(",")]
    cexpr,counts=condition_means(rna,args.dataset,conds)
    if args.baseline not in cexpr: raise ValueError("Baseline absent.")
    print("[INFO] cohort condition sample counts:",counts)
    sym2ent=M.load_symbol_to_entrez_mapping(args.mapping)
    gene_ids=[sym2ent.get(str(s).lower().strip(),str(s)) for s in gene_symbols]

    with open(args.diet_bounds) as f: raw=json.load(f)
    diet={M.canonical_code(k):{rid:[float(v[0]),float(v[1])] for rid,v in d.items()} for k,d in raw.items()}
    controls=json.loads(Path(args.positive_controls).read_text())
    base0=load_json_model(args.model)
    _,trans,intset=M.classify_reactions(base0,transporter_strategy="either")
    model_gene_ids={str(g.id) for g in base0.genes}
    n_model_hits=len(set(gene_ids)&model_gene_ids)
    print(f"[PRECHECK] production-aligned gene IDs present in iMM1415: {n_model_hits}")
    if n_model_hits < 100:
        raise RuntimeError("Too few genes map to iMM1415; aborting sensitivity sweep.")

    qs=[float(x) for x in args.quantiles.split(",")]
    fs=[float(x) for x in args.floors.split(",")]
    caps=[float(x) for x in args.caps.split(",")]
    objs=[x.strip() for x in args.objectives.split(",")]
    grid=list(itertools.product(qs,fs,caps,objs))
    flux_store={}; pc=[]; meta=[]
    t0=time.time()

    for idx,(q,floor,cap,obj) in enumerate(grid,1):
        label=f"q{q:g}_floor{floor:g}_cap{cap:g}_{obj}"
        base,_=M.set_objective_reaction(base0.copy(),objective_id=obj,sense="max")
        condition_flux={}
        for c,vec in cexpr.items():
            mdl=base.copy()
            try: mdl.solver=args.solver
            except Exception: pass
            mdl,_=M.apply_diet_bounds_layer1(mdl,code=c,diet_bounds=diet,diet_units="model")
            mdl,*_=M.apply_expression_constraints_scoped(
                mdl,gene_ids,vec,trans,eflux_quantile=q,eflux_floor=floor,
                eflux_cap=cap,label="L2",symbol_to_entrez=None)
            mdl,*_=M.apply_expression_constraints_scoped(
                mdl,gene_ids,vec,intset,eflux_quantile=q,eflux_floor=floor,
                eflux_cap=cap,label="L3",symbol_to_entrez=None)
            mdl=M.validate_model(mdl)
            sol=M.solve_flux(mdl,solve_mode="pfba",pfba_fraction=1.0)
            if sol.status!="optimal": raise RuntimeError(f"{label}/{c}: {sol.status}")
            condition_flux[c]=sol.fluxes
        flux=pd.DataFrame(condition_flux)
        flux.to_csv(out/f"flux_{label}.csv")
        flux_store[label]=flux
        setting={"setting":label,"quantile":q,"floor":floor,"cap":cap,"objective":obj,
                 "is_primary":q==0.95 and floor==0.1 and cap==1000 and obj=="BIOMASS_mm_1_no_glygln"}
        meta.append(setting)
        qpc=evaluate_controls(condition_flux,controls,args.baseline)
        if len(qpc):
            for k,v in setting.items():qpc[k]=v
            pc.append(qpc)
        print(f"[{idx}/{len(grid)}] {label} ({time.time()-t0:.0f}s)")

    pd.DataFrame(meta).to_csv(out/"settings.csv",index=False)
    if pc:
        pcd=pd.concat(pc,ignore_index=True); pcd.to_csv(out/"positive_control_long.csv",index=False)
        summ=(pcd.groupby(["setting","quantile","floor","cap","objective","is_primary"])
              .agg(controls_passed=("passed","sum"),controls_total=("passed","size"))
              .reset_index())
        summ["pass_rate"]=summ.controls_passed/summ.controls_total
        summ.to_csv(out/"sensitivity_summary.csv",index=False)

    # Flux / contrast stability to primary setting, separated by objective.
    rows=[]
    for obj in objs:
        labels=[m["setting"] for m in meta if m["objective"]==obj]
        primary=[m["setting"] for m in meta if m["objective"]==obj and m["quantile"]==0.95 and m["floor"]==0.1 and m["cap"]==1000]
        ref=primary[0] if primary else labels[0]
        for lab in labels:
            for c in cexpr:
                a=flux_store[ref][c].values; b=flux_store[lab][c].values
                rho=float(spearmanr(a,b).correlation)
                agree=float(np.mean(np.sign(a)==np.sign(b)))
                rows.append({"objective":obj,"reference_setting":ref,"setting":lab,
                             "comparison_type":"condition_flux","condition":c,
                             "spearman":rho,"directional_agreement":agree})
            for c in cexpr:
                if c==args.baseline:continue
                a=(flux_store[ref][c]-flux_store[ref][args.baseline]).values
                b=(flux_store[lab][c]-flux_store[lab][args.baseline]).values
                rho=float(spearmanr(a,b).correlation)
                agree=float(np.mean(np.sign(a)==np.sign(b)))
                rows.append({"objective":obj,"reference_setting":ref,"setting":lab,
                             "comparison_type":"contrast_change","condition":f"{c}_vs_{args.baseline}",
                             "spearman":rho,"directional_agreement":agree})
    pd.DataFrame(rows).to_csv(out/"flux_stability_to_primary.csv",index=False)

    # Objective sensitivity at the shared q/floor/cap primary normalization.
    if len(objs)>=2:
        rows=[]
        ref_obj=objs[0]
        ref=f"q{0.95:g}_floor{0.1:g}_cap{1000:g}_{ref_obj}"
        if ref in flux_store:
            for obj in objs[1:]:
                lab=f"q{0.95:g}_floor{0.1:g}_cap{1000:g}_{obj}"
                if lab not in flux_store:continue
                for c in cexpr:
                    if c==args.baseline:continue
                    a=(flux_store[ref][c]-flux_store[ref][args.baseline]).values
                    b=(flux_store[lab][c]-flux_store[lab][args.baseline]).values
                    rows.append({"reference_objective":ref_obj,"comparison_objective":obj,
                                 "contrast":f"{c}_vs_{args.baseline}",
                                 "spearman_change":float(spearmanr(a,b).correlation),
                                 "directional_agreement_change":float(np.mean(np.sign(a)==np.sign(b)))})
        pd.DataFrame(rows).to_csv(out/"objective_sensitivity.csv",index=False)
    print(f"[DONE] {out}")

if __name__=="__main__":
    main()

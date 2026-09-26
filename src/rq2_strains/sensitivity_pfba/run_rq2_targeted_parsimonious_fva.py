#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from cobra.io import load_json_model

from fva_refinement_common import (
    fva_difference,
    load_modeling_helpers,
    nanmean_expression,
    run_parsimonious_fva,
)


def parse_eps(text: str) -> list[float]:
    vals=[]
    for x in re.split(r"[,;\s]+",str(text).strip()):
        if x:
            v=float(x)
            if v<0: raise ValueError("parsimony epsilon must be >=0")
            vals.append(v)
    return vals


def sample_columns_for_group(rna: pd.DataFrame, group: str) -> list[str]:
    rx=re.compile(rf"(?<![A-Za-z0-9]){re.escape(group)}(?![A-Za-z0-9])",re.I)
    return [c for c in rna.columns if c not in {"Gene_Symbol","Gene_ID","Entrez"} and rx.search(str(c))]


def main():
    ap=argparse.ArgumentParser(description="Targeted parsimonious FVA for the 16-reaction RQ2 universal magnitude core.")
    ap.add_argument("--config",default="run_layered_config_RQ2_pfba_FVA.json")
    ap.add_argument("--modeling_script",default="map_fixv5_multigroupsv8_layered_manuscript_rerun.py")
    ap.add_argument("--rq2_results",default="RERUN_RQ2_PFBA_FVA")
    ap.add_argument("--output",default="RERUN_RQ2_PFBA_FVA/targeted_parsimonious_FVA")
    ap.add_argument("--parsimony_eps",default="0.01,0.05")
    ap.add_argument("--tolerance",type=float,default=1e-7)
    ap.add_argument("--processes",type=int,default=1)
    args=ap.parse_args()

    cfg=json.loads(Path(args.config).read_text(encoding="utf-8"))
    common=cfg.get("common_args",{})
    out=Path(args.output); out.mkdir(parents=True,exist_ok=True)
    (out/"per_strain").mkdir(exist_ok=True)
    helpers=load_modeling_helpers(args.modeling_script)

    core_path=Path(args.rq2_results)/"cross_strain"/"RQ2_universal_magnitude_core.csv"
    stats_path=Path(args.rq2_results)/"cross_strain"/"RQ2_all_strain_reaction_stats_long.csv"
    core=pd.read_csv(core_path)
    targets=core["ReactionID"].astype(str).tolist()
    if len(targets)!=16:
        print(f"[WARN] Expected 16 universal-core reactions, found {len(targets)}")
    stats=pd.read_csv(stats_path)
    stats=stats[stats["ReactionID"].astype(str).isin(targets)].copy()
    stats["PrimaryDirection"]=np.where(stats["HFD_minus_SCD"]>0,"Up",np.where(stats["HFD_minus_SCD"]<0,"Down","Zero"))

    model_file=common.get("model_file","iMM1415.json")
    diet_file=common.get("diet_bounds_json","expanded_diet_bounds_flat.json")
    mapping_file=common.get("mapping_file","mouse_entrez_to_symbol.csv")
    objective_id=common.get("objective_id","BIOMASS_mm_1_no_glygln")
    objective_sense=common.get("objective_sense","max")
    eflux_quantile=float(common.get("eflux_quantile",0.95))
    eflux_floor=float(common.get("eflux_floor",0.1))
    eflux_cap=float(common.get("eflux_cap",1000.0))
    transporter_strategy=common.get("transporter_strategy","either")
    solver=common.get("solver","gurobi")
    fva_fraction=float(common.get("fva_fraction",1.0))

    symbol_map=helpers.load_symbol_to_entrez_mapping(mapping_file)
    raw=json.loads(Path(diet_file).read_text(encoding="utf-8"))
    diet_bounds={helpers.canonical_code(k):{rid:[float(v[0]),float(v[1])] for rid,v in d.items()} for k,d in raw.items()}

    base=load_json_model(model_file)
    base=helpers.configure_solver(base,solver=solver,deterministic=True)
    base,obj_id=helpers.set_objective_reaction(base,objective_id=objective_id,objective_regex=None,sense=objective_sense)
    _,trans_set,int_set=helpers.classify_reactions(base,transporter_strategy=transporter_strategy,
        transporter_regex=None,transporter_subsystem_regex=None,compartments_for_transport=("e",))
    print(f"[MODEL] reactions={len(base.reactions)} targets={len(targets)} objective={obj_id}")

    eps_values=parse_eps(args.parsimony_eps)
    long=[]; manifest=[]; qc_all=[]

    for job in cfg.get("jobs",[]):
        strain=str(job["name"]); inp=Path(job["input_csv"])
        rna=pd.read_csv(inp)
        if "Gene_Symbol" not in rna.columns:
            raise ValueError(f"{inp} lacks Gene_Symbol")
        hcols=sample_columns_for_group(rna,"HFD")
        scols=sample_columns_for_group(rna,"SCD")
        if not hcols or not scols:
            raise ValueError(f"{strain}: HFD={len(hcols)} SCD={len(scols)} columns")
        print(f"\n[STRAIN] {strain}: HFD n={len(hcols)}, SCD n={len(scols)}")
        gene_symbols=rna["Gene_Symbol"].astype(str).str.lower().tolist()
        gene_ids=[symbol_map.get(sym.strip().lower(),sym) for sym in gene_symbols]

        def build(group,cols):
            vec=nanmean_expression(rna,cols)
            mdl=base.copy(); mdl=helpers.configure_solver(mdl,solver=solver,deterministic=True,logger=lambda *_:None)
            mdl,l1=helpers.apply_diet_bounds_layer1(mdl,code=group,diet_bounds=diet_bounds,diet_units="model",mw_map={},gDW=1.0,hours_per_day=24.0)
            mdl,l2,*_=helpers.apply_expression_constraints_scoped(mdl,gene_ids,vec,trans_set,
                eflux_quantile=eflux_quantile,eflux_floor=eflux_floor,eflux_cap=eflux_cap,label="L2",symbol_to_entrez=symbol_map)
            mdl,l3,*_=helpers.apply_expression_constraints_scoped(mdl,gene_ids,vec,int_set,
                eflux_quantile=eflux_quantile,eflux_floor=eflux_floor,eflux_cap=eflux_cap,label="L3",symbol_to_entrez=symbol_map)
            mdl=helpers.validate_model(mdl)
            sol=helpers.solve_flux(mdl,solve_mode="pfba",pfba_fraction=1.0)
            return mdl,{"Strain":strain,"Group":group,"N_replicates":len(cols),"BiologicalObjective":float(sol.objective_value),
                        "pFBATotalFlux":float(getattr(sol,"total_flux",np.abs(sol.fluxes).sum())),"L1":int(l1),"L2":int(l2),"L3":int(l3)}

        hmdl,hinfo=build("HFD",hcols); smdl,sinfo=build("SCD",scols)
        manifest.extend([hinfo,sinfo])
        sstats=stats[stats["Strain"].astype(str)==strain].set_index("ReactionID")

        for eps in eps_values:
            hfva,hqc,hmeta=run_parsimonious_fva(hmdl,targets,fraction=fva_fraction,parsimony_epsilon=eps,processes=args.processes,tolerance=args.tolerance)
            sfva,sqc,smeta=run_parsimonious_fva(smdl,targets,fraction=fva_fraction,parsimony_epsilon=eps,processes=args.processes,tolerance=args.tolerance)
            label=f"eps{eps:g}".replace(".","p")
            hfva.reset_index().to_csv(out/"per_strain"/f"{strain}_HFD_pfva_{label}.csv",index=False)
            sfva.reset_index().to_csv(out/"per_strain"/f"{strain}_SCD_pfva_{label}.csv",index=False)
            hqc["Strain"]=strain; hqc["Group"]="HFD"; hqc["Epsilon"]=eps
            sqc["Strain"]=strain; sqc["Group"]="SCD"; sqc["Epsilon"]=eps
            qc_all.extend([hqc,sqc])
            d=fva_difference(hfva,sfva,tolerance=args.tolerance).reset_index()
            d["Strain"]=strain; d["ParsimonyEpsilon"]=eps
            d["PrimaryDirection"]=d["ReactionID"].map(sstats["PrimaryDirection"])
            d["HFD_minus_SCD_pFBA"]=d["ReactionID"].map(sstats["HFD_minus_SCD"])
            d["MatchesPrimary"]=d["Robust"] & (d["RobustDirection"]==d["PrimaryDirection"])
            d.to_csv(out/"per_strain"/f"{strain}_HFD_vs_SCD_pfva_{label}_difference.csv",index=False)
            long.append(d)

            pd.DataFrame([{**hinfo,"ParsimonyEpsilon":eps,**hmeta},{**sinfo,"ParsimonyEpsilon":eps,**smeta}]).to_csv(
                out/"per_strain"/f"{strain}_pfva_{label}_metadata.csv",index=False)

    pd.DataFrame(manifest).to_csv(out/"representative_model_manifest.csv",index=False)
    if qc_all: pd.concat(qc_all,ignore_index=True).to_csv(out/"FVA_numerical_QC.csv",index=False)
    all_long=pd.concat(long,ignore_index=True); all_long.to_csv(out/"RQ2_targeted_pfva_all_strains_long.csv",index=False)

    rows=[]
    for (eps,rid),rdf in all_long.groupby(["ParsimonyEpsilon","ReactionID"]):
        n=rdf["Strain"].nunique(); nrob=int(rdf["Robust"].sum()); nmatch=int(rdf["MatchesPrimary"].sum())
        rows.append({"ParsimonyEpsilon":eps,"ReactionID":rid,"N_strains":n,"N_robust":nrob,"N_match_primary":nmatch,
                     "AllStrainsRobustSameAsPrimary":bool(n>0 and nmatch==n),
                     "RobustDirectionsObserved":";".join(sorted(set(rdf.loc[rdf.Robust,"RobustDirection"].astype(str))))})
    summary=pd.DataFrame(rows); summary.to_csv(out/"RQ2_targeted_pfva_cross_strain_summary.csv",index=False)
    high=summary.groupby("ParsimonyEpsilon").agg(N_reactions=("ReactionID","nunique"),
        N_all_strains_robust=("AllStrainsRobustSameAsPrimary","sum")).reset_index()
    high.to_csv(out/"RQ2_targeted_pfva_high_level_summary.csv",index=False)
    print("\n[SUMMARY]")
    print(high.to_string(index=False))
    print(f"\n[DONE] {out.resolve()}")

if __name__=="__main__":
    main()

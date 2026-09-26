#!/usr/bin/env python3
"""
rq4_hepatic_integration_final.py
================================
Final RQ4 host-side runner that reuses the biologically corrected v13 constraint
functions but fixes two provenance issues:

1. pFBA fraction is explicit and defaults to 1.0.
2. The biological objective value is measured BEFORE pFBA and stored separately.
   COBRApy pfba().objective_value is the secondary minimum-total-flux objective
   and must not be reported as biomass/ATPM.

It also applies the revised deterministic Gurobi settings.

Required
--------
--source_script  rq4_hepatic_integration_CORRECTED_v13.py
--hepatic_model  iMM1415.json
--expression_data C57BL/6J HFD/SCD expression CSV
--portal_metabolites portal_metabolites_for_hepatic_model.json
--condition ND_SCD or DD_HFD
"""
from __future__ import annotations
import argparse, importlib.util, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
sys.path.insert(0,str(HERE.parent/"common"))
from deterministic_solver_v2 import configure_solver_deterministic, canonicalize_solution

def load_module(path):
    spec=importlib.util.spec_from_file_location("rq4_source",path)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def optimize(mdl,use_pfba,pfba_fraction,solver):
    configure_solver_deterministic(mdl,solver=solver,verbose=False)
    bio=float(mdl.slim_optimize(error_value=None))
    if use_pfba:
        from cobra.flux_analysis import pfba
        sol=pfba(mdl,fraction_of_optimum=float(pfba_fraction))
    else:
        sol=mdl.optimize()
    sol=canonicalize_solution(sol,mdl)
    try:
        sol.biological_objective_value=bio
        sol.secondary_objective_value=float(sol.objective_value)
    except Exception: pass
    return sol,bio

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source_script",required=True)
    ap.add_argument("--hepatic_model",required=True)
    ap.add_argument("--expression_data",required=True)
    ap.add_argument("--portal_metabolites",required=True)
    ap.add_argument("--condition",required=True)
    ap.add_argument("--results_dir",required=True)
    ap.add_argument("--objective_mode",choices=["biomass","atpm","functional"],default="biomass")
    ap.add_argument("--pfba_fraction",type=float,default=1.0)
    ap.add_argument("--no-pfba",dest="use_pfba",action="store_false")
    ap.set_defaults(use_pfba=True)
    ap.add_argument("--solver",default="gurobi")
    ap.add_argument("--eflux_floor",type=float,default=0.1)
    ap.add_argument("--eflux_cap",type=float,default=1000.0)
    ap.add_argument("--gene_mapping")
    ap.add_argument("--diet_bounds")
    ap.add_argument("--portal_scaling",type=float,default=0.1)
    ap.add_argument("--portal_mode",choices=["availability","soft","hard","forced"],default="forced")
    ap.add_argument("--flux_threshold",type=float,default=0.01)
    args=ap.parse_args()
    if not (0 < args.pfba_fraction <= 1):
        raise ValueError("pfba_fraction must be in (0,1].")

    import cobra
    R=load_module(args.source_script)
    out=Path(args.results_dir); out.mkdir(parents=True,exist_ok=True)
    model=cobra.io.load_json_model(args.hepatic_model)
    portal=json.loads(Path(args.portal_metabolites).read_text())
    if args.condition in portal:
        pm=portal[args.condition]
    elif "portal_metabolites" in portal:
        pm=portal["portal_metabolites"].get(args.condition,{})
    else:
        pm={}
    print(f"[INFO] condition={args.condition}; portal metabolites={len(pm)}")

    common=dict(
        expression_file=args.expression_data,
        condition=args.condition,
        diet_file=args.diet_bounds,
        objective_mode=args.objective_mode,
        eflux_floor=args.eflux_floor,
        eflux_cap=args.eflux_cap,
        portal_scaling=args.portal_scaling,
        gene_mapping_file=args.gene_mapping,
        verbose=True
    )
    mb,sb=R.apply_corrected_eflux_with_microbiome(
        model,portal_metabolites=None,**common
    )
    mm,sm=R.apply_corrected_eflux_with_microbiome(
        model,portal_metabolites=pm,portal_mode=args.portal_mode,**common
    )
    solb,biob=optimize(mb,args.use_pfba,args.pfba_fraction,args.solver)
    solm,biom=optimize(mm,args.use_pfba,args.pfba_fraction,args.solver)
    if solb.status!="optimal" or solm.status!="optimal":
        raise RuntimeError(f"Optimization failed baseline={solb.status}, microbiome={solm.status}")

    fc=R.compare_flux_solutions(solb,solm,model,threshold=args.flux_threshold)
    fc.to_csv(out/f"{args.condition}_flux_comparison.csv",index=False)
    pd.DataFrame({
        "ReactionID":solb.fluxes.index,
        "baseline_flux":solb.fluxes.values,
        "microbiome_flux":solm.fluxes.reindex(solb.fluxes.index).values
    }).to_csv(out/f"{args.condition}_flux_vectors.csv",index=False)

    metadata={
        "condition":args.condition,
        "objective_mode":args.objective_mode,
        "use_pfba":args.use_pfba,
        "pfba_fraction":args.pfba_fraction,
        "portal_scaling":args.portal_scaling,
        "portal_mode":args.portal_mode,
        "flux_threshold":args.flux_threshold,
        "biological_objective_baseline":biob,
        "biological_objective_microbiome":biom,
        "biological_objective_delta":biom-biob,
        "secondary_pfba_objective_baseline":float(solb.objective_value),
        "secondary_pfba_objective_microbiome":float(solm.objective_value),
        "n_microbiome_attributable":int(fc["microbiome_attributable"].sum()) if "microbiome_attributable" in fc else None,
        "n_reactions":len(fc),
        "stats_baseline":sb,
        "stats_microbiome":sm,
    }
    # make JSON safe
    def clean(x):
        if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [clean(v) for v in x]
        if isinstance(x,(np.integer,)):return int(x)
        if isinstance(x,(np.floating,)):return float(x)
        if isinstance(x,np.ndarray):return x.tolist()
        return x if isinstance(x,(str,int,float,bool,type(None))) else str(x)
    (out/f"{args.condition}_run_metadata.json").write_text(json.dumps(clean(metadata),indent=2))
    print(f"[OK] {out}")

if __name__=="__main__":
    main()

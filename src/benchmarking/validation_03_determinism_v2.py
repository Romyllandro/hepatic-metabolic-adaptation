#!/usr/bin/env python3
"""
validation_03_determinism_v2.py
===============================
Determinism check aligned with revised pFBA production.

Builds one fully constrained three-layer representative model for a selected
within-study condition and solves it N times with the production-aligned
deterministic Gurobi configuration.
"""
from __future__ import annotations
import argparse, importlib.util, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"common"))
from deterministic_solver_v2 import deterministic_optimize, configure_solver_deterministic

def load_module(path):
    spec=importlib.util.spec_from_file_location("prod",path)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--flux_script",required=True)
    ap.add_argument("--rnaseq",required=True)
    ap.add_argument("--dataset",default="GSE101657")
    ap.add_argument("--condition",default="HFD")
    ap.add_argument("--model",required=True)
    ap.add_argument("--diet_bounds",required=True)
    ap.add_argument("--mapping",required=True)
    ap.add_argument("--objective",default="BIOMASS_mm_1_no_glygln")
    ap.add_argument("--n",type=int,default=10)
    ap.add_argument("--tol",type=float,default=1e-8)
    ap.add_argument("--solver",default="gurobi")
    ap.add_argument("--out",default="determinism_v2.csv")
    args=ap.parse_args()

    from cobra.io import load_json_model
    M=load_module(args.flux_script)
    rna=pd.read_csv(args.rnaseq)
    if "Gene_Symbol" not in rna.columns:
        raise ValueError(f"Gene_Symbol column not found: {list(rna.columns[:8])}")
    gene_symbols=rna["Gene_Symbol"].astype(str).str.lower().tolist()
    sample_cols=list(rna.columns[2:])
    hits=[x for x in sample_cols if x.startswith(f"{args.condition}_{args.dataset}_")]
    if not hits:raise ValueError("No matching sample columns.")
    X=rna[hits].apply(pd.to_numeric,errors="coerce").to_numpy(dtype=float)
    vec=np.nanmean(X,axis=1)
    sym=M.load_symbol_to_entrez_mapping(args.mapping)
    genes=[sym.get(str(s).lower().strip(),str(s)) for s in gene_symbols]
    with open(args.diet_bounds) as f:raw=json.load(f)
    diet={M.canonical_code(k):{rid:[float(v[0]),float(v[1])] for rid,v in d.items()} for k,d in raw.items()}
    base=load_json_model(args.model)
    base,_=M.set_objective_reaction(base,objective_id=args.objective,sense="max")
    _,trans,intset=M.classify_reactions(base,transporter_strategy="either")

    def build():
        m=base.copy()
        m,_=M.apply_diet_bounds_layer1(m,code=args.condition,diet_bounds=diet,diet_units="model")
        m,*_=M.apply_expression_constraints_scoped(m,genes,vec,trans,eflux_quantile=.95,
              eflux_floor=.1,eflux_cap=1000,label="L2",symbol_to_entrez=None)
        m,*_=M.apply_expression_constraints_scoped(m,genes,vec,intset,eflux_quantile=.95,
              eflux_floor=.1,eflux_cap=1000,label="L3",symbol_to_entrez=None)
        return M.validate_model(m)

    sols=[]
    bios=[]
    for i in range(args.n):
        sol=deterministic_optimize(build(),use_pfba=True,fraction_of_optimum=1.0,solver=args.solver)
        sols.append(sol.fluxes.values)
        bios.append(sol.biological_objective_value)
    arr=np.asarray(sols)
    deltas=np.max(np.abs(arr-arr[0]),axis=1)
    out=pd.DataFrame({"repeat":range(args.n),"max_abs_delta_vs_first":deltas,
                      "biological_objective":bios})
    out.to_csv(args.out,index=False)
    mx=float(deltas.max())
    print(out.to_string(index=False))
    print(f"max |delta flux|={mx:.3e} -> {'PASS' if mx<args.tol else 'FAIL'}")
    if mx>=args.tol:raise SystemExit(2)

if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""
Direct same-input, same-objective three-layer FBA versus pFBA test on GSE101657.

v2 fixes:
- temporary expression filename now obeys the production parser requirement:
    *_<GROUPS>_gene_expression.csv
- defaults to SCD,HFD,KD so both HFD-vs-SCD and KD-vs-SCD can be evaluated
  in the same controlled cohort.
"""
from __future__ import annotations
import argparse, subprocess, sys, re
from pathlib import Path
import pandas as pd

def run(c,cwd,dry=False):
    print("\n> "+" ".join(map(str,c)))
    if dry:return
    p=subprocess.run([str(x) for x in c],cwd=str(cwd))
    if p.returncode: raise RuntimeError(f"failed ({p.returncode})")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project_root",default=".")
    ap.add_argument("--modeling_script",required=True)
    ap.add_argument("--expression",default="GSEMERGED_SCD_HFD_KD_WD_gene_expression.csv")
    ap.add_argument("--model",default="iMM1415.json")
    ap.add_argument("--diet_bounds",default="expanded_diet_bounds_flat.json")
    ap.add_argument("--mapping",default="mouse_entrez_to_symbol.csv")
    ap.add_argument("--dataset",default="GSE101657")
    ap.add_argument("--conditions",default="SCD,HFD,KD")
    ap.add_argument("--baseline",default="SCD")
    ap.add_argument("--objective",default="BIOMASS_mm_1_no_glygln")
    ap.add_argument("--solver",default="gurobi")
    ap.add_argument("--quantile",type=float,default=0.95)
    ap.add_argument("--floor",type=float,default=0.1)
    ap.add_argument("--cap",type=float,default=1000)
    ap.add_argument("--output",default="reviewer_results/meeting_closure/FBA_pFBA_GSE101657")
    ap.add_argument("--dry-run",action="store_true")
    a=ap.parse_args()

    root=Path(a.project_root).resolve()
    out=root/a.output; out.mkdir(parents=True,exist_ok=True)
    expr=root/a.expression
    df=pd.read_csv(expr)

    requested=[x.strip().upper() for x in a.conditions.split(",") if x.strip()]
    if a.baseline.upper() not in requested:
        raise ValueError(f"Baseline {a.baseline} must be in --conditions")

    ann=[c for c in ["Gene_ID","Gene_Symbol"] if c in df.columns]
    if "Gene_Symbol" not in ann:
        raise ValueError("Expression matrix must contain Gene_Symbol.")

    sample_cols=[]
    present=[]
    for cond in requested:
        hits=[c for c in df.columns
              if str(c).upper().startswith(cond+"_") and a.dataset.upper() in str(c).upper()]
        if hits:
            present.append(cond); sample_cols.extend(hits)
        print(f"[INPUT] {a.dataset}/{cond}: {len(hits)} sample columns")

    if a.baseline.upper() not in present:
        raise ValueError(f"No {a.baseline} samples found for {a.dataset}.")
    if len(present)<2:
        raise ValueError(f"Need >=2 conditions; found {present}")

    # Critical filename contract used by the production model parser:
    # *_<GROUP1>_<GROUP2>_..._gene_expression.csv
    group_chunk="_".join(present)
    subset=out/f"{a.dataset}_{group_chunk}_gene_expression.csv"
    df[ann+sample_cols].to_csv(subset,index=False)
    print(f"[INPUT] Wrote parser-compatible subset: {subset.name}")

    core=[
        sys.executable, Path(a.modeling_script).resolve(), subset,
        "--model_file", (root/a.model).resolve(),
        "--results_dir", None,
        "--explicit_groups",",".join(present),
        "--baseline_code",a.baseline.upper(),
        "--test_type","t-test",
        "--diet_bounds_json",(root/a.diet_bounds).resolve(),
        "--eflux_quantile",str(a.quantile),"--eflux_floor",str(a.floor),"--eflux_cap",str(a.cap),
        "--objective_id",a.objective,"--objective_sense","max",
        "--transporter_strategy","either",
        "--mapping_file",(root/a.mapping).resolve(),
        "--solver",a.solver,"--fva_mode","none","--write_replicates_long"
    ]

    for mode in ["fba","pfba"]:
        c=list(core)
        rdir=out/mode
        c[c.index(None)]=rdir
        c += ["--solve_mode",mode]
        if mode=="pfba":
            c += ["--pfba_fraction","1.0"]
        run(c,root,a.dry_run)

    if not a.dry_run:
        fba=out/"fba/flux_analysis/reaction_flux_comparison_extended.csv"
        pfba=out/"pfba/flux_analysis/reaction_flux_comparison_extended.csv"
        cmp_script=Path(__file__).resolve().parent/"compare_fba_pfba.py"
        for test in present:
            if test==a.baseline.upper():
                continue
            cdir=out/"comparison"/f"{test}_vs_{a.baseline.upper()}"
            run([sys.executable,cmp_script,
                 "--fba",fba,"--pfba",pfba,
                 "--test",test,"--baseline",a.baseline.upper(),
                 "--output",cdir],root,False)
        print(f"[DONE] {out}")

if __name__=="__main__":
    main()

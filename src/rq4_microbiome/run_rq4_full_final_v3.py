#!/usr/bin/env python3
from __future__ import annotations
import argparse, subprocess, sys
from pathlib import Path

def call(cmd):
    print("\n> "+" ".join(map(str,cmd)))
    p=subprocess.run(cmd)
    if p.returncode:
        raise SystemExit(p.returncode)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--reviewed_map",required=True)
    ap.add_argument("--metatranscriptome",required=True)
    ap.add_argument("--agora_dir",required=True)
    ap.add_argument("--hepatic_model",required=True)
    ap.add_argument("--hepatic_expression",required=True)
    ap.add_argument("--gene_mapping",required=True)
    ap.add_argument("--diet_bounds",required=True)
    ap.add_argument("--output",default="reviewer_results/RQ4_final_full_v3")
    ap.add_argument("--solver",default="gurobi")
    ap.add_argument("--resume",action="store_true")
    ap.add_argument("--dry-run",action="store_true")
    args=ap.parse_args()

    here=Path(__file__).resolve().parent
    out=Path(args.output)
    maps=out/"maps"
    maps.mkdir(parents=True,exist_ok=True)

    call([
        sys.executable,str(here/"prepare_rq4_species_maps_final_v2.py"),
        "--reviewed_map",args.reviewed_map,
        "--agora_dir",args.agora_dir,
        "--output_dir",str(maps)
    ])

    cmd=[
        sys.executable,str(here/"run_rq4_sensitivity_suite_final_v3.py"),
        "--community_script",str(here/"rq4_microbiome_community_modeling_final_v3.py"),
        "--hepatic_runner",str(here/"rq4_hepatic_integration_final_v2.py"),
        "--hepatic_source_script",str(here/"rq4_hepatic_integration_CORRECTED_v14_source.py"),
        "--attribution_script",str(here/"rq4_attribution_analysis_source.py"),
        "--metatranscriptome",args.metatranscriptome,
        "--agora_dir",args.agora_dir,
        "--species_map_primary",str(maps/"rq4_species_map_primary_27.csv"),
        "--species_map_sensitivity",str(maps/"rq4_species_map_mapping_sensitivity_28.csv"),
        "--hepatic_model",args.hepatic_model,
        "--hepatic_expression",args.hepatic_expression,
        "--gene_mapping",args.gene_mapping,
        "--diet_bounds",args.diet_bounds,
        "--tradeoffs","0.3,0.5,0.7",
        "--portal_scales","0.05,0.10,0.20",
        "--pfba_fraction","1.0",
        "--solver",args.solver,
        "--output",str(out)
    ]
    if args.resume: cmd.append("--resume")
    if args.dry_run: cmd.append("--dry-run")
    call(cmd)

if __name__=="__main__":
    main()

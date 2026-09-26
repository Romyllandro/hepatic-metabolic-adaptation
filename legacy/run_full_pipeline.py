#!/usr/bin/env python3
"""
run_full_pipeline.py
====================
Convenience orchestrator for the manuscript reproduction pipeline.

The stage scripts remain independently runnable; this wrapper only resolves paths,
writes RQ1/RQ2 config files, and executes them in the frozen order.

Stages:
  preflight
  rq1
  rq2
  rq1_fva
  rq2_pfva
  benchmark
  ablation
  parameter_sensitivity
  determinism
  rq3
  rq4
  cross_rq
  figures
  verify

Optional:
  benchmark_ecoli   (requires external E. coli benchmark inputs)

Use --dry-run to print commands without executing.
"""
from __future__ import annotations
import argparse, json, subprocess, sys, shlex
from pathlib import Path

def load_json(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))

def cmd_text(cmd):
    return " ".join(shlex.quote(str(x)) for x in cmd)

def run(cmd,cwd,dry=False):
    print("\n$ "+cmd_text(cmd))
    if dry:return
    r=subprocess.run([str(x) for x in cmd],cwd=str(cwd))
    if r.returncode:
        raise RuntimeError(f"Command failed ({r.returncode}): {cmd_text(cmd)}")

def need(p,label):
    p=Path(p)
    if not p.exists(): raise FileNotFoundError(f"{label}: {p}")
    return p.resolve()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",default="pipeline_config.json")
    ap.add_argument("--stages",default="all",
                    help="all or comma-separated stage names")
    ap.add_argument("--dry-run",action="store_true")
    ap.add_argument("--resume",action="store_true")
    args=ap.parse_args()

    cfg=load_json(args.config)
    root=Path(cfg.get("project_root",".")).resolve()
    src=(root/cfg.get("source_root","src")).resolve()
    out=(root/cfg.get("results_root","reviewer_results")).resolve()
    out.mkdir(parents=True,exist_ok=True)
    resolved=out/"resolved_configs"; resolved.mkdir(exist_ok=True)

    py=cfg.get("python_executable",sys.executable)
    inp=cfg["inputs"]
    par=cfg.get("parameters",{})

    stage_order=[
        "preflight","rq1","rq2","rq1_fva","rq2_pfva",
        "benchmark","ablation","parameter_sensitivity","determinism",
        "rq3","rq4","cross_rq","figures","verify"
    ]
    requested=stage_order if args.stages=="all" else [x.strip() for x in args.stages.split(",") if x.strip()]
    bad=[x for x in requested if x not in stage_order+["benchmark_ecoli"]]
    if bad: raise ValueError(f"Unknown stages: {bad}")

    def P(key):
        x=Path(inp[key])
        return x if x.is_absolute() else root/x

    if "preflight" in requested:
        required=["model","diet_bounds","gene_mapping","rq1_expression","rq4_metatranscriptome",
                  "rq4_reviewed_species_map","rq3_sc_data","rq3_metadata"]
        for k in required: need(P(k),k)
        need(Path(inp["agora_dir"]),"agora_dir")
        for strain,path in inp["rq2_strain_expression"].items():
            p=Path(path); p=p if p.is_absolute() else root/p
            need(p,f"RQ2 {strain}")
        print("[PASS] input-file preflight")
        if not args.dry_run:
            # Runtime import preflight
            for pkg in ["pandas","numpy","scipy","statsmodels","sklearn","matplotlib","cobra","micom","gurobipy"]:
                run([py,"-c",f"import {pkg}; print('{pkg}: OK')"],root,False)

    # Resolve common paths
    model=P("model"); diet=P("diet_bounds"); mapping=P("gene_mapping")
    rq1_expr=P("rq1_expression")
    rq1_out=out/"RERUN_RQ1_PFBA"
    rq2_out=out/"RERUN_RQ2_PFBA_FVA"

    # RQ1
    if "rq1" in requested:
        rq1cfg={
            "project_root":str(root),
            "python_executable":py,
            "input_expression":str(rq1_expr),
            "model_file":str(model),
            "diet_bounds_json":str(diet),
            "mapping_file":str(mapping),
            "results_dir":str(rq1_out),
            "scripts":{
                "flux":str(src/"rq1_bulk/map_fixv5_multigroupsv8_layered_manuscript_rerun.py"),
                "batch":str(src/"rq1_bulk/batch_correct_flux_by_dataset_3d_revised.py"),
                "comprehensive":str(src/"rq1_bulk/comprehensive_flux_analysis_revised.py"),
                "combined":str(src/"rq1_bulk/combined_analysis_interactive_generic_revised.py"),
                # Kept for runner compatibility; primary representative FVA is also refined in rq1_fva.
                "fva_summary":str(src/"rq1_bulk/run_rq1_cohort_representative_fva.py"),
            },
            "test_type":"t-test","eflux_quantile":0.95,"eflux_floor":0.1,"eflux_cap":1000,
            "objective_id":"BIOMASS_mm_1_no_glygln","objective_sense":"max",
            "transporter_strategy":"either","edge_abs_diff_threshold":0.2,
            "solve_mode":"pfba","pfba_fraction":1.0,"solver":par.get("solver","gurobi"),
            # Avoid the lightweight legacy representative-FVA summarizer inside the orchestrator.
            # The publication refinement is executed in the rq1_fva stage.
            "fva_mode":"none","fva_fraction":1.0,"fva_processes":1,
            "permutations":999,"fdr":0.05,"effect_size":0.5,
            "expected_group_counts":{"SCD":37,"HFD":38,"KD":12,"WD":12},
        }
        p=resolved/"rq1_config.json"; p.write_text(json.dumps(rq1cfg,indent=2),encoding="utf-8")
        run([py,src/"rq1_bulk/run_flux_pipeline_RQ1_revision.py","--config",p],root,args.dry_run)

    # RQ2
    rq2cfg={
        "runner":{
            "project_root":str(root),"python_executable":py,
            "script":str(src/"rq1_bulk/map_fixv5_multigroupsv8_layered_manuscript_rerun.py"),
            "log_dir":str(rq2_out/"logs_layered_runs"),
            "parallel":1,"stop_on_error":True,
            "success_marker":".runner_success.json",
        },
        "common_args":{
            "model_file":str(model),"test_type":"t-test","diet_bounds_json":str(diet),
            "eflux_quantile":0.95,"eflux_floor":0.1,"eflux_cap":1000,
            "objective_id":"BIOMASS_mm_1_no_glygln","objective_sense":"max",
            "transporter_strategy":"either","edge_abs_diff_threshold":0.2,
            "mapping_file":str(mapping),"write_replicates_long":True,
            "solve_mode":"pfba","pfba_fraction":1.0,"solver":par.get("solver","gurobi"),
            "fva_mode":"representative","fva_fraction":1.0,"fva_processes":1,
        },
        "jobs":[]
    }
    for strain,path in inp["rq2_strain_expression"].items():
        pp=Path(path); pp=pp if pp.is_absolute() else root/pp
        rq2cfg["jobs"].append({
            "name":strain,"input_csv":str(pp),
            "results_dir":str(rq2_out/f"results_{strain}_GSE182668")
        })
    rq2cfg_path=resolved/"rq2_config.json"
    rq2cfg_path.write_text(json.dumps(rq2cfg,indent=2),encoding="utf-8")

    if "rq2" in requested:
        c=[py,src/"rq2_strains/run_layered_executor_revised.py","--config",rq2cfg_path]
        if args.resume:c.append("--resume")
        run(c,root,args.dry_run)
        run([py,src/"rq2_strains/rq2_cross_strain_analysis_revised.py",
             "--config",rq2cfg_path,"--output",rq2_out/"cross_strain_analysis",
             "--theta","0.20","--thresholds","0.05,0.10,0.15,0.20,0.25,0.30,0.50",
             "--fdr","0.05","--effect_size","0.5"],root,args.dry_run)

    if "rq1_fva" in requested:
        run([py,src/"rq1_bulk/run_rq1_cohort_representative_fva.py",
             "--expression",rq1_expr,"--model",model,"--diet_bounds",diet,"--mapping",mapping,
             "--modeling_script",src/"rq1_bulk/map_fixv5_multigroupsv8_layered_manuscript_rerun.py",
             "--rq1_results",rq1_out,
             "--output",out/"cohort_representative_FVA",
             "--fva_fraction","1.0","--parsimony_eps","0.01,0.05",
             "--tolerance","1e-7","--processes","1"],root,args.dry_run)

    if "rq2_pfva" in requested:
        run([py,src/"rq2_strains/run_rq2_targeted_parsimonious_fva.py",
             "--config",rq2cfg_path,
             "--modeling_script",src/"rq1_bulk/map_fixv5_multigroupsv8_layered_manuscript_rerun.py",
             "--rq2_results",rq2_out,
             "--output",out/"targeted_parsimonious_FVA",
             "--parsimony_eps","0.01,0.05","--tolerance","1e-7","--processes","1"],root,args.dry_run)

    # Benchmark
    bench_out=out/"benchmark_liver_GSE101657"
    if "benchmark" in requested:
        run([py,src/"benchmarking/benchmark_liver_methods.py",
             "--flux_script",src/"rq1_bulk/map_fixv5_multigroupsv8_layered_manuscript_rerun.py",
             "--rnaseq",rq1_expr,"--dataset","GSE101657","--conditions","SCD,HFD",
             "--baseline","SCD","--test","HFD","--model",model,"--diet_bounds",diet,
             "--mapping",mapping,"--positive_controls",root/"configs/positive_controls.json",
             "--methods","three_layer,eflux,gimme,pfba",
             "--objectives","BIOMASS_mm_1_no_glygln,ATPM","--solver",par.get("solver","gurobi"),
             "--output",bench_out],root,args.dry_run)
        if not args.dry_run:
            run([py,src/"benchmarking/recompute_benchmark_stats_robust_v1.py",
                 "--benchmark_dir",bench_out,"--test","HFD","--baseline","SCD"],root,False)

    if "ablation" in requested:
        run([py,src/"benchmarking/validation_02_layer_ablation_v2.py",
             "--flux_script",src/"rq1_bulk/map_fixv5_multigroupsv8_layered_manuscript_rerun.py",
             "--rnaseq",rq1_expr,"--dataset","GSE101657","--conditions","SCD,HFD,KD",
             "--baseline","SCD","--model",model,"--diet_bounds",diet,"--mapping",mapping,
             "--objective","BIOMASS_mm_1_no_glygln","--out",out/"layer_ablation_v2.csv"],root,args.dry_run)

    if "parameter_sensitivity" in requested:
        run([py,src/"benchmarking/sweep_threelayer_pipeline_v2.py",
             "--flux_script",src/"rq1_bulk/map_fixv5_multigroupsv8_layered_manuscript_rerun.py",
             "--model",model,"--rnaseq",rq1_expr,"--dataset","GSE101657",
             "--conditions","SCD,HFD,KD","--mapping",mapping,"--diet_bounds",diet,
             "--positive_controls",root/"configs/positive_controls.json",
             "--quantiles","0.90,0.95,0.99","--floors","0.05,0.1","--caps","1000",
             "--objectives","BIOMASS_mm_1_no_glygln,ATPM",
             "--results_dir",out/"threelayer_sensitivity_GSE101657"],root,args.dry_run)

    if "determinism" in requested:
        run([py,src/"benchmarking/validation_03_determinism_v2.py",
             "--flux_script",src/"rq1_bulk/map_fixv5_multigroupsv8_layered_manuscript_rerun.py",
             "--rnaseq",rq1_expr,"--dataset","GSE101657","--condition","HFD",
             "--model",model,"--diet_bounds",diet,"--mapping",mapping,
             "--objective","BIOMASS_mm_1_no_glygln","--n","10","--solver",par.get("solver","gurobi"),
             "--out",out/"determinism_v2.csv"],root,args.dry_run)

    # Optional external E. coli benchmark
    if "benchmark_ecoli" in requested:
        eco=cfg.get("ecoli_benchmark",{})
        for k in ["model","transcriptomes","measured_modules","uptakes"]:
            if not eco.get(k): raise ValueError(f"ecoli_benchmark.{k} is required")
        run([py,src/"benchmarking/benchmark_ecoli_13c_methods.py",
             "--model",eco["model"],"--objective",eco.get("objective","BIOMASS_Ec_iJO1366_core_53p95M"),
             "--transcriptomes",eco["transcriptomes"],"--measured_modules",eco["measured_modules"],
             "--uptakes",eco["uptakes"],"--regime","both","--solver",par.get("solver","gurobi"),
             "--output",out/"benchmark_ecoli_13c"],root,args.dry_run)

    # RQ3
    rq3_out=out/"RQ3_sensitivity"
    if "rq3" in requested:
        run([py,src/"rq3_single_cell/run_rq3_sensitivity_suite.py",
             "--rq3_script",src/"rq3_single_cell/rq3_FINAL_COMPLETE_REVISED.py",
             "--sc_data",P("rq3_sc_data"),"--sc_metadata",P("rq3_metadata"),
             "--diet_bounds_file",diet,
             "--condition_mapping",json.dumps(inp.get("rq3_condition_mapping",{"Chow":"SCD","WesternDiet":"WD"})),
             "--model",model,"--quantiles","0.95,0.99",
             "--objectives","BIOMASS_mm_1_no_glygln,ATPM",
             "--output",rq3_out],root,args.dry_run)

    # RQ4
    rq4_out=out/"RQ4_final_full_v3"
    if "rq4" in requested:
        c=[py,src/"rq4_microbiome/run_rq4_full_final_v3.py",
           "--reviewed_map",P("rq4_reviewed_species_map"),
           "--metatranscriptome",P("rq4_metatranscriptome"),
           "--agora_dir",Path(inp["agora_dir"]),
           "--hepatic_model",model,
           "--hepatic_expression",P("rq4_hepatic_expression"),
           "--gene_mapping",mapping,"--diet_bounds",diet,
           "--solver",par.get("solver","gurobi"),"--output",rq4_out]
        if args.resume:c.append("--resume")
        run(c,root,args.dry_run)

    # Cross-RQ
    cross_out=out/"final_cross_rq"
    if "cross_rq" in requested:
        run([py,src/"cross_rq/run_cross_rq_final.py",
             "--rq1",rq1_out,
             "--rq2",rq2_out,
             "--rq3_results",rq3_out/"q99_BIOMASS_mm_1_no_glygln",
             "--rq4_results",rq4_out,
             "--rq3_comparison","WesternDiet_vs_Chow",
             "--anchor_abs_diff","0.20","--rq4_flux_threshold","0.01",
             "--strict_expected_counts","--output",cross_out],root,args.dry_run)

    # Figures
    if "figures" in requested:
        figout=out/"final_figures"
        run([py,src/"reporting/generate_figures_rq1_rq2.py",
             "--rq1",rq1_out,"--rq2",rq2_out,
             "--rq2-pfva",out/"targeted_parsimonious_FVA",
             "--output",figout/"rq1_rq2"],root,args.dry_run)
        run([py,src/"reporting/generate_final_figures_rq4_crossrq.py",
             "--rq4",rq4_out,"--cross_rq",cross_out,
             "--output",figout/"rq4_crossrq"],root,args.dry_run)

    if "verify" in requested:
        run([py,src/"utils/verify_final_counts.py",
             "--cross-summary",cross_out/"cross_RQ_summary_final.json"],root,args.dry_run)

    print("\n[DONE] requested pipeline stages finished.")

if __name__=="__main__":
    main()

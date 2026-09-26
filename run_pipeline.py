#!/usr/bin/env python3
"""
run_pipeline.py -- end-to-end orchestrator for the multi-scale hepatic metabolic
modeling pipeline (Madadjim, Vechetti & Cui; Briefings in Bioinformatics, BIB-26-1735).

The stage scripts under src/ are the exact scripts used for the manuscript and remain
independently runnable. This wrapper only (1) resolves paths from pipeline_config.json,
(2) writes the job configs the stage scripts expect, and (3) calls them in dependency order.

PRIMARY analyses (numbers reported as primary in the manuscript)
  preflight        check inputs and Python packages
  rq1              bulk multi-cohort standard FBA (99 samples) -> dataset offset removal ->
                   Welch/BH statistics, PCA, PERMANOVA -> signature summary (262 / 10)
  rq2              nine-strain standard FBA, aggregate design (primary tiers 243 / 11)
                   + replicate-preserving design (187 / 15) -> conservation tiers,
                   theta sweep, within-strain and between-strain statistics
  rq3              single-cell pseudobulk pFBA (P95/P99 x biomass/ATPM x |dv| thresholds)
  rq3_integration  abundance-weighted attribution, bulk comparison, hierarchy
  rq4              MICOM community (27-species strict map) -> hepatic coupling -> attribution
  validation       positive controls, liver benchmark (layered / E-Flux / GIMME / diet-only pFBA),
                   layer ablation, parameter sweep, determinism
  closure_tests    same-input FBA vs pFBA, E-Flux cap sensitivity, ATPM FVA, RQ4 pFBA on/off
  figures          Figures 2-7 (Figure 1 is a schematic)
  verify           recompute the manuscript headline numbers and compare with expectations

SENSITIVITY analyses (reported as sensitivity / robustness in the manuscript)
  rq1_pfba, rq2_pfba, rq1_fva, rq2_pfva, cross_rq_pfba

OPTIONAL (needs the external E. coli inputs shipped in inputs/ecoli_benchmark)
  benchmark_ecoli  five-method 13C-MFA benchmark (Supplementary Table SL3)

Examples
  python run_pipeline.py --config pipeline_config.json --stages preflight
  python run_pipeline.py --config pipeline_config.json --stages all --dry-run
  python run_pipeline.py --config pipeline_config.json --stages rq1,rq2,verify
  python run_pipeline.py --config pipeline_config.json --stages verify --use-reference   # no solver needed
"""
from __future__ import annotations
import argparse, json, os, shlex, subprocess, sys
from pathlib import Path

PRIMARY = ["preflight", "rq1", "rq2", "rq3", "rq3_integration", "rq4",
           "validation", "closure_tests", "figures", "verify"]
SENSITIVITY = ["rq1_pfba", "rq2_pfba", "rq1_fva", "rq2_pfva", "cross_rq_pfba"]
OPTIONAL = ["benchmark_ecoli"]
STRAINS = ["129S1SvImJ", "AJ", "C57BL6J", "CASTEiJ", "DBA2J", "NODShiLtJ", "NZOHlLtJ", "PWKPhJ", "WSBEiJ"]


def q(x): return " ".join(shlex.quote(str(i)) for i in x)


class Runner:
    def __init__(self, root, dry, env=None):
        self.root, self.dry, self.env = root, dry, env

    def __call__(self, cmd, cwd=None, extra_env=None):
        print("\n$ " + q(cmd), flush=True)
        if self.dry:
            return
        env = dict(os.environ); env.update(self.env or {}); env.update(extra_env or {})
        r = subprocess.run([str(c) for c in cmd], cwd=str(cwd or self.root), env=env)
        if r.returncode:
            raise SystemExit(f"[FAILED] exit {r.returncode}: {q(cmd)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="pipeline_config.json")
    ap.add_argument("--stages", default="all",
                    help="'all' (primary), 'all+sensitivity', 'everything', or a comma list of stage names")
    ap.add_argument("--dry-run", action="store_true", help="print commands only")
    ap.add_argument("--resume", action="store_true", help="skip finished strain jobs (RQ2) / RQ4 scenarios")
    ap.add_argument("--use-reference", action="store_true",
                    help="verify/figures: read the frozen tables in reference_outputs/ instead of results/")
    a = ap.parse_args()

    cfg = json.loads(Path(a.config).read_text(encoding="utf-8"))
    root = Path(cfg.get("project_root", ".")).resolve()
    src = (root / cfg.get("source_root", "src")).resolve()
    out = (root / cfg.get("results_root", "results")).resolve()
    ref = root / "reference_outputs"
    py = cfg.get("python_executable", sys.executable)
    par = cfg.get("parameters", {}); inp = cfg["inputs"]
    solver = par.get("solver", "gurobi"); obj = par.get("objective", "BIOMASS_mm_1_no_glygln")
    run = Runner(root, a.dry_run, {"PYTHONPATH": os.pathsep.join(
        [str(src / d) for d in ["common", "benchmarking", "rq1_bulk/sensitivity_pfba", "rq4_microbiome"]])})

    if a.stages == "all": stages = PRIMARY
    elif a.stages == "all+sensitivity": stages = PRIMARY[:-2] + SENSITIVITY + PRIMARY[-2:]
    elif a.stages == "everything": stages = PRIMARY[:-2] + SENSITIVITY + OPTIONAL + PRIMARY[-2:]
    else: stages = [s.strip() for s in a.stages.split(",") if s.strip()]
    bad = [s for s in stages if s not in PRIMARY + SENSITIVITY + OPTIONAL]
    if bad: raise SystemExit(f"Unknown stage(s): {bad}")

    def P(key, base=inp):
        p = Path(base[key]); return p if p.is_absolute() else root / p
    model, diet, mapping = P("model"), P("diet_bounds"), P("gene_mapping")
    rq1_expr = P("rq1_expression")
    strain_expr = {s: (Path(v) if Path(v).is_absolute() else root / v) for s, v in inp["rq2_strain_expression"].items()}
    out.mkdir(parents=True, exist_ok=True); cfgdir = out / "resolved_configs"; cfgdir.mkdir(exist_ok=True)
    fba = src / "rq1_bulk/primary_fba/map_fixv5_multigroupsv8_layered_manuscript_run.py"
    pfba = src / "rq1_bulk/sensitivity_pfba/map_fixv5_multigroupsv8_layered_manuscript_rerun.py"
    common_flux = ["--model_file", model, "--diet_bounds_json", diet, "--mapping_file", mapping,
                   "--eflux_quantile", par.get("eflux_quantile", 0.95), "--eflux_floor", par.get("eflux_floor", 0.1),
                   "--eflux_cap", par.get("eflux_cap", 1000), "--objective_id", obj, "--objective_sense", "max",
                   "--transporter_strategy", "either",
                   "--edge_abs_diff_threshold", par.get("edge_abs_diff_threshold", 0.2), "--write_replicates_long"]

    R1 = out / "RQ1_primary_FBA"; R2A = out / "RQ2_primary_FBA_aggregate"; R2R = out / "RQ2_FBA_replicate"
    R3 = out / "RQ3_sensitivity"; R3I = out / "RQ3_integration"; R4 = out / "RQ4_final_full_v3"
    FLAT = out / "RQs_Output_Files"   # flat copy with manuscript file names (production: RQs_Output_Files_final)

    # ------------------------------------------------------------------ preflight
    if "preflight" in stages:
        missing = []
        for k in ["model", "diet_bounds", "gene_mapping", "rq1_expression", "rq4_metatranscriptome",
                  "rq4_reviewed_species_map", "rq3_sc_data", "rq3_metadata"]:
            if not P(k).exists(): missing.append(f"{k}: {P(k)}")
        for s, p in strain_expr.items():
            if not p.exists(): missing.append(f"RQ2 {s}: {p}")
        if not Path(inp["agora_dir"]).exists(): missing.append(f"agora_dir: {inp['agora_dir']}")
        for m in missing: print("[MISSING]", m)
        print("[NOTE] RQ3 needs the single-cell matrix and RQ4 needs AGORA2 (see docs/DATA_SOURCES.md);"
              " other stages run without them." if missing else "[PASS] all inputs present")
        for pkg in ["pandas", "numpy", "scipy", "statsmodels", "sklearn", "matplotlib", "openpyxl", "cobra", "micom", "gurobipy"]:
            run([py, "-c", f"import {pkg}; print('{pkg}', getattr({pkg}, '__version__', 'ok'))"])

    # ------------------------------------------------------------------ RQ1 primary (standard FBA)
    if "rq1" in stages:
        # Production (11 Apr 2026): run_flux_pipeline_v3.py -> flux script, dataset-offset removal,
        # combined analysis, then comprehensive_flux_analysis.py on the offset-removed table.
        run([py, fba, rq1_expr, "--results_dir", R1, "--test_type", "mann-whitney",
             "--rank_product_for", "HFD,KD,WD", "--no_fva", *common_flux])
        flux = R1 / "flux_analysis/reaction_flux_comparison_extended.csv"
        flux_bc = R1 / "flux_analysis/reaction_flux_comparison_extended_batch_corrected.csv"
        run([py, src / "rq1_bulk/primary_fba/batch_correct_flux_by_dataset_3d.py", "--input", flux,
             "--output_csv", flux_bc, "--outdir", R1 / "batch_correction_outputs"])
        run([py, src / "rq1_bulk/primary_fba/combined_analysis_interactive_generic.py", "--input", flux_bc,
             "--output_dir", R1 / "combined_analysis", "--no_interactive"])
        # The manuscript statistics were produced on the dataset-offset-removed table (production log,
        # 11 Apr 2026). Set parameters.rq1_inference_matrix = "uncorrected" to run the same Welch/BH
        # analysis on the original feasible per-sample fluxes instead (see docs/PROVENANCE_AND_KNOWN_ISSUES.md).
        inf = flux if par.get("rq1_inference_matrix", "batch_corrected") == "uncorrected" else flux_bc
        run([py, src / "rq1_bulk/primary_fba/comprehensive_flux_analysis.py", "--input", inf,
             "--output", R1 / "comprehensive_analysis"])
        run([py, src / "rq1_bulk/primary_fba/rq1_signature_summary.py",
             "--stats_dir", R1 / "comprehensive_analysis/csv_outputs", "--output", R1 / "summary"])
        run([py, src / "utils/flatten_outputs.py", "--out", FLAT, "--rq1", R1, "--annotations", P("reaction_annotations")])

    # ------------------------------------------------------------------ RQ2 primary (standard FBA)
    def executor_cfg(name, results, common):
        c = {"runner": {"project_root": str(root), "python_executable": py, "script": str(fba),
                        "log_dir": str(results / "logs_layered_runs"), "parallel": 1, "stop_on_error": True,
                        "success_marker": ".runner_success.json"},
             "common_args": common,
             "jobs": [{"name": s, "input_csv": str(strain_expr[s]), "results_dir": str(results / f"results_{s}_GSE182668")}
                      for s in STRAINS]}
        p = cfgdir / f"{name}.json"; p.write_text(json.dumps(c, indent=2)); return p
    base_args = {"model_file": str(model), "diet_bounds_json": str(diet), "mapping_file": str(mapping),
                 "eflux_quantile": par.get("eflux_quantile", 0.95), "eflux_floor": par.get("eflux_floor", 0.1),
                 "eflux_cap": par.get("eflux_cap", 1000), "objective_id": obj, "objective_sense": "max",
                 "transporter_strategy": "either", "edge_abs_diff_threshold": par.get("edge_abs_diff_threshold", 0.2),
                 "write_replicates_long": True, "rank_product_for": "HFD,KD,WD", "no_fva": False}
    if "rq2" in stages:
        agg = executor_cfg("rq2_fba_aggregate", R2A, {**base_args, "test_type": "mann-whitney", "aggregate": True})
        rep = executor_cfg("rq2_fba_replicate", R2R, {**base_args, "test_type": "t-test", "aggregate": False})
        for c in (agg, rep):
            cmd = [py, src / "rq2_strains/primary_fba/run_layered_executor.py", "--config", c]
            if a.resume: cmd.append("--resume")
            run(cmd)
        run([py, src / "rq2_strains/primary_fba/rq2_conservation_tiers.py", "--rq2_root", R2A,
             "--output", out / "RQ2_tiers_aggregate"])
        run([py, src / "rq2_strains/primary_fba/rq2_conservation_tiers.py", "--rq2_root", R2R,
             "--output", out / "RQ2_tiers_replicate"])
        run([py, src / "rq2_strains/primary_fba/rq2_layer3_strain_specificity.py", "--aggregate_dir", R2A,
             "--replicate_dir", R2R, "--edge_threshold", "0.2", "--output_dir", out / "RQ2_layer3_strain_specificity"])
        run([py, src / "rq2_strains/primary_fba/rq2_build_evidence_table.py", "--aggregate_dir", R2A,
             "--replicate_dir", R2R, "--layer3_dir", out / "RQ2_layer3_strain_specificity",
             "--output_dir", out / "RQ2_evidence_table"])
        run([py, src / "utils/flatten_outputs.py", "--out", FLAT, "--rq2", R2A])

    # ------------------------------------------------------------------ RQ3
    if "rq3" in stages:
        run([py, src / "rq3_single_cell/run_rq3_sensitivity_suite.py",
             "--rq3_script", src / "rq3_single_cell/rq3_FINAL_COMPLETE_REVISED.py",
             "--sc_data", P("rq3_sc_data"), "--sc_metadata", P("rq3_metadata"), "--diet_bounds_file", diet,
             "--condition_mapping", json.dumps(inp.get("rq3_condition_mapping", {"Chow": "SCD", "WesternDiet": "WD"})),
             "--model", model, "--quantiles", "0.95,0.99", "--objectives", f"{obj},ATPM",
             "--abs_thresholds", "0.05,0.10,0.20", "--output", R3])
    if "rq3_integration" in stages:
        rq3ref = R3 / f"q99_{obj}"
        rq1stats = FLAT / "RQ1_flux_pairwise_stats.csv"
        run([py, src / "rq3_single_cell/integration/rq1_rq2_rq3_integration_analysis_REVISED.py",
             "--rq1_pairwise_stats", rq1stats, "--rq3_stats", rq3ref / "statistics/statistical_tests.csv",
             "--rq3_aggregation", rq3ref / "aggregation_summary.csv", "--bulk_comparison", "WD_vs_SCD",
             "--cellular_comparison", "WesternDiet_vs_Chow", "--output_dir", R3I])
        run([py, src / "rq3_single_cell/integration/rq2_rq3_multi_strain_integration_ENHANCED.py",
             "--base_dir", FLAT, "--rq1_stats", rq1stats, "--rq3_stats", rq3ref / "statistics/statistical_tests.csv",
             "--rq3_aggregation", rq3ref / "aggregation_summary.csv", "--output_dir", out / "RQ3_cross_strain_projection"])
        env = {"BULK_FLUX_FILE": str(R1 / "flux_analysis/reaction_flux_comparison_extended_batch_corrected.csv"),
               "SC_STATS_FILE": str(rq3ref / "statistics/statistical_tests.csv"),
               "ABUNDANCE_FILE": str(R3I / "tables/cell_abundance_used.csv"),
               "OUTPUT_DIR": str(out / "RQ3_bulk_validation")}
        run([py, src / "rq3_single_cell/integration/bulk_validation_LOCAL.py"], extra_env=env)
        run([py, src / "rq3_single_cell/integration/hierarchical_attribution_LOCAL.py"],
            extra_env={"CONTRIBUTION_FILE": str(R3I / "tables/phase2_contribution_analysis.csv"),
                       "OUTPUT_DIR": str(out / "RQ3_hierarchical_attribution")})

    # ------------------------------------------------------------------ RQ4
    if "rq4" in stages:
        cmd = [py, src / "rq4_microbiome/run_rq4_full_final_v3.py", "--reviewed_map", P("rq4_reviewed_species_map"),
               "--metatranscriptome", P("rq4_metatranscriptome"), "--agora_dir", Path(inp["agora_dir"]),
               "--hepatic_model", model, "--hepatic_expression", P("rq4_hepatic_expression"),
               "--gene_mapping", mapping, "--diet_bounds", diet, "--solver", solver, "--output", R4]
        if a.resume: cmd.append("--resume")
        run(cmd)

    # ------------------------------------------------------------------ validation (reviewer requests)
    if "validation" in stages:
        pc = P("positive_controls")
        run([py, src / "validation/positive_controls/validation_04_positive_controls.py",
             "--data_dir", ref / "RQ1_primary_FBA" if a.use_reference else FLAT,
             "--out", out / "validation/positive_control_report.csv"])
        bench = out / "validation/benchmark_liver_GSE101657"
        run([py, src / "benchmarking/benchmark_liver_methods.py", "--flux_script", pfba, "--rnaseq", rq1_expr,
             "--dataset", "GSE101657", "--conditions", "SCD,HFD", "--baseline", "SCD", "--test", "HFD",
             "--model", model, "--diet_bounds", diet, "--mapping", mapping, "--positive_controls", pc,
             "--methods", "three_layer,eflux,gimme,pfba", "--objectives", f"{obj},ATPM", "--solver", solver,
             "--output", bench])
        run([py, src / "benchmarking/recompute_benchmark_stats_robust_v1.py", "--benchmark_dir", bench,
             "--test", "HFD", "--baseline", "SCD"])
        run([py, src / "benchmarking/validation_02_layer_ablation_v2.py", "--flux_script", pfba, "--rnaseq", rq1_expr,
             "--dataset", "GSE101657", "--conditions", "SCD,HFD,KD", "--baseline", "SCD", "--model", model,
             "--diet_bounds", diet, "--mapping", mapping, "--objective", obj, "--out", out / "validation/layer_ablation_v2.csv"])
        run([py, src / "benchmarking/sweep_threelayer_pipeline_v2.py", "--flux_script", pfba, "--model", model,
             "--rnaseq", rq1_expr, "--dataset", "GSE101657", "--conditions", "SCD,HFD,KD", "--mapping", mapping,
             "--diet_bounds", diet, "--positive_controls", pc, "--quantiles", "0.90,0.95,0.99", "--floors", "0.05,0.1",
             "--caps", "1000", "--objectives", f"{obj},ATPM", "--results_dir", out / "validation/threelayer_sensitivity_GSE101657"])
        run([py, src / "benchmarking/validation_03_determinism_v2.py", "--flux_script", pfba, "--rnaseq", rq1_expr,
             "--dataset", "GSE101657", "--condition", "HFD", "--model", model, "--diet_bounds", diet, "--mapping", mapping,
             "--objective", obj, "--n", "10", "--solver", solver, "--out", out / "validation/determinism_v2.csv"])

    if "closure_tests" in stages:
        run([py, src / "validation/meeting_closure_tests/run_all_remaining_meeting_tests.py",
             "--project_root", root, "--pipeline_root", root, "--agora_dir", Path(inp["agora_dir"]),
             "--expression", rq1_expr, "--model", model, "--diet_bounds", diet, "--mapping", mapping,
             "--positive_controls", P("positive_controls"), "--frozen_rq1", out / "RQ1_pFBA_sensitivity",
             "--frozen_rq4", R4, "--rq4_metatranscriptome", P("rq4_metatranscriptome"),
             "--rq4_hepatic_expression", P("rq4_hepatic_expression"), "--solver", solver])

    # ------------------------------------------------------------------ sensitivity (pFBA / FVA)
    rq1p = out / "RQ1_pFBA_sensitivity"; rq2p = out / "RQ2_pFBA_sensitivity"
    if "rq1_pfba" in stages:
        c = {"project_root": str(root), "python_executable": py, "input_expression": str(rq1_expr),
             "model_file": str(model), "diet_bounds_json": str(diet), "mapping_file": str(mapping), "results_dir": str(rq1p),
             "scripts": {"flux": str(pfba),
                         "batch": str(src / "rq1_bulk/sensitivity_pfba/batch_correct_flux_by_dataset_3d_revised.py"),
                         "comprehensive": str(src / "rq1_bulk/sensitivity_pfba/comprehensive_flux_analysis_revised.py"),
                         "combined": str(src / "rq1_bulk/sensitivity_pfba/combined_analysis_interactive_generic_revised.py"),
                         "fva_summary": str(src / "rq1_bulk/sensitivity_pfba/run_rq1_cohort_representative_fva.py")},
             "test_type": "t-test", "eflux_quantile": 0.95, "eflux_floor": 0.1, "eflux_cap": 1000, "objective_id": obj,
             "objective_sense": "max", "transporter_strategy": "either", "edge_abs_diff_threshold": 0.2,
             "solve_mode": "pfba", "pfba_fraction": 1.0, "solver": solver, "fva_mode": "none", "fva_fraction": 1.0,
             "fva_processes": 1, "permutations": 999, "fdr": 0.05, "effect_size": 0.5,
             "expected_group_counts": {"SCD": 37, "HFD": 38, "KD": 12, "WD": 12}}
        p = cfgdir / "rq1_pfba.json"; p.write_text(json.dumps(c, indent=2))
        run([py, src / "rq1_bulk/sensitivity_pfba/run_flux_pipeline_RQ1_revision.py", "--config", p])
    rq2pcfg = cfgdir / "rq2_pfba.json"
    if any(s in stages for s in ["rq2_pfba", "rq2_pfva"]):
        c = {"runner": {"project_root": str(root), "python_executable": py, "script": str(pfba),
                        "log_dir": str(rq2p / "logs_layered_runs"), "parallel": 1, "stop_on_error": True,
                        "success_marker": ".runner_success.json"},
             "common_args": {**{k: v for k, v in base_args.items() if k not in ("no_fva", "rank_product_for")},
                             "test_type": "t-test", "solve_mode": "pfba", "pfba_fraction": 1.0, "solver": solver,
                             "fva_mode": "representative", "fva_fraction": 1.0, "fva_processes": 1},
             "jobs": [{"name": s, "input_csv": str(strain_expr[s]), "results_dir": str(rq2p / f"results_{s}_GSE182668")}
                      for s in STRAINS]}
        rq2pcfg.write_text(json.dumps(c, indent=2))
    if "rq2_pfba" in stages:
        cmd = [py, src / "rq2_strains/sensitivity_pfba/run_layered_executor_revised.py", "--config", rq2pcfg]
        if a.resume: cmd.append("--resume")
        run(cmd)
        run([py, src / "rq2_strains/sensitivity_pfba/rq2_cross_strain_analysis_revised.py", "--config", rq2pcfg,
             "--output", rq2p / "cross_strain_analysis", "--theta", "0.20",
             "--thresholds", "0.05,0.10,0.15,0.20,0.25,0.30,0.50", "--fdr", "0.05", "--effect_size", "0.5"])
    if "rq1_fva" in stages:
        run([py, src / "rq1_bulk/sensitivity_pfba/run_rq1_cohort_representative_fva.py", "--expression", rq1_expr,
             "--model", model, "--diet_bounds", diet, "--mapping", mapping, "--modeling_script", pfba,
             "--rq1_results", rq1p, "--output", out / "cohort_representative_FVA", "--fva_fraction", "1.0",
             "--parsimony_eps", "0.01,0.05", "--tolerance", "1e-7", "--processes", "1"])
    if "rq2_pfva" in stages:
        run([py, src / "rq2_strains/sensitivity_pfba/run_rq2_targeted_parsimonious_fva.py", "--config", rq2pcfg,
             "--modeling_script", pfba, "--rq2_results", rq2p, "--output", out / "targeted_parsimonious_FVA",
             "--parsimony_eps", "0.01,0.05", "--tolerance", "1e-7", "--processes", "1"])
    if "cross_rq_pfba" in stages:
        run([py, src / "cross_rq/run_cross_rq_final.py", "--rq1", rq1p, "--rq2", rq2p,
             "--rq3_results", R3 / f"q99_{obj}", "--rq4_results", R4, "--rq3_comparison", "WesternDiet_vs_Chow",
             "--anchor_abs_diff", "0.20", "--rq4_flux_threshold", "0.01", "--strict_expected_counts",
             "--output", out / "final_cross_rq_pfba"])

    # ------------------------------------------------------------------ optional E. coli benchmark
    if "benchmark_ecoli" in stages:
        e = cfg["ecoli_benchmark"]; E = lambda k: P(k, e)
        eo = out / "benchmark_ecoli"
        run([py, src / "benchmarking/benchmark_ecoli_13c_methods.py", "--model", E("model"), "--objective", e["objective"],
             "--transcriptomes", E("transcriptomes"), "--measured_modules", E("measured_modules"), "--uptakes", E("uptakes"),
             "--regime", "both", "--solver", solver, "--methods", "eflux,gimme,pfba,riptide", "--output", eo / "methods"])
        (eo / "reduced_engine").mkdir(parents=True, exist_ok=True)
        run([py, src / "validation/engine_recovery_13c/run_benchmark.py", "--prepared", E("prepared_dir"),
             "--model", E("model"), "--objective", e["objective"],
             "--validation", src / "validation/engine_recovery_13c/validation_01_benchmark_recovery.py",
             "--solver", solver, "--regime", "both", "--out", eo / "reduced_engine/benchmark_results.csv"],
            cwd=eo / "reduced_engine")
        run([py, src / "benchmarking/merge_5method_benchmark.py", "--methods_long", eo / "methods/benchmark_long_format_flux.csv",
             "--reduced_long", eo / "reduced_engine/benchmark_long_format_flux_reduced_engine.csv",
             "--output", eo / "five_method_final"])

    # ------------------------------------------------------------------ figures
    if "figures" in stages:
        F = out / "figures"; F.mkdir(parents=True, exist_ok=True)
        d1 = ref / "RQ1_primary_FBA" if a.use_reference else FLAT
        d2 = ref / "RQ2_primary_FBA_aggregate" if a.use_reference else FLAT
        run([py, src / "reporting/figures/fig2_hybrid_v2.py", "--data-dir", d1, "--out", F / "Figure2.png"])
        run([py, src / "reporting/figures/fig3_hybrid_heatmap.py", d2, F / "Figure3.png"])
        fi = ref / "figure_inputs"
        run([py, src / "reporting/figures/generate_Fig4_5_RQ3_panelv4.py",
             "--supp_xlsx", fi / "RQ3_Section23_Supplementary_Tables.xlsx",
             "--val_file", fi / "Section24_Validation_Supplementary_Tables.xlsx", "--out_dir", F, "--panel"])
        run([py, src / "reporting/figures/generate_Fig6_hybrid.py"],
            extra_env={"RQ4_RESULTS": str(ref / "RQ4_final_full_v3" if a.use_reference else R4),
                       "RQ4_ANNOTATIONS": str(P("reaction_annotations")), "FIG_OUT": str(F / "figure6")})
        run([py, src / "reporting/figures/generate_Fig7_relabel.py"], extra_env={"FIG_OUT": str(F / "figure7")})

    # ------------------------------------------------------------------ verify
    if "verify" in stages:
        cmd = [py, src / "utils/verify_manuscript_numbers.py", "--expected", ref / "EXPECTED_MANUSCRIPT_RESULTS.json"]
        if a.use_reference:
            cmd += ["--rq1_stats", ref / "RQ1_primary_FBA", "--rq1_prefix", "RQ1_reaction_stats_",
                    "--rq2_flat", ref / "RQ2_primary_FBA_aggregate", "--rq2_replicate_flat", ref / "RQ2_FBA_replicate",
                    "--benchmark", ref / "benchmarks/benchmark_5method_summary.csv"]
        else:
            cmd += ["--rq1_stats", R1 / "comprehensive_analysis/csv_outputs", "--rq1_prefix", "reaction_stats_",
                    "--rq2_root", R2A, "--rq2_replicate_root", R2R]
        run(cmd)
    print("\n[DONE] stages:", ", ".join(stages))


if __name__ == "__main__":
    main()

# Stage-to-manuscript map

| Stage | Script / entry point | Main manuscript role |
|---|---|---|
| RQ1 | `src/rq1_bulk/run_flux_pipeline_RQ1_revision.py` | Bulk multi-cohort diet effects |
| RQ1 FVA | `run_rq1_cohort_representative_fva.py` | Alternative-optima robustness |
| RQ2 | `src/rq2_strains/run_layered_executor_revised.py` + `rq2_cross_strain_analysis_revised.py` | Genetic background / conservation |
| RQ2 pFVA | `run_rq2_targeted_parsimonious_fva.py` | Universal-core robustness |
| Direct benchmark | `src/benchmarking/benchmark_liver_methods.py` | E-Flux/GIMME/pFBA comparison |
| Layer ablation | `validation_02_layer_ablation_v2.py` | Three-layer architecture |
| Parameter sensitivity | `sweep_threelayer_pipeline_v2.py` | P90/P95/P99, floor, objective |
| Determinism | `validation_03_determinism_v2.py` | Solver reproducibility |
| RQ3 | `src/rq3_single_cell/run_rq3_sensitivity_suite.py` | Single-cell attribution robustness |
| RQ4 | `src/rq4_microbiome/run_rq4_full_final_v3.py` | MICOM + host integration |
| Cross-RQ | `src/cross_rq/run_cross_rq_final.py` | Final pathway-level convergence |
| Reporting | `src/reporting/` | Figures 2, 3, 6, 7 |

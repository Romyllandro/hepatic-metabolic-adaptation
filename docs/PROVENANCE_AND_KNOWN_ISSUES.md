# Provenance and known issues

This release was assembled from the production working directory (`pn_pipeline`, Windows workstation,
Python 3.9, Gurobi). Every production script below was copied **byte-for-byte**, except the path changes
listed in `CHANGELOG.md`. Provenance was established from:

* production run manifests (`pipeline_manifests/*.yaml`) and logs;
* per-job `.runner_success.json` and `analysis_report.json` files;
* console transcripts;
* MD5 comparison of every released script against the production copies.

## 1. Which run produced each primary number

| Result | Production run (date) | Exact command (paths shortened) |
|---|---|---|
| RQ1 fluxes (99 samples, FBA) | `Final_run_3142026/Step_1_RQ1_{uncorrected,batch_corrected}` (11 Apr 2026) | `run_flux_pipeline_v3.py --results_dir Final_run_3142026/Step_1_RQ1 --no_fva --parallel` → `map_fixv5_multigroupsv8_layered_manuscript_run.py GSEMERGED_… --test_type mann-whitney --eflux_quantile 0.95 --eflux_floor 0.1 --eflux_cap 1000 --objective_id BIOMASS_mm_1_no_glygln --transporter_strategy either --edge_abs_diff_threshold 0.2 --rank_product_for HFD,KD,WD --no_fva` |
| RQ1 statistics, PCA, PERMANOVA | `Final_run_3142026/Step_1_RQ1_visualizations` (11 Apr 2026) | `comprehensive_flux_analysis.py --input …/reaction_flux_comparison_extended_batch_corrected.csv --output Final_run_3142026/Step_1_RQ1_visualizations` |
| RQ2 primary (aggregate design) | `972026_Step_2_RQ2_aggregated` (7 Sep 2026) | `map_fixv5_…_run.py <strain>.csv --test_type mann-whitney --aggregate …` (FVA on; FVA does not affect Diff) |
| RQ2 replicate-preserving design | `972026_Step_2_RQ2` (7 Sep 2026) | `run_layered_executor.py --config run_layered_config.json` (`test_type t-test`, `aggregate false`) |
| RQ1/RQ2 pFBA sensitivity | `RERUN_RQ1_PFBA`, `RERUN_RQ2_PFBA_FVA` (8 Sep 2026) | `run_flux_pipeline_RQ1_revision.py`, `run_layered_executor_revised.py` |
| RQ3 | `reviewer_results/RQ3_sensitivity` (9 Sep 2026) | `run_rq3_sensitivity_suite.py` |
| RQ4 | `reviewer_results/RQ4_final_full_v3` (9 Sep 2026) | `run_rq4_full_final_v3.py` |
| E. coli 5-method benchmark | `reviewer_results/benchmark_ecoli_13c_21reaction_riptide_zerofill` + `benchmark_long_format_flux_reduced_engine.csv` (23–24 Sep 2026) | reproduced exactly by `merge_5method_benchmark.py` |

The flat tables used for the manuscript (`RQs_Output_Files_final`) are byte-identical to the run
directories above. This was verified for RQ1 statistics and for all nine RQ2 edge files.

`src/utils/verify_manuscript_numbers.py` recomputes 34 headline values from these frozen tables, and all
34 match. Regenerating Figures 2, 6 and 7 from the frozen tables gives pixel-identical images.

## 2. Discrepancies between the code path and the manuscript text

These are reported so the text and the code can be reconciled before publication. The release
reproduces what was actually run.

1. **Input to the RQ1 reaction-level statistics.** The Welch/BH statistics were computed on the
   dataset-offset-removed flux table (`…_batch_corrected.csv`); see the console log of 11 Apr 2026.
   These statistics give 95 / 154 / 134 / 32 / 227, the 262 union and the 10-reaction signature.
   Methods 2.3 says the adjusted matrix was used *only* for PCA/visualization and that primary inference
   used unadjusted pooled samples.
   * To run the same statistics on the unadjusted fluxes, set
     `"rq1_inference_matrix": "uncorrected"` in `pipeline_config.json`. The numbers will change.
   * Either the text or the analysis should be updated.
2. **PERMANOVA.** The primary run's PERMANOVA (`comprehensive_flux_analysis.py`) permutes labels
   freely (no restriction within dataset) on z-scored PCA input from the offset-removed matrix.
   * A dataset-restricted PERMANOVA on raw feasible fluxes exists only for the pFBA sensitivity run
     (`comprehensive_flux_analysis_revised.py`: F = 3.83, R² = 0.108, p = 0.001).
   * The manuscript describes a restricted PERMANOVA on the unadjusted FBA matrix.
3. **"450 across all estimable pairwise comparisons".** 450 is the union over all six contrasts,
   *including* KD–WD (321 significant). The manuscript calls KD–WD non-estimable. The union over the five
   estimable contrasts is **360**. `rq1_signature_summary.py` reports both values.
4. **Engine-level test type.** The flux engine was called with `--test_type mann-whitney`. Its own
   `flux_pairwise_stats.csv` is not the source of the manuscript statistics, which come from Welch tests
   in `comprehensive_flux_analysis.py`. The argument is kept as run, for exact reproduction.
5. **Figure 7 / Supplementary Table S9 reference set.** The 15 HFD reference reactions and their
   membership in the genetic, single-cell and microbiome sets come from the April cross-RQ analysis
   (`src/cross_rq/legacy_15_reference/cross_rq_all_tables.py`). `generate_Fig7_relabel.py` embeds the
   approved membership table.
   * A check on 25 Sep 2026 against the final RQ3 (562-reaction) and RQ4 v3 (47-reaction) sets gave
     12/15 single-cell and 1/15 microbiome overlaps, not the 6/15 and 7/15 shown.
   * The later `run_cross_rq_final.py` (pFBA path; 20 anchors) is included as a sensitivity analysis.
   * This should be resolved before the cross-scale numbers are treated as final.
6. **Solver dependence.** Exact reaction sets depend on the vertex the LP solver returns. With GLPK the
   A/J strain reaches the same optimum as Gurobi but gives 64 instead of 96 responsive reactions. Use
   Gurobi to reproduce the published numbers; see `INSTALL.md`.

## 3. Items that are not code

* Figure 1 is a hand-made schematic. The submitted image contains typographical errors
  ("Evaluata GPR rule", mislabelled icons). Fix these before publication.
* The duplicate Figure 2 in the manuscript file is a document issue. The correct Figure 2 is
  `reference_outputs/manuscript_figures/Figure2.png`, generated by `fig2_hybrid_v2.py`.

## 4. Not included

* Superseded development scripts: the numbered `rq4_microbiome_community_modeling_v2…v8` series, early
  benchmark versions whose gene mapping matched 0 genes, the invalid RQ4 v2 handoff, and `pn_ai/`
  (a separate project).
* Word documents of the manuscript and the response to reviewers.
* `__pycache__` folders (excluded by `.gitignore`).

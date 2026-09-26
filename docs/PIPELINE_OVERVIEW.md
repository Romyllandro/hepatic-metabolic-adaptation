# Pipeline overview

## 1. Core modeling engine (all scales)

Model: iMM1415 (3,726 reactions, 2,775 metabolites, 1,375 genes).
Objective: `BIOMASS_mm_1_no_glygln`, used as a generic biosynthetic/metabolic-demand objective; ATP
maintenance (`ATPM`) is the sensitivity objective.

Reactions are split into three mutually exclusive classes, and each class has its own constraint layer:

| Layer | Reactions | Bound |
|---|---|---|
| 1 Environment | exchange reactions (`EX_`) | `l_j^diet ≤ v_j ≤ u_j^diet` from `expanded_diet_bounds_flat.json` |
| 2 Transport | reactions linking the extracellular compartment to intracellular compartments | `max(l_j⁰, −C_j) ≤ v_j ≤ min(u_j⁰, C_j)` |
| 3 Intracellular | all remaining reactions | same E-Flux clamp as layer 2 |

The reaction capacity `C_j` comes from normalised expression `x_g` through the GPR rule: the minimum
over genes for AND relationships and the maximum for OR relationships. Each sample is normalised on its
own: positive values are divided by the sample's 95th percentile (99th for single-cell pseudobulk), then
clipped to a floor of 0.1 and a cap of 1000 (single-cell: floor 0.001, cap 10,000).

The engine is implemented in `map_fixv5_multigroupsv8_layered_manuscript_run.py`, which does standard FBA.
Its sibling `..._rerun.py` adds `--solve_mode pfba` and FVA.

## 2. Stages

| Stage | Entry script(s) | Solver | Key parameters |
|---|---|---|---|
| `rq1` | `primary_fba/map_fixv5_..._run.py` → `batch_correct_flux_by_dataset_3d.py` → `comprehensive_flux_analysis.py` → `rq1_signature_summary.py` | FBA, one LP per sample (99) | P95 / 0.1 / 1000; Welch t-test + BH-FDR < 0.05; PCA; PERMANOVA (999 permutations) |
| `rq2` | `primary_fba/run_layered_executor.py` (×2 designs) → `rq2_conservation_tiers.py` → `rq2_layer3_strain_specificity.py` → `rq2_build_evidence_table.py` | FBA | aggregate design (primary) and replicate design; θ = 0.20 (sweep 0.10–0.30); per-strain BH-FDR; HC3 diet × strain interaction |
| `rq3` | `run_rq3_sensitivity_suite.py` → `rq3_FINAL_COMPLETE_REVISED.py` | pFBA per cell type × condition | P95/P99 × biomass/ATPM; threshold-responsive = 2 of 3 criteria (fold 1.5, relative 0.5, |Δv| 0.05/0.10/0.20) |
| `rq3_integration` | `rq1_rq2_rq3_integration_analysis_REVISED.py`, `rq2_rq3_multi_strain_integration_ENHANCED.py`, `bulk_validation_LOCAL.py`, `hierarchical_attribution_LOCAL.py` | none | abundance weighting; baseline vs Δflux concordance; functional/anatomical/lineage hierarchy |
| `rq4` | `run_rq4_full_final_v3.py` | MICOM cooperative tradeoff + host pFBA | 27-species strict map (28-species sensitivity); tradeoff 0.5 (0.3/0.7); portal scale 0.10 (0.05/0.20); forced / soft / availability-only coupling; biomass / ATPM |
| `validation` | `validation_04_positive_controls.py`, `benchmark_liver_methods.py`, `validation_02_layer_ablation_v2.py`, `sweep_threelayer_pipeline_v2.py`, `validation_03_determinism_v2.py` | FBA/pFBA | GSE101657 liver benchmark: layered / E-Flux / GIMME / diet-only pFBA |
| `closure_tests` | `run_all_remaining_meeting_tests.py` | FBA/pFBA | same-input FBA vs pFBA; caps 100/1,000/10,000; ATPM FVA; RQ4 community/host pFBA on/off |
| `benchmark_ecoli` | `benchmark_ecoli_13c_methods.py` + `run_benchmark.py` → `merge_5method_benchmark.py` | FBA/pFBA | 8 carbon sources × 21 reactions; AC/DC regimes; RIPTiDe zero-filled |
| `rq1_pfba`, `rq2_pfba`, `rq1_fva`, `rq2_pfva`, `cross_rq_pfba` | files in `sensitivity_pfba/`, `cross_rq/run_cross_rq_final.py` | pFBA + (parsimonious) FVA | pFBA fraction 1.0; parsimony ε 0.01/0.05 |
| `figures` | `src/reporting/figures/*` | none | Figures 2–7 |
| `verify` | `src/utils/verify_manuscript_numbers.py` | none | compares with `reference_outputs/EXPECTED_MANUSCRIPT_RESULTS.json` |

## 3. Data flow

```text
inputs/bulk/GSEMERGED_*.csv ──rq1──> results/RQ1_primary_FBA/ ─┐
inputs/strains/*.csv ─────────rq2──> results/RQ2_primary_FBA_aggregate/, RQ2_FBA_replicate/ ─┤
                                                               ├─ flatten ─> results/RQs_Output_Files/ ─┬─> figures
inputs/single_cell/*  ────────rq3──> results/RQ3_sensitivity/q99_BIOMASS_mm_1_no_glygln/ ──┤           └─> verify
                                         └─ rq3_integration ─> results/RQ3_integration/, RQ3_bulk_validation/
inputs/microbiome + AGORA2 ───rq4──> results/RQ4_final_full_v3/
```

`results/RQs_Output_Files/` is a flat copy with the file names used in the manuscript tables
(`RQ1_reaction_stats_HFD_vs_SCD.csv`, `RQ2_AJ_edges_HFD_vs_SCD.csv`, and so on). In the production
directory it was called `RQs_Output_Files_final`.

## 4. Run time (production workstation, Gurobi)

| Stage | Approximate time |
|---|---|
| rq1 (99 LPs, no FVA) | 3–5 min |
| rq2 aggregate + replicate (9 strains each, FVA on) | 1–2 h |
| rq3 sensitivity suite (4 settings) | 1–3 h |
| rq4 full + sensitivity | several hours |
| rq1_pfba with FVA | about 1.5 h |
| verify / figures with `--use-reference` | < 1 min |

## 5. Interpretation rules built into the design

* Magnitude tiers (RQ2) and threshold-responsive classifications (RQ3) are deterministic model
  classifications, not significance tests. Replicate-level statistics are reported separately.
* NZO/HlLtJ has one HFD sample. It stays in the magnitude analysis but is non-estimable for statistics.
* The single-cell atlas compares WD with chow. Overlap with bulk HFD responses is contextual convergence,
  not condition-matched validation.
* The cross-strain single-cell panel reweights one atlas across strains. It is not nine independent
  experiments.
* Only positive community exports are coupled to the liver in RQ4. Negative MICOM exchange means
  microbial uptake, not delivery to the host.
* FBA and pFBA reach the same optimum but can select different flux vectors. Reaction-level membership
  is therefore reported together with pFBA and FVA sensitivity.

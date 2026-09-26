# Running the pipeline on your own data

The engine is not specific to liver or to these diets. You need four things:

* a genome-scale model in COBRApy JSON or SBML;
* an expression matrix;
* a diet (or medium) bounds file;
* an Entrez–symbol mapping, if your model's GPRs use a different gene identifier than your matrix.

## 1. Expression matrix (bulk)

A CSV with genes in rows and samples in columns:

```text
Gene_Symbol,Gene_ID,SCD_GSE101657_Control_SCD_0,...,HFD_GSE101657_HFD_0,...
Tbc1d30,74694,668,...
```

* The first two columns are the gene symbol and the Entrez ID.
* Each sample column name starts with a **group code**, followed by `_`. The engine infers groups from
  this prefix (`--infer_groups_from_filename`), or you can pass `--explicit_groups`.
* For multi-cohort runs, include a dataset tag matching `GSE\d+` or `GSM\d+` in every column name. The
  dataset-offset removal and the restricted permutation tests read this tag
  (`--dataset_regex` changes the pattern).
* Values should be non-negative, on a linear scale (counts, TPM or FPKM). Each sample is normalised
  internally to its own P95.
* The baseline group defaults to SCD, then WD, then the first group found. Override it with
  `--baseline_code`.

## 2. Diet / medium bounds

A JSON file with one entry per group code, mapping exchange reactions to `[lower, upper]` bounds in model
units (mmol gDW⁻¹ h⁻¹). Negative lower bounds allow uptake.

```json
{
  "SCD": {"EX_glc__D_e": [-10, 1000], "EX_ala__L_e": [-2, 1000]},
  "HFD": {"EX_glc__D_e": [-4, 1000], "EX_hdca_e": [-3, 1000]}
}
```

Exchange reactions that are not listed keep their model defaults. If your bounds are in g/day or
mmol/day, set `--diet_bounds_units`, and pass `--mw_table` and `--gDW` for the conversion.

## 3. Run the engine on one dataset

```bash
python src/rq1_bulk/primary_fba/map_fixv5_multigroupsv8_layered_manuscript_run.py my_expression.csv \
  --model_file my_model.json --diet_bounds_json my_diets.json --mapping_file mouse_entrez_to_symbol.csv \
  --objective_id BIOMASS_mm_1_no_glygln --eflux_quantile 0.95 --eflux_floor 0.1 --eflux_cap 1000 \
  --transporter_strategy either --edge_abs_diff_threshold 0.2 --write_replicates_long \
  --results_dir results/my_run --no_fva --solver gurobi
```

Main outputs:

| File | Contents |
|---|---|
| `flux_analysis/reaction_flux_comparison_extended.csv` | per-sample fluxes, group means, `Diff(A-B)` |
| `stats_comparison/flux_pairwise_stats.csv` | engine-level pairwise tests (`--test_type`) |
| `cytoscape_edges/edges_<A>_vs_<B>.csv` | reaction–metabolite edges for reactions with abs(Diff) ≥ threshold |
| `analysis_report.json` | parameters, sample counts, constrained-reaction counts per layer, objective values |

Useful options:

* `--aggregate`: average replicates within each group before solving (the RQ2 primary design).
* `--no_fva`: skip FVA.
* The `_rerun.py` engine in `sensitivity_pfba/` adds `--solve_mode pfba --pfba_fraction 1.0`.

## 4. Downstream statistics and summaries

```bash
python src/rq1_bulk/primary_fba/batch_correct_flux_by_dataset_3d.py --input <flux.csv> --output_csv <flux_bc.csv>
python src/rq1_bulk/primary_fba/comprehensive_flux_analysis.py --input <flux or flux_bc> --output results/my_stats
python src/rq1_bulk/primary_fba/rq1_signature_summary.py --stats_dir results/my_stats/csv_outputs --output results/my_summary
```

`rq1_signature_summary.py` assumes the contrast names used in this study. Edit `SCD_REF` and `ESTIMABLE`
at the top of the script for your own groups.

For several genetic backgrounds or conditions, list one job per expression file in an executor config
(see `results/resolved_configs/rq2_fba_aggregate.json` after a run) and run:

```bash
python src/rq2_strains/primary_fba/run_layered_executor.py --config my_jobs.json --resume
python src/rq2_strains/primary_fba/rq2_conservation_tiers.py --rq2_root <results root> --output <dir> --theta 0.2
```

To use strain names other than the nine founders, edit `STRAINS` in `rq2_conservation_tiers.py` and
`flatten_outputs.py`.

## 5. Single-cell data

Supply a genes × cells CSV and a metadata CSV (`cell_id, cell_type, diet`), then run
`src/rq3_single_cell/run_rq3_sensitivity_suite.py`. Set `--condition_mapping` to map your condition
labels to diet keys in the bounds JSON, for example `{"Chow":"SCD","WesternDiet":"WD"}`.
`rq3_data_preparation.py` converts h5ad or Seurat-exported data to this format.

## 6. Microbiome coupling

`src/rq4_microbiome/run_rq4_full_final_v3.py` needs:

* a species-activity table (columns `Genome_name`, one column per condition);
* a reviewed species-to-AGORA2 map (see `configs/curated_species_to_agora_REVIEWED_v3.csv`; use exact
  mappings, not fuzzy matching);
* a local AGORA2 `.mat` directory.

The runner executes the primary configuration and all sensitivity scenarios. To run a single host scenario with your own portal settings, call `rq4_hepatic_integration_final_v2.py` directly (`--portal_scaling`, `--portal_mode forced|soft|availability`, `--objective_mode`, `--no-pfba`).

## 7. Good practice

* Always run a sensitivity arm: pFBA vs FBA, at least two objectives, a threshold sweep, and FVA on the
  reactions you report. Exact reaction sets depend on the solution chosen by the solver.
* Record solver name, version and settings (`src/utils/capture_environment.py`).

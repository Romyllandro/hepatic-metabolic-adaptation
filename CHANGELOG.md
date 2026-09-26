# Changelog

## v1.0.0 — manuscript release (2026-09-26)

First public release, matching the revised manuscript (BIB-26-1735). The primary analyses use standard
FBA; pFBA and FVA are sensitivity analyses.

### Added (new code; production scripts not modified)

* `run_pipeline.py`: a single orchestrator with primary, sensitivity and optional stages, `--dry-run`,
  `--resume` and `--use-reference`.
* `src/rq1_bulk/primary_fba/rq1_signature_summary.py`: per-contrast counts, unions and the cross-diet
  signature from the `reaction_stats_*` tables.
* `src/rq2_strains/primary_fba/rq2_conservation_tiers.py`: union, universal core, tiers, per-strain counts
  and θ sweep. It uses the same exact `≥ θ` comparison as the engine's edge filter.
* `src/benchmarking/merge_5method_benchmark.py`: rebuilds Supplementary Table SL3. It reproduces the
  submitted table exactly.
* `src/utils/flatten_outputs.py`: builds the flat `RQ1_*` / `RQ2_*` tables used by the figure scripts.
* `src/utils/verify_manuscript_numbers.py` and `reference_outputs/EXPECTED_MANUSCRIPT_RESULTS.json`:
  34 automated checks.
* Documentation under `docs/`, plus `CITATION.cff`, `.zenodo.json`, licences and `tools/make_zenodo_bundle.py`.

### Changed (minimal edits to production scripts)

* `src/rq3_single_cell/integration/bulk_validation_LOCAL.py` and `hierarchical_attribution_LOCAL.py`:
  hard-coded Windows input and output paths replaced by environment variables with relative defaults. No
  analysis code changed.
* `src/reporting/figures/generate_Fig6_hybrid.py` and `generate_Fig7_relabel.py`: data and output
  directories read from environment variables (`RQ4_RESULTS`, `RQ4_ANNOTATIONS`, `FIG_OUT`); the output
  directory is created with `parents=True`.
* `src/validation/meeting_closure_tests/run_all_remaining_meeting_tests.py` and related scripts: paths to
  the pFBA engine updated to its new location in `src/rq1_bulk/sensitivity_pfba/`.
* `src/benchmarking/benchmark_ecoli_13c_methods.py` and `benchmark_methods_core.py`: replaced with the
  final production versions (RIPTiDe added, zero-filled scoring).

### Reorganised

* The pFBA/FVA engine and its runners moved from `src/rq1_bulk/` and `src/rq2_strains/` into
  `sensitivity_pfba/` subfolders. The standard-FBA production scripts were added under `primary_fba/`.
* The earlier pFBA-primary package entry points (`run_full_pipeline.py`, its config and README) and
  their expected counts are kept in `legacy/` for provenance.

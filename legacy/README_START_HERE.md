# Multi-scale constraint-based hepatic metabolism pipeline
## Frozen manuscript reproduction package

This directory packages the **final corrected code path** used for the revised
manuscript. It deliberately excludes superseded/broken intermediate scripts
(e.g., early benchmark gene-mapping versions and the invalid RQ4 v2 community
handoff).

The pipeline covers:

1. **RQ1 — multi-cohort bulk liver**: three-layer E-Flux + replicate-level pFBA,
   dataset-aware inference, PCA/QC, and FVA refinement.
2. **RQ2 — genetic background**: nine strain-specific replicate-preserving runs,
   conservation analysis, threshold sensitivity, and targeted parsimonious FVA.
3. **Reviewer validation**: direct three-layer/E-Flux/GIMME/diet-only pFBA
   benchmark, layer ablation, parameter/objective sensitivity, determinism,
   and robust post-hoc statistics.
4. **RQ3 — single-cell**: cell-type pseudobulk E-Flux/pFBA plus P95/P99,
   biomass/ATPM, and threshold sensitivity.
5. **RQ4 — microbiome**: strict curated 27-species primary map, 28-species
   mapping sensitivity, MICOM tradeoff sensitivity, host objective/portal
   scaling/coupling sensitivity, and hepatic attribution.
6. **Final cross-RQ integration**: reaction tracing and pathway-level
   convergence across RQ1-RQ4.
7. **Figures**: reproducible Figure 2, Figure 3, Figure 6, and Figure 7 scripts.

---

# 1. Install

Recommended:

```bat
conda env create -f environment.yml
conda activate metabolic-modeling
```

or install into the existing environment:

```bat
pip install -r requirements.txt
```

A valid Gurobi installation/license is required for the upstream metabolic
optimization stages.

For an archival release, capture the exact working environment:

```bat
python src\utils\capture_environment.py --output environment_capture
```

---

# 2. Add input data

Read:

```text
inputs/README_INPUTS.md
inputs/INPUT_MANIFEST.csv
```

Copy the example configuration:

```bat
copy pipeline_config.example.json pipeline_config.json
```

Then edit at minimum:
- RQ3 expression/metadata paths;
- `agora_dir`;
- any input paths that differ from the example layout.

The reviewed RQ4 species-map snapshot is already included under `configs/`.

---

# 3. Dry-run / preflight

```bat
python run_full_pipeline.py --config pipeline_config.json --stages preflight
```

To inspect the full command chain without running optimizations:

```bat
python run_full_pipeline.py --config pipeline_config.json --stages all --dry-run
```

---

# 4. Full manuscript analysis

```bat
python run_full_pipeline.py --config pipeline_config.json --stages all
```

For long RQ2/RQ4 jobs, interrupted stages can be resumed:

```bat
python run_full_pipeline.py --config pipeline_config.json --stages rq2,rq4 --resume
```

You can also run stages separately, for example:

```bat
python run_full_pipeline.py --config pipeline_config.json --stages rq1,rq2
python run_full_pipeline.py --config pipeline_config.json --stages rq1_fva,rq2_pfva
python run_full_pipeline.py --config pipeline_config.json --stages benchmark,ablation,parameter_sensitivity,determinism
python run_full_pipeline.py --config pipeline_config.json --stages rq3
python run_full_pipeline.py --config pipeline_config.json --stages rq4
python run_full_pipeline.py --config pipeline_config.json --stages cross_rq,figures,verify
```

The optional external 13C-MFA E. coli benchmark is not part of `all`; run it only
after supplying its external inputs:

```bat
python run_full_pipeline.py --config pipeline_config.json --stages benchmark_ecoli
```

---

# 5. Frozen analysis settings

## RQ1 / RQ2 bulk E-Flux
- percentile: P95
- floor: 0.1
- cap: 1000
- objective: `BIOMASS_mm_1_no_glygln`
- solve: pFBA, 100% of biological optimum
- deterministic Gurobi: Threads=1, Method=1, NumericFocus=3, Seed=1

## RQ3 reference
- percentile: P99
- floor: 0.001
- cap: 10000
- objective: `BIOMASS_mm_1_no_glygln`
- sensitivity: P95/P99 × biomass/ATPM; |delta v|=0.05/0.10/0.20
- RQ3 is WesternDiet-vs-Chow and is treated as orthogonal cellular-context evidence.

## RQ4 reference
- strict primary species map: 27
- mapping sensitivity: 28
- MICOM tradeoff: 0.50 (sensitivity 0.30/0.70)
- MICOM pFBA
- host objective: biomass-associated demand (sensitivity ATPM)
- host pFBA fraction: 1.0
- portal scale: 0.10 (sensitivity 0.05/0.20)
- portal mode: forced (sensitivity soft + availability-only)

---

# 6. Expected final QC

`reference_outputs/EXPECTED_FROZEN_RESULTS.json` stores the revision-level expected
counts.

Final cross-RQ analysis should reproduce:

```text
RQ1 HFD anchors: 20
RQ1 all-three SCD-referenced shared reactions: 27
RQ2 universal magnitude core: 16
RQ3 P99/biomass responsive reactions: 562
RQ4 HFD microbiome-attributable reactions (|delta v| >= 0.01): 47

RQ1 anchor overlap with RQ2 universal: 0
RQ1 anchor overlap with RQ3 WD/Chow: 18
RQ1 anchor overlap with RQ4 thresholded microbiome attribution: 6
RQ1 anchor overlap with RQ4 microbiome-dominant set: 1
pathways represented across all four RQs: 6
```

The final cross-RQ script is run with `--strict_expected_counts`, so it stops if
the headline frozen counts drift.

---

# 7. Important reproducibility/interpretation rules

- **Do not infer statistics from batch-adjusted fluxes.** Batch-adjusted flux
  matrices are for PCA/QC/visualization. Primary RQ1 inference uses original
  feasible per-sample pFBA fluxes with dataset-aware modeling.
- **Do not average RQ2 biological replicates before inference.** Preserve
  per-animal samples; NZO HFD n=1 remains non-estimable statistically.
- **Do not call RQ3 replicate-level differential inference.** It is deterministic
  pseudobulk threshold attribution, with explicit sensitivity analysis.
- **Do not run the superseded RQ4 v2 path.** The valid path is
  `src/rq4_microbiome/run_rq4_full_final_v3.py`.
- **Do not use fuzzy species matching for final RQ4.** The final pipeline uses a
  strict curated whitelist.
- **Do not treat negative MICOM community exchange fluxes as microbial delivery.**
  The v3 pipeline propagates direction-aware positive export; in the primary
  result acetate is the only curated positive portal export used for host coupling.
- **Do not claim a universal reaction-level metabolic backbone.** Final evidence
  supports pathway/transport-system convergence across scales.
- **Do not claim three-layer E-Flux mathematically outperforms conventional
  E-Flux.** Under matched settings the direct benchmark is identical; novelty is
  modular separation/auditability and multi-scale tracing.

---

# 8. Sharing / archival release checklist

Before uploading to Zenodo or a public repository:

1. Add permitted derived input matrices under `inputs/`.
2. Add external-resource links/DOIs for iMM1415, AGORA2, and raw GEO data.
3. Run `src/utils/capture_environment.py`.
4. Run the complete pipeline or at least the final `verify` stage.
5. Include `MANIFEST_SHA256.json`.
6. Add an explicit software/data license chosen by the authors.
7. Add repository DOI/version to the manuscript Data Availability statement.

`LICENSE_NOT_INCLUDED.txt` is intentional: this package does not assume which
license the authors wish to use.

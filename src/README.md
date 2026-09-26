# src/

| Folder | Role in the manuscript | Status |
|---|---|---|
| `rq1_bulk/primary_fba/` | Bulk multi-cohort FBA, statistics, PCA/PERMANOVA, signature (Section 3.2, Fig. 2) | primary |
| `rq1_bulk/sensitivity_pfba/` | pFBA engine (`--solve_mode pfba`), FVA, dataset-adjusted HC3 model, restricted PERMANOVA | sensitivity |
| `rq2_strains/primary_fba/` | Nine-strain FBA executor, conservation tiers, θ sweep, within/between-strain statistics (Section 3.3, Fig. 3) | primary |
| `rq2_strains/sensitivity_pfba/` | pFBA strain runs, cross-strain pFBA analysis, targeted parsimonious FVA | sensitivity |
| `rq3_single_cell/` | Pseudobulk pFBA + P95/P99 × biomass/ATPM × threshold suite (Section 3.4) | primary |
| `rq3_single_cell/integration/` | Abundance-weighted attribution, bulk comparison (Section 3.5), hierarchy, cross-strain projection (Figs. 4–5) | primary |
| `rq4_microbiome/` | MICOM community, species map, hepatic coupling, attribution (Section 3.6, Fig. 6) | primary |
| `cross_rq/` | `run_cross_rq_final.py` (pFBA path) and `legacy_15_reference/` (source of the 15 reference reactions in Fig. 7) | see docs/PROVENANCE_AND_KNOWN_ISSUES.md |
| `benchmarking/` | Liver and E. coli method benchmarks, layer ablation, parameter sweep, determinism (Section 3.1) | validation |
| `validation/` | Positive controls, 13C engine recovery, reviewer closure tests (FBA vs pFBA, caps, ATPM FVA, RQ4 pFBA on/off) | validation |
| `reporting/` | Figure scripts (Figs. 2–7) and the supplementary-table builder | reporting |
| `common/` | Deterministic Gurobi settings | shared |
| `utils/` | Verification, output flattening, environment capture | shared |

All scripts accept `--help`. `run_pipeline.py` in the repository root shows the exact arguments used for
each stage (`--dry-run`).

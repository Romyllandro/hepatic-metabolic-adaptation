# Data sources

## Included in the repository (`inputs/`)

| File | Content | Source |
|---|---|---|
| `models/iMM1415.json` | mouse genome-scale model | Sigurdsson et al., *BMC Syst Biol* 2010, 4:140; BiGG Models (http://bigg.ucsd.edu/models/iMM1415) |
| `models/iJO1366.json` | *E. coli* model (benchmark only) | Orth et al., *Mol Syst Biol* 2011, 7:535; BiGG Models |
| `bulk/GSEMERGED_SCD_HFD_KD_WD_gene_expression.csv` | merged liver expression, 99 samples (SCD 37, HFD 38, KD 12, WD 12) | derived from GEO GSE101657, GSE159090, GSE160646, GSE188344, GSE246221, GSE248297 |
| `bulk/RQ1_sample_metadata.csv` | sample → group → dataset | this study |
| `bulk/expanded_diet_bounds_flat.json` | exchange bounds for SCD, LFD, HFD, KD, WD | this study (see Supplementary Methods) |
| `bulk/mouse_entrez_to_symbol.csv` | Entrez ↔ symbol mapping | MyGene.info query (taxid 10090) |
| `bulk/RQ4_reaction_annotations.csv` | reaction names and subsystems for iMM1415 | extracted from the model |
| `strains/male-<STRAIN>-GSE182668_HFD_SCD_gene_expression.csv` | per-strain liver expression, 9 founder strains | derived from GEO GSE182668 (Bachmann et al., *iScience* 2022) |
| `microbiome/Meta_GSE104913.csv` | species-level metatranscriptomic activity, ND/SCD and DD/HFD | derived from GEO GSE104913 (Chen et al., *Nat Biotechnol* 2020) |
| `single_cell/cell_metadata_fixed.csv` | cell → cell type → diet (Chow / WesternDiet) | derived from GEO GSE218300 (Bendixen et al., *J Hepatol* 2024) |
| `ecoli_benchmark/*` | transcriptomes, uptake rates, 13C-MFA module fluxes, per-condition inputs | Gerosa et al., *Cell Syst* 2015 (Table S2, `mmc2.xlsx`); Bhadra-Lobo et al., *PLOS ONE* 2020 (S1 archive) |

## Not in the git repository

| Resource | Size | How to obtain |
|---|---|---|
| Single-cell expression matrix `expression_matrix.csv` (genes × cells) | 5.7 GB | Zenodo archive (`single_cell_expression_matrix.csv.gz`), or rebuild from GEO GSE218300 with `src/rq3_single_cell/integration/rq3_data_preparation.py`. Place it at `inputs/single_cell/expression_matrix.csv`. |
| AGORA2 reconstructions (`.mat`) | several GB | https://www.vmh.life/#downloadview (Heinken et al., *Nat Biotechnol* 2023). Set `inputs.agora_dir` in `pipeline_config.json`. |
| Gurobi | – | https://www.gurobi.com (free academic licence) |
| Full intermediate results of the production run | ~1.5 GB | Zenodo archive (`production_results.zip`) |

## GEO accessions used

| Accession | Use | Reference |
|---|---|---|
| GSE101657 | bulk: SCD, HFD, KD | Newman et al., *Cell Metab* 2017 |
| GSE159090 | bulk: SCD, HFD, WD | Smati et al., *Gut* 2022 |
| GSE160646 | bulk: SCD, HFD | Geißler et al., *J Nutr Biochem* 2022 |
| GSE188344 | bulk: SCD, HFD | Matsushita et al., *Cell Chem Biol* 2022 |
| GSE246221 | bulk: SCD, HFD | Jeong et al., *Nat Commun* 2024 |
| GSE248297 | bulk: SCD, HFD, KD | Gallop et al., *Sci Adv* 2025 |
| GSE182668 | nine strains, HFD vs SCD | Bachmann et al., *iScience* 2022 |
| GSE218300 | single-cell liver atlas, WD vs chow | Bendixen et al., *J Hepatol* 2024 |
| GSE104913 | gut metatranscriptome | Chen et al., *Nat Biotechnol* 2020 |

## Licences

Derived expression matrices are redistributed under the terms of GEO public data. Models from BiGG and
the published benchmark tables keep their original licences and should be cited as above. Data generated
by this study (diet bounds, curated species map, all results) are released under CC BY 4.0.

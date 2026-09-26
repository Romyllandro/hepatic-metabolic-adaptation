# Input data

The code package intentionally does **not** redistribute large raw/derived biological
datasets, the iMM1415 model, or AGORA2 reconstructions. Put those files under
`inputs/` (or change their paths in `pipeline_config.json`).

## Required manuscript inputs

RQ1 uses six GEO liver cohorts:
- GSE101657
- GSE159090
- GSE160646
- GSE188344
- GSE246221
- GSE248297

The frozen merged RQ1 matrix has 99 samples:
- SCD 37
- HFD 38
- KD 12
- WD 12

RQ2 uses the nine GSE182668 strain-specific HFD/SCD matrices listed in the
example config.

RQ3 requires the exact single-cell expression matrix and metadata used in the
manuscript. Preserve the `Chow -> SCD` and `WesternDiet -> WD` condition mapping.

RQ4 requires:
- `Meta_GSE104913.csv`
- AGORA2 `.mat` directory
- reviewed species map included under `configs/`
- C57BL/6J GSE182668 hepatic expression for host integration

## Data-release recommendation

For a public Zenodo/GitHub release, distribute:
1. all scripts/configs in this package;
2. derived expression matrices that you are permitted to redistribute;
3. the reviewed species map;
4. frozen result tables/figures;
5. a DOI/link for external model/data resources that cannot be redistributed.

Before archiving, run:
`python src/utils/capture_environment.py --output environment_capture`
and include that folder in the release.

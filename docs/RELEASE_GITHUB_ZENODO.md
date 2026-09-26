# Depositing this release on GitHub and Zenodo

Two deposits, linked:

* **GitHub** holds everything in this folder: code, docs, small inputs and frozen reference tables
  (about 60 MB).
* **Zenodo** archives the GitHub release automatically and gives it a DOI. A second Zenodo record (or
  extra files in the same record) holds the large data that cannot go on GitHub.

## Step 1 — check the release locally (Windows, from `pn_pipeline`)

```bat
cd hepatic-metabolic-adaptation_release_v1.0.0
copy pipeline_config.example.json pipeline_config.json
python run_pipeline.py --config pipeline_config.json --stages verify --use-reference
python run_pipeline.py --config pipeline_config.json --stages figures --use-reference
python src\utils\capture_environment.py --output environment_capture
```

The verify command must end with `34/34 checks passed`. Commit the `environment_capture` folder: it
records the exact package and Gurobi versions from your `metabolic-modeling` environment.

## Step 2 — confirm the metadata

* `LICENSE` (MIT) and `LICENSE-DATA` (CC BY 4.0) are proposed defaults. Change them if the lab prefers
  different licences.
* In `CITATION.cff` and `.zenodo.json`, check the author list, ORCIDs, affiliations and version.

## Step 3 — push to GitHub

```bat
cd hepatic-metabolic-adaptation_release_v1.0.0
git init
git add .
git commit -m "Release v1.0.0: code and frozen results for BIB-26-1735"
git branch -M main
git remote add origin https://github.com/sbbi-unl/hepatic-metabolic-adaptation.git
git push -u origin main
```

If the repository already has history, push to a new branch and open a pull request instead of
overwriting it.

## Step 4 — connect Zenodo and tag the release

1. Log in at https://zenodo.org with GitHub, open **GitHub** in the account menu, and switch the
   repository on.
2. On GitHub, create a release with tag `v1.0.0`, titled
   "v1.0.0 – manuscript release (Briefings in Bioinformatics)".
3. Zenodo archives the release and mints a DOI. The metadata comes from `.zenodo.json`.

## Step 5 — upload the large data to Zenodo

From the `pn_pipeline` folder:

```bat
python hepatic-metabolic-adaptation_release_v1.0.0\tools\make_zenodo_bundle.py --pn_pipeline . --out zenodo_upload
```

This writes to `zenodo_upload/`:

| File | Contents |
|---|---|
| `single_cell_expression_matrix.csv.gz` | the 5.7 GB single-cell matrix, compressed |
| `production_results.zip` | full production outputs: RQ1–RQ4 runs, sensitivity runs, benchmarks, closure tests, logs and manifests |
| `SHA256SUMS.txt` | checksums for both files |

Create a new Zenodo upload ("Dataset") with these files. In the upload, link it to the software DOI from
Step 4 ("is supplemented by").

## Step 6 — update the manuscript

In the Data Availability statement, replace "will be deposited in Zenodo with a DOI upon acceptance" with
the two DOIs (software and data) and the GitHub URL.

#!/usr/bin/env python3
"""
make_zenodo_bundle.py
=====================
Build the large-data companion files for the Zenodo data record from the production
working directory (pn_pipeline). Run on the workstation that holds the production outputs.

Creates, in --out:
  single_cell_expression_matrix.csv.gz   gzip of data/scRNA_52w_data/expression_matrix.csv (5.7 GB)
  production_results.zip                 full production outputs (see RESULT_DIRS below)
  SHA256SUMS.txt                         checksums

Usage (Windows, from pn_pipeline):
  python hepatic-metabolic-adaptation_release_v1.0.0\\tools\\make_zenodo_bundle.py --pn_pipeline . --out zenodo_upload
  add --skip_single_cell or --skip_results to build only one part
"""
import argparse, gzip, hashlib, shutil, zipfile
from pathlib import Path

# (source path inside pn_pipeline, archive path inside production_results.zip)
RESULT_DIRS = [
    ("RQs_Output_Files_final", "manuscript_tables/RQs_Output_Files_final"),
    ("Final_run_3142026/Step_1_RQ1_uncorrected", "RQ1_primary_FBA/Step_1_RQ1_uncorrected"),
    ("Final_run_3142026/Step_1_RQ1_batch_corrected", "RQ1_primary_FBA/Step_1_RQ1_batch_corrected"),
    ("Final_run_3142026/Step_1_RQ1_visualizations_manuscript_used", "RQ1_primary_FBA/Step_1_RQ1_visualizations"),
    ("Final_run_3142026/Step_3b_RQ3_Results_manus", "RQ3_primary_run/Step_3b_RQ3_Results"),
    ("972026_Step_2_RQ2_aggregated", "RQ2_primary_FBA_aggregate"),
    ("972026_Step_2_RQ2", "RQ2_FBA_replicate"),
    ("RERUN_RQ1_PFBA", "sensitivity/RQ1_pFBA"),
    ("RERUN_RQ2_PFBA_FVA", "sensitivity/RQ2_pFBA_FVA"),
    ("reviewer_results/rq2_layer3_strain_specificity", "RQ2_statistics/layer3_strain_specificity"),
    ("reviewer_results/RQ3_sensitivity", "RQ3_sensitivity"),
    ("reviewer_results/RQ4_final_full_v3", "RQ4_final_full_v3"),
    ("reviewer_results/final_cross_rq", "sensitivity/final_cross_rq_pfba"),
    ("reviewer_results/benchmark_liver_GSE101657_v4", "validation/benchmark_liver_GSE101657"),
    ("reviewer_results/benchmark_ecoli_13c_21reaction_riptide_zerofill", "validation/benchmark_ecoli_methods"),
    ("reviewer_results/benchmark_ecoli_13c_5method_final", "validation/benchmark_ecoli_5method_final"),
    ("reviewer_results/threelayer_sensitivity_GSE101657", "validation/threelayer_sensitivity_GSE101657"),
    ("reviewer_results/meeting_closure", "validation/closure_tests"),
    ("revision_runs", "validation/revision_runs"),
    ("pipeline_logs", "logs/pipeline_logs"),
    ("pipeline_manifests", "logs/pipeline_manifests"),
]
RESULT_FILES = [
    ("reviewer_results/layer_ablation_v2.csv", "validation/layer_ablation_v2.csv"),
    ("benchmark_long_format_flux_reduced_engine.csv", "validation/benchmark_long_format_flux_reduced_engine.csv"),
    ("reviewer_results/Supplementary_Tables_ready_2026-09-26.xlsx", "Supplementary_Tables.xlsx"),
    ("Final_run_3142026/Console_4112026.txt", "logs/console_RQ1_RQ2_2026-04-11.txt"),
]
SKIP_PARTS = {"__pycache__"}


def sha256(p, buf=1 << 24):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(buf), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pn_pipeline", default=".")
    ap.add_argument("--out", default="zenodo_upload")
    ap.add_argument("--skip_single_cell", action="store_true")
    ap.add_argument("--skip_results", action="store_true")
    a = ap.parse_args()
    root = Path(a.pn_pipeline).resolve(); out = Path(a.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    made = []

    if not a.skip_single_cell:
        src = root / "data/scRNA_52w_data/expression_matrix.csv"
        dst = out / "single_cell_expression_matrix.csv.gz"
        print(f"[gzip] {src} -> {dst} (this takes a while)")
        with open(src, "rb") as fi, gzip.open(dst, "wb", compresslevel=6) as fo:
            shutil.copyfileobj(fi, fo, length=1 << 24)
        made.append(dst)

    if not a.skip_results:
        dst = out / "production_results.zip"
        print(f"[zip] {dst}")
        missing = []
        with zipfile.ZipFile(dst, "w", compression=zipfile.ZIP_DEFLATED, allowZip64=True) as z:
            for s, arc in RESULT_DIRS:
                d = root / s
                if not d.exists():
                    missing.append(s); continue
                for f in d.rglob("*"):
                    if f.is_file() and not (SKIP_PARTS & set(f.parts)):
                        z.write(f, f"{arc}/{f.relative_to(d).as_posix()}")
            for s, arc in RESULT_FILES:
                f = root / s
                if f.exists():
                    z.write(f, arc)
                else:
                    missing.append(s)
        for m in missing:
            print(f"[warn] not found, skipped: {m}")
        made.append(dst)

    with open(out / "SHA256SUMS.txt", "w") as fh:
        for f in made:
            fh.write(f"{sha256(f)}  {f.name}\n")
    print("[done]", *[f"{f.name} ({f.stat().st_size / 1e9:.2f} GB)" for f in made], sep="\n  ")


if __name__ == "__main__":
    main()

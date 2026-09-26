#!/usr/bin/env python3
"""
build_benchmark_inputs.py
=========================
Build the tidy input CSVs for Prong-1 benchmark recovery (validation_01) from the
Bhadra-Lobo, Kim & Lun (2020) S1 Dataset (E. coli arm).

Reference: Bhadra-Lobo S, Kim MK, Lun DS. PLoS ONE 15(9):e0238689 (2020).
13C measured fluxes & microarray expression originate from Gerosa et al. 2015,
Cell Systems 1(4):270-282, DOI 10.1016/j.cels.2015.09.008.

WHAT IT BUILDS
--------------
1. transcriptomes.csv        organism, condition, gene_id, expression          (linear scale)
2. uptakes.csv               organism, condition, exchange_id, lower_bound, upper_bound, note
3. reference_<method>_<regime>.csv
                             organism, condition, reaction_id, reference_flux   (authors' own predictions)

WHY THESE THREE, AND THE ENGINEERING BENEFIT
--------------------------------------------
* transcriptomes.csv  -> the model input. The archive ships expression already on a LINEAR
  scale (the authors exponentiated any log data back). Re-emitting it in one tidy long table
  means your pipeline reads one schema for every condition instead of 8 ad-hoc files: fewer
  code paths, fewer chances for a per-file bug.
* uptakes.csv         -> documents the carbon-source regime. The DC models encode the carbon
  source structurally (active exchange open, the other seven closed), NOT as a numeric rate.
  Capturing that diff makes the constraint regime auditable and reviewer-checkable rather than
  hidden inside 8 SBML files.
* reference_*.csv     -> the authors' OWN published predictions, so you can validate your engine
  prediction-vs-prediction (implementation faithfulness) with ZERO external downloads. This is
  the tightest test of *your code*: if your E-Flux2 output matches theirs reaction-by-reaction,
  any later disagreement with 13C data is about the method, not your implementation.

WHAT IT DOES NOT BUILD
----------------------
measured_fluxes.csv (the 13C ground truth) cannot be built from this archive: the file
`ecoli2_measured_flux9.txt` referenced by the MATLAB correlation script is NOT included.
Reconstruct it from Gerosa 2015, applying the authors' aggregation (read from
cal_correl_ecoli2.m): OR -> sum, AND -> min, leading '-' -> negate. See README.md.

USAGE
-----
    pip install cobra                       # needs python-libsbml (pulled in by cobra)
    python build_benchmark_inputs.py --archive /path/to/extracted/s1 --out ./ecoli_inputs
    # options:
    #   --which all|transcriptomes|uptakes|reference   (default all)
    #   --method E-Flux2|SPOT|Lee|pFBA                  (reference block; default E-Flux2)
    #   --regime AC|DC|FullAC                           (reference file; default AC)

--archive may point either at the extracted root (the dir that CONTAINS
Ecoli_model_and_data/) or directly at Ecoli_model_and_data/. Both are auto-detected.
"""
import argparse, csv, os, sys, contextlib

# Canonical carbon-source order used throughout the dataset (expression files AND the
# 8-column blocks of the predicted-flux files are in exactly this order).
CARBONS = ["acetate", "fructose", "galactose", "gluconate",
           "glucose", "glycerol", "pyruvate", "succinate"]

# The 8 carbon-source exchange reactions in iJO1366 (MOST's _LPAREN_e_RPAREN_ naming).
CARBON_EX = {
    "acetate":   "EX_ac_LPAREN_e_RPAREN_",
    "fructose":  "EX_fru_LPAREN_e_RPAREN_",
    "galactose": "EX_gal_LPAREN_e_RPAREN_",
    "gluconate": "EX_glcn_LPAREN_e_RPAREN_",
    "glucose":   "EX_glc_LPAREN_e_RPAREN_",
    "glycerol":  "EX_glyc_LPAREN_e_RPAREN_",
    "pyruvate":  "EX_pyr_LPAREN_e_RPAREN_",
    "succinate": "EX_succ_LPAREN_e_RPAREN_",
}

# Column-block layout of *_ecoli_predicted_flux.txt: col0 = Rxns, then four 8-wide blocks.
METHOD_BLOCK = {"E-Flux2": 0, "SPOT": 1, "Lee": 2, "pFBA": 3}

ORGANISM = "E. coli"


def find_data_root(archive: str) -> str:
    """Accept either the extracted root or the Ecoli_model_and_data dir itself."""
    cand = os.path.join(archive, "Ecoli_model_and_data")
    if os.path.isdir(cand):
        return cand
    if os.path.isdir(os.path.join(archive, "expression")):
        return archive
    sys.exit(f"[error] could not find Ecoli_model_and_data/ under {archive!r}")


def build_transcriptomes(root: str, out: str) -> None:
    path = os.path.join(out, "transcriptomes.csv")
    n = 0
    with open(path, "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["organism", "condition", "gene_id", "expression"])
        for c in CARBONS:
            src = os.path.join(root, "expression", f"expo_{c}.csv")
            with open(src) as fi:
                for line in fi:
                    line = line.strip()
                    if not line:
                        continue
                    gid, val = line.split(",")          # format: b-number,value (no header)
                    w.writerow([ORGANISM, c, gid, val])
                    n += 1
    print(f"[ok] transcriptomes.csv : {n} rows ({n // len(CARBONS)} genes x {len(CARBONS)} conditions)")


def build_uptakes(root: str, out: str) -> None:
    try:
        from cobra.io import read_sbml_model
    except ImportError:
        sys.exit("[error] cobra not installed. Run: pip install cobra")

    devnull = open(os.devnull, "w")

    def load(p):
        # libsbml prints hundreds of harmless 'discouraged encoding' warnings to stderr;
        # silence them at the fd level so the console stays readable.
        with contextlib.redirect_stderr(devnull):
            return read_sbml_model(p)

    rows = []
    for c in CARBONS:
        dc = load(os.path.join(root, "models", f"DC_{c}_iJO1366.xml"))
        dcb = {r.id: (r.lower_bound, r.upper_bound) for r in dc.reactions}
        for name, exid in CARBON_EX.items():
            if exid not in dcb:
                sys.exit(f"[error] {exid} missing from DC_{c} model")
            lb, ub = dcb[exid]
            note = "OPEN (active carbon source)" if lb < 0 else "closed"
            rows.append([ORGANISM, c, exid, lb, ub, note])

    path = os.path.join(out, "uptakes.csv")
    with open(path, "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["organism", "condition", "exchange_id", "lower_bound", "upper_bound", "note"])
        w.writerows(rows)

    # sanity: each condition must have exactly one OPEN carbon, and it must match the name
    for c in CARBONS:
        opens = [r[2] for r in rows if r[1] == c and "OPEN" in r[5]]
        assert opens == [CARBON_EX[c]], f"DC_{c}: expected only {CARBON_EX[c]} open, got {opens}"
    print(f"[ok] uptakes.csv : {len(rows)} rows (8 carbon exchanges x {len(CARBONS)} conditions); "
          f"active-carbon sanity check passed")


def build_reference(root: str, out: str, method: str, regime: str) -> None:
    if method not in METHOD_BLOCK:
        sys.exit(f"[error] --method must be one of {list(METHOD_BLOCK)}")
    src = os.path.join(root, "predicted_fluxes_and_correlation_scripts",
                       f"{regime}_ecoli_predicted_flux.txt")
    if not os.path.isfile(src):
        sys.exit(f"[error] predicted-flux file not found: {src}")

    block = METHOD_BLOCK[method]
    cols = [1 + block * 8 + i for i in range(8)]       # 8 columns for this method, in CARBON order
    out_name = f"reference_{method.replace('-', '').lower()}_{regime}.csv"
    path = os.path.join(out, out_name)
    n = 0
    with open(src) as fi:
        fi.readline()                                  # skip header row
        with open(path, "w", newline="") as fo:
            w = csv.writer(fo)
            w.writerow(["organism", "condition", "reaction_id", "reference_flux"])
            for line in fi:
                p = line.rstrip("\n").split("\t")
                if len(p) <= cols[-1]:
                    continue
                rid = p[0]
                for ci, c in zip(cols, CARBONS):
                    w.writerow([ORGANISM, c, rid, p[ci]])
                    n += 1
    print(f"[ok] {out_name} : {n} rows ({n // len(CARBONS)} reactions x {len(CARBONS)} conditions) "
          f"[{method} / {regime} — authors' own predictions]")


def main():
    ap = argparse.ArgumentParser(description="Build Bhadra-Lobo 2020 E. coli benchmark input CSVs.")
    ap.add_argument("--archive", required=True,
                    help="Path to the extracted S1 Dataset (root or Ecoli_model_and_data/).")
    ap.add_argument("--out", default="./ecoli_inputs", help="Output directory.")
    ap.add_argument("--which", default="all",
                    choices=["all", "transcriptomes", "uptakes", "reference"])
    ap.add_argument("--method", default="E-Flux2", choices=list(METHOD_BLOCK))
    ap.add_argument("--regime", default="AC", choices=["AC", "DC", "FullAC"])
    args = ap.parse_args()

    root = find_data_root(args.archive)
    os.makedirs(args.out, exist_ok=True)

    if args.which in ("all", "transcriptomes"):
        build_transcriptomes(root, args.out)
    if args.which in ("all", "uptakes"):
        build_uptakes(root, args.out)
    if args.which in ("all", "reference"):
        build_reference(root, args.out, args.method, args.regime)

    print(f"[done] outputs in {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
build_measured_fluxes.py
========================
Extract the 13C measured intracellular net fluxes (the Prong-1 ground truth) from
Gerosa et al. 2015 Data S1 (Elsevier supplement file mmc2.xlsx) into a tidy long table.

Reference: Gerosa L, et al. Cell Systems 1(4):270-282 (2015), DOI 10.1016/j.cels.2015.09.008.
This is the SAME 13C dataset Bhadra-Lobo 2020 used for the E. coli benchmark.

WHY USE GEROSA DATA S1 DIRECTLY (instead of the Bhadra-Lobo bridge file)
-----------------------------------------------------------------------
Bhadra-Lobo's S1 Dataset references `ecoli2_measured_flux9.txt`, which is NOT in their
archive. Gerosa Data S1 is the PRIMARY source and is self-contained: its 'Metabolic fluxes'
sheet carries both the measured module fluxes (8 carbon sources) AND a 'Rnx-Flux relationship'
column mapping each module to model reaction IDs. Using the primary source removes the missing
dependency and is more defensible for review.

WHAT IT EMITS
-------------
measured_flux_modules.csv :
    organism, condition, module_short, module_name, measured_flux, sd, rnx_relationship
  - one row per (module, carbon source); conditions in pipeline-lowercase names
  - measured_flux units: mmol * gCDW^-1 * h^-1
  - rnx_relationship is Gerosa's mapping string (BiGG-style ids; ',' = same-flux chain,
    '+'/'-' = net combination). Aggregation to compare with model predictions:
        comma-separated within parens -> reactions share the module's net flux (use representative/mean)
        '+'  -> add reaction fluxes ;  '-' -> subtract
    A handful of ids need iJO1366 aliases (e.g. ACt2r->ACt2rpp, GLCpts->GLCptspp,
    ZWF->G6PDH2r); resolved in the correlation step, not here.

NOTE ON CARBON ORDER: Gerosa's column order (…Glucose, Glycerol, Gluconate…) differs from the
expression files (…Gluconate, Glucose, Glycerol…). This script keys by NAME, not position.

USAGE
-----
    pip install openpyxl
    python build_measured_fluxes.py --gerosa /path/to/mmc2.xlsx --out ./ecoli_inputs
"""
import argparse, csv, os, sys

# Gerosa 'Metabolic fluxes' sheet layout (0-indexed columns):
#   0        = module short name (flux data rows start at row index 2)
#   1..8     = flux for Acetate,Fructose,Galactose,Glucose,Glycerol,Gluconate,Pyruvate,Succinate
#   9        = blank
#   10..17   = standard deviation, same carbon order
#   18       = blank
#   19       = module number
#   20       = module long name      } these three are offset +1 row vs the flux columns
#   21       = module short name     } (module k's metadata sits one row above flux-row k)
#   22       = Rnx-Flux relationship }
FLUX_COLS = list(range(1, 9))
SD_COLS   = list(range(10, 18))
# carbon header order in the sheet -> pipeline-lowercase condition names
CARBON_ORDER = ["acetate", "fructose", "galactose", "glucose",
                "glycerol", "gluconate", "pyruvate", "succinate"]
ORGANISM = "E. coli"


def main():
    ap = argparse.ArgumentParser(description="Extract 13C measured fluxes from Gerosa 2015 Data S1.")
    ap.add_argument("--gerosa", required=True, help="Path to mmc2.xlsx (Gerosa Data S1).")
    ap.add_argument("--out", default="./ecoli_inputs", help="Output directory.")
    ap.add_argument("--sheet", default="Metabolic fluxes", help="Sheet name (default: 'Metabolic fluxes').")
    args = ap.parse_args()

    try:
        from openpyxl import load_workbook
    except ImportError:
        sys.exit("[error] openpyxl not installed. Run: pip install openpyxl")

    os.makedirs(args.out, exist_ok=True)
    wb = load_workbook(args.gerosa, read_only=True, data_only=True)
    if args.sheet not in wb.sheetnames:
        sys.exit(f"[error] sheet {args.sheet!r} not found. Sheets: {wb.sheetnames}")
    rows = list(wb[args.sheet].iter_rows(values_only=True))

    # module metadata: short -> (long_name, relationship), from the offset block (cols 20-22)
    meta = {}
    for r in rows[1:]:
        if len(r) > 22 and r[21] is not None:
            meta[str(r[21]).strip()] = (
                str(r[20]).strip() if r[20] is not None else "",
                str(r[22]).strip() if r[22] is not None else "",
            )

    out_path = os.path.join(args.out, "measured_flux_modules.csv")
    n = 0
    n_modules = 0
    with open(out_path, "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["organism", "condition", "module_short", "module_name",
                    "measured_flux", "sd", "rnx_relationship"])
        for r in rows[2:]:
            if not r or r[0] is None:
                continue
            short = str(r[0]).strip()
            long_name, rel = meta.get(short, ("", ""))
            n_modules += 1
            for fcol, scol, cond in zip(FLUX_COLS, SD_COLS, CARBON_ORDER):
                flux = r[fcol] if fcol < len(r) else None
                sd = r[scol] if scol < len(r) else None
                w.writerow([ORGANISM, cond, short, long_name, flux, sd, rel])
                n += 1

    print(f"[ok] measured_flux_modules.csv : {n} rows "
          f"({n_modules} modules x {len(CARBON_ORDER)} conditions)")
    print(f"[done] {os.path.abspath(out_path)}")


if __name__ == "__main__":
    main()

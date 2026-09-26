#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_measured_reaction_fluxes.py
=================================
Turn the module-level Gerosa measured fluxes (measured_flux_modules.csv) into the
reaction-level tidy table that prepare_benchmark_inputs.py / validation_01 expect:

    organism, condition, reaction_id, measured_flux

WHY THIS STEP EXISTS
--------------------
Gerosa reports fluxes per "flux module" (e.g. PGI, Middle_EMP) and maps each module to
model reactions in its 'Rnx-Flux relationship' column. validation_01 compares ONE predicted
reaction flux per measured value, so we must reduce each module to a single representative
iJO1366 reaction:

  * single reaction, e.g. "( PGI )"                  -> that reaction (aliased to iJO1366 id)
  * comma-chain,    e.g. "(GAPD, PGK)" / "(EDD,EDA)" -> all members carry the SAME net flux at
                                                        steady state, so any one is a valid
                                                        representative; we take the first that
                                                        exists in iJO1366.
  * composite net,  e.g. "(PFK - FBP)", "(MDH + MQO)" -> a signed COMBINATION of reactions; it
                                                        cannot be a single reaction_id, so it is
                                                        DROPPED from the per-reaction file and
                                                        listed in the mapping report. (~5 modules;
                                                        ~29 of 35 mapped modules remain — ample
                                                        for the benchmark correlation.)

Gerosa's reaction names are BiGG-style but predate iJO1366's periplasm suffixes; the ALIASES
table below resolves the handful that differ (verified to exist in iJO1366). Using Gerosa's
own mapping (the primary 13C source) is more defensible than the unavailable Bhadra-Lobo bridge
file; report r against their published band rather than expecting a bit-identical match.

USAGE
-----
    python build_measured_reaction_fluxes.py \
        --modules ecoli_inputs/measured_flux_modules.csv \
        --model iJO1366.json \
        --out ecoli_inputs/measured_fluxes.csv
    # then feed --out to prepare_benchmark_inputs.py --measured_fluxes
"""
import argparse, os, re, sys
import pandas as pd

# Gerosa token -> iJO1366 reaction id (only where the names differ). All targets verified present.
ALIASES = {
    "ACt2r": "ACt2rpp", "D_LACt2": "D_LACt2pp", "FRUpts": "FRUpts2pp", "GALabc": "GALabcpp",
    "GLCpts": "GLCptspp", "GlcnUpt": "GLCNt2rpp", "GlycUpt": "GLYCtpp", "PYRt2r": "PYRt2rpp",
    "SUCCt": "SUCCt2_2pp", "ZWF": "G6PDH2r", "ACONT": "ACONTa", "ACONT2": "ACONTb",
    "SUCDH3": "SUCDi", "MQO": "MDH2",
    # FUM_SEC is a Gerosa secretion pseudo-reaction with no internal iJO1366 analogue -> unresolved.
}


def load_model_ids(model_path):
    from cobra.io import load_json_model, read_sbml_model
    if model_path.endswith(".json"):
        m = load_json_model(model_path)
        return set(r.id for r in m.reactions)
    else:
        import contextlib
        with contextlib.redirect_stderr(open(os.devnull, "w")):
            m = read_sbml_model(model_path)
        # SBML encodes ids as 'R_<bigg>'; strip to match the BiGG/json namespace
        return set(r.id[2:] if r.id.startswith("R_") else r.id for r in m.reactions)


def is_composite(rel):
    # '+' or '-' anywhere = signed combination of reactions (ids themselves use '_' not '-')
    return ("+" in rel) or ("-" in rel)


def tokens(rel):
    return [t for t in re.split(r"[()\,\s]+", rel) if t]


def resolve_module(rel, model_ids):
    """Return (reaction_id, status). status in {ok, composite, unresolved, nomap}."""
    if not rel or rel.strip() == "":
        return None, "nomap"
    if is_composite(rel):
        return None, "composite"
    for tok in tokens(rel):                      # first member that maps into the model wins
        rid = ALIASES.get(tok, tok)
        if rid in model_ids:
            return rid, "ok"
    return None, "unresolved"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modules", required=True, help="measured_flux_modules.csv from build_measured_fluxes.py")
    ap.add_argument("--model", required=True, help="iJO1366.json (or SBML)")
    ap.add_argument("--out", default="ecoli_inputs/measured_fluxes.csv")
    ap.add_argument("--keep_uptake", action="store_true",
                    help="keep *_Ex uptake/secretion modules (default: drop, intracellular only)")
    args = ap.parse_args()

    df = pd.read_csv(args.modules)
    model_ids = load_model_ids(args.model)

    # resolve each unique module once
    rep, status = {}, {}
    for short, rel in df[["module_short", "rnx_relationship"]].drop_duplicates().itertuples(index=False):
        rid, st = resolve_module(str(rel) if pd.notna(rel) else "", model_ids)
        if (not args.keep_uptake) and str(short).endswith("_Ex"):
            rid, st = None, "uptake_skipped"
        rep[short], status[short] = rid, st

    # emit reaction-level rows for resolved modules
    rows = []
    for r in df.itertuples(index=False):
        rid = rep.get(r.module_short)
        if rid:
            rows.append([r.organism, r.condition, rid, r.measured_flux])
    out = pd.DataFrame(rows, columns=["organism", "condition", "reaction_id", "measured_flux"])
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    out.to_csv(args.out, index=False)

    # mapping report (one row per module) for transparency / methods
    rep_rows = [{"module_short": s, "rnx_relationship": rel, "reaction_id": rep[s], "status": status[s]}
                for s, rel in df[["module_short", "rnx_relationship"]].drop_duplicates().itertuples(index=False)]
    rep_path = args.out.replace(".csv", "_mapping_report.csv")
    pd.DataFrame(rep_rows).to_csv(rep_path, index=False)

    n_ok = sum(1 for v in status.values() if v == "ok")
    n_mod = len(status)
    print(f"[ok] {args.out}: {len(out)} rows "
          f"({n_ok}/{n_mod} modules mapped to a reaction x {out['condition'].nunique()} conditions)")
    drops = {s: status[s] for s in status if status[s] != "ok"}
    if drops:
        print("[info] modules not used (status):")
        for s, st in drops.items():
            print(f"        {s:<14} {st:<14} {dict(zip(df.module_short, df.rnx_relationship)).get(s,'')}")
    print(f"[ok] wrote mapping report -> {rep_path}")


if __name__ == "__main__":
    sys.exit(main())

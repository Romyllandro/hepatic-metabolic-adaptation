#!/usr/bin/env python3
"""
benchmark_ecoli_13c_methods.py
==============================
13C-MFA benchmark of established transcriptome-constrained methods.

Purpose
-------
This is the reviewer-facing *external benchmark*. It evaluates:
- standard E-Flux + pFBA
- GIMME + pFBA representative solution
- unconstrained pFBA baseline

against Gerosa et al. measured module fluxes. The liver-specific three-layer
architecture is not forced onto E. coli, because its Layer-1 diet exchanges and
liver-specific transporter/enzyme separation are not the same experimental
construct. The three-layer-vs-established-method comparison is performed
separately by benchmark_liver_methods.py.

Inputs
------
--model                 iJO1366.json
--objective             BIOMASS_Ec_iJO1366_core_53p95M
--transcriptomes        tidy CSV: organism,condition,gene_id,expression
--measured_modules      output from build_measured_fluxes.py:
                        organism,condition,module_short,module_name,
                        measured_flux,sd,rnx_relationship
--uptakes               optional tidy uptake CSV for DC regime:
                        organism,condition,exchange_id,lower_bound,upper_bound

Gerosa's rnx_relationship is evaluated at module level. Comma-separated
reaction groups are treated as same-flux chains and summarized by their median;
+/- terms are added/subtracted. A mapping audit is written so unresolved modules
are visible rather than silently discarded.
"""
from __future__ import annotations
import argparse, json, re, sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
sys.path.insert(0,str(HERE.parent/"common"))
from benchmark_methods_core import (
    apply_eflux,gimme_solve,riptide_solve,uncentered_pearson,centered_pearson,
    safe_spearman,scaled_nrmse
)
from deterministic_solver_v2 import deterministic_optimize, configure_solver_deterministic

DEFAULT_ALIASES = {
    "ACt2r":"ACt2rpp",
    "GLCpts":"GLCptspp",
    "ZWF":"G6PDH2r",
}

def resolve_exchange(model, rid):
    # BUGFIX 2026-09-23: the previous version only expanded __D_e/__L_e stereo
    # variants off the ORIGINAL token (checking r.endswith("_e")), not off each
    # generated candidate -- so an SBML-escaped id like "EX_glc_LPAREN_e_RPAREN_"
    # (which ends in "_RPAREN_", not "_e") never got its "EX_glc_e" candidate
    # expanded to "EX_glc__D_e", the model's actual stereo-qualified glucose
    # exchange id. That silently dropped glucose's DC-regime uptake bound for
    # every condition (left open/closed at whatever the model's default was,
    # identical to the AC regime) instead of raising a resolution error.
    # validation_01_benchmark_recovery.py's _resolve_exchange_id() already does
    # this correctly (expands every "_e"-ending candidate, not just the raw
    # token) -- this mirrors that fix.
    candidates=[str(rid)]
    r=str(rid)
    if "_LPAREN_e_RPAREN_" in r:
        candidates.append(r.replace("_LPAREN_e_RPAREN_","_e"))
    if r.startswith("R_"): candidates.append(r[2:])
    for c in list(candidates):
        if c.endswith("_e"):
            b=c[:-2]; candidates.extend([b+"__D_e",b+"__L_e"])
    for c in candidates:
        try:return model.reactions.get_by_id(c)
        except KeyError: pass
    return None

def apply_uptake(model, sub):
    n=0; missing=[]
    for _,r in sub.iterrows():
        rxn=resolve_exchange(model,r["exchange_id"])
        if rxn is None:
            missing.append(str(r["exchange_id"])); continue
        rxn.lower_bound=float(r["lower_bound"])
        rxn.upper_bound=float(r["upper_bound"])
        n+=1
    return n,missing

def _normalize_rid(x):
    return re.sub(r"[^A-Za-z0-9_]","",str(x).strip().strip("'\""))

def build_resolver(model, aliases):
    ids={r.id for r in model.reactions}
    norm={re.sub(r"[^a-z0-9]","",r.id.lower()):r.id for r in model.reactions}
    def resolve(token):
        t=_normalize_rid(token)
        cand=[t,aliases.get(t,t)]
        if t.startswith("R_"): cand.append(t[2:])
        for c in list(cand):
            if c.endswith("pp"): cand.append(c[:-2])
            else: cand.append(c+"pp")
        for c in cand:
            if c in ids:return c
        key=re.sub(r"[^a-z0-9]","",t.lower())
        return norm.get(key)
    return resolve

def split_top_level(expr):
    """
    Split an expression into signed top-level terms while preserving parentheses.
    Example: -(A,B)+C-D -> [(-1,'(A,B)'),(+1,'C'),(-1,'D')].
    """
    s=str(expr).replace(" ","")
    out=[]; depth=0; start=0; sign=1
    if s.startswith("+"): start=1
    elif s.startswith("-"): sign=-1; start=1
    i=start
    while i<len(s):
        ch=s[i]
        if ch=="(": depth+=1
        elif ch==")": depth=max(0,depth-1)
        elif depth==0 and ch in "+-":
            term=s[start:i]
            if term: out.append((sign,term))
            sign=1 if ch=="+" else -1
            start=i+1
        i+=1
    term=s[start:]
    if term: out.append((sign,term))
    return out

def group_tokens(term):
    t=term.strip()
    if t.startswith("(") and t.endswith(")"):
        t=t[1:-1]
    return [x for x in re.split(r"[,;]",t) if x]

def predicted_module_flux(fluxes, relationship, resolve, zero_fill_pruned=True):
    if relationship is None or str(relationship).strip()=="" or str(relationship).lower()=="nan":
        return np.nan,[],["missing relationship"]
    resolved=[]; issues=[]; total=0.0
    for sign,term in split_top_level(relationship):
        vals=[]; ids=[]
        for tok in group_tokens(term):
            rid=resolve(tok)
            if rid is None:
                # Genuinely unresolvable token (bad mapping/typo) -- a real data
                # problem, always excluded from scoring regardless of zero_fill_pruned.
                issues.append(f"unresolved:{tok}")
                continue
            if rid not in fluxes.index:
                # DESIGN 2026-09-23 (Roland's decision): rid IS a real model reaction
                # (resolve() succeeded) but is absent from THIS solution's flux index
                # -- i.e. the method (currently only RIPTiDe) pruned it out of its
                # own contextualized model. That is a modeling decision the method
                # made, not a mapping failure, so by default we score it as that
                # method's own implicit prediction of zero flux rather than quietly
                # shrinking this method's denominator relative to every other method
                # in the comparison. Tagged "pruned:" (not "unresolved:") in the
                # mapping audit so it stays traceable and distinguishable from a
                # genuine unresolved-token issue.
                if zero_fill_pruned:
                    vals.append(0.0); ids.append(rid)
                    issues.append(f"pruned:{tok}")
                else:
                    issues.append(f"unresolved:{tok}")
                continue
            vals.append(float(fluxes.loc[rid])); ids.append(rid)
        if not vals:
            continue
        # Comma groups represent same module flux in Gerosa's relationship table.
        total += sign*float(np.median(vals))
        resolved.extend(ids)
    if not resolved:
        return np.nan,resolved,issues
    return total,resolved,issues

def solve_method(base,expr,method,args):
    m=base.copy()
    configure_solver_deterministic(m,solver=args.solver,verbose=False)
    if method=="pfba":
        return deterministic_optimize(m,use_pfba=True,fraction_of_optimum=1.0),{}
    if method=="eflux":
        scope=[r.id for r in m.reactions if not r.boundary and not r.id.startswith("EX_")]
        changed=apply_eflux(m,expr,scope,quantile=args.eflux_quantile,
                            floor=args.eflux_floor,cap=args.eflux_cap,
                            scale=args.eflux_scale)
        sol=deterministic_optimize(m,use_pfba=True,fraction_of_optimum=1.0)
        return sol,{"reactions_constrained":changed}
    if method=="gimme":
        return gimme_solve(m,expr,args.objective,
                           threshold_quantile=args.gimme_threshold_quantile,
                           objective_fraction=args.gimme_objective_fraction,
                           solver=args.solver)
    if method=="riptide":
        return riptide_solve(m,expr,args.objective,
                             fraction=args.riptide_fraction,
                             solver=args.solver)
    raise ValueError(method)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--model",required=True)
    ap.add_argument("--objective",required=True)
    ap.add_argument("--transcriptomes",required=True)
    ap.add_argument("--measured_modules",required=True)
    ap.add_argument("--uptakes")
    ap.add_argument("--organism",default="E. coli")
    ap.add_argument("--regime",choices=["AC","DC","both"],default="both")
    ap.add_argument("--solver",default="gurobi")
    ap.add_argument("--methods",default="eflux,gimme,pfba")
    ap.add_argument("--eflux_quantile",type=float,default=0.95)
    ap.add_argument("--eflux_floor",type=float,default=0.1)
    ap.add_argument("--eflux_cap",type=float,default=1000.0)
    ap.add_argument("--eflux_scale",type=float,default=1000.0,
                    help="Capacity scale for the external E. coli E-Flux benchmark.")
    ap.add_argument("--gimme_threshold_quantile",type=float,default=0.25)
    ap.add_argument("--gimme_objective_fraction",type=float,default=0.90)
    ap.add_argument("--riptide_fraction",type=float,default=0.8,
                    help="Fraction of the optimal objective RIPTiDe must retain "
                         "when pruning (riptide.contextualize's own fraction=).")
    ap.add_argument("--no_zero_fill_pruned",action="store_true",
                    help="By default, a reaction a method prunes out of its own "
                         "contextualized model (e.g. RIPTiDe) is scored as that "
                         "method's own implicit prediction of zero flux, so every "
                         "method is compared on the same fixed reaction panel. Pass "
                         "this flag to instead drop pruned reactions from the metric "
                         "(shrinks that method's own denominator).")
    ap.add_argument("--alias_json")
    ap.add_argument("--output",default="benchmark_ecoli_13c")
    args=ap.parse_args()

    import cobra
    from cobra.io import load_json_model
    out=Path(args.output); out.mkdir(parents=True,exist_ok=True)
    tx=pd.read_csv(args.transcriptomes)
    mm=pd.read_csv(args.measured_modules)
    up=pd.read_csv(args.uptakes) if args.uptakes else None
    org=args.organism.strip().lower()
    tx=tx[tx.organism.astype(str).str.strip().str.lower()==org]
    mm=mm[mm.organism.astype(str).str.strip().str.lower()==org]
    if up is not None:
        up=up[up.organism.astype(str).str.strip().str.lower()==org]
    conditions=sorted(set(tx.condition.astype(str)) & set(mm.condition.astype(str)))
    if not conditions: raise SystemExit("No shared benchmark conditions.")

    aliases=dict(DEFAULT_ALIASES)
    if args.alias_json:
        aliases.update(json.loads(Path(args.alias_json).read_text()))
    methods=[x.strip() for x in args.methods.split(",") if x.strip()]
    regimes=["AC","DC"] if args.regime=="both" else [args.regime]

    rows=[]; map_rows=[]; long_rows=[]
    for cond in conditions:
        e=tx[tx.condition.astype(str)==cond]
        expr=dict(zip(e.gene_id.astype(str),pd.to_numeric(e.expression,errors="coerce")))
        meas=mm[mm.condition.astype(str)==cond].copy()
        for regime in regimes:
            base=load_json_model(args.model); base.objective=args.objective
            if regime=="DC":
                if up is None: continue
                n,missing=apply_uptake(base,up[up.condition.astype(str)==cond])
                if n==0:
                    print(f"[WARN] {cond} DC: no uptake bounds applied; skipped")
                    continue
            resolver=build_resolver(base,aliases)
            for method in methods:
                try:
                    sol,meta=solve_method(base,expr,method,args)
                except Exception as exc:
                    rows.append({"condition":cond,"regime":regime,"method":method,
                                 "status":"ERROR","error":str(exc)})
                    continue
                pred=[]; obs=[]
                for _,r in meas.iterrows():
                    val,resolved,issues=predicted_module_flux(
                        sol.fluxes,r.rnx_relationship,resolver,
                        zero_fill_pruned=not args.no_zero_fill_pruned)
                    map_rows.append({
                        "condition":cond,"module_short":r.module_short,
                        "relationship":r.rnx_relationship,
                        "resolved_reactions":";".join(resolved),
                        "issues":";".join(issues)
                    })
                    status="OK" if np.isfinite(val) else (
                        "pruned_zero_filled" if any(i.startswith("pruned:") for i in issues)
                        else "unresolved")
                    long_rows.append({
                        "condition":cond,"regime":regime,"method":method,
                        "reaction_id":r.module_short,
                        "measured_flux":float(r.measured_flux) if pd.notna(r.measured_flux) else np.nan,
                        "predicted_flux":float(val) if np.isfinite(val) else np.nan,
                        "status":status,
                    })
                    if np.isfinite(val) and pd.notna(r.measured_flux):
                        pred.append(abs(float(val))); obs.append(abs(float(r.measured_flux)))
                pred=np.asarray(pred); obs=np.asarray(obs)
                nrmse,alpha=scaled_nrmse(pred,obs)
                rows.append({
                    "condition":cond,"regime":regime,"method":method,"status":"OK",
                    "n_modules":len(pred),
                    "uncentered_pearson":uncentered_pearson(pred,obs),
                    "centered_pearson":centered_pearson(pred,obs),
                    "spearman":safe_spearman(pred,obs),
                    "scaled_nrmse":nrmse,"least_squares_scale":alpha,
                    **meta
                })
                print(f"{cond:10s} {regime} {method:8s} n={len(pred):2d} "
                      f"r_unc={rows[-1]['uncentered_pearson']:.3f}")

    df=pd.DataFrame(rows); df.to_csv(out/"benchmark_results_by_condition.csv",index=False)
    pd.DataFrame(map_rows).drop_duplicates().to_csv(out/"gerosa_module_mapping_audit.csv",index=False)
    pd.DataFrame(long_rows).to_csv(out/"benchmark_long_format_flux.csv",index=False)
    ok=df[df.status=="OK"].copy()
    if len(ok):
        summary=(ok.groupby(["regime","method"])
                 .agg(n_conditions=("condition","nunique"),
                      mean_uncentered_pearson=("uncentered_pearson","mean"),
                      mean_centered_pearson=("centered_pearson","mean"),
                      mean_spearman=("spearman","mean"),
                      mean_scaled_nrmse=("scaled_nrmse","mean"))
                 .reset_index())
        summary.to_csv(out/"benchmark_method_summary.csv",index=False)
        print("\n",summary.to_string(index=False))
    print(f"[OK] outputs: {out}")

if __name__=="__main__":
    main()

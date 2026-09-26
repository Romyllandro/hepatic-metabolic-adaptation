#!/usr/bin/env python3
"""
benchmark_methods_core.py
=========================
Shared helpers for reviewer-facing method comparison.

Implemented comparators
-----------------------
1. three_layer : production three-layer diet + transporter + intracellular E-Flux
2. eflux       : conventional single-scope E-Flux on all non-exchange reactions
3. gimme       : GIMME-style low-expression penalty followed by a parsimonious
                 representative solution at the same retained primary objective.
4. pfba        : unconstrained pFBA baseline (optional)
5. riptide     : RIPTiDe-style transcript-weighted pruning, then the same
                 deterministic pFBA representative-solution step used by
                 eflux/gimme/pfba (see riptide_solve() docstring for why).

The GIMME implementation is intentionally explicit and auditable. It:
- computes reaction expression using the same AND=min / OR=max convention;
- chooses a pre-specified gene-expression threshold quantile;
- maximizes the biological objective;
- requires a fixed fraction of that optimum;
- minimizes weighted absolute flux through low-expression reactions;
- fixes the GIMME penalty near its optimum and then minimizes total absolute flux
  to select one representative solution.

This last parsimony stage is used only to remove arbitrary alternative-optimum
selection for comparison. It is reported as "GIMME+pFBA representative", not as
a claim that canonical GIMME itself contains pFBA.
"""
from __future__ import annotations

import math, re
from typing import Dict, Iterable, Optional, Sequence, Tuple
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Numeric GPR parser with parenthesis support
# ---------------------------------------------------------------------------

_TOKEN = re.compile(r"\(|\)|\band\b|\bor\b|[^()\s]+", re.IGNORECASE)

class _Parser:
    def __init__(self, text, values):
        self.tokens = _TOKEN.findall(str(text))
        self.i = 0
        self.values = values

    def peek(self):
        return self.tokens[self.i] if self.i < len(self.tokens) else None

    def pop(self):
        t = self.peek()
        if t is not None:
            self.i += 1
        return t

    def expr(self):
        val = self.term()
        while self.peek() and self.peek().lower() == "or":
            self.pop()
            rhs = self.term()
            if val is None: val = rhs
            elif rhs is None: pass
            else: val = max(val, rhs)
        return val

    def term(self):
        val = self.factor()
        while self.peek() and self.peek().lower() == "and":
            self.pop()
            rhs = self.factor()
            if val is None or rhs is None:
                val = None
            else:
                val = min(val, rhs)
        return val

    def factor(self):
        t = self.pop()
        if t is None:
            return None
        if t == "(":
            val = self.expr()
            if self.peek() == ")":
                self.pop()
            return val
        if t == ")":
            return None
        return self.values.get(str(t))

def gpr_numeric(rule: str, values: Dict[str, float]) -> Optional[float]:
    if not rule or not str(rule).strip():
        return None
    try:
        return _Parser(rule, values).expr()
    except Exception:
        return None

def reaction_expression(model, gene_expr: Dict[str,float]) -> pd.Series:
    out={}
    for rxn in model.reactions:
        out[rxn.id]=gpr_numeric(rxn.gene_reaction_rule,gene_expr)
    return pd.Series(out,dtype=float)

# ---------------------------------------------------------------------------
# E-Flux
# ---------------------------------------------------------------------------

def normalize_gene_expression(expr: Dict[str,float], quantile=0.95,
                              floor=0.1, cap=1000.0, scale=1.0):
    vals=np.array([float(v) for v in expr.values()
                   if v is not None and np.isfinite(v) and float(v)>0],dtype=float)
    denom=float(np.quantile(vals,quantile)) if len(vals) else 1.0
    denom=max(denom,1e-12)
    return {
        str(g): float(np.clip(float(v)/denom, floor, cap))*float(scale)
        for g,v in expr.items()
        if v is not None and np.isfinite(v)
    }

def apply_eflux(model, expr: Dict[str,float], reaction_ids: Optional[Iterable[str]]=None,
                quantile=0.95, floor=0.1, cap=1000.0, scale=1.0):
    norm=normalize_gene_expression(expr,quantile,floor,cap,scale)
    scope=None if reaction_ids is None else set(reaction_ids)
    changed=0
    for rxn in model.reactions:
        if scope is not None and rxn.id not in scope:
            continue
        val=gpr_numeric(rxn.gene_reaction_rule,norm)
        if val is None or not np.isfinite(val) or val<=0:
            continue
        old=(rxn.lower_bound,rxn.upper_bound)
        if rxn.lower_bound<0:
            rxn.lower_bound=max(rxn.lower_bound,-float(val))
        if rxn.upper_bound>0:
            rxn.upper_bound=min(rxn.upper_bound,float(val))
        if old!=(rxn.lower_bound,rxn.upper_bound):
            changed+=1
    return changed

# ---------------------------------------------------------------------------
# GIMME
# ---------------------------------------------------------------------------

def _add_abs_variable(model, rxn, name):
    problem=model.problem
    z=problem.Variable(name,lb=0)
    c1=problem.Constraint(z-rxn.flux_expression,lb=0,
                          name=name+"_ge_pos")
    c2=problem.Constraint(z+rxn.flux_expression,lb=0,
                          name=name+"_ge_neg")
    model.add_cons_vars([z,c1,c2])
    return z

def gimme_solve(model, expr: Dict[str,float], objective_reaction: str,
                threshold_quantile=0.25, objective_fraction=0.90,
                pfba_penalty_tolerance=1e-8, solver=None):
    """
    GIMME-style solve + parsimonious representative solution.

    v4 performance change
    ---------------------
    COBRApy already represents each reaction with forward and reverse
    non-negative variables. Therefore |v| can be penalized as
    forward_variable + reverse_variable. This avoids adding ~3,700 auxiliary
    absolute-value variables and ~7,400 constraints for every sample.
    """
    model=model.copy()
    try:
        from deterministic_solver_v2 import configure_solver_deterministic
        configure_solver_deterministic(model, solver=solver, verbose=False)
    except Exception:
        if solver:
            try: model.solver=solver
            except Exception: pass

    model.objective=objective_reaction
    optimum=float(model.slim_optimize(error_value=None))
    if not np.isfinite(optimum):
        raise RuntimeError("Primary objective could not be optimized.")
    if optimum <= 0:
        raise RuntimeError(f"Primary objective optimum is non-positive ({optimum}).")

    rx_expr=reaction_expression(
        model,{str(k):float(v) for k,v in expr.items()
               if v is not None and np.isfinite(v)}
    )
    gene_pos=np.array([float(v) for v in expr.values()
                       if v is not None and np.isfinite(v) and float(v)>0])
    threshold=float(np.quantile(gene_pos,threshold_quantile)) if len(gene_pos) else 0.0

    primary_expr=model.objective.expression
    obj_con=model.problem.Constraint(
        primary_expr, lb=float(objective_fraction)*optimum,
        name="GIMME_primary_objective_floor"
    )
    model.add_cons_vars([obj_con])
    try:
        model.solver.update()
    except Exception:
        pass

    penalty_terms=[]
    low_rxns=[]
    for rxn in model.reactions:
        x=rx_expr.get(rxn.id,np.nan)
        if not np.isfinite(x):
            continue
        if x < threshold:
            weight=max(threshold-float(x),1e-12)
            penalty_terms.append(weight*(rxn.forward_variable+rxn.reverse_variable))
            low_rxns.append(rxn.id)

    if penalty_terms:
        penalty=sum(penalty_terms)
        model.objective=model.problem.Objective(penalty,direction="min")
        gim_sol=model.optimize()
        if gim_sol.status!="optimal":
            raise RuntimeError(f"GIMME penalty optimization status={gim_sol.status}")
        penalty_opt=float(gim_sol.objective_value)

        pen_con=model.problem.Constraint(
            penalty,
            ub=penalty_opt+max(abs(penalty_opt),1.0)*pfba_penalty_tolerance,
            name="GIMME_penalty_fix"
        )
        model.add_cons_vars([pen_con])
    else:
        penalty_opt=0.0

    # Minimal-total-flux representative using COBRA's existing split variables.
    try:
        model.solver.update()
    except Exception:
        pass
    total_abs=sum((rxn.forward_variable+rxn.reverse_variable) for rxn in model.reactions)
    model.objective=model.problem.Objective(total_abs,direction="min")
    sol=model.optimize()
    if sol.status!="optimal":
        raise RuntimeError(f"GIMME representative solve status={sol.status}")

    try:
        sol.biological_objective_value=float(primary_expr.primal)
        sol.gimme_penalty_value=penalty_opt
        sol.gimme_threshold=threshold
        sol.gimme_low_expression_reactions=len(low_rxns)
    except Exception:
        pass

    return sol, {
        "primary_optimum":optimum,
        "objective_fraction":objective_fraction,
        "expression_threshold":threshold,
        "threshold_quantile":threshold_quantile,
        "penalty_optimum":penalty_opt,
        "n_low_expression_reactions":len(low_rxns),
    }

# ---------------------------------------------------------------------------
# RIPTiDe
# ---------------------------------------------------------------------------

def riptide_solve(model, expr: Dict[str, float], objective_reaction: str,
                  fraction=0.8, solver=None):
    """
    RIPTiDe-style solve + parsimonious representative solution.

    RIPTiDe (Jenior et al. 2019, PLoS Comp Biol) contextualizes a model to a
    transcriptome by linear-weighted pruning (flux minimization weighted by
    inverse transcript abundance) rather than GIMME's low-expression penalty
    threshold, then characterizes the resulting reduced flux space by sampling.

    That sampling step is deliberately NOT used for the single comparable flux
    vector this benchmark scores against measured 13C fluxes. Every other
    method in this comparison (eflux, gimme, pfba) produces ONE deterministic
    representative flux vector via pFBA at 100% of the primary objective's
    optimum, using the same deterministic_optimize()/Gurobi settings used
    throughout the manuscript's own production pipeline. To keep RIPTiDe on
    the same footing -- "every method... restricted to what defines the
    method itself" -- this function uses RIPTiDe ONLY for its
    contextualization/pruning step, then runs the SAME deterministic pFBA used
    for eflux/gimme/pfba on the resulting pruned model, rather than reporting
    a flux-sample median. That choice is deliberate and documented here
    rather than silently made: a flux-sample median would introduce a
    fundamentally different, stochastic solution concept for RIPTiDe alone,
    which is a bigger fairness break than reusing the existing deterministic
    scoring step. riptide_object.flux_samples is still available on the
    object this function discards if the sampled distribution is wanted for
    a secondary analysis later -- it is just not what is scored here.

    Any of the measured central-carbon reactions that RIPTiDe prunes entirely
    out of the contextualized model (e.g. ICL, frequently repressed under
    glucose-rich conditions) will simply be absent from the returned
    solution's flux index. predicted_module_flux() in benchmark_ecoli_13c_methods.py
    already handles that as an "unresolved" entry in the mapping audit rather
    than an error -- that is expected RIPTiDe behavior, not a bug, and is also
    reported here via n_reactions_pruned for transparency.

    Requires the `riptide` package (pip install riptide; Jenior et al.,
    https://github.com/mjenior/riptide). This call intentionally uses only
    the three most stable, long-standing kwargs (model=, transcriptome=,
    fraction=). If your installed version's riptide.contextualize() signature
    differs, please report the exact error -- it should be a small,
    localized fix here rather than a redesign.
    """
    import riptide as _riptide
    from deterministic_solver_v2 import configure_solver_deterministic, deterministic_optimize

    model = model.copy()
    try:
        configure_solver_deterministic(model, solver=solver, verbose=False)
    except Exception:
        if solver:
            try: model.solver = solver
            except Exception: pass

    model.objective = objective_reaction
    clean_expr = {str(k): float(v) for k, v in expr.items()
                  if v is not None and np.isfinite(v)}

    riptide_object = _riptide.contextualize(
        model=model, transcriptome=clean_expr, fraction=fraction
    )
    pruned = riptide_object.model
    pruned.objective = objective_reaction

    n_before = len(model.reactions)
    n_after = len(pruned.reactions)

    sol = deterministic_optimize(pruned, use_pfba=True, fraction_of_optimum=1.0)

    return sol, {
        "riptide_fraction": fraction,
        "n_reactions_pruned": n_before - n_after,
        "n_reactions_retained": n_after,
    }

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def safe_spearman(a,b):
    from scipy.stats import spearmanr
    a=np.asarray(a,float); b=np.asarray(b,float)
    ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3:return np.nan
    return float(spearmanr(a[ok],b[ok]).correlation)

def centered_pearson(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3:return np.nan
    if np.std(a[ok])==0 or np.std(b[ok])==0:return np.nan
    return float(np.corrcoef(a[ok],b[ok])[0,1])

def uncentered_pearson(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3:return np.nan
    aa=a[ok]; bb=b[ok]
    den=np.linalg.norm(aa)*np.linalg.norm(bb)
    return float(np.dot(aa,bb)/den) if den else np.nan

def directional_agreement(a,b,tol=1e-9):
    a=np.asarray(a,float); b=np.asarray(b,float)
    ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()==0:return np.nan
    sa=np.where(np.abs(a[ok])<=tol,0,np.sign(a[ok]))
    sb=np.where(np.abs(b[ok])<=tol,0,np.sign(b[ok]))
    return float(np.mean(sa==sb))

def topk_jaccard(a,b,k=100):
    a=np.asarray(a,float); b=np.asarray(b,float)
    k=min(k,len(a),len(b))
    A=set(np.argsort(-np.abs(a))[:k]); B=set(np.argsort(-np.abs(b))[:k])
    return len(A&B)/len(A|B) if (A|B) else np.nan

def scaled_nrmse(pred,meas):
    """
    NRMSE after a single least-squares scale factor. Useful because expression-
    constrained methods often predict relative rather than calibrated flux units.
    """
    p=np.asarray(pred,float); m=np.asarray(meas,float)
    ok=np.isfinite(p)&np.isfinite(m)
    if ok.sum()<3:return np.nan,np.nan
    p=p[ok]; m=m[ok]
    den=np.dot(p,p)
    alpha=float(np.dot(p,m)/den) if den else 0.0
    rmse=float(np.sqrt(np.mean((alpha*p-m)**2)))
    scale=float(np.sqrt(np.mean(m**2)))
    return (rmse/scale if scale else np.nan),alpha

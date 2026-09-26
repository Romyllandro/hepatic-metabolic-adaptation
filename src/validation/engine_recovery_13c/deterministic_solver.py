#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
deterministic_solver.py
=======================
Cross-platform deterministic solver configuration for COBRApy pipelines.
Combines strict algorithmic constraints for Gurobi with aggressive 
floating-point truncation to ensure identical CSV outputs on Windows/Linux.
"""

from __future__ import annotations

import logging
from typing import Optional
import numpy as np
from cobra.flux_analysis import pfba

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration Constants
# ---------------------------------------------------------------------------

DEFAULT_ZERO_CUTOFF = 1e-10
DEFAULT_ROUND_DECIMALS = 12

GUROBI_DETERMINISTIC_PARAMS: dict = {
    "Method":          1,      # 1 = dual simplex (avoids AVX/hardware specific barrier math)
    "Threads":         1,      # Eliminates race conditions in threading schedulers
    "Crossover":       0,      # Disable barrier-to-basis crossover
    "NumericFocus":    3,      # Max precision over speed
    "ScaleFlag":       0,      # Disable auto-scaling heuristics
    "FeasibilityTol":  1e-9,   # Strict feasibility
    "OptimalityTol":   1e-9,   # Strict optimality
    "Seed":            42,     # Fixed tie-breaking
}

GLPK_DETERMINISTIC_PARAMS: dict = {}

# ---------------------------------------------------------------------------
# Core Solver Configuration
# ---------------------------------------------------------------------------

def configure_solver_deterministic(model, solver: Optional[str] = None, verbose: bool = True) -> None:
    """Apply deterministic solver settings to a COBRApy model in-place."""
    _log = logger.info if verbose else logger.debug

    if solver is not None:
        try:
            model.solver = solver
            _log(f"[DET_SOLVER] Solver set to: {solver}")
        except Exception as exc:
            raise RuntimeError(f"Cannot assign solver '{solver}' to model: {exc}") from exc

    solver_name = _detect_solver_name(model)
    _log(f"[DET_SOLVER] Detected solver interface: {solver_name}")

    if "gurobi" in solver_name.lower():
        _apply_gurobi_params(model, _log)
    elif "glpk" in solver_name.lower():
        _apply_glpk_params(model, _log)
    elif "cplex" in solver_name.lower():
        _apply_cplex_params(model, _log)
    else:
        logger.warning(f"[DET_SOLVER] Unknown solver '{solver_name}'; results may not be reproducible.")

def _detect_solver_name(model) -> str:
    try:
        return model.solver.interface.__name__.lower()
    except AttributeError:
        pass
    try:
        return type(model.solver).__name__.lower()
    except AttributeError:
        return "unknown"

def _apply_gurobi_params(model, log_fn) -> None:
    try:
        grb_model = model.solver.problem
    except AttributeError:
        logger.warning("[DET_SOLVER] Cannot access gurobi model; parameters not applied.")
        return

    applied, failed = [], []
    for param, value in GUROBI_DETERMINISTIC_PARAMS.items():
        try:
            setattr(grb_model.Params, param, value)
            applied.append(f"{param}={value}")
        except Exception as exc:
            failed.append(f"{param} ({exc})")

    if applied:
        log_fn(f"[DET_SOLVER] Gurobi params applied: {', '.join(applied)}")
    if failed:
        logger.warning(f"[DET_SOLVER] Gurobi params FAILED to set: {', '.join(failed)}")

def _apply_glpk_params(model, log_fn) -> None:
    log_fn("[DET_SOLVER] GLPK detected — inherently deterministic, no params needed.")

def _apply_cplex_params(model, log_fn) -> None:
    try:
        cpx = model.solver.problem
        cpx.parameters.parallel.set(1)
        cpx.parameters.threads.set(1)
        cpx.parameters.randomseed.set(42)
        log_fn("[DET_SOLVER] CPLEX params applied: parallel=1, threads=1, randomseed=42")
    except Exception as exc:
        logger.warning(f"[DET_SOLVER] CPLEX determinism params could not be set: {exc}")

# ---------------------------------------------------------------------------
# Output Canonicalization & Execution
# ---------------------------------------------------------------------------

def canonicalize_solution(solution, model, zero_cutoff: float = DEFAULT_ZERO_CUTOFF, decimals: int = DEFAULT_ROUND_DECIMALS):
    """Zero out tiny numerical noise and force strict reaction indexing order."""
    if solution is None or not hasattr(solution, "fluxes"):
        return solution
    
    try:
        # Reindex to guarantee exact same row order across OS outputs
        order = [r.id for r in model.reactions]
        fluxes = solution.fluxes.reindex(order)
        
        # Aggressive truncation of microscopic floating point artifacts
        fluxes = fluxes.where(np.abs(fluxes) >= zero_cutoff, 0.0).round(decimals)
        solution.fluxes = fluxes
    except Exception as e:
        logger.warning(f"[DET_SOLVER] Canonicalization partial failure: {e}")
        try:
            solution.fluxes.loc[np.abs(solution.fluxes) < zero_cutoff] = 0.0
            solution.fluxes = solution.fluxes.round(decimals)
        except Exception:
            pass
            
    return solution

def deterministic_optimize(model, use_pfba: bool = True, fraction_of_optimum: float = 1.0):
    """
    The master execution function. 
    Applies strict solver parameters, runs the optimization, and sanitizes the output.
    """
    # 1. Lock down the solver math
    configure_solver_deterministic(model, verbose=False)

    biological_objective = None
    try:
        biological_objective = float(model.slim_optimize(error_value=None))
    except Exception:
        pass

    # 2. Execute
    if use_pfba:
        solution = pfba(model, fraction_of_optimum=fraction_of_optimum)
        if biological_objective is not None:
            try:
                solution.biological_objective_value = biological_objective
            except Exception:
                pass
    else:
        solution = model.optimize()
        try:
            solution.biological_objective_value = float(solution.objective_value)
        except Exception:
            pass

    # 3. Sanitize floating point artifacts
    return canonicalize_solution(solution, model)
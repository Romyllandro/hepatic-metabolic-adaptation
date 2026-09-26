#!/usr/bin/env python3
"""
deterministic_solver_v2.py
==========================
Deterministic COBRApy/Gurobi helper aligned with the revised RQ1/RQ2 production reruns.

Production-aligned Gurobi settings
----------------------------------
Threads=1
Method=1
NumericFocus=3
Seed=1

Crossover is not forced because Method=1 is dual simplex (no barrier crossover step).
ScaleFlag is left at the solver default because disabling scaling can worsen numerical
conditioning; the earlier helper changed it without a scientific need.

The helper also separates the biological objective value from the pFBA secondary
minimum-total-flux objective. COBRApy's pfba() Solution.objective_value is the
secondary objective and must not be reported as biomass/ATPM.
"""
from __future__ import annotations

from typing import Optional
import numpy as np

DEFAULT_ZERO_CUTOFF = 1e-10
DEFAULT_ROUND_DECIMALS = 12

GUROBI_DETERMINISTIC_PARAMS = {
    "Threads": 1,
    "Method": 1,
    "NumericFocus": 3,
    "Seed": 1,
}

def _solver_name(model) -> str:
    try:
        return model.solver.interface.__name__.lower()
    except Exception:
        return type(model.solver).__name__.lower()

def configure_solver_deterministic(model, solver: Optional[str] = None, verbose: bool = False):
    if solver:
        model.solver = solver
    name = _solver_name(model)
    applied = {}
    if "gurobi" in name:
        problem = model.solver.problem
        for key, val in GUROBI_DETERMINISTIC_PARAMS.items():
            try:
                problem.setParam(key, val)
                applied[key] = val
            except Exception:
                try:
                    setattr(problem.Params, key, val)
                    applied[key] = val
                except Exception:
                    pass
    if verbose:
        print(f"[DET] backend={name}; applied={applied}")
    return applied

def canonicalize_solution(solution, model, zero_cutoff=DEFAULT_ZERO_CUTOFF,
                          decimals=DEFAULT_ROUND_DECIMALS):
    if solution is None or not hasattr(solution, "fluxes"):
        return solution
    order = [r.id for r in model.reactions]
    x = solution.fluxes.reindex(order).copy()
    x = x.mask(x.abs() < zero_cutoff, 0.0).round(decimals)
    solution.fluxes = x
    return solution

def deterministic_optimize(model, use_pfba: bool = True,
                           fraction_of_optimum: float = 1.0,
                           solver: Optional[str] = None):
    """
    Optimize and return a canonicalized Solution.

    IMPORTANT: use_pfba is a boolean. Call with keyword syntax:
        deterministic_optimize(model, use_pfba=True)
    Never pass "fba"/"pfba" as the positional second argument.
    """
    if not isinstance(use_pfba, bool):
        raise TypeError(
            "use_pfba must be bool. Passing the string 'fba' is truthy and silently "
            "runs pFBA in older helpers; this v2 helper rejects that mistake."
        )
    configure_solver_deterministic(model, solver=solver, verbose=False)

    biological_optimum = float(model.slim_optimize(error_value=None))

    if use_pfba:
        from cobra.flux_analysis import pfba
        solution = pfba(model, fraction_of_optimum=float(fraction_of_optimum))
    else:
        solution = model.optimize()

    try:
        solution.biological_objective_value = biological_optimum
        solution.secondary_objective_value = float(solution.objective_value)
    except Exception:
        pass
    return canonicalize_solution(solution, model)

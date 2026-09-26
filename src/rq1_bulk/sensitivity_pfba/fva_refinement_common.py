#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import math
import re
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis


def load_modeling_helpers(script_path: str | Path):
    path = Path(script_path)
    if not path.exists():
        raise FileNotFoundError(f"Modeling script not found: {path}")
    spec = importlib.util.spec_from_file_location("layered_modeling", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import modeling script: {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def nanmean_expression(df: pd.DataFrame, columns: Sequence[str]) -> np.ndarray:
    if not columns:
        raise ValueError("No replicate columns supplied for representative expression profile.")
    sub = df.loc[:, list(columns)].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    counts = np.isfinite(sub).sum(axis=1)
    sums = np.nansum(sub, axis=1)
    out = np.full(sub.shape[0], np.nan, dtype=float)
    np.divide(sums, counts, out=out, where=counts > 0)
    return out


def clean_fva_bounds(df: pd.DataFrame, tolerance: float = 1e-7) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Clean solver-scale FVA numerical noise without changing material intervals.

    Values with |bound| < tolerance are set to zero. Tiny min>max inversions whose
    gap is <= tolerance are collapsed to the midpoint. Material inversions are not
    silently repaired; they are recorded as QC failures and set to NaN.
    """
    out = df.copy()
    if "ReactionID" in out.columns:
        out = out.set_index("ReactionID")
    mn = pd.to_numeric(out["minimum"], errors="coerce")
    mx = pd.to_numeric(out["maximum"], errors="coerce")
    mn = mn.mask(mn.abs() < tolerance, 0.0)
    mx = mx.mask(mx.abs() < tolerance, 0.0)

    inv = mn.notna() & mx.notna() & (mn > mx)
    gap = mn - mx
    tiny = inv & (gap <= tolerance)
    material = inv & ~tiny

    mid = (mn + mx) / 2.0
    mid = mid.mask(mid.abs() < tolerance, 0.0)
    mn = mn.mask(tiny, mid)
    mx = mx.mask(tiny, mid)
    mn = mn.mask(material, np.nan)
    mx = mx.mask(material, np.nan)

    out["minimum"] = mn
    out["maximum"] = mx
    out["range"] = mx - mn
    qc = pd.DataFrame({
        "ReactionID": out.index.astype(str),
        "TinyInversionCollapsed": tiny.to_numpy(dtype=bool),
        "MaterialInversionExcluded": material.to_numpy(dtype=bool),
    })
    return out, qc


def run_standard_fva(model, reaction_ids: Sequence[str], fraction: float, processes: int = 1,
                     tolerance: float = 1e-7) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = flux_variability_analysis(
        model,
        reaction_list=list(reaction_ids),
        fraction_of_optimum=float(fraction),
        processes=max(1, int(processes)),
    )
    raw = raw.copy()
    raw.index.name = "ReactionID"
    return clean_fva_bounds(raw, tolerance=tolerance)


def _original_objective_optimum(model) -> float:
    value = model.slim_optimize(error_value=np.nan)
    if not np.isfinite(value):
        raise RuntimeError("Primary biological objective could not be optimized.")
    return float(value)


def run_parsimonious_fva(model, reaction_ids: Sequence[str], fraction: float = 1.0,
                         parsimony_epsilon: float = 0.01, processes: int = 1,
                         tolerance: float = 1e-7) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """FVA under both primary-objective and near-pFBA total-flux constraints.

    The model first retains `fraction` of the biological objective. Total absolute
    reaction flux (forward + reverse variables) is then minimized. FVA is finally
    performed while total flux is constrained to <= (1+epsilon) times that minimum.
    This asks how variable target reactions remain among solutions that are both
    objective-equivalent and approximately as parsimonious as pFBA.
    """
    eps = float(parsimony_epsilon)
    if eps < 0:
        raise ValueError("parsimony_epsilon must be >= 0")
    frac = float(fraction)
    if not (0 < frac <= 1.0):
        raise ValueError("fraction must satisfy 0 < fraction <= 1")

    with model as m:
        primary_expr = m.objective.expression
        primary_direction = str(m.objective.direction)
        primary_opt = _original_objective_optimum(m)

        if primary_direction.lower().startswith("max"):
            primary_bound = frac * primary_opt
            c_primary = m.problem.Constraint(
                primary_expr, lb=primary_bound, name="pfva_primary_objective_constraint"
            )
        else:
            # For a minimization objective, allow up to 1/fraction of the optimum.
            primary_bound = primary_opt / frac
            c_primary = m.problem.Constraint(
                primary_expr, ub=primary_bound, name="pfva_primary_objective_constraint"
            )
        m.add_cons_vars(c_primary)

        total_flux_expr = sum(
            (rxn.forward_variable + rxn.reverse_variable) for rxn in m.reactions
        )
        m.objective = m.problem.Objective(total_flux_expr, direction="min")
        min_total_flux = m.slim_optimize(error_value=np.nan)
        if not np.isfinite(min_total_flux):
            raise RuntimeError("Failed to obtain the pFBA total-flux minimum.")
        min_total_flux = float(min_total_flux)

        total_flux_upper = min_total_flux * (1.0 + eps) + float(tolerance)
        c_parsimony = m.problem.Constraint(
            total_flux_expr,
            ub=total_flux_upper,
            name=f"pfva_total_flux_eps_{str(eps).replace('.', '_')}",
        )
        m.add_cons_vars(c_parsimony)

        # Restore the biological objective before invoking COBRApy FVA. The explicit
        # primary-objective constraint above already preserves the required optimum;
        # fraction_of_optimum is repeated inside FVA as a conservative redundant check.
        m.objective = m.problem.Objective(primary_expr, direction=primary_direction)
        raw = flux_variability_analysis(
            m,
            reaction_list=list(reaction_ids),
            fraction_of_optimum=frac,
            processes=max(1, int(processes)),
        )
        raw = raw.copy()
        raw.index.name = "ReactionID"
        cleaned, qc = clean_fva_bounds(raw, tolerance=tolerance)
        meta = {
            "primary_optimum": primary_opt,
            "primary_fraction": frac,
            "pFBA_min_total_flux": min_total_flux,
            "parsimony_epsilon": eps,
            "total_flux_upper_bound": total_flux_upper,
        }
        return cleaned, qc, meta


def fva_difference(treatment: pd.DataFrame, control: pd.DataFrame,
                   tolerance: float = 1e-7) -> pd.DataFrame:
    idx = treatment.index.intersection(control.index)
    t = treatment.loc[idx]
    c = control.loc[idx]
    out = pd.DataFrame(index=idx)
    out.index.name = "ReactionID"
    out["TreatmentMin"] = t["minimum"]
    out["TreatmentMax"] = t["maximum"]
    out["ControlMin"] = c["minimum"]
    out["ControlMax"] = c["maximum"]
    out["DiffMin"] = out["TreatmentMin"] - out["ControlMax"]
    out["DiffMax"] = out["TreatmentMax"] - out["ControlMin"]
    out["RobustDirection"] = np.where(
        out["DiffMin"] > tolerance,
        "Up",
        np.where(out["DiffMax"] < -tolerance, "Down", "Ambiguous"),
    )
    out["Robust"] = out["RobustDirection"] != "Ambiguous"
    out["Tolerance"] = float(tolerance)
    return out


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")

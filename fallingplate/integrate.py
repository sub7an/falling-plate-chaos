"""integrate.py — thin wrapper over scipy.integrate.solve_ivp.

DOP853 (production/Lyapunov) or RK45 (coarse); tolerances from config, never
hard-coded. Provides event helpers with a grazing guard. See PLAN §2.
"""
from __future__ import annotations

from typing import Any, Callable

import numpy as np
from scipy.integrate import solve_ivp


def integrate_trajectory(rhs: Callable, s0: np.ndarray, t_span: tuple, cfg: dict) -> Any:
    """Integrate one trajectory with scipy.integrate.solve_ivp.

    `cfg` supplies the numerics (no hard-coded magic numbers; PLAN §2):
      - method: "DOP853" (production/Lyapunov) or "RK45" (coarse). Default DOP853.
      - rtol, atol: tolerances (required; fixed by the §4.3 convergence study).
      - dense_output: enable dense output for event/Δx(T) localization (default True).
      - t_eval, max_step, args: optional, forwarded to solve_ivp.

    Returns a scipy OdeResult. Physics enters only through `rhs`.
    """
    method = cfg.get("method", "DOP853")
    kwargs = dict(
        method=method,
        rtol=cfg["rtol"],
        atol=cfg["atol"],
        dense_output=cfg.get("dense_output", True),
    )
    if "t_eval" in cfg:
        kwargs["t_eval"] = cfg["t_eval"]
    if "max_step" in cfg:
        kwargs["max_step"] = cfg["max_step"]
    if "args" in cfg:
        kwargs["args"] = cfg["args"]
    return solve_ivp(rhs, t_span, np.asarray(s0, dtype=float), **kwargs)


def depth_event(y0: float, depth: float, graze_tol: float) -> Callable:
    """Event for the first transversal downward crossing of y = y0 - depth.

    direction = -1 (moving down). Crossings with |y_dot| < graze_tol are flagged
    as grazing and skipped by the caller, which reports the flagged fraction per
    depth (a non-negligible fraction means the depth is too shallow). PLAN §1, §2.
    """
    raise NotImplementedError("Scaffold only; see PLAN §1, §2.")


def poincare_event(section: dict) -> Callable:
    """Directional Poincare-section crossing event (default theta_dot = 0)."""
    raise NotImplementedError("Scaffold only; see PLAN §1, §2.")

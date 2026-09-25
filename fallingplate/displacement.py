"""displacement.py — landing displacement and sensitivity.

Fixed-depth measure (with grazing guard) plus a fixed-time cross-check. E5
sensitivity: displacement is an unbounded running sum, so expect an exponential
(Lyapunov) phase then diffusive ~sqrt(D) (chaotic) or ~linear (regular) growth,
NOT saturation; fit all three models and compare. See PLAN §1, §5-E5.

Trajectories are sampled arrays: t [N], y6 [6,N] = (x, y, vx', vy', theta, omega)
from model.rhs_full. Depth is measured downward from the release height y0.
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np

from fallingplate import model


class Landing(NamedTuple):
    dx: float          # horizontal displacement at the landing crossing (nan if none)
    grazing_flag: bool  # a grazing (sub-threshold |ydot|) crossing was skipped
    found: bool


def landing_displacement(t: np.ndarray, y6: np.ndarray, depth: float,
                         y0: float = 0.0, graze_tol: float = 1e-3) -> Landing:
    """Δx at the FIRST transversal downward crossing of y = y0 − depth.

    Only downward crossings (ydot < 0) count. A crossing whose |ydot| < graze_tol
    is flagged as grazing and skipped in favour of the next genuine one (PLAN §1).
    """
    x, y = y6[0], y6[1]
    vx, vy, th = y6[2], y6[3], y6[4]
    target = y0 - depth
    grazed = False
    for i in range(y.size - 1):
        if y[i] > target >= y[i + 1]:                 # y decreasing through target
            f = (y[i] - target) / (y[i] - y[i + 1])
            # ydot at the crossing (lab-frame), interpolated state
            vxc = vx[i] + f * (vx[i + 1] - vx[i])
            vyc = vy[i] + f * (vy[i + 1] - vy[i])
            thc = th[i] + f * (th[i + 1] - th[i])
            _, ydot = model.lab_velocity(vxc, vyc, thc)
            if ydot >= 0:
                continue                              # not a downward crossing
            if abs(ydot) < graze_tol:
                grazed = True
                continue                              # grazing: take the next
            dx = x[i] + f * (x[i + 1] - x[i])
            return Landing(float(dx), grazed, True)
    return Landing(float("nan"), grazed, False)


def landing_displacement_fixed_time(t: np.ndarray, y6: np.ndarray, t_land: float) -> float:
    """Cross-check: horizontal displacement x at fixed time t_land (interp)."""
    return float(np.interp(t_land, t, y6[0]))


def spread_vs_depth(dx_ensemble: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """sigma(D) = std of Δx across the IC ensemble, per depth. `dx_ensemble` is
    [n_members, n_depths]; `valid` the same-shape finite/landed mask. Depths with
    < 2 valid members give nan. PLAN §5-E5."""
    n_depths = dx_ensemble.shape[1]
    sigma = np.full(n_depths, np.nan)
    for j in range(n_depths):
        col = dx_ensemble[valid[:, j], j]
        if col.size >= 2:
            sigma[j] = np.std(col)
    return sigma


def _r2(y, yhat):
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0


def fit_growth_models(depths: np.ndarray, sigma: np.ndarray) -> dict:
    """Fit exponential (log σ = a + b·D), linear (σ = m·D) and sqrt (σ = k·√D) to
    σ(D); return coefficients, R² per model, the dominant model, and the
    exponential slope b (to compare with lambda_max/U from E2). PLAN §5-E5."""
    m = np.isfinite(sigma) & (sigma > 0) & np.isfinite(depths)
    D, s = np.asarray(depths)[m], np.asarray(sigma)[m]
    out: dict = {"n": int(D.size)}
    if D.size < 3:
        return out | {"dominant": None}

    # exponential: linear fit of log(sigma) vs D
    b, a = np.polyfit(D, np.log(s), 1)
    out["exp"] = {"slope": float(b), "intercept": float(a),
                  "r2": _r2(np.log(s), a + b * D)}
    # linear through origin: sigma = m D
    mlin = float(np.sum(D * s) / np.sum(D * D))
    out["linear"] = {"slope": mlin, "r2": _r2(s, mlin * D)}
    # sqrt through origin: sigma = k sqrt(D)
    ksq = float(np.sum(np.sqrt(D) * s) / np.sum(D))
    out["sqrt"] = {"coeff": ksq, "r2": _r2(s, ksq * np.sqrt(D))}
    out["dominant"] = max(("exp", "linear", "sqrt"), key=lambda k_: out[k_]["r2"])
    return out

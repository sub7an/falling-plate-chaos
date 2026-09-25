"""Validation ladder rung 3 — tolerance convergence (PLAN §4.3).

Checks that the integrated trajectory is Cauchy-convergent as rtol/atol tighten
and that DOP853 and RK45 agree at tight tolerance. Run in a REGULAR regime
(I* = 1.1 fluttering, spec §7) over a modest horizon so exponential separation
does not masquerade as non-convergence; convergence of lambda_max and the
observed-order study at the |w|w kink (spec §9) belong to the Lyapunov milestone.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from fallingplate import integrate, model

CONFIG = Path(__file__).resolve().parents[1] / "config" / "model_params.yaml"
IC = [0.0, 0.01, 0.5, 0.0]     # spec §7 tilted IC
T = 30.0


def _final_state(method: str, rtol: float, atol: float) -> np.ndarray:
    p = model.Params.from_yaml(str(CONFIG), I_star=1.1)
    sol = integrate.integrate_trajectory(
        lambda t, s: model.rhs(t, s, p),
        IC,
        (0.0, T),
        {"method": method, "rtol": rtol, "atol": atol, "t_eval": [T]},
    )
    return sol.y[:, -1]


def test_tolerance_cauchy_convergence():
    # Successive tolerance refinements must bring the final state closer together.
    tols = [(1e-6, 1e-8), (1e-8, 1e-10), (1e-10, 1e-12), (1e-12, 1e-14)]
    states = [_final_state("DOP853", rt, at) for rt, at in tols]
    gaps = [np.linalg.norm(states[i + 1] - states[i]) for i in range(len(states) - 1)]
    # Monotone shrinking gaps -> Cauchy sequence.
    assert gaps[0] > gaps[1] > gaps[2]
    # And the finest refinement has effectively converged.
    assert gaps[-1] < 1e-6


def test_integrator_cross_check():
    # Independent method (RK45) must agree with DOP853 at tight tolerance.
    s_dop = _final_state("DOP853", 1e-12, 1e-14)
    s_rk = _final_state("RK45", 1e-12, 1e-14)
    assert np.linalg.norm(s_dop - s_rk) < 1e-6

"""Validation ladder rung 5 — Lyapunov code verified on KNOWN systems (PLAN §4).

No plate physics here. The Lyapunov machinery in fallingplate.lyapunov must
reproduce textbook exponents on the logistic map and the Lorenz system BEFORE it
is ever trusted on the plate model (PLAN §2, §4 rung 5; M2 gate).

Targets:
  * Logistic x_{n+1} = r x_n (1-x_n):
      r=4.0  -> lambda = ln 2 = 0.6931  (fully chaotic)
      r=3.2  -> lambda = 0.5 ln|f'(x1) f'(x2)| = -0.9163  (stable period-2, <0)
      r=3.5  -> lambda < 0                                 (stable period-4)
  * Lorenz (sigma=10, rho=28, beta=8/3):
      lambda_max ~ 0.906  (variational AND independent Benettin agree)
      middle exponent ~ 0  (flow direction)
      sum of exponents = -(sigma + 1 + beta) = -13.6667  (constant divergence)
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from fallingplate import lyapunov as L

SIGMA, RHO, BETA = 10.0, 28.0, 8.0 / 3.0


# --------------------------- Logistic map ---------------------------

def _logistic(r):
    step = lambda x: r * x * (1.0 - x)
    jac = lambda x: np.array([[r * (1.0 - 2.0 * x[0])]])
    return step, jac


def _logistic_lambda(r):
    step, jac = _logistic(r)
    spec = L.map_lyapunov_spectrum(
        step, jac, [0.1234567], {"n_iter": 200_000, "n_transient": 5_000}
    )
    return spec[0]


def test_logistic_chaotic_r4():
    # Fully chaotic: lambda = ln 2 exactly.
    assert _logistic_lambda(4.0) == pytest.approx(math.log(2.0), abs=1e-3)


def test_logistic_period2_r32():
    # Stable period-2 window: negative, equal to the analytic period-2 value.
    r = 3.2
    x1 = (r + 1 + math.sqrt((r + 1) * (r - 3))) / (2 * r)
    x2 = (r + 1 - math.sqrt((r + 1) * (r - 3))) / (2 * r)
    analytic = 0.5 * math.log(abs(r * (1 - 2 * x1) * r * (1 - 2 * x2)))
    lam = _logistic_lambda(r)
    assert lam < 0
    assert lam == pytest.approx(analytic, abs=1e-2)   # analytic ~ -0.9163


def test_logistic_period4_r35():
    # Stable period-4 window: negative.
    assert _logistic_lambda(3.5) < 0


# --------------------------- Lorenz system ---------------------------

def _lorenz_rhs(t, s):
    x, y, z = s
    return [SIGMA * (y - x), x * (RHO - z) - y, x * y - BETA * z]


def _lorenz_jac(t, s):
    x, y, z = s
    return np.array(
        [[-SIGMA, SIGMA, 0.0], [RHO - z, -1.0, -x], [y, x, -BETA]]
    )


_LORENZ_CFG = dict(
    t_total=1000.0, dt_renorm=0.5, t_transient=50.0,
    method="DOP853", rtol=1e-9, atol=1e-9,
)


@pytest.fixture(scope="module")
def lorenz_spectrum():
    return L.variational_spectrum(_lorenz_rhs, _lorenz_jac, [1.0, 1.0, 1.0], _LORENZ_CFG)


def test_lorenz_variational_lambda_max(lorenz_spectrum):
    assert lorenz_spectrum[0] == pytest.approx(0.906, abs=1e-2)


def test_lorenz_middle_exponent_is_zero(lorenz_spectrum):
    # A bounded flow has a zero exponent along the trajectory direction.
    assert abs(lorenz_spectrum[1]) < 1e-2


def test_lorenz_spectrum_sum_is_divergence(lorenz_spectrum):
    # Constant phase-space contraction: sum lambda_i = trace(J) = -(sigma+1+beta).
    assert lorenz_spectrum.sum() == pytest.approx(-(SIGMA + 1.0 + BETA), abs=1e-2)


def test_lorenz_benettin_agrees_with_variational(lorenz_spectrum):
    # Independent (Jacobian-free) cross-check of lambda_max.
    lmax_b = L.benettin_lmax(_lorenz_rhs, [1.0, 1.0, 1.0], 1e-8, dict(_LORENZ_CFG, seed=1))
    assert lmax_b == pytest.approx(0.906, abs=2e-2)
    assert lmax_b == pytest.approx(lorenz_spectrum[0], abs=2e-2)

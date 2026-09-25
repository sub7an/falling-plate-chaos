"""Unit tests for the E5 displacement machinery (M7): lab-frame quadrature,
landing displacement with grazing guard, and growth-model fitting. See PLAN §5-E5.
"""
from __future__ import annotations

import math

import numpy as np

from fallingplate import displacement, model


def test_lab_velocity_rotates_body_to_lab():
    # theta=0: identity.
    assert model.lab_velocity(2.0, 1.0, 0.0) == (2.0, 1.0)
    # theta=pi/2: (vx,vy) -> (-vy, vx).
    dx, dy = model.lab_velocity(2.0, 1.0, math.pi / 2)
    assert dx == abs(dx) * np.sign(dx) and abs(dx + 1.0) < 1e-12 and abs(dy - 2.0) < 1e-12


def test_rhs_full_prepends_position_quadrature():
    p = model.Params(C_T=1.2, C_R=math.pi, A=1.4, B=1.0, mu1=0.2, mu2=0.2, I_star=2.2)
    s6 = np.array([3.0, -4.0, 0.3, -0.7, 1.1, 0.4])   # x,y,vx,vy,theta,omega
    d = model.rhs_full(0.0, s6, p)
    assert d.shape == (6,)
    assert (d[0], d[1]) == model.lab_velocity(0.3, -0.7, 1.1)
    assert np.allclose(d[2:], model.rhs(0.0, np.array([0.3, -0.7, 1.1, 0.4]), p))


def test_landing_displacement_downward_crossing():
    t = np.array([0.0, 1.0, 2.0])
    # y decreases 0 -> -3 -> -6; x = 0,5,10; theta=0 so ydot = vy = -3 (downward).
    y6 = np.array([[0.0, 5.0, 10.0],      # x
                   [0.0, -3.0, -6.0],     # y
                   [0.0, 0.0, 0.0],       # vx
                   [-3.0, -3.0, -3.0],    # vy
                   [0.0, 0.0, 0.0],       # theta
                   [0.0, 0.0, 0.0]])      # omega
    L = displacement.landing_displacement(t, y6, depth=4.0, graze_tol=1e-3)
    assert L.found and not L.grazing_flag
    assert abs(L.dx - (5.0 + (1.0 / 3.0) * 5.0)) < 1e-9   # interp at target y=-4


def test_landing_displacement_grazing_is_skipped():
    t = np.array([0.0, 1.0, 2.0])
    # crosses -4 but |ydot| ~ 1e-5 < graze_tol -> flagged and skipped; no other crossing.
    y6 = np.array([[0.0, 5.0, 10.0],
                   [0.0, -3.0, -6.0],
                   [0.0, 0.0, 0.0],
                   [-1e-5, -1e-5, -1e-5],   # ydot ~ -1e-5 (grazing)
                   [0.0, 0.0, 0.0],
                   [0.0, 0.0, 0.0]])
    L = displacement.landing_displacement(t, y6, depth=4.0, graze_tol=1e-3)
    assert L.grazing_flag and not L.found


def test_fit_growth_models_identifies_dominant():
    D = np.arange(1.0, 25.0)
    assert displacement.fit_growth_models(D, 1e-3 * np.exp(0.3 * D))["dominant"] == "exp"
    assert displacement.fit_growth_models(D, 0.5 * D)["dominant"] == "linear"
    assert displacement.fit_growth_models(D, 0.7 * np.sqrt(D))["dominant"] == "sqrt"


def test_spread_vs_depth_needs_two_valid():
    dx = np.array([[1.0, np.nan], [3.0, np.nan]])
    valid = np.isfinite(dx)
    sig = displacement.spread_vs_depth(dx, valid)
    assert sig[0] == np.std([1.0, 3.0]) and math.isnan(sig[1])

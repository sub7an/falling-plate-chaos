"""Validation ladder rungs 1-2 for fallingplate.model (PLAN §4).

Every assertion traces to docs/model_spec.md: the governing equations
(5.1)-(5.3), closures (4.7)-(4.9), fixed points (§6), and the numeric anchors
(§7). No physics value is invented here; anchors are quoted from the spec.

Rung 1 — unit tests: RHS shape/finiteness, closure signs & zero-speed limits,
         translation/time/angle invariance, analytic-vs-FD Jacobian.
Rung 2 — limiting cases: zero-coefficient ballistic free-fall (analytic),
         fixed-point residuals, edge-on saddle eigenvalues, flutter/tumble
         either side of the heteroclinic I*_C.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from fallingplate import integrate, model

CONFIG = Path(__file__).resolve().parents[1] / "config" / "model_params.yaml"


@pytest.fixture(scope="module")
def p() -> model.Params:
    """Spec parameters (C_T, C_R, A, B, mu1, mu2) at the heteroclinic I*_C."""
    return model.Params.from_yaml(str(CONFIG))


def _zero_coeffs(I_star: float) -> model.Params:
    """All fluid coefficients set to zero: ballistic limit (PLAN §4 rung 2)."""
    return model.Params(C_T=0.0, C_R=0.0, A=0.0, B=0.0, mu1=0.0, mu2=0.0, I_star=I_star)


# ---------------------------------------------------------------------------
# Rung 1 — unit tests
# ---------------------------------------------------------------------------

def test_params_from_yaml_matches_spec(p):
    # Spec §5 table.
    assert p.C_T == 1.2
    assert p.C_R == pytest.approx(math.pi)
    assert p.A == 1.4
    assert p.B == 1.0
    assert p.mu1 == 0.2
    assert p.mu2 == 0.2
    # I_star override leaves coefficients untouched.
    q = model.Params.from_yaml(str(CONFIG), I_star=2.2)
    assert q.I_star == 2.2
    assert (q.C_T, q.A, q.B) == (p.C_T, p.A, p.B)


def test_rhs_shape_and_finite(p):
    d = model.rhs(0.0, np.array([0.3, -0.7, 1.1, 0.4]), p)
    assert d.shape == (4,)
    assert np.all(np.isfinite(d))


def test_dtheta_equals_omega(p):
    for omega in (-1.3, 0.0, 2.5):
        d = model.rhs(0.0, np.array([0.2, 0.5, 0.9, omega]), p)
        assert d[2] == omega  # dtheta/dt = omega, exactly (spec §4)


def test_translation_and_time_invariance(p):
    # RHS carries no x, y and no explicit t (PLAN A1): value must not depend on t.
    s = np.array([0.4, -0.2, 0.7, 0.3])
    assert np.allclose(model.rhs(0.0, s, p), model.rhs(123.4, s, p))


def test_theta_is_cyclic(p):
    # theta is the only periodic coordinate (period 2*pi, spec §2).
    s0 = np.array([0.4, -0.2, 0.7, 0.3])
    s1 = np.array([0.4, -0.2, 0.7 + 2 * math.pi, 0.3])
    assert np.allclose(model.rhs(0.0, s0, p), model.rhs(0.0, s1, p))


def test_torque_is_dissipative_and_odd(p):
    # tau = (mu1 + mu2|w|) w : same sign as w (opposes rotation via -tau in 5.3),
    # and odd in w.
    for omega in (0.1, 1.0, 3.0):
        assert model.torque(omega, p) > 0
        assert model.torque(-omega, p) == pytest.approx(-model.torque(omega, p))
    assert model.torque(0.0, p) == 0.0


def test_fluid_force_opposes_translation(p):
    # Pure edge-on motion (vy'=0): force is along vx' and opposes it.
    fx, fy = model.fluid_force(0.9, 0.0, p)
    assert fx > 0 and fy == 0.0          # +vx' -> +Fx' -> -Fx' decelerates in (5.1)
    fx2, _ = model.fluid_force(-0.9, 0.0, p)
    assert fx2 < 0
    # Pure broadside motion (vx'=0): force is along vy' and opposes it.
    gx, gy = model.fluid_force(0.0, 0.8, p)
    assert gy > 0 and gx == 0.0


def test_zero_speed_limits_are_finite(p):
    # Spec §9: velocity-drag / fluid-force terms -> 0 at s = 0; only C_R*w survives.
    assert model.fluid_force(0.0, 0.0, p) == (0.0, 0.0)
    assert model.gamma(0.0, 0.0, 0.0, p) == 0.0
    assert model.gamma(0.0, 0.0, 1.5, p) == pytest.approx((2 / math.pi) * p.C_R * 1.5)
    d = model.rhs(0.0, np.array([0.0, 0.0, 0.3, 0.7]), p)
    assert np.all(np.isfinite(d))


def test_jacobian_matches_finite_difference(p):
    rng = np.random.default_rng(20260922)
    for _ in range(25):
        # Random states with clearly non-zero speed (Jacobian is C^1 only there).
        vx, vy = rng.uniform(-2, 2, size=2)
        if math.hypot(vx, vy) < 0.3:
            vx += 0.5
        theta = rng.uniform(-math.pi, math.pi)
        omega = rng.uniform(-2, 2)
        s = np.array([vx, vy, theta, omega])

        J = model.jacobian(0.0, s, p)
        Jfd = np.empty((4, 4))
        for k in range(4):
            h = 1e-6 * max(1.0, abs(s[k]))
            sp, sm = s.copy(), s.copy()
            sp[k] += h
            sm[k] -= h
            Jfd[:, k] = (model.rhs(0.0, sp, p) - model.rhs(0.0, sm, p)) / (2 * h)
        assert np.allclose(J, Jfd, rtol=1e-5, atol=1e-7)


# ---------------------------------------------------------------------------
# Rung 2 — limiting cases and spec anchors
# ---------------------------------------------------------------------------

def test_ballistic_free_fall_broadside():
    # Zero fluid coefficients, aligned broadside IC (theta=0, vx'=0, w=0):
    # the card must not rotate and vy'(t) = vy0 - t/(I*+1) exactly (added-mass
    # inertia I*+1 from 5.2). Lab-frame fall is straight and vertical.
    #
    # RESTRICTION: the closed form holds ONLY for this decoupled IC. The -vx'vy'
    # term in (5.3) drives omega (and hence theta) whenever vx'vy' != 0, so a
    # general IC does NOT free-fall in a straight line even with zero fluid
    # coefficients. Here vx'=0 keeps vx'vy'=0, so omega and theta stay pinned.
    I_star, vy0, T = 1.7, 0.5, 3.0
    pz = _zero_coeffs(I_star)
    sol = integrate.integrate_trajectory(
        lambda t, s: model.rhs(t, s, pz),
        [0.0, vy0, 0.0, 0.0],
        (0.0, T),
        {"method": "DOP853", "rtol": 1e-11, "atol": 1e-13, "t_eval": [T]},
    )
    vx, vy, theta, omega = sol.y[:, -1]
    assert abs(vx) < 1e-10 and abs(theta) < 1e-10 and abs(omega) < 1e-10
    assert vy == pytest.approx(vy0 - T / (I_star + 1.0), abs=1e-9)


def test_ballistic_free_fall_edge_on():
    # Aligned edge-on IC (theta=pi/2, vy'=0, w=0): vx'(t) = vx0 - t/I* (5.1).
    # RESTRICTION (as above): the closed form holds only for this decoupled IC.
    # vy'=0 keeps the -vx'vy' driving term in (5.3) zero, so omega/theta stay
    # pinned; a general IC couples theta and omega and does not free-fall straight.
    I_star, vx0, T = 1.7, 0.5, 3.0
    pz = _zero_coeffs(I_star)
    sol = integrate.integrate_trajectory(
        lambda t, s: model.rhs(t, s, pz),
        [vx0, 0.0, math.pi / 2, 0.0],
        (0.0, T),
        {"method": "DOP853", "rtol": 1e-11, "atol": 1e-13, "t_eval": [T]},
    )
    vx, vy, theta, omega = sol.y[:, -1]
    assert abs(vy) < 1e-10 and abs(omega) < 1e-10
    assert theta == pytest.approx(math.pi / 2, abs=1e-12)
    assert vx == pytest.approx(vx0 - T / I_star, abs=1e-9)


def test_fixed_point_residuals(p):
    # Spec §6 / (6.1),(6.2): edge-on V and broadside W, with the paper's values.
    V = math.sqrt(math.pi / (p.A - p.B))
    W = math.sqrt(math.pi / (p.A + p.B))
    assert V == pytest.approx(2.8025, abs=1e-3)   # spec §6
    assert W == pytest.approx(1.1441, abs=1e-3)   # spec §6

    fixed_points = [
        (-V, 0.0, math.pi / 2, 0.0),        # edge-on   (verified in spec §7)
        (+V, 0.0, 3 * math.pi / 2, 0.0),    # edge-on   (paired)
        (0.0, -W, 0.0, 0.0),                # broadside (verified in spec §7)
        (0.0, +W, math.pi, 0.0),            # broadside (paired)
    ]
    for fp in fixed_points:
        res = model.rhs(0.0, np.array(fp), p)
        assert np.max(np.abs(res)) < 1e-13, fp


def test_edge_on_saddle_eigenvalues():
    # Spec §7 anchor 4, at the heteroclinic I*_C = 1.2191.
    p = model.Params.from_yaml(str(CONFIG), I_star=1.2191)
    V = math.sqrt(math.pi / (p.A - p.B))
    J = model.jacobian(0.0, np.array([-V, 0.0, math.pi / 2, 0.0]), p)
    eig = np.linalg.eigvals(J)

    # Analytic decoupled stable eigenvalue: lambda_s = -2 sqrt((A-B)/pi) / I*.
    lam_s_analytic = -2.0 * math.sqrt((p.A - p.B) / math.pi) / p.I_star
    assert lam_s_analytic == pytest.approx(-0.58539, abs=1e-4)

    real = np.sort(eig[np.abs(eig.imag) < 1e-9].real)
    assert real[0] == pytest.approx(-0.5854, abs=2e-4)   # lambda_s
    assert real[-1] == pytest.approx(+0.3813, abs=2e-4)  # lambda_u (unstable)

    complex_pair = eig[np.abs(eig.imag) >= 1e-9]
    assert complex_pair.real[0] == pytest.approx(-0.9861, abs=2e-4)
    assert np.max(np.abs(complex_pair.imag)) == pytest.approx(2.5949, abs=2e-3)


def _net_rotation(I_star: float, T: float = 500.0) -> float:
    """Net unwrapped rotation over [100, T] for the spec's flutter/tumble IC."""
    p = model.Params.from_yaml(str(CONFIG), I_star=I_star)
    t_eval = np.linspace(0.0, T, int(T / 0.05) + 1)
    sol = integrate.integrate_trajectory(
        lambda t, s: model.rhs(t, s, p),
        [0.0, 0.01, 0.5, 0.0],            # spec §7: tilted IC used for the checks
        (0.0, T),
        {"method": "DOP853", "rtol": 1e-11, "atol": 1e-13, "t_eval": t_eval},
    )
    theta_unwrapped = np.unwrap(sol.y[2])
    post = t_eval >= 100.0                 # discard the edge-on transient
    return theta_unwrapped[post][-1] - theta_unwrapped[post][0]


def test_flutter_tumble_across_heteroclinic():
    # Spec §7 anchor 3: I* = 1.2190 flutters (bounded, ~zero net rotation);
    # I* = 1.2192 tumbles (secular rotation, several revolutions).
    net_flutter = _net_rotation(1.2190)
    net_tumble = _net_rotation(1.2192)
    assert abs(net_flutter) < 2 * math.pi        # flutter: no net revolutions
    assert abs(net_tumble) > 4 * math.pi         # tumble: monotone winding

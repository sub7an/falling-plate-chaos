"""lyapunov.py — Lyapunov exponents.

PRIMARY: variational (tangent-linear) spectrum for a continuous flow, integrated
as one augmented system with the state (shared adaptive steps) + periodic QR
reorthonormalization.
CROSS-CHECK: Benettin lambda_max as an augmented reference+perturbed system
(shared steps, no separately-integrated orbits), independent of the Jacobian.
DISCRETE: map_lyapunov_spectrum for iterated maps (logistic verification, and
later Poincare return maps).

The reduced-state metric treats theta periodically (`angular_metric`).

Verify on logistic + Lorenz FIRST (PLAN §4 rung 5). No physics here — all model
knowledge enters through the caller-supplied `rhs`/`jac`. See PLAN §2.
"""
from __future__ import annotations

from typing import Callable, Optional

import numpy as np
from scipy.integrate import solve_ivp


def _fd_jacobian(rhs: Callable, t: float, y: np.ndarray, eps: float = 1e-7) -> np.ndarray:
    """Central finite-difference Jacobian of `rhs` at (t, y). Fallback only."""
    n = y.size
    J = np.empty((n, n))
    f = np.asarray(rhs(t, y), dtype=float)
    for k in range(n):
        h = eps * max(1.0, abs(y[k]))
        yp, ym = y.copy(), y.copy()
        yp[k] += h
        ym[k] -= h
        J[:, k] = (np.asarray(rhs(t, yp)) - np.asarray(rhs(t, ym))) / (2 * h)
    return J


def variational_spectrum(
    rhs: Callable,
    jac: Optional[Callable],
    s0: np.ndarray,
    cfg: dict,
) -> np.ndarray:
    """Full Lyapunov spectrum of a continuous flow via tangent-linear + QR.

    Primary method. The state and a set of `k` tangent vectors Q are advanced by
    ONE augmented ODE (shared adaptive steps): d(state)/dt = rhs, dQ/dt = J·Q.
    After each segment of length `dt_renorm`, Q is reorthonormalized by QR and
    ln|diag(R)| accumulated. lambda_i = mean growth rate / time.

    cfg keys: t_total, dt_renorm, t_transient (default 0), method (DOP853),
    rtol, atol, k (# exponents, default state dim). If `jac` is None a
    finite-difference Jacobian is used. Returns exponents sorted descending.
    """
    s0 = np.asarray(s0, dtype=float)
    n = s0.size
    k = int(cfg.get("k", n))
    dt = float(cfg["dt_renorm"])
    n_seg = int(round(float(cfg["t_total"]) / dt))
    n_tr = int(round(float(cfg.get("t_transient", 0.0)) / dt))
    method = cfg.get("method", "DOP853")
    rtol, atol = cfg["rtol"], cfg["atol"]
    jfun = jac if jac is not None else (lambda t, y: _fd_jacobian(rhs, t, y))

    def augmented(t, Y):
        state = Y[:n]
        Q = Y[n:].reshape(n, k)
        dstate = np.asarray(rhs(t, state), dtype=float)
        J = np.asarray(jfun(t, state), dtype=float)
        return np.concatenate([dstate, (J @ Q).ravel()])

    state = s0.copy()
    Q = np.eye(n, k)  # orthonormal columns (standard basis)
    sum_log = np.zeros(k)
    count = 0
    t0 = 0.0
    for i in range(n_seg):
        Y0 = np.concatenate([state, Q.ravel()])
        sol = solve_ivp(
            augmented, (t0, t0 + dt), Y0,
            method=method, rtol=rtol, atol=atol, dense_output=False,
        )
        Yf = sol.y[:, -1]
        state = Yf[:n]
        Z = Yf[n:].reshape(n, k)
        Q, R = np.linalg.qr(Z)
        if i >= n_tr:
            sum_log += np.log(np.abs(np.diag(R)))
            count += 1
        t0 += dt

    return np.sort(sum_log / (count * dt))[::-1]


def benettin_lmax(rhs: Callable, s0: np.ndarray, delta0: float, cfg: dict) -> float:
    """lambda_max of a continuous flow via Benettin (augmented, shared steps).

    Independent of the Jacobian: a reference and a perturbed trajectory are
    integrated together in one augmented ODE (identical adaptive steps), and the
    separation is periodically measured and renormalized to `delta0`. This is
    the standard Benettin renormalization method and cross-checks the variational
    lambda_max.

    cfg keys: t_total, dt_renorm, t_transient (default 0), method, rtol, atol,
    seed (default 0) for the initial perturbation direction.
    """
    s0 = np.asarray(s0, dtype=float)
    n = s0.size
    dt = float(cfg["dt_renorm"])
    n_seg = int(round(float(cfg["t_total"]) / dt))
    n_tr = int(round(float(cfg.get("t_transient", 0.0)) / dt))
    method = cfg.get("method", "DOP853")
    rtol, atol = cfg["rtol"], cfg["atol"]

    rng = np.random.default_rng(cfg.get("seed", 0))
    e = rng.standard_normal(n)
    e /= np.linalg.norm(e)
    ref = s0.copy()
    pert = ref + delta0 * e

    def augmented(t, Y):
        r = Y[:n]
        p = Y[n:]
        return np.concatenate([np.asarray(rhs(t, r)), np.asarray(rhs(t, p))])

    sum_log = 0.0
    count = 0
    t0 = 0.0
    for i in range(n_seg):
        Y0 = np.concatenate([ref, pert])
        sol = solve_ivp(
            augmented, (t0, t0 + dt), Y0,
            method=method, rtol=rtol, atol=atol, dense_output=False,
        )
        Yf = sol.y[:, -1]
        ref = Yf[:n]
        pert = Yf[n:]
        d = pert - ref
        dist = np.linalg.norm(d)
        if i >= n_tr:
            sum_log += np.log(dist / delta0)
            count += 1
        pert = ref + (delta0 / dist) * d  # renormalize back to delta0
        t0 += dt

    return sum_log / (count * dt)


def map_lyapunov_spectrum(
    step: Callable,
    jac: Callable,
    x0: np.ndarray,
    cfg: dict,
) -> np.ndarray:
    """Lyapunov spectrum of an iterated map x -> step(x) via QR.

    Discards `n_transient` iterations to land on the attractor, then accumulates
    ln|diag(R)| of the QR factorization of J·Q over `n_iter` iterations.
    lambda_i = mean growth rate per iteration. For a 1-D map this reduces to
    the standard mean of ln|f'(x)|. Returns exponents sorted descending.

    cfg keys: n_iter, n_transient (default 0), k (# exponents, default dim).
    """
    x = np.atleast_1d(np.asarray(x0, dtype=float))
    n = x.size
    k = int(cfg.get("k", n))
    n_iter = int(cfg["n_iter"])
    n_tr = int(cfg.get("n_transient", 0))

    for _ in range(n_tr):
        x = np.atleast_1d(np.asarray(step(x), dtype=float))

    Q = np.eye(n, k)
    sum_log = np.zeros(k)
    for _ in range(n_iter):
        J = np.asarray(jac(x), dtype=float).reshape(n, n)
        Q, R = np.linalg.qr(J @ Q)
        sum_log += np.log(np.abs(np.diag(R)))
        x = np.atleast_1d(np.asarray(step(x), dtype=float))

    return np.sort(sum_log / n_iter)[::-1]


def angular_metric(u_a: np.ndarray, u_b: np.ndarray, angle_idx: int) -> float:
    """Reduced-state distance with the theta component wrapped to (-pi, pi]."""
    u_a = np.asarray(u_a, dtype=float)
    u_b = np.asarray(u_b, dtype=float)
    d = u_a - u_b
    d[angle_idx] = (d[angle_idx] + np.pi) % (2 * np.pi) - np.pi
    return float(np.linalg.norm(d))

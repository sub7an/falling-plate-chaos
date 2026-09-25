"""regimes.py — regime classification and Poincare maps.

Classify steady / fluttering / tumbling (period-1, period-2, ...) / chaotic from
a post-transient trajectory using the mean rotation rate Omega, omega sign
changes, a Poincare return-point count, and lambda_max (from lyapunov.py). No
physics here — this operates on an already-integrated trajectory. See PLAN §1.

State convention (spec §2): u = (vx', vy', theta, omega). We integrate dtheta/dt
= omega WITHOUT wrapping, so `theta` arrives already unwrapped and its secular
growth is the net rotation.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np


class Regime(Enum):
    STEADY = "steady"
    FLUTTERING = "fluttering"
    TUMBLING = "tumbling"
    MIXTURE = "periodic mixture of fluttering and tumbling"
    CHAOTIC = "chaotic"
    UNKNOWN = "unknown"


@dataclass
class RegimeReport:
    regime: Regime
    label: str          # human-readable, includes period / sub-type
    Omega: float        # mean rotation rate <omega>
    net_rev: float      # net revolutions over the window
    gross_rev: float    # integral of |omega| / 2pi over the window
    directedness: float # |net_rev| / gross_rev in [0,1]: ~1 tumbling, ~0 fluttering
    theta_amp: float    # peak-to-peak of theta (full swing; ~net_rev*2pi if tumbling)
    n_sign_changes: int # omega zero-crossings
    poincare_period: int  # tumbling return-point count (0 if section not crossed)
    lam_max: float | None


def mean_rotation_rate(t: np.ndarray, omega: np.ndarray) -> float:
    """Omega = <omega> over the window (time-averaged)."""
    return float(np.trapezoid(omega, t) / (t[-1] - t[0]))


def _sign_changes(omega: np.ndarray) -> int:
    s = np.sign(omega)
    s = s[s != 0]
    if s.size < 2:
        return 0
    return int(np.count_nonzero(np.diff(s) != 0))


def poincare_map(t: np.ndarray, y: np.ndarray, phase: float = np.pi):
    """Tumbling Poincare section: reduced state (vx', vy', omega) at upward
    crossings of theta = phase + 2*pi*k. `theta` is already unwrapped, so a
    fluttering orbit stays in one 2*pi band and yields no crossings (empty).

    Section choice logged in docs/decisions.md: the spec's default omega=0
    section is not crossed by a one-signed tumbling orbit, so we section on the
    cyclic angle instead. Returns (states [m,3], crossing_times [m]).
    """
    vx, vy, theta, omega = y
    m = (theta - phase) / (2.0 * np.pi)
    states, times = [], []
    for i in range(len(m) - 1):
        # integers of m crossed on this step, either rotation direction
        lo, hi = (m[i], m[i + 1]) if m[i] <= m[i + 1] else (m[i + 1], m[i])
        for k in range(int(np.ceil(lo)), int(np.floor(hi)) + 1):
            denom = m[i + 1] - m[i]
            if denom == 0.0:
                continue
            f = (k - m[i]) / denom
            states.append([vx[i] + f * (vx[i + 1] - vx[i]),
                           vy[i] + f * (vy[i + 1] - vy[i]),
                           omega[i] + f * (omega[i + 1] - omega[i])])
            times.append(t[i] + f * (t[i + 1] - t[i]))
    return np.asarray(states), np.asarray(times)


def poincare_period(states: np.ndarray, tol_rel: float = 1e-3, keep: int = 24) -> int:
    """Return-point count from the section-state sequence: 1, 2, 4, or -1
    (higher / aperiodic). 0 if the section was not crossed (not tumbling)."""
    v = np.asarray(states)
    if v.ndim != 2 or v.shape[0] < 4:
        return 0
    v = v[-keep:]
    scale = max(float(np.mean(np.linalg.norm(v, axis=1))), 1e-12)
    tol = tol_rel * scale
    if np.max(np.ptp(v, axis=0)) < tol:
        return 1
    for period in (2, 4):
        branches = [v[j::period] for j in range(period)]
        if all(np.max(np.ptp(b, axis=0)) < tol for b in branches if len(b) > 1):
            centres = np.array([b.mean(axis=0) for b in branches])
            if np.max(np.ptp(centres, axis=0)) > tol:
                return period
    return -1


def classify_regime(
    t: np.ndarray,
    y: np.ndarray,
    lam_max: float | None = None,
    *,
    lam_pos_thresh: float = 2e-2,
    tumbling_directedness: float = 0.8,
    flutter_directedness: float = 0.2,
    min_tumble_rev: float = 3.0,
    small_amp_thresh: float = np.pi / 2,
) -> RegimeReport:
    """Label a post-transient trajectory y = (vx', vy', theta, omega) over t.

    Discriminator (PLAN §1) built on rotation *directedness*
    R = |<omega>| / <|omega|> (fraction of rotation that does not cancel):
      * lambda_max > 0                       -> chaotic
      * R >= tumbling_directedness (persistent one-way spin, many revolutions)
                                             -> tumbling; period from Poincare count
      * R <= flutter_directedness (rotation cancels; bounded swing)
                                             -> fluttering (small-amplitude flagged)
      * in between (both tumble and flutter episodes)
                                             -> periodic mixture
    Thresholds are arguments, not magic numbers.
    """
    vx, vy, theta, omega = y
    Omega = mean_rotation_rate(t, omega)
    net_rev = float((theta[-1] - theta[0]) / (2.0 * np.pi))
    gross = float(np.trapezoid(np.abs(omega), t) / (2.0 * np.pi))
    R = abs(net_rev) / gross if gross > 1e-9 else 0.0
    theta_amp = float(np.ptp(theta))
    n_sc = _sign_changes(omega)

    states, _ = poincare_map(t, y)
    period = poincare_period(states)

    if lam_max is not None and lam_max > lam_pos_thresh:
        sub = "tumbling" if R >= tumbling_directedness else (
            "fluttering" if R <= flutter_directedness else "mixed")
        return RegimeReport(Regime.CHAOTIC, f"chaotic ({sub})", Omega, net_rev,
                            gross, R, theta_amp, n_sc, period, lam_max)

    if R >= tumbling_directedness and abs(net_rev) >= min_tumble_rev:
        pname = {1: "period-one", 2: "period-two", 4: "period-four"}.get(period, "periodic")
        if period == -1:
            pname = "aperiodic"
        return RegimeReport(Regime.TUMBLING, f"{pname} tumbling", Omega, net_rev,
                            gross, R, theta_amp, n_sc, period, lam_max)

    # Low directedness: rotation cancels. Pure fluttering never completes a full
    # turn (theta_amp < 2*pi); if theta swings past a full revolution yet nets to
    # ~zero, it is doing tumble flips that reverse -> flutter/tumble mixture.
    if R <= flutter_directedness and theta_amp < 2.0 * np.pi:
        amp = "small-amplitude " if theta_amp < small_amp_thresh else ""
        return RegimeReport(Regime.FLUTTERING, f"{amp}fluttering", Omega, net_rev,
                            gross, R, theta_amp, n_sc, period, lam_max)

    return RegimeReport(Regime.MIXTURE, "periodic mixture of fluttering and tumbling",
                        Omega, net_rev, gross, R, theta_amp, n_sc, period, lam_max)

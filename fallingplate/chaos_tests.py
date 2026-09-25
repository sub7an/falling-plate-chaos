"""chaos_tests.py — independent chaos diagnostics.

0-1 test for chaos (Gottwald-Melbourne): sample the observable stroboscopically
on the Poincare section or at ~ the dominant period (never the fine integrator
step, which biases K high); report the MEDIAN K over many random c values.
K ~ 0 regular, K ~ 1 chaotic. See PLAN §2.
"""
from __future__ import annotations

import numpy as np


def _K_for_c(phi: np.ndarray, c: float) -> float:
    """Correlation-method K for one frequency c (Gottwald-Melbourne, with the
    Vosc mean-square-displacement correction)."""
    n = phi.size
    j = np.arange(1, n + 1)
    p = np.cumsum(phi * np.cos(j * c))
    q = np.cumsum(phi * np.sin(j * c))
    ncut = n // 10
    nn = np.arange(1, ncut + 1)
    phi_mean = phi.mean()
    # mean-square displacement of (p, q), averaged over start points
    M = np.empty(ncut)
    for k in nn:
        dp = p[k:] - p[:-k]
        dq = q[k:] - q[:-k]
        M[k - 1] = np.mean(dp * dp + dq * dq)
    # Vosc correction removes the oscillatory ballistic part
    Vosc = phi_mean ** 2 * (1.0 - np.cos(nn * c)) / (1.0 - np.cos(c))
    D = M - Vosc
    # correlation between n and D(n)
    dn = nn - nn.mean()
    dD = D - D.mean()
    denom = np.sqrt(np.sum(dn * dn) * np.sum(dD * dD))
    return float(np.sum(dn * dD) / denom) if denom > 0 else 0.0


def zero_one_test(series: np.ndarray, n_c: int = 100, seed: int = 0) -> float:
    """Median K over n_c random frequencies c in (pi/5, 4pi/5).

    Input MUST already be sub-sampled stroboscopically (Poincare section or ~
    dominant period) — never the fine integrator step, which biases K high.
    K ~ 0 regular, K ~ 1 chaotic. PLAN §2.
    """
    phi = np.asarray(series, dtype=float)
    phi = phi[np.isfinite(phi)]
    if phi.size < 100 or np.ptp(phi) < 1e-12:
        return 0.0  # too short or constant (a fixed point) -> regular
    rng = np.random.default_rng(seed)
    cs = rng.uniform(np.pi / 5.0, 4.0 * np.pi / 5.0, size=n_c)
    return float(np.median([_K_for_c(phi, c) for c in cs]))


def return_map(section_points: np.ndarray) -> np.ndarray:
    """Successive-crossing return map (x_n, x_{n+1}) for qualitative structure.
    PLAN §2, §5-E3."""
    x = np.asarray(section_points, dtype=float)
    return np.column_stack([x[:-1], x[1:]])

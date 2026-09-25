"""Verify the 0-1 test (Gottwald-Melbourne) on known series before it gates the
plate (E2 cross-check). K ~ 1 chaotic, K ~ 0 regular. See PLAN §2, §5-E2.
"""
from __future__ import annotations

import numpy as np

from fallingplate import chaos_tests


def _logistic(r, n=3000, x0=0.1, burn=1000):
    x = x0
    for _ in range(burn):
        x = r * x * (1 - x)
    xs = np.empty(n)
    for i in range(n):
        x = r * x * (1 - x)
        xs[i] = x
    return xs


def test_zero_one_chaotic():
    assert chaos_tests.zero_one_test(_logistic(4.0), n_c=100, seed=1) > 0.9
    assert chaos_tests.zero_one_test(_logistic(3.8), n_c=100, seed=1) > 0.9


def test_zero_one_regular():
    assert chaos_tests.zero_one_test(_logistic(3.2), n_c=100, seed=1) < 0.1   # period-2
    assert chaos_tests.zero_one_test(_logistic(3.5), n_c=100, seed=1) < 0.1   # period-4


def test_zero_one_degenerate_series_is_regular():
    assert chaos_tests.zero_one_test(np.full(500, 0.7)) == 0.0   # constant
    assert chaos_tests.zero_one_test(np.zeros(10)) == 0.0        # too short


def test_return_map_shape():
    x = np.arange(5.0)
    rm = chaos_tests.return_map(x)
    assert rm.shape == (4, 2)
    assert np.allclose(rm[:, 0], x[:-1]) and np.allclose(rm[:, 1], x[1:])

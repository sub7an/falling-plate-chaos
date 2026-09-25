"""Smoke tests for the scaffold (M0 acceptance: pytest collects and these pass).

No physics is asserted here. Real validation-ladder tests (PLAN §4) arrive with
the spec and implementation.
"""
import importlib

import pytest

MODULES = [
    "fallingplate",
    "fallingplate.model",
    "fallingplate.integrate",
    "fallingplate.regimes",
    "fallingplate.lyapunov",
    "fallingplate.chaos_tests",
    "fallingplate.displacement",
    "fallingplate.periodic_orbits",
    "fallingplate.io_cache",
    "fallingplate.sweeps",
    "fallingplate.plotting",
]


@pytest.mark.parametrize("name", MODULES)
def test_modules_import(name):
    assert importlib.import_module(name) is not None


def test_physics_is_implemented():
    """model.rhs is implemented (M1); it returns a length-4 derivative vector."""
    import numpy as np

    from fallingplate import model

    p = model.Params(C_T=1.2, C_R=np.pi, A=1.4, B=1.0, mu1=0.2, mu2=0.2, I_star=1.5)
    d = model.rhs(0.0, np.array([0.0, 0.1, 0.5, 0.0]), p)
    assert d.shape == (4,)
    assert np.all(np.isfinite(d))

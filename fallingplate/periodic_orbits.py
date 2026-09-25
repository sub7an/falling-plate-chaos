"""periodic_orbits.py — periodic orbits, Floquet analysis, continuation.

Used by E3/E4 to locate period-doublings precisely (Floquet multiplier crossing
-1) rather than by eye from a coarse sweep. Tumbling orbits: use the WRAPPED angle
in the shooting residual, record the winding number, and fix the phase with the
Poincare section. See PLAN §3, §5-E4.
"""
from __future__ import annotations

from typing import Callable

import numpy as np


def find_orbit_shooting(rhs: Callable, guess: np.ndarray, section: dict, cfg: dict) -> dict:
    """Newton shooting for a periodic orbit. Returns orbit state, period, and
    winding number. Residual uses the wrapped angle; phase fixed by `section`."""
    raise NotImplementedError("Scaffold only; see PLAN §5-E4.")


def monodromy(rhs: Callable, orbit: dict, cfg: dict) -> np.ndarray:
    """Monodromy matrix via the variational integrator (lyapunov.py). PLAN §5-E4."""
    raise NotImplementedError("Scaffold only; see PLAN §5-E4.")


def floquet_multipliers(monodromy_matrix: np.ndarray) -> np.ndarray:
    """Eigenvalues of the monodromy matrix; period-doubling at multiplier = -1."""
    raise NotImplementedError("Scaffold only; see PLAN §5-E4.")


def continue_in_Istar(rhs_family: Callable, orbit0: dict, i_star_grid: np.ndarray, cfg: dict) -> list:
    """Pseudo-arclength / natural continuation of an orbit in I_star, tracking
    Floquet multipliers to locate bifurcations. PLAN §5-E4."""
    raise NotImplementedError("Scaffold only; see PLAN §5-E4.")

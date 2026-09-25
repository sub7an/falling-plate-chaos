"""sweeps.py — parallel I* sweep driver.

Coarse-before-fine with adaptive refinement near transitions; joblib over the
available cores. Physics-free: `task` is a caller-supplied callable of one I*.
See PLAN §3, §7.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from joblib import Parallel, delayed


def sweep_istar(task: Callable, i_star_grid: np.ndarray, cfg: dict) -> list:
    """Run `task(I_star)` across the grid in parallel. Results are returned in
    grid order. `cfg['n_jobs']` selects cores (default -1 = all). PLAN §3, §7."""
    n_jobs = cfg.get("n_jobs", -1)
    grid = np.asarray(i_star_grid, dtype=float)
    return Parallel(n_jobs=n_jobs)(delayed(task)(float(I)) for I in grid)


def refine_near_transitions(i_star_grid: np.ndarray, labels: list, cfg: dict) -> np.ndarray:
    """Return I* points to ADD: midpoints of adjacent cells whose regime label
    differs (coarse-before-fine refinement near transitions). PLAN §7."""
    grid = np.asarray(i_star_grid, dtype=float)
    order = np.argsort(grid)
    grid = grid[order]
    labels = [labels[i] for i in order]
    extra = [0.5 * (grid[i] + grid[i + 1])
             for i in range(len(grid) - 1) if labels[i] != labels[i + 1]]
    return np.array(sorted(extra))

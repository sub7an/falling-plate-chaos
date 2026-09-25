"""Fast unit tests for the E1 plumbing added at M4 (no plate integration):
content-addressed cache round-trip, provenance sidecar, transition refinement,
and that the bifurcation figure is written. See PLAN §3, §4 rung 1, §5-E1.
"""
from __future__ import annotations

import json

import numpy as np

from experiments.e3_route import _n_branches, _successive_doublings
from fallingplate import io_cache, plotting
from fallingplate.sweeps import refine_near_transitions


def test_cache_key_depends_on_params_and_sources():
    k1 = io_cache.cache_key({"a": 1, "b": 2})
    assert k1 == io_cache.cache_key({"b": 2, "a": 1})   # order-independent
    assert k1 != io_cache.cache_key({"a": 1, "b": 3})   # value-sensitive


def test_cached_round_trip_and_sidecar(tmp_path, monkeypatch):
    monkeypatch.setattr(io_cache, "DATA_DIR", tmp_path)
    calls = {"n": 0}

    def compute():
        calls["n"] += 1
        return {"x": np.arange(5), "y": np.linspace(0, 1, 5)}

    params = {"foo": 1}
    key = io_cache.cache_key(params)

    r1 = io_cache.cached(key, compute, params, seed=7)
    r2 = io_cache.cached(key, compute, params, seed=7)  # served from disk
    assert calls["n"] == 1                               # computed once
    assert np.array_equal(r1["x"], r2["x"])
    assert np.allclose(r1["y"], r2["y"])

    sidecar = json.loads((tmp_path / f"{key}.json").read_text())
    assert sidecar["seed"] == 7
    assert set(sidecar["provenance"]) >= {"git_hash", "python", "numpy", "scipy"}


def test_refine_near_transitions_adds_midpoints():
    grid = np.array([1.0, 1.1, 1.2, 1.3])
    labels = ["flutter", "flutter", "tumble", "tumble"]
    extra = refine_near_transitions(grid, labels, {})
    assert np.allclose(extra, [1.15])   # only across the flutter->tumble boundary


def test_n_branches_counts_section_clusters():
    assert _n_branches(np.full(20, 0.6), tol=5e-4, cap=16) == 1          # period-1
    assert _n_branches(np.array([0.6, -0.79] * 10), tol=5e-4, cap=16) == 2  # period-2
    assert _n_branches(np.linspace(-1, 1, 200), tol=5e-4, cap=16) == -1  # smear -> aperiodic
    assert _n_branches(np.empty(0), tol=5e-4, cap=16) == 0


def test_successive_doublings_needs_clean_cascade():
    istar = np.array([1.45, 1.50, 1.55, 1.60, 1.65])
    # 1 -> 2 -> 4 -> 8 clean cascade => three doubling points (2,4,8)
    assert _successive_doublings(istar, np.array([1, 2, 4, 8, -1])) == [1.50, 1.55, 1.60]
    # 1 -> 2 then straight to chaos (no 4) => only one doubling recorded
    assert _successive_doublings(istar, np.array([1, 2, -1, -1, -1])) == [1.50]


def test_bifurcation_diagram_writes_file(tmp_path):
    data = {
        "istar": np.array([1.0, 1.5, 2.0, 2.5]),
        "category": np.array(["fluttering", "period-1 tumbling",
                              "aperiodic (chaos-candidate)", "fluttering"]),
        "scatter_istar": np.array([1.0, 1.0, 1.5, 2.0, 2.5]),
        "scatter_omega": np.array([0.5, -0.5, 0.7, 0.1, 1.0]),
    }
    out = tmp_path / "bif.png"
    path = plotting.bifurcation_diagram(data, str(out))
    assert out.exists() and out.stat().st_size > 0
    assert path == str(out)


def test_hysteresis_diagram_writes_file(tmp_path):
    I = np.linspace(1.1, 3.0, 10)
    data = {"istar_up": I, "up": np.linspace(0.2, 0.6, 10),
            "istar_down": I, "down": np.linspace(0.2, 0.6, 10),
            "basin_istar": np.repeat([1.4, 2.2], 3),
            "basin_absomega": np.array([0.49, 0.49, 0.49, 0.59, 0.59, 0.60])}
    out = tmp_path / "hyst.png"
    assert plotting.hysteresis_diagram(data, str(out)) == str(out)
    assert out.exists() and out.stat().st_size > 0

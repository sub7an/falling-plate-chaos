"""io_cache.py — content-addressed cache + provenance sidecars.

Cache key = hash(params) + hash of the physics/numerics source files
(model.py, integrate.py, lyapunov.py, ...). The git hash is NOT part of the key;
it is recorded in the sidecar only, so the cache invalidates on physics/numerics
edits but not on unrelated commits. Sidecar records: git hash, Python version,
NumPy/SciPy version, params, seed. See PLAN §3, §7, §9.
"""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np
import scipy

PKG_DIR = Path(__file__).resolve().parent
DATA_DIR = PKG_DIR.parent / "data"

# Files whose contents affect numeric output; hashed into the cache key.
PHYSICS_NUMERICS_SOURCES = ("model.py", "integrate.py", "lyapunov.py", "regimes.py")


def _hash_sources(sources: Iterable[str]) -> str:
    h = hashlib.sha256()
    for name in sorted(sources):
        h.update(name.encode())
        h.update((PKG_DIR / name).read_bytes())
    return h.hexdigest()


def cache_key(params: dict, sources: Iterable[str] = PHYSICS_NUMERICS_SOURCES) -> str:
    """Deterministic key from params + physics/numerics source-file hashes."""
    h = hashlib.sha256()
    h.update(json.dumps(params, sort_keys=True, default=str).encode())
    h.update(_hash_sources(sources).encode())
    return h.hexdigest()[:16]


def provenance() -> dict:
    """Collect git hash, Python/NumPy/SciPy versions for the sidecar. PLAN §9."""
    try:
        git = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=PKG_DIR, stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        git = "unknown"
    return {
        "git_hash": git,
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
    }


def cached(key: str, compute: Callable[[], dict], params: dict, seed: int) -> dict:
    """Return cached arrays for `key`; else compute, store as .npz, and write a
    .json sidecar (params, seed, provenance). `compute` returns a dict of arrays."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    npz = DATA_DIR / f"{key}.npz"
    if npz.exists():
        with np.load(npz, allow_pickle=False) as d:
            return {k: d[k] for k in d.files}
    result = compute()
    np.savez(npz, **result)
    sidecar = {"key": key, "params": params, "seed": seed, "provenance": provenance()}
    (DATA_DIR / f"{key}.json").write_text(json.dumps(sidecar, indent=2, default=str))
    return result

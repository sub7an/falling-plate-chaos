"""E6 — multistability / hysteresis. See PLAN §5-E6.

Two probes for coexisting attractors: (1) up vs down I* continuation sweeps
(seed each I* from the previous converged state) — a difference between branches
is hysteresis; (2) basin sampling — many random ICs at diagnostic I*, clustered
by a convergent order parameter <|omega|> (gross rotation rate; ergodic per
attractor). |<omega>| separately flags the CW/CCW symmetry pair. Null results
(single attractor, no hysteresis) are valid.
"""
from __future__ import annotations

from dataclasses import asdict
from functools import partial
from pathlib import Path

import numpy as np
import yaml
from joblib import Parallel, delayed
from scipy.integrate import solve_ivp

from fallingplate import io_cache, model, plotting

ROOT = Path(__file__).resolve().parents[1]


def _order_params(I_star: float, ic, cfg: dict):
    """Return (<|omega|>, <omega>, final_state) for a trajectory from `ic`."""
    p = model.Params.from_yaml(str(ROOT / cfg["model_config"]), I_star=I_star)
    tr, ig = cfg["trajectory"], cfg["integrator"]
    te = np.arange(0.0, tr["T"] + tr["dt"], tr["dt"])
    sol = solve_ivp(lambda t, s: model.rhs(t, s, p), (0.0, tr["T"]), ic,
                    method=ig["method"], rtol=ig["rtol"], atol=ig["atol"], t_eval=te)
    m = sol.t >= tr["transient"]
    omega = sol.y[3, m]
    return float(np.mean(np.abs(omega))), float(np.mean(omega)), sol.y[:, -1]


def _continuation(grid, ic0, cfg):
    ic = np.array(ic0, float)
    absw = []
    for I in grid:
        aw, _, ic = _order_params(float(I), ic, cfg)
        absw.append(aw)
    return np.array(absw)


def _basin_point(args):
    I_star, ic, cfg = args
    aw, mw, _ = _order_params(I_star, ic, cfg)
    return I_star, aw, mw


def main(config_path: str = "config/e6.yaml") -> None:
    cfg = yaml.safe_load((ROOT / config_path).read_text())
    sr = cfg["sweep_range"]
    grid = np.round(np.arange(sr["start"], sr["stop"] + sr["step"] / 2, sr["step"]), 3)
    ic0 = cfg["trajectory"]["ic"]

    key_params = {"sweep_range": sr, "basin": cfg["basin"], "trajectory": cfg["trajectory"],
                  "integrator": cfg["integrator"], "seed": cfg["seed"],
                  "coeffs": asdict(model.Params.from_yaml(str(ROOT / cfg["model_config"])))}
    key = io_cache.cache_key(key_params)

    def compute():
        up = _continuation(grid, ic0, cfg)
        _, _, ic_top = _order_params(float(grid[-1]), ic0, cfg)   # seed down-sweep
        down = _continuation(grid[::-1], ic_top, cfg)[::-1]
        # basin sampling
        rng = np.random.default_rng(cfg["seed"])
        b = cfg["basin"]
        tasks = []
        for I in b["istar"]:
            for _ in range(b["n_ic"]):
                ic = [rng.uniform(-b["vel_scale"], b["vel_scale"]),
                      rng.uniform(-b["vel_scale"], b["vel_scale"]),
                      rng.uniform(0, 2 * np.pi), rng.uniform(-1, 1)]
                tasks.append((float(I), ic, cfg))
        res = Parallel(n_jobs=cfg["sweep"]["n_jobs"])(delayed(_basin_point)(t) for t in tasks)
        bi = np.array([r[0] for r in res]); ba = np.array([r[1] for r in res])
        bm = np.array([r[2] for r in res])
        return {"istar": grid, "up": up, "down": down,
                "basin_istar": bi, "basin_absomega": ba, "basin_meanomega": bm}

    data = io_cache.cached(key, compute, key_params, cfg["seed"])

    fig_data = {"istar_up": data["istar"], "up": data["up"],
                "istar_down": data["istar"], "down": data["down"],
                "basin_istar": data["basin_istar"], "basin_absomega": data["basin_absomega"]}
    fig_path = ROOT / cfg["outputs"]["figure"]
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    plotting.hysteresis_diagram(fig_data, str(fig_path))

    # --- hysteresis: where do up and down branches differ? ---
    scale = max(float(np.max(data["up"])), 1e-9)
    hyst = np.abs(data["up"] - data["down"]) > cfg["cluster"]["tol"] * scale
    print(f"E6 multistability (cache {key})")
    if hyst.any():
        print(f"  HYSTERESIS: up != down at {int(hyst.sum())} I* points: "
              f"{np.array2string(data['istar'][hyst], precision=3)}")
    else:
        print("  No hysteresis: up and down continuation coincide across [1.1, 3.0].")

    # --- basin: distinct attractors per I* (cluster <|omega|>, symmetry-invariant) ---
    print("  Basin sampling (random ICs) — distinct attractors by <|omega|>:")
    b = cfg["basin"]
    tol = cfg["cluster"]["tol"] * scale
    for I in b["istar"]:
        sel = np.isclose(data["basin_istar"], I)
        vals = np.sort(data["basin_absomega"][sel])
        clusters = 1 + int(np.sum(np.diff(vals) > tol))
        mw = data["basin_meanomega"][sel]
        signs = np.unique(np.sign(np.round(mw, 2)))
        cw_ccw = "±ω pair" if (-1 in signs and 1 in signs) else "one direction"
        print(f"    I*={I:<6}: {clusters} <|ω|>-cluster(s) over {sel.sum()} ICs "
              f"(range {vals.min():.3f}-{vals.max():.3f}); net rotation: {cw_ccw}")
    print(f"figure -> {fig_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

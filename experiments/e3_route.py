"""E3 — route to chaos. See PLAN §5-E3.

Runs only because E2 (M5) confirmed a chaotic band. Characterises the transition
from the period-2 orbit (born ~1.474) to chaos (~1.61): a fine theta=pi
bifurcation zoom + lambda_max, counting SUCCESSIVE period-doublings. The E4
Feigenbaum test requires >= 3 successive doublings; this experiment decides
whether that precondition is met (it must be established, not assumed).
"""
from __future__ import annotations

from dataclasses import asdict
from functools import partial
from pathlib import Path

import numpy as np
import yaml
from scipy.integrate import solve_ivp

from fallingplate import io_cache, model, plotting
from fallingplate import lyapunov as L
from fallingplate.sweeps import sweep_istar

ROOT = Path(__file__).resolve().parents[1]


def _n_branches(omega: np.ndarray, tol: float, cap: int) -> int:
    """Distinct theta=pi section branches (period counter). > cap => aperiodic (-1)."""
    v = np.sort(np.asarray(omega))
    if v.size == 0:
        return 0
    n = 1
    for i in range(1, v.size):
        if v[i] - v[i - 1] > tol:
            n += 1
        if n > cap:
            return -1
    return n


def e3_task(I_star: float, cfg: dict) -> dict:
    p = model.Params.from_yaml(str(ROOT / cfg["model_config"]), I_star=I_star)
    ic = cfg["trajectory"]["ic"]
    ig = cfg["integrator"]

    def g(t, s):  # zeros at theta = phase (mod 2pi); smooth (no wrap)
        return np.sin((s[2] - cfg["section"]["phase"]) / 2.0)
    g.direction = 0

    sol = solve_ivp(lambda t, s: model.rhs(t, s, p), (0.0, ig["T"]), ic,
                    method=ig["method"], rtol=ig["rtol"], atol=ig["atol"],
                    events=g, max_step=ig["max_step"])
    te, ye = sol.t_events[0], sol.y_events[0]
    omega = ye[te >= ig["transient"]][:, 3] if ye.size else np.empty(0)
    omega = omega[-ig["keep"]:]
    period = _n_branches(omega, cfg["period"]["cluster_tol"], cfg["period"]["cap"])

    lc = cfg["lyapunov"]
    lam = float(L.variational_spectrum(
        lambda t, s: model.rhs(t, s, p), lambda t, s: model.jacobian(t, s, p), ic,
        dict(t_total=lc["t_total"], dt_renorm=lc["dt_renorm"], t_transient=lc["t_transient"],
             method=lc["method"], rtol=lc["rtol"], atol=lc["atol"], k=1))[0])

    return {"I_star": I_star, "period": period, "lam_max": lam, "omega": omega}


def _assemble(results: list) -> dict:
    results = sorted(results, key=lambda r: r["I_star"])
    si, sw = [], []
    for r in results:
        si.extend([r["I_star"]] * r["omega"].size)
        sw.extend(r["omega"].tolist())
    return {
        "istar": np.array([r["I_star"] for r in results]),
        "period": np.array([r["period"] for r in results]),
        "lam_max": np.array([r["lam_max"] for r in results]),
        "scatter_istar": np.array(si),
        "scatter_omega": np.array(sw),
    }


def _successive_doublings(istar: np.ndarray, period: np.ndarray) -> list:
    """I* where the branch count first reaches 2, then 4, then 8, ... in order
    (a clean cascade). Stops at the first missing power of two."""
    order = np.argsort(istar)
    I, P = istar[order], period[order]
    doublings, target = [], 2
    for i in range(len(I)):
        if P[i] == target:
            doublings.append(float(I[i]))
            target *= 2
    return doublings


def main(config_path: str = "config/e3.yaml") -> None:
    cfg = yaml.safe_load((ROOT / config_path).read_text())
    s = cfg["istar_sweep"]
    grid = np.linspace(s["start"], s["stop"], s["num"])
    task = partial(e3_task, cfg=cfg)

    key_params = {
        "istar_sweep": s, "trajectory": cfg["trajectory"], "section": cfg["section"],
        "integrator": cfg["integrator"], "lyapunov": cfg["lyapunov"], "period": cfg["period"],
        "coeffs": asdict(model.Params.from_yaml(str(ROOT / cfg["model_config"]))),
    }
    key = io_cache.cache_key(key_params)
    data = io_cache.cached(key, lambda: _assemble(sweep_istar(task, grid, cfg["sweep"])),
                           key_params, cfg["seed"])

    doublings = _successive_doublings(data["istar"], data["period"])
    data = dict(data)
    data["doublings"] = np.array(doublings)

    fig_path = ROOT / cfg["outputs"]["figure"]
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    plotting.route_diagram(data, str(fig_path))

    # --- route summary ---
    I, P = data["istar"], data["period"]
    order = np.argsort(I)
    print(f"E3 route to chaos: {len(I)} I* points in [{I.min():.3f}, {I.max():.3f}] (cache {key})")
    # period bands
    start = order[0]
    print("Period bands (I* range -> branch count; -1 = aperiodic/chaotic):")
    seg0 = 0
    Is, Ps = I[order], P[order]
    for i in range(1, len(Ps) + 1):
        if i == len(Ps) or Ps[i] != Ps[seg0]:
            tag = "chaotic" if Ps[seg0] == -1 else f"{Ps[seg0]} branch(es)"
            print(f"  [{Is[seg0]:.3f}, {Is[i-1]:.3f}]  {tag}")
            seg0 = i
    print(f"successive period-doublings (clean cascade): {doublings}  (n={len(doublings)})")
    if len(doublings) >= 3:
        print("  => >=3 successive doublings: E4 Feigenbaum test IS warranted.")
    else:
        print("  => < 3 successive doublings: E4 Feigenbaum test NOT warranted (null result). "
              "Route is not a clean period-doubling cascade; do not force Feigenbaum.")
    print(f"figure -> {fig_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

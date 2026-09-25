"""E2 — Lyapunov vs I*. THE CHAOS GATE. See PLAN §5-E2.

For each I*: lambda_max via the verified variational spectrum (M2, analytic
Jacobian), with uncertainty from two averaging times; and median K from the 0-1
test on the theta=0 Poincare-section omega series (independent cross-check).
GATE: if lambda_max <= 0 everywhere within uncertainty AND K ~ 0, STOP and
propose alternatives — do not force a chaos narrative.
"""
from __future__ import annotations

from dataclasses import asdict
from functools import partial
from pathlib import Path

import numpy as np
import yaml
from scipy.integrate import solve_ivp

from fallingplate import chaos_tests, io_cache, model, plotting, regimes
from fallingplate import lyapunov as L
from fallingplate.sweeps import sweep_istar

ROOT = Path(__file__).resolve().parents[1]


def _lam_max(p: model.Params, ic, lc: dict, t_total: float) -> float:
    cfg = dict(t_total=t_total, dt_renorm=lc["dt_renorm"], t_transient=lc["t_transient"],
               method=lc["method"], rtol=lc["rtol"], atol=lc["atol"], k=1)
    return float(L.variational_spectrum(
        lambda t, s: model.rhs(t, s, p), lambda t, s: model.jacobian(t, s, p), ic, cfg)[0])


def _section_omega(p: model.Params, ic, zc: dict) -> np.ndarray:
    def g(t, s):
        return np.sin(s[2] / 2.0)   # theta = 0 (mod 2pi)
    g.direction = 0
    sol = solve_ivp(lambda t, s: model.rhs(t, s, p), (0.0, zc["T"]), ic,
                    method=zc["method"], rtol=zc["rtol"], atol=zc["atol"],
                    events=g, max_step=zc["max_step"])
    te, ye = sol.t_events[0], sol.y_events[0]
    return ye[te >= zc["transient"]][:, 3] if ye.size else np.empty(0)


def e2_task(I_star: float, cfg: dict) -> dict:
    p = model.Params.from_yaml(str(ROOT / cfg["model_config"]), I_star=I_star)
    ic = cfg["trajectory"]["ic"]
    lc = cfg["lyapunov"]

    lam_long = _lam_max(p, ic, lc, lc["t_long"])
    lam_short = _lam_max(p, ic, lc, lc["t_short"])
    series = _section_omega(p, ic, cfg["zero_one"])
    K = chaos_tests.zero_one_test(series, n_c=cfg["zero_one"]["n_c"], seed=cfg["seed"])

    return {
        "I_star": I_star,
        "lam_max": lam_long,
        "lam_unc": abs(lam_long - lam_short),
        "K": K,
        "n_section": series.size,
    }


def _assemble(results: list) -> dict:
    results = sorted(results, key=lambda r: r["I_star"])
    return {k: np.array([r[k] for r in results])
            for k in ("I_star", "lam_max", "lam_unc", "K", "n_section")} | {
        "istar": np.array([r["I_star"] for r in results])}


def _load_e1(cfg) -> dict | None:
    """Load the cached E1 bundle (bifurcation backdrop), if present."""
    try:
        e1cfg = yaml.safe_load((ROOT / cfg["e1_config"]).read_text())
        s = e1cfg["istar_sweep"]
        key_params = {
            "istar_sweep": s, "integrator": e1cfg["integrator"],
            "trajectory": e1cfg["trajectory"], "section": e1cfg["section"],
            "refine": e1cfg["sweep"].get("refine", False),
            "coeffs": asdict(model.Params.from_yaml(str(ROOT / e1cfg["model_config"]))),
        }
        npz = io_cache.DATA_DIR / f"{io_cache.cache_key(key_params)}.npz"
        if npz.exists():
            with np.load(npz, allow_pickle=False) as d:
                return {k: d[k] for k in d.files}
    except Exception:
        pass
    return None


def main(config_path: str = "config/e2.yaml") -> None:
    cfg = yaml.safe_load((ROOT / config_path).read_text())
    s = cfg["istar_sweep"]
    grid = np.linspace(s["start"], s["stop"], s["num"])
    task = partial(e2_task, cfg=cfg)

    key_params = {
        "istar_sweep": s, "trajectory": cfg["trajectory"],
        "lyapunov": cfg["lyapunov"], "zero_one": cfg["zero_one"], "seed": cfg["seed"],
        "coeffs": asdict(model.Params.from_yaml(str(ROOT / cfg["model_config"]))),
    }
    key = io_cache.cache_key(key_params)
    data = io_cache.cached(key, lambda: _assemble(sweep_istar(task, grid, cfg["sweep"])),
                           key_params, cfg["seed"])

    # figure
    fig_path = ROOT / cfg["outputs"]["figure"]
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    plotting.lyapunov_overlay(data, str(fig_path), bif=_load_e1(cfg))

    # --- GATE decision ---
    I = data["istar"]
    lam, unc, K = data["lam_max"], data["lam_unc"], data["K"]
    lam_pos, chaotic_K = cfg["gate"]["lam_pos"], cfg["zero_one"]["chaotic_K"]
    confirmed = (lam - unc > lam_pos) & (K > chaotic_K)

    print(f"E2 Lyapunov gate: {len(I)} I* points (cache {key})")
    print(f"  lambda_max peak = {lam.max():+.3f} at I*={I[np.argmax(lam)]:.3f}; "
          f"max median K = {K.max():.3f}")
    if confirmed.any():
        lo, hi = I[confirmed].min(), I[confirmed].max()
        print(f"  CHAOS CONFIRMED (lambda_max-unc>0 AND K>{chaotic_K}) at "
              f"{int(confirmed.sum())} points in I* ~ [{lo:.3f}, {hi:.3f}]")
        print("  GATE: PASS -> chaotic band exists; proceed to E3/E4.")
    else:
        print("  GATE: FAIL -> no confirmed chaos (lambda_max<=0 within unc or K~0). "
              "STOP and propose alternatives (PLAN §5-E2).")
    # points that are positive-lambda OR high-K but not both (for review)
    disagree = (confirmed != ((lam - unc > lam_pos) | (K > chaotic_K)))
    if disagree.any():
        pts = ", ".join(f"{i:.3f}(λ={l:+.2f}±{u:.2f},K={k:.2f})"
                        for i, l, u, k in zip(I[disagree], lam[disagree], unc[disagree], K[disagree]))
        print(f"  method disagreement (one signal only): {pts}")
    print(f"figure -> {fig_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

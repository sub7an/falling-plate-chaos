"""E1 — regime map vs I*. See PLAN §5-E1.

Hypothesis: ordered regime bands exist; at least one aperiodic band (chaos
candidate, to be confirmed by E2's lambda_max). Reads config/e1.yaml; writes a
cached array bundle + provenance sidecar to data/ and a bifurcation-style figure
to figures/. Uses the verified model.py (M1) and regimes.py (M3); NO Lyapunov
here (that is E2/M5).
"""
from __future__ import annotations

from dataclasses import asdict
from functools import partial
from pathlib import Path

import numpy as np
import yaml
from scipy.integrate import solve_ivp

from fallingplate import io_cache, model, plotting, regimes
from fallingplate.sweeps import refine_near_transitions, sweep_istar

ROOT = Path(__file__).resolve().parents[1]


def _category(regime: str, period: int, label: str) -> str:
    """Display category for the regime strip. period == -1 (aperiodic section)
    is the E1 chaos candidate, to be confirmed by E2."""
    if period == -1:
        return "aperiodic (chaos-candidate)"
    if regime == "tumbling":
        name = {1: "period-1", 2: "period-2", 4: "period-4"}.get(period, "periodic")
        return f"{name} tumbling"
    if regime == "fluttering":
        return label  # "fluttering" or "small-amplitude fluttering"
    if regime == regimes.Regime.MIXTURE.value:
        return "flutter/tumble mixture"
    return regime


def regime_task(I_star: float, cfg: dict) -> dict:
    """Integrate one I*, classify the regime, and return diagnostics + the
    theta=0 Poincare-section omega values (event-located, integrator precision)."""
    p = model.Params.from_yaml(str(ROOT / cfg["model_config"]), I_star=I_star)
    tr = cfg["trajectory"]
    itg = cfg["integrator"]
    T, transient, dt = tr["T"], tr["transient"], tr["dt"]

    def sec_event(t, s):
        return np.sin((s[2] - cfg["section"]["phase"]) / 2.0)
    sec_event.direction = 0

    t_eval = np.arange(0.0, T + dt, dt)
    sol = solve_ivp(
        lambda t, s: model.rhs(t, s, p), (0.0, T), tr["ic"],
        method=itg["method"], rtol=itg["rtol"], atol=itg["atol"],
        t_eval=t_eval, events=sec_event, max_step=itg["max_step"],
    )
    mask = sol.t >= transient
    rep = regimes.classify_regime(sol.t[mask], sol.y[:, mask])

    te, ye = sol.t_events[0], sol.y_events[0]
    keep_mask = te >= transient
    sec_omega = ye[keep_mask][:, 3] if ye.size else np.empty(0)
    if sec_omega.size >= 4:
        period = regimes.poincare_period(ye[keep_mask][:, [0, 1, 3]])
    else:
        period = rep.poincare_period

    return {
        "I_star": I_star,
        "regime": rep.regime.value,
        "category": _category(rep.regime.value, period, rep.label),
        "period": period,
        "Omega": rep.Omega, "net_rev": rep.net_rev, "R": rep.directedness,
        "theta_amp": rep.theta_amp,
        "sec_omega": np.asarray(sec_omega[-tr["keep"]:], dtype=float),
    }


def _assemble(results: list) -> dict:
    results = sorted(results, key=lambda r: r["I_star"])
    istar = np.array([r["I_star"] for r in results])
    scatter_i, scatter_w = [], []
    for r in results:
        scatter_i.extend([r["I_star"]] * r["sec_omega"].size)
        scatter_w.extend(r["sec_omega"].tolist())
    return {
        "istar": istar,
        "category": np.array([r["category"] for r in results]),
        "regime": np.array([r["regime"] for r in results]),
        "period": np.array([r["period"] for r in results]),
        "Omega": np.array([r["Omega"] for r in results]),
        "net_rev": np.array([r["net_rev"] for r in results]),
        "R": np.array([r["R"] for r in results]),
        "theta_amp": np.array([r["theta_amp"] for r in results]),
        "scatter_istar": np.array(scatter_i),
        "scatter_omega": np.array(scatter_w),
    }


def main(config_path: str = "config/e1.yaml") -> None:
    cfg = yaml.safe_load((ROOT / config_path).read_text())
    s = cfg["istar_sweep"]
    coarse = np.linspace(s["start"], s["stop"], s["num"])
    task = partial(regime_task, cfg=cfg)

    key_params = {
        "istar_sweep": s,
        "integrator": cfg["integrator"],
        "trajectory": cfg["trajectory"],
        "section": cfg["section"],
        "refine": cfg["sweep"].get("refine", False),
        "coeffs": asdict(model.Params.from_yaml(str(ROOT / cfg["model_config"]))),
    }
    key = io_cache.cache_key(key_params)

    def compute() -> dict:
        results = sweep_istar(task, coarse, cfg["sweep"])
        if cfg["sweep"].get("refine", False):
            cats = [r["category"] for r in results]
            extra = refine_near_transitions(coarse, cats, cfg["sweep"])
            if extra.size:
                results += sweep_istar(task, extra, cfg["sweep"])
        return _assemble(results)

    data = io_cache.cached(key, compute, key_params, cfg["seed"])

    # figure
    fig_path = ROOT / cfg["outputs"]["figure"]
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    plotting.bifurcation_diagram(data, str(fig_path))

    # compact band summary (never dump arrays)
    istar = data["istar"]
    order = np.argsort(istar)
    cats = data["category"][order]
    iss = istar[order]
    print(f"E1 regime map: {len(iss)} I* points in [{iss.min():.3f}, {iss.max():.3f}]  "
          f"(cache {key})")
    print("Bands (I* range -> regime):")
    start = 0
    for i in range(1, len(cats) + 1):
        if i == len(cats) or cats[i] != cats[start]:
            print(f"  [{iss[start]:.3f}, {iss[i-1]:.3f}]  {cats[start]}")
            start = i
    print(f"figure -> {fig_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

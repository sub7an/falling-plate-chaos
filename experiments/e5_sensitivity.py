"""E5 — landing-displacement sensitivity. See PLAN §5-E5.

For each I*: an IC-perturbation ensemble of falling cards; sigma(D) = std of the
horizontal landing displacement at fixed depth D. Displacement is an unbounded
running sum, so we expect (chaotic band) an exponential Lyapunov phase
log sigma ~ log sigma0 + (lambda_max/U) D, then diffusive ~sqrt(D); (regular)
bounded / no growth. We fit exp/linear/sqrt, compare the exp slope to
lambda_max/U (lambda_max from the verified variational method, U = mean descent
speed), and report which model dominates the late phase.
"""
from __future__ import annotations

from dataclasses import asdict
from functools import partial
from pathlib import Path

import numpy as np
import yaml
from scipy.integrate import solve_ivp

from fallingplate import displacement, io_cache, model, plotting
from fallingplate import lyapunov as L
from fallingplate.sweeps import sweep_istar

ROOT = Path(__file__).resolve().parents[1]


def _member(I_star: float, dstate: np.ndarray, cfg: dict):
    p = model.Params.from_yaml(str(ROOT / cfg["model_config"]), I_star=I_star)
    ig = cfg["integrator"]
    ic6 = np.concatenate([[0.0, 0.0], np.asarray(cfg["ensemble"]["ic"]) + dstate])
    t_eval = np.arange(0.0, ig["T"] + ig["dt"], ig["dt"])
    sol = solve_ivp(lambda t, s: model.rhs_full(t, s, p), (0.0, ig["T"]), ic6,
                    method=ig["method"], rtol=ig["rtol"], atol=ig["atol"], t_eval=t_eval)
    return sol.t, sol.y


def e5_task(I_star: float, cfg: dict) -> dict:
    en = cfg["ensemble"]
    rng = np.random.default_rng(cfg["seed"])
    perts = en["delta0"] * rng.standard_normal((en["n_members"], 4))
    perts[0] = 0.0
    d = cfg["depths"]
    depths = np.linspace(d["start"], d["stop"], d["num"])

    trajs = [_member(I_star, perts[k], cfg) for k in range(en["n_members"])]
    t0, y0 = trajs[0]
    U = float(-(y0[1][-1] - y0[1][0]) / (t0[-1] - t0[0]))   # mean descent speed (chords/time)

    dx = np.full((en["n_members"], depths.size), np.nan)
    graze = 0
    for k, (t, y6) in enumerate(trajs):
        for j, D in enumerate(depths):
            Ld = displacement.landing_displacement(t, y6, float(D), graze_tol=en["graze_tol"])
            if Ld.found:
                dx[k, j] = Ld.dx
            graze += int(Ld.grazing_flag)
    sigma = displacement.spread_vs_depth(dx, np.isfinite(dx))

    # sigma vs TIME (removes the depth<->time U conversion): d log sigma_x / dt -> lambda_max
    t_grid = trajs[0][0]
    sigma_time = np.vstack([y6[0] for _, y6 in trajs]).std(axis=0)

    lc = cfg["lyapunov"]
    p = model.Params.from_yaml(str(ROOT / cfg["model_config"]), I_star=I_star)
    lam = float(L.variational_spectrum(
        lambda t, s: model.rhs(t, s, p), lambda t, s: model.jacobian(t, s, p),
        cfg["ensemble"]["ic"],
        dict(t_total=lc["t_total"], dt_renorm=lc["dt_renorm"], t_transient=lc["t_transient"],
             method=lc["method"], rtol=lc["rtol"], atol=lc["atol"], k=1))[0])

    return {"I_star": I_star, "depths": depths, "sigma": sigma, "U": U, "lam_max": lam,
            "t_grid": t_grid, "sigma_time": sigma_time,
            "graze_frac": graze / (en["n_members"] * depths.size)}


def _exp_phase_slope(x: np.ndarray, sig: np.ndarray, sat_frac: float):
    """Slope of log(sig) vs x over the exponential phase (x = depth or time), and
    the phase boundary. Returns (slope, x_max_of_phase)."""
    m = np.isfinite(sig) & (sig > 0) & (x > 0)
    if m.sum() < 3:
        return np.nan, float(x[-1])
    smax = sig[m].max()
    below = x[m][sig[m] < sat_frac * smax]
    xmax = float(below.max()) if below.size >= 3 else float(x[len(x) // 3])
    em = m & (x <= xmax)
    slope = float(np.polyfit(x[em], np.log(sig[em]), 1)[0]) if em.sum() >= 3 else np.nan
    return slope, xmax


def _analyse(r: dict, fitcfg: dict) -> dict:
    D, s = r["depths"], r["sigma"]
    m = np.isfinite(s) & (s > 0)
    smax, smin = np.nanmax(s[m]), np.nanmin(s[m])
    regular = (smax / smin) < fitcfg["min_growth_ratio"]
    exp_slope, exp_dmax, late_dom = np.nan, float(D[-1]), None
    exp_slope_time = np.nan
    if not regular:
        exp_slope, exp_dmax = _exp_phase_slope(D, s, fitcfg["sat_frac"])
        # late-phase: sqrt vs linear (log-log slope over D > exp_dmax)
        lm = m & (D > exp_dmax)
        if lm.sum() >= 3:
            p = float(np.polyfit(np.log(D[lm]), np.log(s[lm]), 1)[0])
            late_dom = f"~D^{p:.2f} ({'sqrt-like' if p < 0.75 else 'linear-like' if p < 1.4 else 'super-linear'})"
        # TIME-based slope (should equal lambda_max; removes the U conversion)
        exp_slope_time, _ = _exp_phase_slope(r["t_grid"], r["sigma_time"], fitcfg["sat_frac"])
    lam = r["lam_max"]
    return {"regular": bool(regular), "exp_slope": exp_slope, "exp_dmax": exp_dmax,
            "lam_over_U": lam / r["U"] if r["U"] > 0 else np.nan,
            "exp_slope_time": exp_slope_time,
            "slope_time_over_lam": exp_slope_time / lam if (lam and np.isfinite(exp_slope_time)) else np.nan,
            "late_dominant": late_dom, "growth_ratio": float(smax / smin)}


def main(config_path: str = "config/e5.yaml") -> None:
    cfg = yaml.safe_load((ROOT / config_path).read_text())
    istar = cfg["istar"]
    task = partial(e5_task, cfg=cfg)

    key_params = {"istar": istar, "ensemble": cfg["ensemble"], "integrator": cfg["integrator"],
                  "depths": cfg["depths"], "lyapunov": cfg["lyapunov"], "seed": cfg["seed"],
                  "analysis_version": 2,  # v2: adds sigma-vs-time (dlogsigma/dt -> lambda_max)
                  "coeffs": asdict(model.Params.from_yaml(str(ROOT / cfg["model_config"])))}
    key = io_cache.cache_key(key_params)

    def compute():
        results = sorted(sweep_istar(task, np.array(istar), cfg["sweep"]),
                         key=lambda r: r["I_star"])
        out = {"istar": np.array([r["I_star"] for r in results]),
               "depths": results[0]["depths"],
               "sigma": np.vstack([r["sigma"] for r in results]),
               "U": np.array([r["U"] for r in results]),
               "lam_max": np.array([r["lam_max"] for r in results]),
               "t_grid": results[0]["t_grid"],
               "sigma_time": np.vstack([r["sigma_time"] for r in results]),
               "graze_frac": np.array([r["graze_frac"] for r in results])}
        return out

    data = io_cache.cached(key, compute, key_params, cfg["seed"])

    # per-I* analysis
    analyses = []
    for i, I in enumerate(data["istar"]):
        r = {"depths": data["depths"], "sigma": data["sigma"][i],
             "U": float(data["U"][i]), "lam_max": float(data["lam_max"][i]),
             "t_grid": data["t_grid"], "sigma_time": data["sigma_time"][i]}
        analyses.append(_analyse(r, cfg["fit"]))

    fig_data = dict(data)
    fig_data["regular"] = np.array([a["regular"] for a in analyses])
    fig_data["exp_slope"] = np.array([a["exp_slope"] for a in analyses])
    fig_data["exp_dmax"] = np.array([a["exp_dmax"] for a in analyses])
    fig_path = ROOT / cfg["outputs"]["figure"]
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    plotting.sensitivity_fits(fig_data, str(fig_path))

    band = cfg["chaotic_band"]
    print(f"E5 landing-displacement sensitivity (cache {key}); chaotic band {band}")
    print(f"{'I*':>5} {'regime':>9} {'growth x':>9} {'slopeD':>8} {'lam/U':>7} "
          f"{'slope_t':>8} {'lam_max':>8} {'slope_t/lam':>11} {'graze':>6}")
    for I, a, U, lam, gz in zip(data["istar"], analyses, data["U"], data["lam_max"],
                                data["graze_frac"]):
        reg = "regular" if a["regular"] else "chaotic"
        es = f"{a['exp_slope']:.4f}" if np.isfinite(a["exp_slope"]) else "  --"
        lu = f"{a['lam_over_U']:.4f}" if np.isfinite(a["lam_over_U"]) else "--"
        st = f"{a['exp_slope_time']:.4f}" if np.isfinite(a["exp_slope_time"]) else "  --"
        rr = f"{a['slope_time_over_lam']:.3f}" if np.isfinite(a["slope_time_over_lam"]) else " --"
        print(f"{I:>5.2f} {reg:>9} {a['growth_ratio']:>9.1e} {es:>8} {lu:>7} "
              f"{st:>8} {lam:>8.4f} {rr:>11} {gz:>6.2f}")
    print("  slopeD vs lam/U mixes in the depth<->time (U) conversion; slope_t = dlogSigma/dt is the")
    print("  direct test: slope_t should equal lambda_max. Regular bands = bounded, no growth.")
    print(f"figure -> {fig_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

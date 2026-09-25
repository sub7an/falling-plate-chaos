"""plotting.py — figure builders.

Never prints arrays or long logs; consumes cached arrays and writes figures to
figures/. See PLAN §3, §7.
"""
from __future__ import annotations

from typing import Any

import matplotlib

matplotlib.use("Agg")  # headless; no display
import matplotlib.pyplot as plt
import numpy as np


def _cell_edges(x: np.ndarray) -> np.ndarray:
    """Edges midway between sorted grid points (for a colored strip)."""
    x = np.sort(np.asarray(x, dtype=float))
    mids = 0.5 * (x[1:] + x[:-1])
    return np.concatenate([[x[0] - (mids[0] - x[0]) if mids.size else x[0] - 0.5],
                           mids,
                           [x[-1] + (x[-1] - mids[-1]) if mids.size else x[-1] + 0.5]])


def bifurcation_diagram(data: Any, out_path: str) -> str:
    """E1 figure: Poincare section (theta=0) omega vs I* + a regime strip.

    `data` keys: scatter_istar, scatter_omega (flattened section points);
    istar, category (per-grid regime category string). Returns the written path.
    PLAN §5-E1.
    """
    istar = np.asarray(data["istar"], dtype=float)
    categories = np.asarray(data["category"]).astype(str)
    order = np.argsort(istar)
    istar_s, cat_s = istar[order], categories[order]

    cats = sorted(set(cat_s.tolist()))
    cmap = plt.get_cmap("tab10" if len(cats) <= 10 else "tab20")
    color = {c: cmap(i % cmap.N) for i, c in enumerate(cats)}

    fig, (ax, strip) = plt.subplots(
        2, 1, figsize=(9, 6), sharex=True,
        gridspec_kw={"height_ratios": [5, 1], "hspace": 0.05},
    )

    ax.scatter(data["scatter_istar"], data["scatter_omega"], s=2, c="k", alpha=0.6, linewidths=0)
    ax.set_ylabel(r"$\omega$ at section ($\theta=0$)")
    ax.set_title(r"E1 - regime map: Poincare section $\omega$ vs $I^*$")
    ax.grid(alpha=0.2)

    edges = _cell_edges(istar_s)
    for i, c in enumerate(cat_s):
        strip.axvspan(edges[i], edges[i + 1], color=color[c])
    strip.set_yticks([])
    strip.set_xlabel(r"$I^*$")
    strip.set_ylabel("regime", rotation=0, ha="right", va="center")

    handles = [plt.Rectangle((0, 0), 1, 1, color=color[c]) for c in cats]
    strip.legend(handles, cats, loc="upper center", bbox_to_anchor=(0.5, -0.6),
                 ncol=min(3, len(cats)), fontsize=8, frameon=False)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def hysteresis_diagram(data: Any, out_path: str) -> str:
    """E6 figure: up vs down continuation of <|omega|>(I*) (hysteresis) plus
    basin-sampling points (coexisting attractors). `data` keys: istar_up, up,
    istar_down, down, basin_istar, basin_absomega. PLAN §5-E6.
    """
    fig, (ax, axb) = plt.subplots(2, 1, figsize=(9, 6), sharex=True,
                                  gridspec_kw={"hspace": 0.08})
    ax.plot(data["istar_up"], data["up"], "-o", ms=3, color="C0", label="up sweep")
    ax.plot(data["istar_down"], data["down"], "-s", ms=3, color="C3",
            mfc="none", label="down sweep")
    ax.set_ylabel(r"$\langle|\omega|\rangle$")
    ax.set_title(r"E6 - up/down continuation (hysteresis) and basin sampling")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2)

    axb.scatter(data["basin_istar"], data["basin_absomega"], s=14, c="k", alpha=0.6)
    axb.set_ylabel(r"$\langle|\omega|\rangle$ (random ICs)")
    axb.set_xlabel(r"$I^*$")
    axb.grid(alpha=0.2)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def sensitivity_fits(data: Any, out_path: str) -> str:
    """E5 figure: sigma(D) per I* (semilog, with exponential-phase fits) plus a
    log-log panel of the chaotic curves against sqrt(D) and linear references.
    `data` keys: istar, depths, sigma [n_istar,n_depths], regular (bool),
    exp_slope, lam_over_U, exp_dmax. PLAN §5-E5.
    """
    istar = np.asarray(data["istar"], dtype=float)
    D = np.asarray(data["depths"], dtype=float)
    sig = np.asarray(data["sigma"], dtype=float)
    regular = np.asarray(data["regular"], dtype=bool)

    fig, (ax, axll) = plt.subplots(1, 2, figsize=(12, 5))
    cmap = plt.get_cmap("viridis")
    for i, I in enumerate(istar):
        c = cmap(i / max(1, len(istar) - 1))
        ls = "--" if regular[i] else "-"
        lbl = f"I*={I:.2f}" + ("  (regular)" if regular[i] else "  (chaotic)")
        ax.semilogy(D, sig[i], ls, color=c, lw=1.4, label=lbl)
        if not regular[i] and np.isfinite(data["exp_slope"][i]):
            b = data["exp_slope"][i]
            dmax = data["exp_dmax"][i]
            mask = D <= dmax
            a = np.log(sig[i][mask][0]) - b * D[mask][0]
            ax.semilogy(D[mask], np.exp(a + b * D[mask]), ":", color=c, lw=2)
    ax.set_xlabel(r"depth $D$ (chords)")
    ax.set_ylabel(r"$\sigma(D) = \mathrm{std}(\Delta x)$")
    ax.set_title(r"E5 - landing-displacement spread (dotted = exp-phase fit)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2, which="both")

    for i, I in enumerate(istar):
        if regular[i]:
            continue
        c = cmap(i / max(1, len(istar) - 1))
        axll.loglog(D, sig[i], "-", color=c, lw=1.4, label=f"I*={I:.2f}")
    # reference slopes anchored to the chaotic set
    chaotic = ~regular
    if chaotic.any():
        s0 = np.nanmax(sig[chaotic][:, len(D) // 2])
        axll.loglog(D, s0 * (D / D[len(D) // 2]) ** 0.5, "k--", lw=1, label=r"$\sqrt{D}$")
        axll.loglog(D, s0 * (D / D[len(D) // 2]) ** 1.0, "k:", lw=1, label=r"linear $D$")
    axll.set_xlabel(r"depth $D$ (chords)")
    axll.set_ylabel(r"$\sigma(D)$")
    axll.set_title("chaotic curves (log-log) vs power-law references")
    axll.legend(fontsize=8)
    axll.grid(alpha=0.2, which="both")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def route_diagram(data: Any, out_path: str) -> str:
    """E3 figure: fine bifurcation zoom (theta=pi section omega vs I*) over the
    route-to-chaos window, with lambda_max below and detected period-doublings
    marked. `data` keys: scatter_istar, scatter_omega, istar, lam_max, period,
    doublings (I* values). PLAN §5-E3.
    """
    fig, (ax, axl) = plt.subplots(2, 1, figsize=(9, 6), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 2], "hspace": 0.08})
    ax.scatter(data["scatter_istar"], data["scatter_omega"], s=2, c="k", alpha=0.5, linewidths=0)
    ax.set_ylabel(r"$\omega$ at section ($\theta=\pi$)")
    ax.set_title(r"E3 - route to chaos: bifurcation zoom + $\lambda_{\max}$")
    ax.grid(alpha=0.2)
    for d in np.atleast_1d(data.get("doublings", [])):
        ax.axvline(float(d), color="C0", ls="--", lw=1)
        ax.text(float(d), ax.get_ylim()[1], f" {float(d):.3f}", color="C0",
                fontsize=8, va="top", ha="left")

    I = np.asarray(data["istar"], dtype=float)
    order = np.argsort(I)
    axl.axhline(0.0, color="0.5", lw=1)
    axl.plot(I[order], np.asarray(data["lam_max"])[order], "o-", ms=3, lw=1, color="C3")
    axl.set_ylabel(r"$\lambda_{\max}$")
    axl.set_xlabel(r"$I^*$")
    axl.grid(alpha=0.2)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def lyapunov_overlay(data: Any, out_path: str, bif: Any = None) -> str:
    """E2 figure: lambda_max(I*)±unc and median K(I*) over the bifurcation diagram.

    `data` keys: istar, lam_max, lam_unc, K. Optional `bif` (E1 data dict with
    scatter_istar/scatter_omega) draws the bifurcation backdrop. PLAN §5-E2.
    """
    istar = np.asarray(data["istar"], dtype=float)
    order = np.argsort(istar)
    I = istar[order]
    lam = np.asarray(data["lam_max"])[order]
    unc = np.asarray(data["lam_unc"])[order]
    K = np.asarray(data["K"])[order]

    n_panels = 3 if bif is not None else 2
    fig, axes = plt.subplots(n_panels, 1, figsize=(9, 2.4 * n_panels), sharex=True,
                             gridspec_kw={"hspace": 0.08})
    axes = np.atleast_1d(axes)
    row = 0

    if bif is not None:
        axes[row].scatter(bif["scatter_istar"], bif["scatter_omega"],
                          s=2, c="k", alpha=0.5, linewidths=0)
        axes[row].set_ylabel(r"$\omega$ at $\theta=0$")
        axes[row].set_title(r"E2 - $\lambda_{\max}$ and 0-1 test $K$ vs $I^*$ (chaos GATE)")
        axes[row].grid(alpha=0.2)
        row += 1

    ax = axes[row]
    ax.axhline(0.0, color="0.5", lw=1)
    chaotic = (lam - unc) > 0
    ax.errorbar(I, lam, yerr=unc, fmt="o-", ms=3, lw=1, color="C3", capsize=2)
    ax.fill_between(I, 0, lam, where=chaotic, color="C3", alpha=0.15)
    ax.set_ylabel(r"$\lambda_{\max}$")
    ax.grid(alpha=0.2)
    if bif is None:
        ax.set_title(r"E2 - $\lambda_{\max}$ and 0-1 test $K$ vs $I^*$ (chaos GATE)")
    row += 1

    ax = axes[row]
    ax.axhline(0.5, color="0.5", lw=1, ls="--")
    ax.plot(I, K, "s-", ms=3, lw=1, color="C0")
    ax.set_ylabel(r"median $K$")
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel(r"$I^*$")
    ax.grid(alpha=0.2)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path

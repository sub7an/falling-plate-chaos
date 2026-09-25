# Findings — chaos in the 2D quasi-steady falling-plate model

**AI-generated synthesis for author review.** This consolidates the results of
milestones M1–M8 into a thematic answer to the research question. The authoritative, chronological
record (with every non-obvious choice and correction) is `docs/decisions.md`; every number here is
reproducible via `python run_all.py`. **Findings are properties of the phenomenological model
(Andersen–Pesavento–Wang 2005, thin-card limit), not of real leaves.** Physics comes only from
`docs/model_spec.md`.

## Research question and answer

*Does this model have a chaotic regime as the dimensionless moment of inertia I\* varies; if so,
where, by what route, and how sensitive is landing displacement to release conditions?*

**Yes.** There is a chaotic band at **I\* ≈ [1.63, 2.83]**, confirmed independently by a positive
maximal Lyapunov exponent and the 0–1 test. The route into it is **not** a Feigenbaum
period-doubling cascade (only two doublings occur, then windowed chaos; no torus). Inside the
chaotic band, horizontal landing displacement is **exponentially sensitive** to release conditions,
growing at a rate equal (within estimation uncertainty) to the Lyapunov exponent; outside it the
displacement is insensitive. The model shows **no hysteresis or non-trivial multistability**.

## Regime sequence vs I\* (E1)

As I\* increases over [1.0, 3.0]: **fluttering** → (heteroclinic ≈1.2191) → **period-1 tumbling** →
(doubling ≈1.474) → **period-2** → (≈1.61) **chaos** → (≈2.85) **small-amplitude broadside
fluttering**. Figure: `figures/e1_regime_map.png`.

## Chaos, confirmed by two independent diagnostics (E2 — the gate)

- Maximal Lyapunov exponent (variational tangent-linear method, verified on the logistic map and
  Lorenz before use) is ≈0 in the periodic bands and rises to a plateau **0.12–0.14, peaking
  λ_max = +0.140 at I\* = 2.2** — consistent with the paper's stated 0.13 ± 0.01.
- The 0–1 test median K jumps 0 → ~1 over the **same** band and back.
- Both signals together (AND) confirm chaos in **I\* ≈ [1.63, 2.83]**; each vetoes the other's
  false positives (K blips at bifurcation points; sub-1e-3 numerical λ_max in periodic bands).
- Figure: `figures/e2_lyapunov.png`.

## Route to chaos — NOT Feigenbaum (E3, null for E4)

Only **two** successive period-doublings occur (period-1→2 at I\*≈1.474, 2→4 at I\*≈1.599); there is
no period-8. The period-2 band ends in a narrow (~0.01-wide) transition threaded with periodic
windows (period-3, -5, -6) and chaos. Power spectra show these window points have clean line spectra
and the "many-branch" states have a broadband floor — **no torus / no Ruelle–Takens quasiperiodic
route**. Because a Feigenbaum δ requires ≥3 successive doublings, **E4 (δ → 4.669) was deliberately
not run — a δ from two doublings would be overclaiming.** Figure: `figures/e3_route.png`.

## Landing-displacement sensitivity (E5)

Ensembles of falling cards from slightly perturbed release conditions (N=200, δ₀=1e-6):
- **Chaotic band:** σ(D) = std(Δx) grows exponentially (by 10⁶–10⁷ over D ≤ 200 chords).
  Measured in **time**, its growth rate d(log σ)/dt equals λ_max within ~10% (mixed sign) — the
  Lyapunov exponent sets the displacement-sensitivity growth rate. (The looser depth-based
  comparison to λ_max/U carries an extra depth↔time conversion via a single mean fall speed U; see
  the 2026-09-25 resolution in `decisions.md`.) Beyond the exponential phase σ keeps growing
  sub-linearly (pre-asymptotic; a clean √D diffusive asymptote is not isolated within D ≤ 200).
- **Regular bands:** σ(D) is flat within noise — displacement is insensitive (these are stable
  limit cycles; perturbations decay). This corrects the plan's "linear growth" expectation, which
  would require neutral dynamics.
- Figure: `figures/e5_sensitivity.png`.

## Multistability / hysteresis — none (E6, null)

Up and down I\* continuation sweeps coincide across [1.0, 3.0] (no hysteresis), and basin sampling
(random ICs) finds a single attractor at every diagnostic I\*, including inside periodic windows and
the chaotic band. The only multistability is the trivial **CW/CCW reflection-symmetry pair** (mirror
spin directions with identical statistics). This clears the "coexisting attractors read as chaos"
risk: the chaos is a single attractor. Figure: `figures/e6_multistability.png`.

## Discrepancy with the source paper (flagged, not forced)

At **I\* = 1.45** this model is **period-one** tumbling, not the paper's figure-read period-two; the
first doubling sits at **I\*≈1.474** here. Verified robust (constant section state to ~1e-12; three
ICs; no bistability). Reported as a property of the model, not reconciled by tuning.

## Validation and reproducibility

- Validation ladder (PLAN §4) passes: unit tests, limiting/ballistic cases, tolerance convergence,
  regime reproduction (5/6 anchors), and Lyapunov code verified on logistic (ln2, −0.916) + Lorenz
  (0.906) before any plate number was trusted. `pytest`: 59 passing.
- Every output is content-addressed cached with a provenance sidecar (git hash, Python/NumPy/SciPy
  versions, params, seed); the cache key includes a hash of the physics/numerics source files, so
  editing `model.py` invalidates dependent results.
- One command regenerates everything: `python run_all.py`.

## Open threads (optional, not required to answer the question)

- E5 late phase: confirm the √D diffusive asymptote at larger depth D (currently pre-asymptotic).
- E3 transition zone [1.598, 1.61]: finer characterization of the window structure (intermittency
  vs crisis) — the route is already established as non-Feigenbaum, non-torus.

# Decisions log

Running log of non-obvious choices. Newest first. Each entry: date, decision, rationale.

## 2026-09-25 — E5 slope-vs-λ_max gap RESOLVED (time-domain analysis) — post-M8 follow-up
Autonomous follow-up on the M7 open thread (the 6–16%, inconsistent-sign gap between the fitted
σ(D) exponential slope and λ_max/U). Decisive test: measure σ growth in TIME, not depth —
log σ_x(t) ≈ log σ₀ + λ_max·t, so d(log σ_x)/dt should equal λ_max DIRECTLY (no U conversion).
Added this to `experiments/e5_sensitivity.py` (stores σ_x(t) = std of ensemble x(t); reports
`slope_t` and `slope_t/λ_max`; cache bumped to analysis_version 2; E5 λ_max averaging time raised
to T=4000 for better convergence). Findings:
- **The depth-based gap conflated TWO removable estimation artifacts:** (i) the depth↔time
  conversion via a single mean U (U is depth/time-dependent, differs per I*), and (ii) λ_max
  under-convergence at T=2000 (λ_max(2.0) drifts 0.124→0.117→0.121→0.121 at T=2000/4000/6000/10000,
  i.e. ~5% finite-time scatter — the plate's variational λ_max is not tightly converged at these T).
- **Removing U (time-domain) makes the sign CONSISTENT and shrinks the gap:** slope_t vs λ_max ≈
  1.11 / 0.92 / 0.90 at I*=2.0/2.2/2.5 (≈10%, mixed sign) — versus the old depth-based +6/−13/−16%.
  The relationship is slope_t ≲ λ_max (an ensemble spread cannot exceed the maximal rate), as
  expected physically.
- **RESOLUTION:** the residual ~10% (mixed sign) is within the combined finite-time uncertainty of
  the two estimates being compared — slope_t is a finite-time (t≲90) fit and λ_max a finite-T
  (≈4000) average, both fluctuating in a chaotic flow. **The exponential displacement-sensitivity
  mechanism is confirmed and the rate equals λ_max within estimation uncertainty; the earlier
  apparent "gap" was dominated by the U-conversion and λ_max convergence, not a real physical
  offset.** This UPGRADES the M7 headline from "6–16% gap, not exact" to "rate = λ_max, confirmed
  directly in time; depth-based comparison is looser only because of the U conversion." (The M7
  entry's depth-based numbers stand as recorded; this entry is the resolution.)

## 2026-09-22 — M2 Lyapunov code verified on logistic + Lorenz (PASSED)
Rung 5 gate: the Lyapunov machinery reproduces textbook exponents on known systems
BEFORE touching the plate. Numbers vs targets:
- Logistic r=4.0: λ = **+0.6931** (target ln2 = 0.6931, exact to 4 dp).
- Logistic r=3.2: λ = **−0.9163** (period-2, <0). **CORRECTION (researcher-confirmed):** an
  earlier target of ~−0.65 was WRONG. The correct value is −0.9163, verified independently by
  the researcher and derived analytically here: the stable period-2 points of the logistic map
  are x₁,₂ = [(r+1) ± √((r+1)(r−3))]/(2r); at r=3.2 that is x₁≈0.7995, x₂≈0.5130, giving
  f'(x)=r(1−2x) → f'(x₁)≈−1.917, f'(x₂)≈−0.0832, so
  λ = ½·ln|f'(x₁)·f'(x₂)| = ½·ln(0.1595) = −0.9163. The test asserts this exact value.
- Logistic r=3.5: λ = **−0.8725** (period-4, <0).
- Lorenz variational spectrum = **(+0.9072, −0.0003, −14.574)**: λ_max≈0.906 ✓, zero
  exponent ✓, and **Σλ = −13.6667 = −(σ+1+β)** (constant divergence) — a structural check.
- Lorenz Benettin λ_max = **+0.9038** (independent, Jacobian-free) — agrees with variational.
Non-obvious choices:
- **`benettin_lmax` = two-trajectory augmented (reference+perturbed in ONE ODE, shared
  adaptive steps), Jacobian-free.** Rationale: a genuinely independent cross-check of the
  variational λ_max (no reuse of `jac`), while shared steps avoid the step-mismatch artefact
  of separately-integrated orbits (PLAN §2).
- **Added `map_lyapunov_spectrum` (discrete QR) to `lyapunov.py`** for the logistic check and
  future Poincaré return-map analysis; reduces to mean ln|f'| in 1-D.
- Lorenz run parameters (T=1000, dt_renorm=0.5, transient=50, rtol/atol=1e-9) give λ_max to
  ~1e-3 in ~2–4 s; kept as the verification config.

## 2026-09-22 — M8 / E6 multistability & hysteresis: NO non-trivial multistability (null)
Up/down I* continuation (seed each I* from the previous converged state) over [1.1, 3.0] step 0.05,
order parameter ⟨|ω|⟩ (gross rotation rate, ergodic per attractor); plus basin sampling (12 random
ICs) at diagnostic I* = 1.22 (flutter/tumble edge), 1.40 (tumble), 1.606 & 1.67 (periodic windows),
2.20 & 2.50 (chaos). Cached + sidecar; figure `figures/e6_multistability.png`. Results:
- **No hysteresis:** the up and down ⟨|ω|⟩(I*) branches coincide across the ENTIRE range [1.1, 3.0]
  (small differences only in the chaotic band = finite-time chaotic fluctuations of one attractor).
  The flutter/tumble transition (~1.2191) is NOT hysteretic (flutter ≤1.215, tumble ≥1.220 both ways).
- **No non-trivial multistability:** at every diagnostic I*, all 12 random ICs converge to a SINGLE
  ⟨|ω|⟩ cluster (ranges ≤0.02 wide) — one attractor, including inside the periodic windows and the
  chaotic band. No coexisting distinct attractors found.
- **Only the trivial CW/CCW symmetry pair:** ⟨ω⟩ takes both signs across ICs at every tumbling/chaotic
  I* (same |⟨ω⟩| magnitude, opposite sign) — the expected reflection-symmetry mirror pair (a card has
  no preferred spin direction), not a distinct dynamical state.
**This completes the milestone chain and is consistent throughout:** the ~1.474 doubling had no
bistability (M5-pre), and E6 finds no hysteresis or coexistence anywhere. The model is effectively
monostable (up to the ± spin symmetry) across 1.1–3.0. Valid null result (PLAN §5-E6, §8 guard against
"coexisting attractors read as chaos" — cleared: the chaos is a single attractor, confirmed by E2 λ_max
+ K and now E6 basin sampling).
- New code (M8): `plotting.hysteresis_diagram`, `experiments/e6_multistability.py`, `config/e6.yaml`.

## 2026-09-22 — M7 / E5 landing-displacement sensitivity: Lyapunov→displacement link CONFIRMED
IC-perturbation ensembles (N=200, δ0=1e-6) of falling cards via the 6-D `rhs_full` (position
quadrature, spec §2); σ(D)=std(Δx) at fixed depth D∈[2,200] chords; chaotic band I*∈{2.0,2.2,2.5}
vs regular controls {1.1,1.4}. λ_max via the verified variational method; U = mean descent speed.
Cached + sidecar; figure `figures/e5_sensitivity.png`. Grazing fraction = 0 at all D (depth not
too shallow). Results:
| I*  | regime  | σ growth ×   | exp slope [95% boot CI] | λ_max/U_ind | agree | late phase |
|-----|---------|--------------|-------------------------|-------------|-------|------------|
| 1.10| regular | 2.9e1 (noise)| — (flat)                | ~0          | —     | flat       |
| 1.40| regular | 4.9e1 (noise)| — (flat)                | ~0          | —     | flat       |
| 2.00| chaotic | 3.3e6        | **0.238** [0.227,0.247] | 0.225       | 5.7%  | sub-linear |
| 2.20| chaotic | 4.6e6        | **0.223** [0.216,0.224] | 0.257       | 13.2% | sub-linear |
| 2.50| chaotic | 1.2e7        | **0.173** [0.170,0.174] | 0.207       | 16.4% | sub-linear |
- **KEY RESULT (stated precisely after the bootstrap): σ(D) grows exponentially in the chaotic band
  at a rate of the same magnitude as λ_max/U, confirming the Lyapunov→displacement MECHANISM
  log σ ≈ log σ₀ + (λ_max/U)·D — but the fitted slope and λ_max/U are NOT equal, and the difference
  is real (beyond sampling noise) with an INCONSISTENT sign.** N=200 bootstrap (300 resamples): the
  slope CIs are tight (width ~0.01–0.02) and λ_max/U_ind falls OUTSIDE every CI, so this is not
  estimation noise. But the sign is not systematic in one direction:
    * I*=2.0: slope 0.238 [0.227,0.247] vs λ/U 0.225 → slope **ABOVE** by ~6%
    * I*=2.2: slope 0.223 [0.216,0.224] vs λ/U 0.257 → slope **BELOW** by ~13%
    * I*=2.5: slope 0.173 [0.170,0.174] vs λ/U 0.207 → slope **BELOW** by ~16%
  Slope is high at 2.0, low at 2.2/2.5, and the gap GROWS with I*. Because the direction flips, this
  is NOT a single-direction physical bias; it reflects the leading-order nature of the relation plus
  estimation subtleties, most plausibly: (a) **U is one whole-trajectory mean descent speed, but the
  exponential phase occupies early/shallow depths where the descent speed differs** (the time↔depth
  conversion is only approximate and depth-dependent); (b) **λ_max is a T=2000 estimate**, near but
  not perfectly converged (e.g. λ_max(2.2)≈0.141 at T=2000 vs ≈0.140 at T=6000); (c) the **finite-
  depth exp-fit window measures a finite-time Lyapunov rate**, not the asymptotic exponent. **Honest
  conclusion: the exponential-growth mechanism is confirmed and the rate is O(λ_max/U), but a precise
  slope = λ_max/U equality is NOT established — the residual is 6–16%, real, and not unidirectional.**
  Do NOT round this to "closely matches."
- **U provenance (Check 2):** the experiment's U came from the ensemble's own reference member
  (T=350) — NOT independent. Cross-checked here against an independent long run (T=2000): U differs
  by 0.6% / 2.0% / 6.5% (I*=2.0/2.2/2.5), so U is robust and the comparison is not materially
  circular. The table's λ_max/U uses the INDEPENDENT U.
- **Regular bands are FLAT within noise, NOT linear (Check 1 — CORRECTS the PLAN prediction).**
  σ(D) at I*=1.1, 1.4 has linear slope ~2.5e-8/4.7e-8 and corr(D,σ)≈0.04 — no growth trend. PLAN
  §5-E5 predicted ~linear growth in regular bands, but that assumes NEUTRAL dynamics; these regular
  regimes are **stable (attracting) limit cycles**, so IC perturbations DECAY onto the cycle and
  displacement is insensitive (bounded), not linearly growing. So: chaotic = exponential
  sensitivity; regular = flat/insensitive. (Growth-ratio split: ~30 noise vs ~1e6–1e7; threshold
  1e3, a fit param not in the cache key.)
- **Late phase (D beyond the exponential):** σ keeps growing SUB-linearly — NOT saturating
  (displacement is an unbounded running sum, as PLAN says). The automated log-log slope (0.9–1.5)
  is contaminated by the exponential tail near the phase cut; the large-D (D>100) behaviour is
  gentler, a √D–linear crossover consistent with the expected diffusive phase, but a clean √D
  asymptote is NOT isolated within D≤200 (would need deeper D). Reported honestly as pre-asymptotic.
- New code (M7): `model.lab_velocity`/`model.rhs_full` (spec §2 quadrature), `displacement.py`
  (landing_displacement + grazing guard, spread_vs_depth, fit_growth_models),
  `plotting.sensitivity_fits`, `experiments/e5_sensitivity.py`, `config/e5.yaml`,
  `tests/test_displacement.py`.

## 2026-09-22 — E4 DECISION: STATED NULL RESULT (do not run Feigenbaum)
**Explicit null result.** The Feigenbaum-scaling hypothesis (PLAN E4; δ_n → 4.669) requires ≥3
successive period-doublings. This model exhibits only **2** (period-1→2 at I*≈1.474, period-2→4 at
I*≈1.5985); there is no period-8, and the period-4 band gives way directly to windowed chaos.
**E4 is therefore not run.** Reporting a Feigenbaum δ from 2 doublings (a single ratio, no
convergence) would be overclaiming; per the project operating contract a null result is the correct, valid outcome. This
supersedes the PLAN's conditional E4/M6 Feigenbaum task.

## 2026-09-22 — M6 / E3 refinement (researcher-requested): windowed chaos, NOT a torus
Two rigor checks on the [1.598, 1.61] transition zone.
1. **Scan-density check.** The E3 grid (step 0.002) placed exactly **6** points in [1.598, 1.608]
   — borderline. A fine local scan (step **0.0005**, T=10000, discard 5000) confirms the structure
   is REAL and coherent, not sparsity: a clean **period-4 band [1.5985, 1.6005]** (λ_max≈0), a mixed
   chaotic strip [1.601, 1.604], a clean **period-3 window [1.6045, 1.607]** (λ_max≈0), period-5 at
   1.6075, then chaos from ~1.608.
2. **Torus check via power spectrum (CORRECTS earlier "quasiperiodic/torus" wording).** PSD of ω(t):
   the window points show clean **line spectra** — 1.55 period-2, 1.600 period-4, 1.606 period-3 all
   have sharp harmonics of a single fundamental on a ~1e-10 floor (commensurate; peak ratios ≈ 2,3,4).
   The many-branch "torus candidate" at 1.608 shows a **RAISED broadband floor** matching the 1.615
   chaos control — NOT two incommensurate sharp lines. **So there is NO torus / no Ruelle–Takens
   quasiperiodic route.** The "13–15 section points" states were weakly-chaotic (or unsettled), not
   quasiperiodic. Corrected reading: the route is **truncated period-doubling (1→2→4) followed by a
   chaotic transition threaded with periodic windows (period-3, -5, -6)** — the period-3 window is
   itself a marker that chaos is fully present (Li–Yorke/Sharkovskii).
[The word "torus-like" in the entry below is superseded by this check: those states are chaotic, not tori.]

## 2026-09-22 — M6 / E3 route to chaos: NOT a Feigenbaum cascade → E4 NOT warranted (null)
Fine θ=π event-section bifurcation zoom over I*∈[1.45,1.62] (86 pts) + λ_max, counting SUCCESSIVE
period-doublings. Cached + sidecar; figure `figures/e3_route.png`. Period bands:
| I* band          | section branches | reading                         |
|------------------|------------------|---------------------------------|
| [1.450, 1.474]   | 1                | period-1 tumbling               |
| [1.476, 1.596]   | 2                | period-2 (WIDE band, one doubling) |
| [1.598, 1.608]   | 15 / 4 / 6 / 3 / smear | intricate transition zone  |
| [1.610, 1.620]   | aperiodic        | chaos                           |
- **Successive clean doublings = 2 only: 1→2 at I*≈1.476, 2→4 at I*≈1.600.** There is NO clean
  4→8→16 continuation; instead the period-2 band ends in a NARROW (~0.012-wide) intricate zone with
  a period-4 point (1.600), periodic windows (period-3 at 1.606, period-6 at 1.602), torus-like
  many-branch states (≈13–15 distinct section points at 1.598/1.608 with λ_max≈0 and K≈0 — the
  quasiperiodic signature), and chaos — before the fully chaotic band from ~1.61.
- **λ_max** stays ≈0 across [1.45, 1.59] (periodic), then shows the spike-and-dip structure of an
  intricate transition (spikes at 1.598/1.604, dips at the 1.600/1.606 windows) climbing to ~0.08
  by 1.62 — NOT the smooth onset of a period-doubling accumulation.
- **DECISION (E4 gate): the Feigenbaum precondition (≥3 successive period-doublings) is NOT met
  (only 2). E4 Feigenbaum-δ work is NOT warranted — this is a valid NULL result** (operating contract:
  "Feigenbaum scaling is a hypothesis to TEST, not assume; null results are valid"). The route is a
  short/incomplete doubling (1→2→4) followed by a mixed torus/periodic-window/intermittency
  breakdown to chaos, consistent with the M5-pre finding that the first doubling is non-cascade-like.
  **Do NOT run E4 as a δ→4.669 measurement.** If the route is pursued further it should be
  characterized as (candidate) quasiperiodic/Ruelle-Takens + windows, via return maps + power
  spectra in the narrow [1.598,1.61] zone — but that is optional, not the planned E4.
- Minor feature noted for completeness: a small branch rearrangement / λ_max bump near I*≈1.557
  inside the period-2 band (branch count stays 2); not pursued.

## 2026-09-22 — Pre-M6 checks: chaos onset + K-blip audit (refines M5)
Two researcher-requested checks before E3, reusing the verified `e2_task` (λ_max + K).
**Check 1 — lower chaos-band edge (fine scan, step 0.01):** λ_max≈0 & K≈0 for I*≤1.60; chaos
turns on SHARPLY at **I*≈1.61** (λ_max +0.032, K 0.91), a bit below the E2 coarse-grid 1.63
(grid resolution). Confirmed-chaos is continuous from 1.61 EXCEPT a **periodic window at I*≈1.67**
(λ_max +0.004, K −0.03) embedded in the chaotic band — a feature E3 must characterize.
**Check 2 — K-blip audit across periodic bands (step 0.02):** elevated-K (K>0.3) points found at
**I*=1.46** (near the 1.474 doubling) AND **I*=1.28** — the latter is NOT at a known bifurcation.
**This CORRECTS the M5 note** that blips sat "only at the two bifurcation points (~1.22, ~1.46)":
the lower blip is at 1.28, not the heteroclinic. Root cause established: these are **0-1 test
transient/settling artifacts**, not intrinsic. Verified at I*=1.28 (period-1 tumbling): K = +0.445
(T=6000, discard 1000) → **−0.018 (T=12000, discard 3000) → −0.309 (T=20000, discard 6000)**;
neighbours 1.26/1.30 are negative throughout. The 0-1 test is biased high by incompletely-decayed
transients / long-period orbits near the slow dynamics of the heteroclinic (1.2191) and the doubling.
**All such blips are sub-threshold (K<0.5) with λ_max≈0, so the E2 gate (λ_max>0 AND K>0.5) was
never at risk; the chaotic band is unchanged.** Action for E3/E5: discard a LONGER transient for the
0-1 test in the near-heteroclinic region, or lean on λ_max there.

## 2026-09-22 — M5 / E2 Lyapunov chaos GATE: PASS
36-point I* grid [1.0, 3.0]. Per I*: λ_max via the verified variational spectrum (M2, analytic
Jacobian, k=1) with uncertainty = |λ(T=2500) − λ(T=1250)|; and independent median K from the
0-1 test (Gottwald–Melbourne, Vosc-corrected, 100 random c) on the θ=0 Poincaré-section ω series
(stroboscopic, not the fine step). ~64 s wall (parallel). Cached + sidecar; figure
`figures/e2_lyapunov.png`. The 0-1 test itself was first verified on the logistic map
(K≈1.0 at r=4.0/3.8; K≈0 at r=3.2/3.5) — `tests/test_chaos_tests.py`.
**Result — chaos CONFIRMED by BOTH methods in I* ≈ [1.63, 2.83]:**
- λ_max ≈ 0 (|·|<1e-3) for I* < 1.6 and > 2.85 (periodic bands); rises to a plateau ~0.12–0.14
  across [1.63, 2.83], **peak λ_max = +0.140 at I* = 2.2** (matches M3 and the paper's 0.13±0.01
  upper edge; earlier variational/Benettin agreement established it as physical, not numerical).
- median K jumps 0 → ~1.0 over the SAME band and back — independent confirmation.
- **GATE PASS:** 22 grid points satisfy (λ_max − unc > 0) AND (K > 0.5), band I* ≈ [1.63, 2.83].
- **The AND of two methods earned its keep:** median K showed spurious ~0.4–0.5 blips at I*≈1.22
  (flutter/tumble heteroclinic) and ≈1.46 (period-doubling) — bifurcation points where 0-1 test is
  unreliable — but λ_max≈0 there vetoed them. Likewise a couple of periodic points had tiny
  numerical λ_max~1e-3 that K≈0 vetoed. Neither alone would be trusted; together they are clean.
- **Consistency with E1:** the confirmed chaotic band [1.63, 2.83] essentially equals E1's coarse
  "aperiodic (chaos-candidate)" band [1.62, 2.82] — so that flag was genuine chaos, and it extends
  slightly below the paper's ~1.8 (real early chaos just above the period-2/mixture band), while
  the upper edge ~2.83 matches the paper's ~2.8. Return to (broadside) flutter above ~2.85.
**Branch decision (PLAN §5-E2):** GATE PASS → a chaotic regime exists → proceed to E3 (route to
chaos) / E4 (period-doubling + Feigenbaum TEST). Carry forward: (i) the ~1.474 first doubling is
supercritical-consistent, no bistability (use event-based section + Floquet in E4); (ii) do NOT
assume a Feigenbaum cascade — E3 must first establish ≥3 successive doublings before E4 measures δ.

## 2026-09-22 — Coexistence/bistability test at the I*≈1.4743 transition (pre-M5)
Researcher-requested check: is the "−0.79 branch" a real coexisting attractor (bistability ⇒
subcritical) or a continuation-tracking artifact? Method: extract a full state ON the −0.79 branch
of the 1.475 period-2 orbit, then seed integrations from it (and from the +0.6-branch state) at I*
BELOW the transition, with long settling (T=6000, discard 5000), and compare to the spec-IC orbit.
Results:
- At I*=1.475 the −0.79-branch seed reproduces the **same** period-2 orbit (branches +0.601/−0.790).
  So the −0.79 branch is **real**, not a tracking artifact.
- BELOW the transition (I*=1.474, 1.470, 1.465, 1.460): spec IC, −0.79 seed, and +0.6 seed **all
  converge to the SAME unique period-1 tumbling orbit** (single section branch +0.599…+0.624,
  net_rev 44–68, R 0.89–0.95). **No coexisting period-2. No bistability. No hysteresis** (window
  < 0.001; up- and down-transition both at ~1.4743).
**Conclusions:**
- **CORRECTS the earlier "subcritical/coexistence" read (M3 probe).** The transition is a
  single-attractor period-1 → period-2 bifurcation with NO hysteresis — consistent with a
  **supercritical** doubling (Floquet → −1, to be confirmed in E4), not subcritical.
- The abrupt "−0.79 branch" is a **section-observable artifact**: the period-2 (net-≈0, mixture)
  orbit crosses θ=π twice per cycle — once forward (+0.6), once backward (−0.79) — so a
  full-amplitude second crossing appears by topology, not a √-law amplitude split.
- **E6 does NOT need I*≈1.4743 as a multistability case study** — no coexisting attractors there.
  (E6 should still probe the flutter↔tumble and chaotic bands for multistability.)

## 2026-09-22 — M4 / E1 regime map (first sweep): ordered bands, hypothesis CONFIRMED
41-point coarse I* grid in [1.0, 3.0] + one adaptive refinement pass near label changes
(→45 points), parallel (joblib, all cores), ~19 s wall. Per I*: spec IC (0,0.01,0.5,0),
DOP853 rtol 1e-10, T=1200, discard 300; regime via `classify_regime`; bifurcation diagram
from the θ=0 Poincaré section (event g=sin(θ/2), crossed by BOTH flutter and tumble).
Cached to `data/<key>.npz` + provenance sidecar; figure `figures/e1_regime_map.png`. **No
Lyapunov here** (E2/M5). Bands found:
| I* band        | regime                          |
|----------------|---------------------------------|
| [1.00, 1.20]   | fluttering                      |
| [1.22, 1.45]   | period-1 tumbling               |
| [1.48, 1.60]   | flutter/tumble mixture (period-2 orbit) |
| [1.62, 2.82]   | aperiodic (chaos-candidate)     |
| [2.85, 3.00]   | fluttering (small-amp broadside)|
- **E1 hypothesis (ordered regime bands + ≥1 aperiodic band) is CONFIRMED**, and the ordering
  matches the APW picture and the M3 anchors: flutter → tumble (onset ~heteroclinic 1.2191) →
  period-2 → chaos → broadside flutter.
- **Nuance 1 — the period-2 band (1.48–1.60) is labelled "flutter/tumble mixture," not
  "period-2 tumbling."** Reason: the period-2 orbit born at the ~1.474 transition alternates a
  forward tumble (ω≈+0.6 at section) with a BACKWARD swing (ω≈−0.79), so rotation directedness
  R = |⟨ω⟩|/⟨|ω|⟩ falls below the tumbling threshold. This is the SAME mixed forward/backward
  orbit seen in the M3 probe, and is consistent with that transition being subcritical/coexistence
  rather than a clean supercritical doubling. The paper's "period-two tumbling @1.45" is not
  reproduced as pure tumbling here.
- **Nuance 2 — chaos-candidate band is wider on the low end (1.62) than the paper's ~1.8–2.8.**
  E1's "aperiodic" is only a COARSE flag (θ=0 section period = −1 from the return map); part of
  1.6–1.8 may be high-period or quasiperiodic rather than chaotic. **E2 (λ_max + 0–1 test) must
  confirm the true chaotic sub-band** before any chaos claim — this is exactly the E2 GATE.
- New code (M4): `io_cache` (content-addressed cache + sidecar), `sweeps` (joblib driver +
  transition refinement), `plotting.bifurcation_diagram` (diagram + regime strip), and
  `experiments/e1_regime_map.py`. Fast plumbing tests in `tests/test_pipeline.py`.

## 2026-09-22 — M3 regime reproduction (rung 4): 5/6 anchors, ONE discrepancy
Standard thin-card config (C_T=1.2, C_R=π, A=1.4, B=1.0, μ1=μ2=0.2), spec §7 tilted IC
(0, 0.01, 0.5, 0), DOP853 rtol 1e-11. λ_max via variational (k=1, analytic Jacobian),
cross-checked by Benettin. Result vs spec §7:
| I*   | spec anchor                       | model result                     | λ_max        |
|------|-----------------------------------|----------------------------------|--------------|
| 1.1  | fluttering                        | fluttering ✓                     | +0.0006 (~0) |
| 1.4  | period-one tumbling               | period-one tumbling ✓            | −0.0000 (~0) |
| 1.45 | period-two tumbling               | **period-ONE tumbling** ✗        | +0.0002 (~0) |
| 1.6  | periodic flutter/tumble mixture   | periodic mixture ✓               | +0.0004 (~0) |
| 2.2  | chaotic, λ_max=0.13±0.01          | chaotic ✓                        | **+0.140 var / +0.137 Benettin** |
| 3.0  | small-amplitude broadside flutter | fluttering, smallest amplitude ✓ | +0.0001 (~0) |
- **λ_max**: all five periodic regimes give |λ_max|<1e-3 (limit cycles, as expected). At
  I*=2.2 variational gives 0.1414/0.1409/0.1402 (T=2000/4000/6000) and Benettin 0.1367 — two
  INDEPENDENT methods agree → physical (not numerical) chaos. Value sits at the upper edge of
  the paper's 0.13±0.01, **matching the spec's own throwaway reproduction (0.138–0.139, §7)**.
- **DISCREPANCY at I*=1.45**: this model gives period-ONE tumbling, not the paper's Fig-3
  read-off of period-two. Verified robust: the section state (vx',vy',ω) and revolution time
  (13.3587) are constant to 5 dp across revolutions, and **three different ICs all give
  period-one** (no coexisting period-2 branch found — not multistability). So the
  period-doubling lies ABOVE 1.45 in this model/config, not between 1.4 and 1.45. Consistent
  with the operating contract ("Feigenbaum scaling is a hypothesis to TEST, null results valid"). Precisely
  locating the first doubling (Floquet multiplier → −1) is **E3/E4 work, deliberately NOT done
  here** (no sweep at M3). Flagged for those milestones.
  - **Probe (researcher-requested), event-based locator (E4-grade):** using `solve_ivp` events
    on g(θ)=cos(θ/2) (zero exactly at θ≡π mod 2π), section values are resolved to ~1e-12. Fine
    scan: period-ONE (single value, ~1e-12 scatter) through I*=1.47…**1.4742**; clean period-TWO
    from **I*=1.4744** onward (two branches each constant to ~2e-12), persisting cleanly through
    1.5. **Provisional first period-1→period-2 transition: I* ≈ 1.4743 (bracketed in
    (1.4742, 1.4744)).** The transition is SINGLE and CLEAN — no period-4 sliver, no chaos, no
    non-monotonicity in the gap (checked at 1.471–1.475 and subdivided 1.4742–1.475).
  - **[CORRECTED — see the "coexistence/bistability test" entry above; NOT subcritical.]** I first
    read the abrupt ω-branch jump (0 → 1.39 across ΔI*=0.0002) as a subcritical doubling / jump to a
    coexisting period-2 (expecting hysteresis). The bistability test **refuted** that: no coexistence,
    no hysteresis. The abrupt "−0.79 branch" is a **section-observable artifact** — the period-2
    orbit crosses the θ=π section twice per cycle (forward ω≈+0.6, backward ω≈−0.79), so the doubled
    orbit gains a full-amplitude second crossing by topology, not a √-law amplitude split. Lesson:
    do not infer bifurcation criticality from a single scalar section observable. E4 still owns the
    Floquet classification (multiplier → −1) and whether a full cascade (→ chaos ~1.6) follows.
  - The ~1e-6 scatter in the 1.45 section values is a **linear-interpolation artifact of the
    section locator** (t_eval-grid sampling), NOT an unresolved period-2 and NOT dynamics: it has
    no period-2 alternation, does not decay with longer settling, and scales as **O(dt²)** —
    5.1e-6 → 1.1e-6 → 2.3e-7 → 2.6e-8 as dt = 0.02 → 0.01 → 0.005 → 0.0025 (period-one to ~1e-8
    at dt=0.0025). **Classifier tolerance (1e-3 rel) is fine; no change needed.** For E4, replace
    the linear-interp section locator with **event-based root-finding** (solve_ivp events / dense
    output + Brent) so the return map is accurate to integrator precision, not O(dt²).
- **Poincaré section choice (deviation from spec §10):** spec §10's default section is upward
  crossings of ω=0. A one-signed tumbling orbit never crosses ω=0, so for period counting we
  section instead on the **cyclic angle θ = π (mod 2π)**, recording (vx',vy',ω) at each
  crossing. `poincare_map` detects crossings in EITHER rotation direction (a first pass only
  caught increasing-θ, giving spurious "no crossings" for clockwise tumbling — fixed).
- **Classifier discriminator = rotation directedness R = |⟨ω⟩| / ⟨|ω|⟩** (fraction of rotation
  that does not cancel): R≈1 tumbling, R≈0 fluttering, intermediate → mixture; a low-R orbit
  whose θ swings past 2π (full flips that cancel) is a mixture, not flutter. Chaotic set by
  λ_max>0. Thresholds are function arguments, not magic numbers.

## 2026-09-22 — Injected message during M2 (noted, refused)
During M2 an out-of-band message attempting to disable tooling safeguards was injected as if from the
author. The author confirmed they did NOT send it. No action was taken on it (correct refusal).
Recorded here for provenance; nothing further required.

## 2026-09-22 — M1 model implementation + validation rungs 1–3 (PASSED)
Anchors 3 (flutter 1.2190 / tumble 1.2192) and 4 (edge-on saddle eigenvalues) reproduced
and independently confirmed by the researcher.
- **Analytic Jacobian in `model.py`** (not finite-difference-only). Hand-derived from
  (5.1)–(5.3)/(4.7)–(4.9) and **cross-checked against a central finite difference** over 25
  random non-zero-speed states (`test_jacobian_matches_finite_difference`). Rationale:
  de-risks the variational Lyapunov integrator + monodromy (M2/E4); FD cross-check catches
  derivation errors, and eigenvalues match anchor 4 at I*_C.
- **Rung-3 tolerance convergence uses a regular regime (I*=1.1 flutter) over a short horizon
  (T=30).** Rationale: in chaotic/heteroclinic regimes exponential separation would masquerade
  as non-convergence; a regular orbit isolates the integrator's own convergence. The observed-
  order study at the non-smooth `|ω|ω` kink (spec §9) is **deferred to the Lyapunov milestone**,
  where step-size artefacts in λ_max are the concern.
- **Ballistic free-fall test (rung 2) uses ONLY the decoupled ICs**, not general ICs:
  broadside `(vx'=0, ω=0, θ=0)` fixed → `vy'(t)=vy0−t/(I*+1)`; edge-on `(vy'=0, ω=0, θ=π/2)`
  fixed → `vx'(t)=vx0−t/I*`. Rationale: the `−vx'vy'` term in (5.3) drives ω (hence θ) whenever
  `vx'vy'≠0`, so the simple closed form holds only when `vx'vy'=0` is maintained (one velocity
  component zero). A general zero-coefficient IC does NOT fall in a straight line. Restriction
  is noted in both ballistic test docstrings.

## 2026-09-22 — Spec second-reader verification (PASSED)
Cross-checked `docs/model_spec.md` v0.1 against the [P91] PDF (rendered pages 98–102)
per the spec's own §11 checklist. **Result: every checked item matches; no discrepancies.**
- Eqs (5.1)–(5.3), closures (4.7)–(4.9), general-β (4.4)–(4.6): match exactly (journal p.98).
- Body-frame **primes** on vx′, vy′ (and x′, y′): confirmed present in the paper — resolves the
  §2/§11 "primes lost in extraction" flag.
- **Exact form of (4.8)** confirmed: (1/π)[A − B(vx′²−vy′²)/(vx′²+vy′²)]√(vx′²+vy′²)(vx′,vy′).
- Parameters (§5.2) confirmed: C_T=1.2, C_R=π, A=1.4, B=1.0, μ1=μ2=0.2.
- Fixed points (6.1)/(6.2) confirmed; V≈2.8025, W≈1.1441. ω≡θ̇ confirmed.
- Validation anchors 1–7 confirmed: regimes (Fig 3), λ_max=0.13±0.01 at I*=2.2, chaos for
  I*≈1.8–2.8 (Fig 4b), period-doubling between 1.4 and 1.45, I*_C≈1.2191, ~25 card widths,
  flutter 1.2190 / tumble 1.2192, saddle eigenvalues (−0.5854, +0.3813, −0.9861±2.5949i),
  period law (6.4) with T0=4.3629, noise-sensitivity near I*_C.
**Items for the researcher (spec is read-only for the AI tooling; these are suggestions only):**
1. Hopf curve (6.3), left blank in the spec, reads in the PDF:
   μ1 = (1/4)·√(3/(5π))·(2I*+1)/(I*−1). Not needed for the RHS; optional to add.
2. `references/` PDF is named `..._b.pdf` but its content is [P91] (2005a, "Analysis of
   transitions", pp.91–104) — the correct primary source, just a misleading suffix. The
   distinct [P65]/2005b ("Unsteady aerodynamics", pp.65–90) is NOT in the repo (spec §0 says
   it is not used for the equations). Suggest renaming the file to avoid confusion.
3. Period law (6.4): the PDF glyph is "log"; the spec writes "ln". Natural log is the correct
   reading (required by the 2/λu prefactor theory, Gaspard 1990). No change needed; noted.
**Gate status:** physics in the spec is verified → M1 (implementing model.py) is unblocked.

## 2026-09-22 — Scaffold (M0)
- **Physics source of truth = `docs/model_spec.md` (pending).** During planning it is a
  pending dependency; once implementation starts, a missing/ambiguous spec is a STOP.
  Rationale: researcher must defend every equation; nothing from memory.
- **Integrator: `solve_ivp` DOP853 (production), RK45 (coarse).** Smooth non-stiff ODE;
  high order for the accuracy Lyapunov needs. Tolerances are starting values, to be fixed
  by the §4.3 convergence study.
- **Lyapunov: variational (tangent-linear) PRIMARY; Benettin (augmented, shared steps) +
  0–1 test as independent cross-checks.** Shared-step augmentation avoids step-mismatch
  artefacts; two independent methods triangulate physical vs numerical chaos.
- **0–1 test sampled on the Poincaré section / at ~dominant period; report median K over
  many random c.** Avoids oversampling autocorrelated data (biases K high).
- **Cache key = param hash + hash of physics/numerics source files; git hash in sidecar
  only.** Cache invalidates on physics/numerics edits, not unrelated commits.
- **Provenance sidecar records git hash, Python version, SciPy version, params, seed.**
- **E5 growth model: no saturation expected.** Displacement is an unbounded running sum →
  exponential (Lyapunov) phase, then diffusive ~√D (chaotic) or ~linear (regular). Fit all
  three and compare; perturbation chosen so the exponential phase spans ≥1 decade.
- **E4 bifurcations located by Newton shooting + Floquet multiplier = −1**, not by eye.
  Tumbling orbits: wrapped angle in the residual, record winding number, phase fixed by the
  Poincaré section.
- **Venv + pinned `requirements.txt` from the actual install** (Python 3.13.5): numpy 2.5.3,
  scipy 1.18.1, matplotlib 3.11.2, PyYAML 6.0.3, joblib 1.6.0, pandas 3.0.6, pyarrow 25.0.1,
  pytest 9.1.1.

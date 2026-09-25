# PLAN — Chaos in the 2D quasi-steady falling-plate model

Model: Andersen, Pesavento & Wang (2005a, *J. Fluid Mech.* **541**, 91–104;
2005b, **541**, 65–90). **Equations & coefficients come ONLY from
`docs/model_spec.md`** (pending; a separate deliverable). Nothing here fills
physics from memory. Findings are properties of *the model*, not real leaves.

Research question: does this model have a chaotic regime as the dimensionless
moment of inertia **I\*** varies? If so, where, by what route, and how
sensitive is horizontal displacement at fixed fall depth to tiny changes in
release conditions? Hypothesis to **test, not assume**: transitions show
Feigenbaum-type scaling — applies only if a period-doubling cascade exists.
Null results are valid.

## Environment
macOS 26.0 (arm64), Apple M4, 10 cores, 16 GB RAM. Python 3.13.5 in project
`.venv`. Parallel I\* sweeps over the 10 cores.

## Working assumptions (correct any that are wrong)
- **A1.** RHS is translation-invariant (no explicit x, y); (x, y) are quadratures,
  not part of the attractor. Spec confirms.
- **A2.** Length unit = chord `c`; consistent time/mass scales; **I\*** is a single
  scalar control; other groups (Re, aspect ratio, density ratio) fixed per experiment.
- **A3.** RHS is C¹ over the operating range (drag `|·|` terms continuous); verified
  numerically before Lyapunov work.
- **A4.** "Displacement at fixed depth D" = horizontal travel `Δx` from release to the
  first transversal downward crossing of vertical drop `D`, `D` in chords.
- **A5.** "Chaos" claims require positive Lyapunov exponent **and** an independent
  confirmation.

## 1. Definitions

**State & reduced space.** Full state `s = (x, y, θ, v…, θ̇)` (exact ordering/velocity
frame per spec). Reduced dynamical state (drops x, y):
`u = (velocity components, θ mod 2π, θ̇)` — dim 4 (2 velocity comps). **θ is cyclic**:
all distances use the wrapped angular difference in (−π, π]; tangent space treats θ on
a circle.

**Regime classifier** (post-transient):
- **Steady/gliding:** `θ̇ → const ≈ 0`, θ asymptotes to fixed orientation; no θ̇ sign changes.
- **Fluttering:** bounded oscillation of θ about a mean, mean rotation `Ω = ⟨θ̇⟩ ≈ 0`,
  θ̇ changes sign; net revolutions = 0.
- **Tumbling:** unwrapped θ grows monotonically on average, `|Ω| > 0`, θ̇ keeps one sign
  on average.
- **Chaotic:** bounded, aperiodic, `λ_max > 0`, plus independent confirmation; may be
  chaotic-fluttering or chaotic-tumbling.
- **Discriminators:** `Ω`; Poincaré return-point count; sign of `λ_max`; 0–1 test `K`.

**Poincaré section:** crossings of a spec-defined surface (default `θ̇ = 0` with fixed
direction, or body chord horizontal). Return-point count → periodicity.

**Landing displacement.** Event `y = y₀ − D`, taken only on a **transversal downward**
crossing (`ẏ < 0`). **Grazing guard:** if `|ẏ|` at crossing is below threshold, flag as
grazing, take the next genuine crossing; **report the flagged fraction per D — a
non-negligible fraction means D is too shallow.** **Cross-check:** also report
**displacement at fixed time** `Δx(T)` to confirm the fixed-depth measure is not an
event-localization artefact. Report `Δx` in chords for several D (final set chosen after
E2, see §5-E5).

**Transition criteria vs I\*:** change in classifier label; doubling of Poincaré points;
`λ_max` crossing zero.

## 2. Numerics
- **Integrator:** `scipy.integrate.solve_ivp`, **DOP853** (dense output) for
  production/Lyapunov; **RK45** for coarse scans. Smooth non-stiff ODE → adaptive
  high-order RK. Stiffness spot-check (Radau) in validation.
- **Tolerances (starting values, fixed by §4.3 convergence):** production
  `rtol=1e-10, atol=1e-12`; coarse `rtol=1e-8, atol=1e-10`.
- **Perturbation size:** Benettin `δ₀ ≈ 1e-7`–`1e-6` in normalized reduced state (3–5
  orders above numerical noise, linear regime); verified by a `δ₀`-plateau check.
- **Event detection:** depth event = direction −1 **and** `|ẏ|` grazing threshold;
  Poincaré event directional; dense output for sub-step localization; fixed-time `Δx(T)`
  from dense output at prescribed T.
- **Lyapunov — variational is PRIMARY.** `λ_max` and full spectrum from the
  **tangent-linear equations integrated as one augmented system with the state** (state
  and tangent share identical adaptive steps); periodic Gram–Schmidt/QR reorthonormalization.
  Jacobian analytic if spec/autodiff supplies it, else finite-difference.
- **Independent cross-checks:**
  - **Benettin** as a single **augmented reference+perturbed system with shared steps**
    (no separately-integrated orbits); must agree with the variational `λ_max`.
  - **0–1 test for chaos** (Gottwald–Melbourne) on a scalar observable, **sampled
    stroboscopically on the Poincaré section** or at Δt ≈ dominant period (never the fine
    step, which biases K high); **report median K over many random c values**.
- **Transient discard:** until running statistics converge; start generous, then verify
  insensitivity to discard length.
- **Uncertainty:** `λ_max ± spread` from varying (a) averaging time T, (b) `δ₀`,
  (c) an IC ensemble. Convergence plots for each.

## 3. Architecture
Package `fallingplate/`:
- `model.py` — **only file with physics.** `rhs(t, state, p) -> dstate`, optional
  `jacobian(t, state, p)`; `Params` dataclass from YAML. **Pending spec. No magic numbers.**
- `integrate.py` — `integrate_trajectory(rhs, s0, t_span, cfg)`, `depth_event` (with
  grazing guard), `poincare_event`.
- `regimes.py` — `classify_regime(traj)`, `poincare_map(...)`, rotation/amplitude stats.
- `lyapunov.py` — `variational_spectrum(...)` (primary), `benettin_lmax(...)` (cross-check),
  periodic-angle metric, reusable variational integrator.
- `chaos_tests.py` — `zero_one_test(series, n_c=...)` returning median K; return-map tools.
- `displacement.py` — `landing_displacement(sol, D)` (+ `grazing_flag`),
  `landing_displacement_fixed_time(sol, T)`, sensitivity ensembles.
- `periodic_orbits.py` — `find_orbit_shooting(...)` (Newton shooting; **wrapped angle in
  the residual, record winding number; phase fixed by the Poincaré section**),
  `monodromy(...)` (via the variational integrator), `floquet_multipliers(...)`,
  `continue_in_Istar(...)`.
- `sweeps.py` — parallel I\* driver (`joblib`, 10 cores), adaptive refinement.
- `io_cache.py` — content-addressed cache in `data/`, **key = hash(params) + hash of
  physics/numerics source files** (`model.py`, `integrate.py`, `lyapunov.py`, …).
  Provenance **sidecar** per output records **git hash, Python version, SciPy version,
  params, seed** (git hash lives here only, not in the key).
- `plotting.py` — figure builders; **never prints arrays**.
- `experiments/e1…e6.py` — one script per experiment; each reads a YAML config.
- `config/*.yaml` — ranges, tolerances, D-values, seeds. No magic numbers in code.

**Data flow:** YAML → `sweeps` → `model.rhs` via `integrate` → analysis
(`regimes`/`lyapunov`/`chaos_tests`/`displacement`/`periodic_orbits`) → cached arrays in
`data/` (+ provenance sidecar) → `plotting` → `figures/`.

## 4. Validation ladder — all pass BEFORE any science
1. **Unit tests:** RHS at simple states; force/torque signs & dimensional consistency;
   event correctness; angle-wrap metric; cache round-trip.
2. **Limiting cases:** zero fluid coefficients → analytic ballistic free-fall; any
   spec-supported symmetry / small-oscillation frequency. Exact cases pinned by spec.
3. **Tolerance convergence:** trajectory & `λ_max` vs `rtol/atol` and vs integrator; fix
   the tolerance where results are stable to target precision.
4. **Regime reproduction:** recover steady/fluttering/tumbling and their ordering vs I\*.
   Where the papers report **numeric** values (glide angle, flutter frequency, tumbling
   rate at stated parameters), reproduce **quantitatively within the papers' stated
   precision**; where none is available, **record in `docs/decisions.md` that validation
   was qualitative only** (no silent qualitative passes).
5. **Lyapunov code verified first** on the **logistic map** (`λ=ln2` at r=4; `λ=0` at
   bifurcations) and **Lorenz** (`λ_max≈0.906`). No plate Lyapunov numbers trusted until
   these reproduce.

**Gate:** all five rungs pass → experiments begin.

## 5. Experiments (in order)
- **E1 — Regime map vs I\*.** *Hyp:* ordered regime bands; ≥1 aperiodic band (open).
  Coarse grid (~20–40 I\*), refine near transitions. ~minutes/10 cores.
  Out: `data/e1_regimes.*`. Fig: bifurcation-style diagram + regime strip.
- **E2 — Lyapunov vs I\*.** *Hyp:* `λ_max>0` in candidate band(s). Same grid, refined.
  Out: `λ_max(I*)±unc`, median `K(I*)` cross-check. Fig: `λ_max` & `K` over the
  bifurcation diagram. **GATE:** if `λ_max ≤ 0` everywhere (within uncertainty) and
  `K≈0` → **no chaos; STOP** and propose alternatives (quasiperiodicity/torus study,
  transient chaos, sensitivity without chaos, or a different control parameter).
- **E3 — Route to chaos** *(only if chaos found).* Distinguish period-doubling /
  quasiperiodic (Ruelle–Takens) / intermittency / crisis via bifurcation diagram,
  Poincaré sections, power spectra, return maps near onset.
- **E4 — Feigenbaum test** *(ONLY if ≥3 successive period-doublings resolve).* Locate each
  doubling precisely by **Newton shooting** for the periodic orbit + **continuation in I\***,
  detecting the bifurcation where a **Floquet multiplier crosses −1** (not by eye from a
  coarse sweep). *Hyp to TEST:* `δₙ → 4.669`; report `δₙ ± error` with the limited-ratio
  caveat.
- **E5 — Landing-displacement sensitivity.** IC-ensemble spread `σ(D) = std(Δx)`.
  Displacement is an **unbounded running sum**, so **no saturation is expected**: after
  decorrelation expect **diffusive `~√D`** growth in chaotic bands and **`~linear`** growth
  in regular bands, preceded by an **exponential** (Lyapunov) phase
  `log σ ≈ log σ₀ + (λ_max/U)·D`. **Fit exponential, linear and √ models and compare** which
  dominates in which band; compare the exponential slope to `λ_max/U` from E2. **Choose the
  ensemble perturbation so the exponential phase spans ≥1 decade of σ.** Final `D` set
  chosen after E2 to cover the growth phases. Fig: `σ(D)` with the three fits per band;
  heatmap over (I\*, D).
- **E6 — Multistability / hysteresis.** Multiple ICs + up/down I\* sweeps → coexisting
  attractors / hysteresis. Fig: up-vs-down overlay; basin sampling.

## 6. Milestones (acceptance + stop-and-review after each)
- **M0** scaffold (dirs, PLAN.md, the operating contract, importable empty stubs). ✔ when `pytest` collects.
- **M1** model interface + integrator + ladder rungs 1–3 pass. → review.
- **M2** Lyapunov verified on logistic + Lorenz. → review.
- **M3** regime reproduction (rung 4). → review.
- **M4** E1. **M5** E2 + GATE decision (branch point). **M6** E3/E4. **M7** E5. **M8** E6.
  Stop-and-review after each.

## 7. Cost control & in-context discipline
- Coarse before fine; adaptive refinement only near transitions.
- Cache everything in `data/` (key = param + physics/numerics source-file hash); **never
  recompute**. Provenance sidecar (git hash, Python, SciPy, params, seed) per output.
- **Ask before any run projected >~5 min.**
- Never print arrays/long logs — save to disk, show summaries + figures. the operating contract < 60
  lines. Recommend `/clear` after each milestone.

## 8. Risks & guards
- **Numerical vs physical chaos** → tolerance-convergence of `λ_max`; variational-vs-Benettin
  agreement; independent 0–1 test; pre-verified Lyapunov code.
- **Step-size artefacts** → adaptive stepping + convergence study; events not stepped over.
- **Transient contamination** → discard-length insensitivity study.
- **Coexisting attractors read as chaos** → E6 multistability check.
- **Model-validity limits** → quasi-steady/2D approximation; flag if I\* leaves the spec's
  valid range; report as model properties, not real leaves.
- **Undefined α at zero velocity** → handled per spec item (11); guarded in code & tests.

## 9. Reproducibility & provenance
- `requirements.txt` pinned to the **actual install** (see file). Fixed seeds in config.
- Provenance sidecar per output (git hash, Python, SciPy, params, seed).
- **One command:** `python run_all.py` regenerates every figure from cache or scratch;
  thin `Makefile` wrapper for convenience (macOS).
- `docs/decisions.md` — running decision log. README records **AI-generated vs
  human-done-and-understood**.
- Repo initialized with `git init` at scaffold time.

## Spec checklist — what `docs/model_spec.md` must contain
1. **State variables** — exact list & ordering; velocity frame (lab vs body).
2. **Cyclic variable** — confirm θ is the only periodic coord; its range.
3. **Parameters** — definition of **I\*** (formula); Re, aspect ratio, density ratio,
   gravity/buoyancy, chord `c`; which are held fixed.
4. **Nondimensionalization** — length, time, mass scales.
5. **Force terms** — added-mass/inertial, circulatory/lift, translational & rotational drag;
   coefficients (A, B, C_T, C_R, …) with numeric values.
6. **Torque terms** — same, with coefficients.
7. **Any `|·|`/sign functions** — exact form (smoothness/grazing guards).
8. **Initial conditions** — release convention (orientation, velocity) and defaults.
9. **Poincaré section** — the physically natural section surface.
10. **Validity ranges** — I\*/Re span where the quasi-steady closure is valid.
11. **Angle of attack at zero velocity** — how α (and forces/torques using it) is
    defined/regularized as speed → 0 (the undefined-α limit).
12. **Coefficient traceability** — **every** coefficient tagged with its **paper +
    equation number**, so each value in `model.py` traces to a spec line.

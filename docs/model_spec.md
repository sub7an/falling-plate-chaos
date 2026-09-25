# Model specification: Andersen–Pesavento–Wang (APW) quasi-steady falling-card model

**Status: DRAFT v0.1.** Transcribed from the PDF text of the primary source. It must be cross-checked against the PDF pages by an independent second reader before any physics code depends on it.

## 0. Sources

- **[P91]** Andersen, Pesavento & Wang (2005), "Analysis of transitions between fluttering, tumbling and steady descent of falling cards", J. Fluid Mech. 541, 91–104. **Primary source of this spec.** Equation and section numbers below refer to this paper. Local file: `references/2005_JFM_pp91-104_transitions.pdf`.
- **[P65]** Andersen, Pesavento & Wang (2005), "Unsteady aerodynamics of fluttering and tumbling plates", J. Fluid Mech. 541, 65–90. Dimensional formulation and experiments. Not used for the equations below.
- **[Xu21]** Xu, Li, Li & Liao (2021), arXiv:2102.10386. Secondary, cross-reference only. Some of its numbers disagree with [P91] (see §9). [P91] wins.
- Wang, Birch & Dickinson (2004), J. Exp. Biol. 207, 449–460: cited by [P91] as the source of C_T, C_R, A, B. Not read.

## 1. Formulation used

Dimensionless, thin-card limit (β = b/a → 0): equations (5.1)–(5.3) with closures (4.7)–(4.9) of [P91]. The general-β form is recorded in §8 but is **not** used.

## 2. State and kinematics

Reduced state **u = (vx′, vy′, θ, ω)**, with ω ≡ dθ/dt.

- vx′, vy′: centre-of-mass velocity components in the body-fixed (co-rotating) frame.
- θ: card angle. It is the only cyclic coordinate (period 2π).
- Lab-frame position (x, y) follows by quadrature and does **not** feed back into the dynamics:

```
dx/dt = vx′ cosθ − vy′ sinθ
dy/dt = vx′ sinθ + vy′ cosθ
```

(In the PDF text the primes on the body-frame components were lost in extraction. Check the page.)

## 3. Scales

L = a (semi-major axis); U = (ρs/ρf − 1) g b; time scale L/U.
I* = ρs b / (ρf a); β = b/a; Re = aU/ν.
Re does not appear explicitly. It enters only through the values of the coefficients [P91 §3.3].

## 4. Governing equations (thin-card limit)

```
I*        dvx′/dt = (I*+1) ω vy′ − Γ vy′ − sinθ − Fx′      (5.1)
(I*+1)    dvy′/dt = −I* ω vx′ + Γ vx′ − cosθ − Fy′         (5.2)
¼(I*+½)   dω/dt   = −vx′ vy′ − τ                            (5.3)
          dθ/dt   = ω
```

Closures:

```
Γ = (2/π) [ −C_T vx′ vy′ / √(vx′²+vy′²) + C_R ω ]                          (4.7)

(Fx′, Fy′) = (1/π) [ A − B (vx′²−vy′²)/(vx′²+vy′²) ] √(vx′²+vy′²) (vx′, vy′)   (4.8)

τ = (μ1 + μ2 |ω|) ω                                                         (4.9)
```

## 5. Parameters

| Symbol | Value | Provenance |
|---|---|---|
| C_T | 1.2 | [P91 §5.2], attributed to Wang et al. 2004 |
| C_R | π | same |
| A | 1.4 | same |
| B | 1.0 | same |
| μ1 | 0.2 | **Chosen by the authors** so the time scale and periodic flutter/tumble come out qualitatively right [P91 §3.3, §5.2]. Not fitted to data. |
| μ2 | 0.2 | same |
| I* | control parameter | varied |

The authors conjecture μ1 and μ2 are the coefficients most sensitive to Re [P91 §3.3].

## 6. Fixed points [P91 (6.1), (6.2)]

```
edge-on:    (vx′, vy′, θ, ω) = (∓V, 0, π/2 or 3π/2, 0),  V = √(π/(A−B)) ≈ 2.8025
broadside:  (vx′, vy′, θ, ω) = (0, ∓W, 0 or π, 0),        W = √(π/(A+B)) ≈ 1.1441
```

Sign pairing (derived by hand, not stated in the paper): (−V, π/2), (+V, 3π/2), (−W, 0), (+W, π). The first and third were verified numerically (§7).

Stability [P91 §6.1]: the edge-on points are unstable for all I* and μ1. The broadside-on points lose stability through a supercritical Hopf bifurcation. The Hopf curve formula (6.3) was garbled in text extraction and is **not** transcribed here. Take it from the PDF if needed.

Hopf curve [P91 eq. (6.3), p.101] — *recovered from OCR, not independently re-derived* (added 2026-09-22 by the researcher as an addition, not a correction):

```
μ1 = (1/4) √(3/(5π)) · (2 I* + 1)/(I* − 1)
```

## 7. Validation anchors

All from [P91] (§5.2, §6.2), for μ1 = μ2 = 0.2, C_T = 1.2, C_R = π, A = 1.4, B = 1.0, thin-card limit. Time is dimensionless.

1. **Regimes (Fig. 3):** I* = 1.1 fluttering; 1.4 period-one tumbling; 1.45 period-two tumbling; 1.6 periodic mixture of fluttering and tumbling; 2.2 chaotic; 3.0 small-amplitude broadside-on fluttering. A period-doubling bifurcation lies between I* = 1.4 and 1.45. Chaos occurs for I* roughly 1.8–2.8 (Fig. 4b).
2. **Lyapunov exponent:** λ_max = 0.13 ± 0.01 at I* = 2.2. The method is not stated.
3. **Heteroclinic bifurcation:** I*_C ≈ 1.2191. I* = 1.2190 flutters; I* = 1.2192 tumbles. The card falls edge-on for about 25 card widths before reaching maximum speed (Fig. 5a).
4. **Edge-on saddle eigenvalues at I*_C:** λs = −0.5854 (the vx′ direction, decoupled), λu = +0.3813, λ± = −0.9861 ± 2.5949i.
5. **Period law:** T = T0 + (2/λu) ln(1/|I* − I*_C|), with T0 = 4.3629 (a fit). The factor 2 arises because the saddle is passed twice per period. Fig. 5d shows agreement on both sides of I*_C.
6. Solutions near I*_C are noise-sensitive; the integration must be done with high accuracy there [P91 §6.2].
7. The periodic solutions shown are independent of initial conditions after short transients [P91 §5.2].

### Independent checks already run (scratch script, **not project code**)

- Fixed-point residuals ≈ 1e-16 for the edge-on (−V, 0, π/2, 0) and broadside (0, −W, 0, 0) points.
- Finite-difference Jacobian at the edge-on point, I* = 1.2191: eigenvalues −0.9861 ± 2.5950i, −0.5854, +0.3813. **Match anchor 4.** Analytic check: λs = −2√((A−B)/π)/I* = −0.58539.
- DOP853, rtol 1e-12, atol 1e-14, IC (vx′, vy′, θ, ω) = (0, 0.01, 0.5, 0): I* = 1.2190 flutters and 1.2192 tumbles (**matches anchor 3**); 1.1 flutters; 1.3, 1.4, 1.45 tumble.
- Benettin λ_max at I* = 2.2 (T = 6000, two different ICs and perturbation sizes): 0.138 and 0.139. Anchor 2 says 0.13 ± 0.01, so this is consistent at the upper edge. It is **not** a converged estimate.
- **Not yet checked:** the period law (anchor 5), the period-doubling location, the chaotic range, and the regimes at I* = 1.6 and 3.0.

## 8. General-β form (recorded, not used)

```
(I*+β²)   dvx′/dt = (I*+1) ω vy′ − Γ vy′ − sinθ − Fx′     (4.4)
(I*+1)    dvy′/dt = −(I*+β²) ω vx′ + Γ vx′ − cosθ − Fy′    (4.5)
¼[I*(1+β²) + ½(1−β²)²] dω/dt = (β²−1) vx′ vy′ − τ          (4.6)
```

## 9. Cautions and known discrepancies

- **Singularities at zero speed.** The direction of (vx′, vy′) is undefined at zero speed, so Γ and F need a safe limit there. Γ·v and F are degree-2 homogeneous in the velocity, so they tend to 0 and are C¹, but not C². Initial conditions must have non-zero speed.
- **Invariant subspace.** An initial condition with vx′ = 0, ω = 0, θ = 0 never rotates (observed numerically). Use a tilted initial condition.
- **Non-smooth torque.** The |ω|ω term is C¹ but not C² at ω = 0. Adaptive high-order integrators may lose accuracy at ω = 0 crossings. The convergence study must check the observed order.
- **[Xu21] inconsistency.** [Xu21] gives flutter at 1.2190 and tumble at 1.2191, and separately says 1.2191466312021015 flutters. These cannot both hold, and [P91] gives tumble at 1.2192. Do not use [Xu21]'s 1e-11 claim as a benchmark until reproduced.
- **Scope.** μ1 and μ2 were tuned for qualitatively correct behaviour. Everything computed is a property of this phenomenological model, not of real cards.
- **Relation to [P65]'s I\*.** [P65] defines I* for an elliptical section as ρs h (l²+h²) / (2 ρf l³). With l = 2a, h = 2b this gives I*_[P65] = ½(1+β²) I*_[P91]. In the thin limit, I*_[P65] ≈ I*_[P91]/2, so I*_C ≈ 1.2191 in [P91]'s convention is about 0.61 in [P65]'s convention. [P65] quotes experimental flutter–tumble transitions of 0.2–0.3 (Smith) and 0.4 (Belmonte). Derived by hand; verify.

## 10. Not specified by [P91] (our choices, to be logged in docs/decisions.md)

- Initial conditions used for Figs. 3 and 5.
- The integrator and tolerances used by the authors.
- How λ_max and T were computed.
- The Poincaré section. **Our choice:** upward crossings of ω = 0, unless a better section emerges. It is not from the paper.

## 11. Checklist for the second reader

Compare against the PDF page images (journal page 91 = PDF page 1):

- Equations (4.4)–(4.9) and (5.1)–(5.3): journal p. 98, PDF p. 8.
- Fixed points (6.1), (6.2): journal p. 100, PDF p. 10.
- §6.2 numbers (I*_C, eigenvalues, T0) and eq. (6.4): journal pp. 101–102, PDF pp. 11–12.
- Parameter values in §5.2: journal p. 98, PDF p. 8.
- Confirm the primes on the body-frame velocities and the exact form of (4.8).

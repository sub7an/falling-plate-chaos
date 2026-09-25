# Chaos, its route, and displacement sensitivity in the quasi-steady falling-plate model

**Readable mirror.** The canonical, publication-ready manuscript is `docs/paper.tex` (build with
`make paper` → `docs/paper.pdf`); it is more complete than this file (summary tables, full
front/back matter, and the finalized bibliography). This Markdown version is kept for quick reading
and may lag. Author: M. Subhan Rao (Independent researcher). Computational work and drafting were
AI-assisted under the author's direction; every result is reproducible via `python run_all.py` and
every methodological choice is logged in `docs/decisions.md`. All claims are properties of the
phenomenological model, not of physical falling objects.

---

## Abstract

We study the two-dimensional quasi-steady "falling-card" model of Andersen, Pesavento & Wang
(2005) in its thin-card limit, treating the dimensionless moment of inertia $I^\*$ as a single
control parameter and asking whether, where, and by what route the model becomes chaotic, and how
sensitive the horizontal landing position is to release conditions. Using a validated numerical
pipeline — a variational (tangent-linear) Lyapunov solver cross-checked against Benettin's method
and the Gottwald–Melbourne 0–1 test, all verified first on the logistic map and the Lorenz system —
we confirm a chaotic band at $I^\* \approx [1.63, 2.83]$, with maximal Lyapunov exponent
$\lambda_{\max}$ peaking at $+0.140$ near $I^\*=2.2$, consistent with the value $0.13\pm0.01$ reported
in the source paper. The route into chaos is **not** a Feigenbaum period-doubling cascade: only two
successive doublings occur before the period-2 orbit gives way, through a narrow band threaded with
periodic windows, to chaos; power spectra rule out a quasiperiodic (torus) route. Within the chaotic
band the standard deviation of horizontal landing displacement across an ensemble of perturbed
release conditions grows exponentially, at a rate equal to $\lambda_{\max}$ within estimation
uncertainty; outside the band it is insensitive. Up/down continuation and basin sampling reveal no
hysteresis and no non-trivial multistability. We report one discrepancy with the source paper's
figure: at $I^\*=1.45$ this model tumbles with period one, not period two.

---

## 1. Introduction

Falling cards, leaves, and plates exhibit a rich set of descent styles — steady gliding,
side-to-side fluttering, end-over-end tumbling, and irregular motion — whose transitions have been
mapped experimentally and modelled with quasi-steady aerodynamic force laws. Andersen, Pesavento &
Wang (2005; hereafter APW) [1,2] introduced a low-dimensional ordinary-differential-equation model
in which the fluid forces enter through translational and rotational drag and a circulatory
("lift") term with coefficients fixed from direct force measurements [3]. APW showed that as the
dimensionless moment of inertia $I^\*$ is varied the model reproduces the observed flutter–tumble
sequence and reported a chaotic regime with a positive Lyapunov exponent.

This paper revisits the model as a dynamical system in its own right. Our aims are (i) to confirm,
with independent diagnostics and pre-verified code, whether a chaotic regime exists and to delimit
it; (ii) to characterise the *route* by which chaos appears as $I^\*$ increases, testing — rather
than assuming — the hypothesis of Feigenbaum period-doubling scaling; and (iii) to quantify the
sensitivity of the horizontal landing position to release conditions, connecting it to the Lyapunov
exponent. Throughout we treat null and negative results as first-class outcomes: a chaos claim
requires a positive Lyapunov exponent *and* an independent confirmation, and a Feigenbaum analysis is
undertaken only if a genuine doubling cascade is resolved.

All findings are properties of this phenomenological model. Its coefficients were tuned by APW for
qualitatively correct behaviour and are not fitted to any specific card; the quasi-steady, purely
two-dimensional closure has a limited range of validity. We therefore make no claim about real
falling objects.

## 2. Model

We use the dimensionless thin-card limit ($\beta = b/a \to 0$, with $a,b$ the semi-axes) of APW,
equations (5.1)–(5.3) with closures (4.7)–(4.9). Length is scaled by the semi-major axis $a$,
velocity by $U=\sqrt{(\rho_s/\rho_f-1)\,g\,b}$, and time by $a/U$; the control parameter is
$I^\* = \rho_s b/(\rho_f a)$. The reduced state is $\mathbf{u}=(v_x', v_y', \theta, \omega)$, where
$v_x',v_y'$ are the centre-of-mass velocity components in the body-fixed (co-rotating) frame,
$\theta$ is the card angle (the only cyclic coordinate, period $2\pi$), and $\omega=\dot\theta$. The
lab-frame position $(x,y)$ follows by quadrature and does not feed back into the dynamics
(translation invariance),

$$\dot x = v_x'\cos\theta - v_y'\sin\theta,\qquad \dot y = v_x'\sin\theta + v_y'\cos\theta .$$

The governing equations are

$$
\begin{aligned}
I^\*\,\dot v_x' &= (I^\*{+}1)\,\omega v_y' - \Gamma v_y' - \sin\theta - F_x', \\
(I^\*{+}1)\,\dot v_y' &= -I^\*\,\omega v_x' + \Gamma v_x' - \cos\theta - F_y', \\
\tfrac14(I^\*{+}\tfrac12)\,\dot\omega &= -v_x' v_y' - \tau, \qquad \dot\theta = \omega,
\end{aligned}
$$

with closures (writing $s=\sqrt{v_x'^2+v_y'^2}$)

$$
\Gamma = \tfrac{2}{\pi}\!\left[-C_T\,\frac{v_x' v_y'}{s} + C_R\,\omega\right],\quad
(F_x',F_y') = \tfrac{1}{\pi}\!\left[A - B\,\frac{v_x'^2-v_y'^2}{s^2}\right] s\,(v_x',v_y'),\quad
\tau = (\mu_1+\mu_2|\omega|)\,\omega .
$$

The coefficients, fixed throughout, are $C_T=1.2$, $C_R=\pi$, $A=1.4$, $B=1.0$, and
$\mu_1=\mu_2=0.2$ [1,3]. The direction of $(v_x',v_y')$ is undefined at zero speed; $\Gamma v$ and
$\mathbf{F}$ are degree-two homogeneous and tend to zero there ($C^1$ but not $C^2$), so admissible
initial conditions carry non-zero speed. The torque term $|\omega|\omega$ is likewise $C^1$ but not
$C^2$ at $\omega=0$.

The model has two families of fixed points: **edge-on** $(\mp V,0,\pi/2\text{ or }3\pi/2,0)$ with
$V=\sqrt{\pi/(A-B)}\approx 2.8025$, and **broadside** $(0,\mp W,0\text{ or }\pi,0)$ with
$W=\sqrt{\pi/(A+B)}\approx 1.1441$.

## 3. Numerical methods

**Integration.** Trajectories are integrated with `scipy.integrate.solve_ivp` using the explicit
DOP853 method (adaptive high-order Runge–Kutta) for production runs and RK45 for coarse scans.
Production tolerances are `rtol` $=10^{-10}$, `atol` $=10^{-12}$; a tolerance-convergence study fixed
these (Section 4). Event detection uses the dense output for sub-step localisation.

**Lyapunov exponents.** The maximal exponent $\lambda_{\max}$ is computed by the variational
(tangent-linear) method: the state and a tangent vector are advanced as one augmented system sharing
identical adaptive steps, with periodic QR reorthonormalisation; $\lambda_{\max}$ is the
time-averaged log growth rate. As an independent cross-check we use Benettin's two-trajectory
renormalisation method, implemented as a single augmented reference+perturbed system (shared steps)
so that it does not reuse the analytic Jacobian.

**Independent chaos test.** We apply the Gottwald–Melbourne 0–1 test [4] to a scalar observable
(the angular velocity sampled on a Poincaré section, i.e. stroboscopically rather than at the fine
integrator step), reporting the median indicator $K$ over 100 random frequencies; $K\to 1$ indicates
chaos and $K\to 0$ regularity.

**Poincaré section and period counting.** Because a one-signed tumbling orbit never crosses
$\omega=0$, we section on the cyclic angle, recording the reduced state at event-located crossings of
$\theta \equiv \text{const}\ (\mathrm{mod}\ 2\pi)$; distinct clustered return points give the period.

**Reproducibility.** Every output is content-addressed cached with a provenance sidecar recording the
git hash, Python/NumPy/SciPy versions, parameters, and seed. The cache key includes a hash of the
physics/numerics source files, so any change to the model or solvers invalidates dependent results.
A single command, `python run_all.py`, regenerates the entire figure set.

## 4. Validation

No plate result was trusted until a validation ladder passed. **Limiting cases:** with all fluid
coefficients set to zero the model reduces to constant-acceleration free fall along a principal
axis, reproduced to $\sim10^{-9}$; the fixed-point residuals vanish to $<10^{-13}$. **Linearisation:**
the analytic Jacobian agrees with finite differences and, at the heteroclinic value
$I^\*_C=1.2191$, gives the edge-on saddle eigenvalues $-0.5854$, $+0.3813$, $-0.9861\pm2.5949\,i$,
matching APW; the decoupled stable eigenvalue equals the closed form
$-2\sqrt{(A-B)/\pi}/I^\* = -0.58539$. **Regimes across $I^\*_C$:** integration reproduces fluttering
at $I^\*=1.2190$ and tumbling at $1.2192$. **Tolerance convergence:** the trajectory is
Cauchy-convergent as tolerances tighten and DOP853 and RK45 agree at tight tolerance. **Lyapunov
code:** verified on the logistic map ($\lambda=\ln 2=0.6931$ at $r=4$; the analytic period-2 value
$-0.9163$ at $r=3.2$) and the Lorenz system ($\lambda_{\max}=0.907$, a zero exponent, and the exact
trace sum $\sum\lambda_i=-(\sigma+1+\beta)$); the 0–1 test returns $K\approx1$ for chaotic and
$K\approx0$ for periodic logistic series.

## 5. Results

### 5.1 Regime sequence (Fig. 1)

A parallel sweep over $I^\*\in[1.0,3.0]$ recovers an ordered sequence of regimes:
**fluttering** for $I^\* \lesssim 1.2$; **period-one tumbling** from the heteroclinic transition at
$I^\*_C\approx1.2191$; a **period-two** band; an **aperiodic (chaotic)** band; and a return to
**small-amplitude broadside fluttering** near $I^\*=3.0$. The bifurcation diagram of the
section angular velocity displays the expected structure — paired flutter branches, a single
tumbling branch, a period-2 split, a chaotic smear, and flutter again.

### 5.2 Chaos, confirmed by two diagnostics (Fig. 2)

The maximal Lyapunov exponent is zero within numerical error ($|\lambda_{\max}|<10^{-3}$) throughout
the periodic bands and rises to a plateau of $0.12$–$0.14$ across $I^\*\approx[1.63,2.83]$, peaking at
$\lambda_{\max}=+0.140$ at $I^\*=2.2$; the independent Benettin estimate there is $0.137$, and the
0–1 test median $K$ jumps from $0$ to $\approx1$ over the same band. Requiring *both* a robustly
positive $\lambda_{\max}$ and $K>0.5$ confirms chaos in $I^\*\approx[1.63,2.83]$. The conjunction is
essential: the 0–1 test alone produces spurious elevated $K$ at bifurcation points (transient
artifacts that decay with longer sampling), and finite-time numerical noise produces sub-$10^{-3}$
positive $\lambda_{\max}$ in periodic bands; each diagnostic vetoes the other's false positives. The
peak value is consistent with APW's reported $0.13\pm0.01$.

### 5.3 Route to chaos is not Feigenbaum (Fig. 3)

Resolving the transition finely (event-located sections, integrator-precision return points) shows
only **two** successive period-doublings: period-1 $\to$ 2 at $I^\*\approx1.474$ and 2 $\to$ 4 at
$I^\*\approx1.599$. There is no period-8; the period-4 state gives way, across a narrow window
($\sim0.01$ in $I^\*$), to chaos interspersed with periodic windows (period-3, -5, -6). Power spectra
distinguish these: the window states show clean line spectra (commensurate harmonics of a single
fundamental) whereas the "many-branch" states show a broadband floor identical to the chaotic
control — i.e. weak chaos, not two incommensurate frequencies. There is thus **no quasiperiodic
(torus / Ruelle–Takens) route**. Because a Feigenbaum ratio $\delta$ requires at least three
successive doublings, and only two are present, we do **not** report a scaling exponent: doing so
from two doublings would be unsupported. The route is best described as a truncated doubling
(1$\to$2$\to$4) followed by a window-threaded transition to chaos.

### 5.4 Landing-displacement sensitivity (Fig. 4)

For ensembles of $N=200$ trajectories released from initial conditions perturbed by
$\delta_0=10^{-6}$, we measure the standard deviation $\sigma$ of the horizontal displacement
$\Delta x$ at fixed descent depth $D$ (with a grazing guard on the landing event). In the chaotic
band $\sigma$ grows exponentially, by six to seven orders of magnitude over $D\le200$ chords; in the
regular bands it remains flat within noise, i.e. displacement is insensitive there (these regimes are
attracting limit cycles, so perturbations decay — linear growth would require neutral dynamics).

Measured directly in time, the growth rate $\mathrm{d}(\log\sigma_x)/\mathrm{d}t$ equals the Lyapunov
exponent to within $\sim10\%$ (ratios $1.11,0.92,0.90$ at $I^\*=2.0,2.2,2.5$; mixed sign), consistent
with the relation $\log\sigma_x \approx \log\sigma_0 + \lambda_{\max}\,t$ to within the combined
finite-time uncertainty of both estimates. The corresponding fixed-depth comparison to
$\lambda_{\max}/U$ (with $U$ the mean descent speed) is looser because it folds in a depth-dependent
time-to-depth conversion. Beyond the exponential phase $\sigma$ continues to grow sub-linearly with
$D$; a clean diffusive $\sqrt{D}$ asymptote is not isolated within the depths studied.

### 5.5 No non-trivial multistability (Fig. 5)

Up and down continuation sweeps of $\langle|\omega|\rangle$ coincide across $I^\*\in[1.0,3.0]$: there
is no hysteresis, including at the flutter–tumble transition. Basin sampling with random initial
conditions finds a single attractor at every diagnostic $I^\*$, inside the periodic windows and the
chaotic band alike. The only multistability is the trivial clockwise/counter-clockwise
reflection-symmetry pair (mirror spin directions with identical statistics). This rules out the
possibility that the chaotic diagnostics reflect coexisting attractors rather than a single chaotic
one.

## 6. Discussion

The model possesses a genuine chaotic band, confirmed by independent methods and delimited more
sharply than by regime inspection alone. Two aspects refine or depart from expectations. First, the
**route** is not the period-doubling cascade one might anticipate: a single robust doubling and a
second, closely followed by a window-threaded transition, is a qualitatively different scenario for
which Feigenbaum universality does not apply. Second, we find a clean **operational link** between a
dynamical invariant and an observable of practical interest — the exponential rate of
landing-position spread equals the Lyapunov exponent — while emphasising that the regular regimes are
*insensitive*, so predictability of landing position depends discontinuously on $I^\*$.

We note one **discrepancy** with the source paper. APW's Fig. 3 reads period-two tumbling at
$I^\*=1.45$; in our integration this value is period-*one* (the section state and revolution time are
constant to $\sim10^{-12}$ across revolutions, and three distinct initial conditions converge to the
same orbit with no bistability). The first doubling occurs here at $I^\*\approx1.474$. We report this
as a property of the model under our conventions rather than reconciling it by parameter adjustment;
the initial conditions used for APW's figures are not specified, and small differences in the
bifurcation location are plausible.

**Limitations.** The coefficients are phenomenological and the closure quasi-steady and
two-dimensional; results outside the parameter range where that closure is credible should not be
over-interpreted. The plate Lyapunov exponent is only moderately converged at attainable averaging
times (a few-percent finite-time scatter), which sets the precision of the $\lambda_{\max}$–
sensitivity comparison. The diffusive tail of $\sigma(D)$ and the fine structure of the
$[1.598,1.61]$ transition (intermittency versus crisis) remain open to deeper characterisation but do
not affect the conclusions above.

## 7. Conclusions

For the thin-card APW model with the standard coefficients, as $I^\*$ increases the descent passes
from fluttering through tumbling to a **chaotic band at $I^\*\approx[1.63,2.83]$**
($\lambda_{\max}$ peak $+0.140$), confirmed by variational and Benettin Lyapunov exponents and the
0–1 test. The route is a **truncated period-doubling (two doublings) followed by windowed chaos, not
a Feigenbaum cascade and not a torus**. Horizontal landing displacement is **exponentially sensitive
to release conditions in the chaotic band, at the Lyapunov rate**, and insensitive elsewhere. The
model exhibits **no hysteresis and no non-trivial multistability**. We additionally flag a period
discrepancy with the source paper at $I^\*=1.45$.

## Reproducibility

Code, configurations, and the decision log accompany this manuscript. The validation suite
(`pytest`, 59 tests) covers the physics, the Lyapunov and 0–1 diagnostics, and the analysis pipeline.
`python run_all.py` regenerates Figs. 1–5 (E1, E2, E3, E5; E6) from cache or from scratch; the
Feigenbaum stage (E4) is intentionally skipped and reported as the null result of Section 5.3.

## Figures

- **Fig. 1** `figures/e1_regime_map.png` — bifurcation diagram and regime strip vs $I^\*$.
- **Fig. 2** `figures/e2_lyapunov.png` — $\lambda_{\max}(I^\*)$ and median $K(I^\*)$ over the
  bifurcation diagram; the chaos gate.
- **Fig. 3** `figures/e3_route.png` — route-to-chaos zoom: section vs $I^\*$ with $\lambda_{\max}$.
- **Fig. 4** `figures/e5_sensitivity.png` — $\sigma(D)$ for chaotic and regular $I^\*$, with
  exponential-phase fits and power-law references.
- **Fig. 5** `figures/e6_multistability.png` — up/down continuation and basin sampling.

## References

[1] A. Andersen, U. Pesavento, Z. J. Wang, *Analysis of transitions between fluttering, tumbling and
steady descent of falling cards*, J. Fluid Mech. **541**, 91–104 (2005).
[2] A. Andersen, U. Pesavento, Z. J. Wang, *Unsteady aerodynamics of fluttering and tumbling
plates*, J. Fluid Mech. **541**, 65–90 (2005).
[3] Z. J. Wang, J. M. Birch, M. H. Dickinson, *Unsteady forces and flows in low Reynolds number
hovering flight*, J. Exp. Biol. **207**, 449–460 (2004).
[4] G. A. Gottwald, I. Melbourne, *On the implementation of the 0–1 test for chaos*, SIAM J. Appl.
Dyn. Syst. **8**, 129–145 (2009).
[5] G. Benettin, L. Galgani, A. Giorgilli, J.-M. Strelcyn, *Lyapunov characteristic exponents for
smooth dynamical systems …*, Meccanica **15**, 9–20 (1980).
[6] E. N. Lorenz, *Deterministic nonperiodic flow*, J. Atmos. Sci. **20**, 130–141 (1963).
[7] M. J. Feigenbaum, *Quantitative universality for a class of nonlinear transformations*, J. Stat.
Phys. **19**, 25–52 (1978).

*References [1–3] are the model sources as recorded in `docs/model_spec.md`; [4–7] are the standard
method/verification references. Bibliographic details to be confirmed against originals before
submission.*

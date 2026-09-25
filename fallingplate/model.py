"""model.py — the ONLY file that encodes physics.

Equations and coefficients come ONLY from docs/model_spec.md (Andersen,
Pesavento & Wang 2005, J. Fluid Mech. 541, 91-104 = [P91]). Nothing here is
filled from memory. No magic numbers: every physical coefficient arrives via
`Params`, loaded from config/model_params.yaml. The literals 2/pi, 1/pi, 1/4,
1/2 are structural parts of the equation forms (5.1)-(5.3), (4.7)-(4.9), not
tunable parameters.

Formulation: dimensionless thin-card limit (beta = b/a -> 0), reduced state

    u = (vx', vy', theta, omega),   omega = d theta / dt,

with vx', vy' the centre-of-mass velocity components in the body-fixed
(co-rotating) frame and theta the (cyclic) card angle. Lab position (x, y) is a
quadrature and does NOT feed back (translation invariance, PLAN A1), so it is
not part of the state carried here.

Governing equations (spec §4):

    I*        dvx'/dt = (I*+1) w vy' - Gamma vy' - sin(theta) - Fx'   (5.1)
    (I*+1)    dvy'/dt = -I* w vx' + Gamma vx' - cos(theta) - Fy'      (5.2)
    (1/4)(I*+1/2) dw/dt = -vx' vy' - tau                             (5.3)
              dtheta/dt = w

Closures (spec §4):

    Gamma = (2/pi) [ -C_T vx' vy' / s + C_R w ]                       (4.7)
    (Fx', Fy') = (1/pi) [ A - B (vx'^2 - vy'^2)/s^2 ] s (vx', vy')    (4.8)
    tau = (mu1 + mu2 |w|) w                                           (4.9)

where s = sqrt(vx'^2 + vy'^2).

Zero-speed limit (spec §9): the direction of (vx', vy') is undefined at s = 0,
but Gamma*v and F are degree-2 homogeneous in velocity, so the velocity-drag /
fluid-force terms tend to 0 and are C^1 (not C^2). We return that limit when
s == 0; valid initial conditions must have non-zero speed.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional

import numpy as np
import yaml


@dataclass(frozen=True)
class Params:
    """Physical coefficients for the thin-card model (spec §5).

    Every field traces to docs/model_spec.md §5 -> [P91]. I_star is the control
    parameter (spec §3), varied per experiment. There are no default numeric
    values on purpose: coefficients must come from the spec-driven YAML, never
    from a code default.
    """

    C_T: float
    C_R: float
    A: float
    B: float
    mu1: float
    mu2: float
    I_star: float

    @classmethod
    def from_yaml(cls, path: str, I_star: Optional[float] = None) -> "Params":
        """Load coefficients from a spec-driven YAML file.

        If `I_star` is given it overrides the file's value (experiments/tests
        sweep the control parameter while holding the coefficients fixed).
        """
        with open(path, "r") as fh:
            data = yaml.safe_load(fh)
        p = cls(
            C_T=float(data["C_T"]),
            C_R=float(data["C_R"]),
            A=float(data["A"]),
            B=float(data["B"]),
            mu1=float(data["mu1"]),
            mu2=float(data["mu2"]),
            I_star=float(data["I_star"]),
        )
        if I_star is not None:
            p = replace(p, I_star=float(I_star))
        return p


# --- Closures (4.7)-(4.9). Exposed for unit tests; reused by rhs/jacobian. ---

_TWO_OVER_PI = 2.0 / np.pi
_ONE_OVER_PI = 1.0 / np.pi


def gamma(vx: float, vy: float, omega: float, p: Params) -> float:
    """Circulation-like term Gamma, closure (4.7). Zero-speed safe."""
    s = np.hypot(vx, vy)
    drag = 0.0 if s == 0.0 else -p.C_T * vx * vy / s
    return _TWO_OVER_PI * (drag + p.C_R * omega)


def fluid_force(vx: float, vy: float, p: Params) -> tuple[float, float]:
    """Fluid force (Fx', Fy'), closure (4.8). Zero-speed safe (limit is 0)."""
    s2 = vx * vx + vy * vy
    if s2 == 0.0:
        return 0.0, 0.0
    s = np.sqrt(s2)
    g = p.A - p.B * (vx * vx - vy * vy) / s2
    coeff = _ONE_OVER_PI * g * s
    return coeff * vx, coeff * vy


def torque(omega: float, p: Params) -> float:
    """Dissipative rotational torque tau, closure (4.9)."""
    return (p.mu1 + p.mu2 * abs(omega)) * omega


def _inertia_coeffs(p: Params) -> tuple[float, float, float]:
    """LHS inertia factors of (5.1), (5.2), (5.3) in the thin-card limit."""
    return p.I_star, p.I_star + 1.0, 0.25 * (p.I_star + 0.5)


def rhs(t: float, state: np.ndarray, p: Params) -> np.ndarray:
    """Right-hand side d(state)/dt, state = (vx', vy', theta, omega).

    Translation-invariant (PLAN A1): no explicit x, y, and no explicit t
    dependence. Equations (5.1)-(5.3) plus dtheta/dt = omega.
    """
    vx, vy, theta, omega = state
    m_x, m_y, m_w = _inertia_coeffs(p)

    Gamma = gamma(vx, vy, omega, p)
    Fx, Fy = fluid_force(vx, vy, p)
    tau = torque(omega, p)

    dvx = ((p.I_star + 1.0) * omega * vy - Gamma * vy - np.sin(theta) - Fx) / m_x
    dvy = (-p.I_star * omega * vx + Gamma * vx - np.cos(theta) - Fy) / m_y
    dtheta = omega
    domega = (-vx * vy - tau) / m_w

    return np.array([dvx, dvy, dtheta, domega])


def lab_velocity(vx: float, vy: float, theta: float) -> tuple[float, float]:
    """Lab-frame centre-of-mass velocity (dx/dt, dy/dt) from body-frame
    components (spec §2): dx/dt = vx' cosθ − vy' sinθ, dy/dt = vx' sinθ + vy' cosθ.
    """
    c, s = np.cos(theta), np.sin(theta)
    return vx * c - vy * s, vx * s + vy * c


def rhs_full(t: float, state: np.ndarray, p: Params) -> np.ndarray:
    """6-D RHS for (x, y, vx', vy', theta, omega): the lab-frame position
    quadrature (spec §2) prepended to `rhs`. Position does NOT feed back into the
    reduced dynamics (translation invariance, PLAN A1); it is carried only to
    measure landing displacement (E5).
    """
    x, y, vx, vy, theta, omega = state
    dx, dy = lab_velocity(vx, vy, theta)
    d = rhs(t, np.array([vx, vy, theta, omega]), p)
    return np.array([dx, dy, d[0], d[1], d[2], d[3]])


def jacobian(t: float, state: np.ndarray, p: Params) -> np.ndarray:
    """Analytic Jacobian d(rhs)/d(state), 4x4.

    Derived by hand from (5.1)-(5.3), (4.7)-(4.9); cross-checked against a
    central finite difference of `rhs` in the test suite. Used later by the
    variational Lyapunov integrator / monodromy (PLAN §2). At zero speed the
    velocity-direction derivatives are not well-defined (C^1 not C^2); we return
    the smooth non-velocity part there, since valid ICs have non-zero speed.
    """
    vx, vy, theta, omega = state
    m_x, m_y, m_w = _inertia_coeffs(p)
    Istar = p.I_star
    C_T, C_R, A, B = p.C_T, p.C_R, p.A, p.B

    s2 = vx * vx + vy * vy
    s = np.sqrt(s2)

    # Gamma partials, closure (4.7).
    dGamma_domega = _TWO_OVER_PI * C_R
    if s == 0.0:
        dGamma_dvx = 0.0
        dGamma_dvy = 0.0
        dFx_dvx = dFx_dvy = dFy_dvx = dFy_dvy = 0.0
        G = A  # (vx^2 - vy^2)/s^2 term -> 0 contribution to force below anyway
    else:
        s3 = s2 * s
        s4 = s2 * s2
        dGamma_dvx = _TWO_OVER_PI * (-C_T * vy**3 / s3)
        dGamma_dvy = _TWO_OVER_PI * (-C_T * vx**3 / s3)

        # Fluid-force partials, closure (4.8). G = A - B (vx^2 - vy^2)/s^2.
        G = A - B * (vx * vx - vy * vy) / s2
        dG_dvx = -4.0 * B * vx * vy * vy / s4
        dG_dvy = 4.0 * B * vx * vx * vy / s4

        dFx_dvx = _ONE_OVER_PI * (dG_dvx * s * vx + G * vx * vx / s + G * s)
        dFx_dvy = _ONE_OVER_PI * (dG_dvy * s * vx + G * vx * vy / s)
        dFy_dvx = _ONE_OVER_PI * (dG_dvx * s * vy + G * vx * vy / s)
        dFy_dvy = _ONE_OVER_PI * (dG_dvy * s * vy + G * vy * vy / s + G * s)

    dtau_domega = p.mu1 + 2.0 * p.mu2 * abs(omega)
    Gamma = gamma(vx, vy, omega, p)  # value enters rows 1-2 via the -Gamma vy / +Gamma vx terms

    # Row 1: d(dvx)/d(state), numerator N1 = (I*+1) w vy - Gamma vy - sin(theta) - Fx.
    dN1_dvx = -dGamma_dvx * vy - dFx_dvx
    dN1_dvy = (Istar + 1.0) * omega - dGamma_dvy * vy - Gamma - dFx_dvy
    dN1_dth = -np.cos(theta)
    dN1_dw = (Istar + 1.0) * vy - dGamma_domega * vy

    # Row 2: numerator N2 = -I* w vx + Gamma vx - cos(theta) - Fy.
    dN2_dvx = -Istar * omega + dGamma_dvx * vx + Gamma - dFy_dvx
    dN2_dvy = dGamma_dvy * vx - dFy_dvy
    dN2_dth = np.sin(theta)
    dN2_dw = -Istar * vx + dGamma_domega * vx

    # Row 4: numerator N3 = -vx vy - tau.
    dN3_dvx = -vy
    dN3_dvy = -vx
    dN3_dth = 0.0
    dN3_dw = -dtau_domega

    return np.array(
        [
            [dN1_dvx / m_x, dN1_dvy / m_x, dN1_dth / m_x, dN1_dw / m_x],
            [dN2_dvx / m_y, dN2_dvy / m_y, dN2_dth / m_y, dN2_dw / m_y],
            [0.0, 0.0, 0.0, 1.0],
            [dN3_dvx / m_w, dN3_dvy / m_w, dN3_dth / m_w, dN3_dw / m_w],
        ]
    )

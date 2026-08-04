"""
transport.py
============
Checklist sections 8-9, 12-14, 18-21: dispersion, advection, retardation,
and the 1-D advection-dispersion-reaction equation itself.

SOLUTION METHOD
---------------
This module solves the 1-D ADE

    R dC/dt = D_L d2C/dx2 - v dC/dx

for a semi-infinite domain with a constant-concentration source at x=0
(C(0,t)=C0) and C(x,0)=0, using the classic retarded Ogata-Banks closed-form
solution (Ogata & Banks 1961; the retardation extension via v_eff=v/R,
D_eff=D_L/R is a standard textbook result, e.g. Fetter, "Contaminant
Hydrogeology"):

    C(x,t)/C0 = 0.5 * [ erfc(a) + exp(v*x/D_L) * erfc(b) ]
    a = (R*x - v*t) / (2*sqrt(D_L*R*t))
    b = (R*x + v*t) / (2*sqrt(D_L*R*t))

This substitution was verified symbolically with sympy before being trusted
here (see tests/test_core.py) rather than just hand-checked once. An
analytical closed form was deliberately chosen over a hand-rolled finite-
difference PDE solver for this first working version: it is exact (no grid/
time-step tuning, no numerical dispersion of its own to confuse with the
physical dispersion being modeled), and it is the same solution EPA's own
BIOSCREEN screening model uses. Numerical (finite-difference / 2-D)
transport is the explicitly-deferred extension -- see README "What this
model does NOT do".

NUMERICAL STABILITY
--------------------
The naive formula above overflows in double precision whenever v*x/D_L is
large (which happens routinely at realistic field distances): exp(v*x/D_L)
explodes while erfc(b) vanishes, and inf*0 = nan even though their true
product is small and finite. The fix used in USGS/EPA screening codes is to
rewrite the second term with the scaled complementary error function
erfcx(z) = exp(z^2)*erfc(z), which SciPy computes directly without overflow:

    exp(v*x/D_L) * erfc(b)  =  exp(v*x/D_L - b**2) * erfcx(b)

tests/test_core.py cross-checks the stabilized and naive computations at
parameter values small enough that the naive version does not overflow, so a
broken rewrite would be caught immediately rather than silently shipping
wrong numbers.

DECAY
-----
NOT implemented. Every metal in this project is modeled with lambda=0, per
the checklist's own repeated warning that metals do not undergo first-order
decay the way an organic contaminant might, and this project's own prior
research explicitly concluded "lambda=0 is still correct... Arsenic does not
biodegrade into 'no arsenic'. Lead does not disappear." A retardation+decay
closed form exists in the literature (Bear 1979 / van Genuchten & Alves
1982), but conventions differ across sources on whether the decay constant
should apply to dissolved mass only or to dissolved+sorbed mass together,
and getting that wrong would silently ship incorrect physics. Rather than
guess, `retarded_transport_concentration` raises a clear error if a nonzero
decay is requested -- an honest gap, in the same spirit as the HONEST_GAP
data-quality flags used elsewhere in this project, not a silent limitation.
"""
import numpy as np
from scipy.special import erfc, erfcx

from units import SECONDS_PER_YEAR
from reference_data import MOLECULAR_DIFFUSION_M2_S, ALPHA_T_OVER_ALPHA_L, ALPHA_V_OVER_ALPHA_L


def xu_eckstein_dispersivity_m(transport_distance_m: float) -> float:
    """alpha_L = 0.83 * [log10(L)]^2.414, L = transport distance in meters.
    Xu, M. & Eckstein, Y. (1995), Ground Water 33(6), 905-908. This is the
    field-scale dispersivity estimator EPA's own online screening tool
    (CEAM) recommends; used here instead of a single assumed dispersivity
    constant. Only valid for L > 1 m; floors at L=1 m to avoid a negative/
    undefined result right at the source."""
    l = max(transport_distance_m, 1.0)
    return 0.83 * (np.log10(l)) ** 2.414


def rule_of_thumb_dispersivity_m(transport_distance_m: float, fraction: float = 0.1) -> float:
    """alpha_L = fraction * L (Gelhar 1993 simple rule of thumb, default
    10%), kept as an alternative/sensitivity option alongside Xu-Eckstein."""
    return fraction * max(transport_distance_m, 1.0)


def dispersion_coefficients(velocity_m_s: float, alpha_l_m: float,
                             molecular_diffusion_m2_s: float = MOLECULAR_DIFFUSION_M2_S) -> dict:
    """D_L = alpha_L * v + D_m ; D_T, D_V via the standard dispersivity-ratio
    simplification (checklist #8)."""
    d_l = alpha_l_m * velocity_m_s + molecular_diffusion_m2_s
    alpha_t = ALPHA_T_OVER_ALPHA_L * alpha_l_m
    alpha_v = ALPHA_V_OVER_ALPHA_L * alpha_l_m
    d_t = alpha_t * velocity_m_s + molecular_diffusion_m2_s
    d_v = alpha_v * velocity_m_s + molecular_diffusion_m2_s
    return dict(D_L=d_l, D_T=d_t, D_V=d_v, alpha_L=alpha_l_m, alpha_T=alpha_t, alpha_V=alpha_v)


def retarded_transport_concentration(
    x_m,
    t_years,
    velocity_m_s: float,
    d_l_m2_s: float,
    retardation_factor: float,
    c0_mg_l: float,
    decay_per_year: float = 0.0,
) -> np.ndarray:
    """
    Retarded Ogata-Banks solution, numerically stabilized (see module
    docstring). `x_m` and `t_years` may each be scalars, or one may be an
    array while the other is scalar (broadcasts normally); to get a full
    concentration-vs-distance-and-time grid, pass x_m as a column vector and
    t_years as a row vector (or use scenarios.concentration_grid, which does
    this for you).

    Returns C(x,t) in mg/L, broadcast-shaped from (x_m, t_years).
    """
    if decay_per_year != 0.0:
        raise NotImplementedError(
            "Decay (decay_per_year != 0) is an intentional HONEST GAP in this "
            "module -- see the transport.py module docstring's 'DECAY' section "
            "for why. Metals in this project are always run with decay_per_year=0."
        )

    x = np.asarray(x_m, dtype=float)
    t_s = np.asarray(t_years, dtype=float) * SECONDS_PER_YEAR
    x2, t2 = np.broadcast_arrays(x, t_s)
    x2 = np.array(x2, dtype=float)
    t2 = np.array(t2, dtype=float)

    R = retardation_factor
    v = velocity_m_s
    D_L = d_l_m2_s

    c_over_c0 = np.zeros_like(x2, dtype=float)
    valid = t2 > 0
    xv, tv = x2[valid], t2[valid]

    sqrt_term = np.sqrt(D_L * R * tv)

    a = (R * xv - v * tv) / (2.0 * sqrt_term)
    b = (R * xv + v * tv) / (2.0 * sqrt_term)
    term1 = erfc(a)
    # term2 = exp(v*x/D_L) * erfc(b), stabilized as exp(v*x/D_L - b**2) * erfcx(b)
    exponent = np.clip((v * xv / D_L) - b ** 2, -700, 700)
    term2 = np.exp(exponent) * erfcx(b)

    c_over_c0[valid] = 0.5 * (term1 + term2)
    # t=0 stays 0 by the initial condition C(x,0)=0 (already the array default)

    return c0_mg_l * c_over_c0


def concentration_time_series(x_m: float, t_years: np.ndarray, velocity_m_s: float,
                               d_l_m2_s: float, retardation_factor: float, c0_mg_l: float) -> np.ndarray:
    """Convenience wrapper: concentration at one fixed distance, over a
    vector of times. Returns an array the same length as t_years."""
    t_years = np.asarray(t_years, dtype=float)
    return retarded_transport_concentration(
        np.full_like(t_years, x_m), t_years, velocity_m_s, d_l_m2_s, retardation_factor, c0_mg_l
    )


def concentration_distance_profile(x_m: np.ndarray, t_years: float, velocity_m_s: float,
                                    d_l_m2_s: float, retardation_factor: float, c0_mg_l: float) -> np.ndarray:
    """Convenience wrapper: concentration vs. distance, at one fixed time.
    Returns an array the same length as x_m."""
    x_m = np.asarray(x_m, dtype=float)
    return retarded_transport_concentration(
        x_m, np.full_like(x_m, t_years), velocity_m_s, d_l_m2_s, retardation_factor, c0_mg_l
    )


def peak_arrival(x_m: float, velocity_m_s: float, d_l_m2_s: float,
                  retardation_factor: float, c0_mg_l: float,
                  t_search_years: np.ndarray) -> dict:
    """Checklist #39: 'maximum concentration' and 'time of maximum
    concentration' at a fixed receptor distance, found by direct search over
    a supplied time grid (the retarded Ogata-Banks profile at fixed x rises
    monotonically toward C0 as t->inf for a constant-source boundary, so the
    practical 'peak' within a finite planning horizon is just its value at
    the last time searched; this function still reports it explicitly rather
    than assuming the caller knows that)."""
    c = concentration_time_series(x_m, t_search_years, velocity_m_s, d_l_m2_s, retardation_factor, c0_mg_l)
    i_max = int(np.argmax(c))
    return dict(max_concentration_mg_l=float(c[i_max]), time_of_max_years=float(t_search_years[i_max]),
                concentration_series_mg_l=c)

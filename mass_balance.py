"""
mass_balance.py
================
Checklist sections 17, 19-21, 49: the mass-balance check the spec calls "one
of the best ways to catch mistakes." Every scenario run in this project
computes this at every reported time step; a run that fails it is flagged in
the output tables rather than silently reported as if it were trustworthy.

For metals (lambda=0, no transformation), conservation is:
    M_initial ~= M_remaining_in_ash + M_dissolved + M_sorbed + M_exported

`M_transformed` is kept as an explicit input (always 0.0 for the metals runs
in this project) so the same function works unmodified if this package is
later extended to a contaminant where transformation is real.
"""
from dataclasses import dataclass

MASS_BALANCE_ERROR_THRESHOLD_PCT = 2.0  # checklist doesn't specify a number;
# 2% is a conventional, documented engineering-screening tolerance for a
# closed-form/analytical mass balance (finite-difference numerical models
# often use a looser default, e.g. 5%, because grid discretization error is
# unavoidable there -- ours is looser only because floating point summation
# of many small time steps can accumulate to a few 0.1%s, not because the
# physics is approximate).


@dataclass
class MassBalanceResult:
    initial_mg: float
    accounted_mg: float
    percent_error: float
    passes: bool
    components: dict


def mass_balance_check(
    initial_mg: float,
    remaining_in_ash_mg: float,
    dissolved_mg: float,
    sorbed_mg: float,
    exported_mg: float,
    transformed_mg: float = 0.0,
    threshold_pct: float = MASS_BALANCE_ERROR_THRESHOLD_PCT,
) -> MassBalanceResult:
    """checklist #49: %MBE = |M_initial - M_accounted| / M_initial * 100."""
    accounted = remaining_in_ash_mg + dissolved_mg + sorbed_mg + exported_mg + transformed_mg
    if initial_mg == 0:
        pct_error = 0.0 if accounted == 0 else float("inf")
    else:
        pct_error = abs(initial_mg - accounted) / initial_mg * 100.0
    return MassBalanceResult(
        initial_mg=initial_mg,
        accounted_mg=accounted,
        percent_error=pct_error,
        passes=pct_error <= threshold_pct,
        components=dict(remaining_in_ash_mg=remaining_in_ash_mg, dissolved_mg=dissolved_mg,
                         sorbed_mg=sorbed_mg, exported_mg=exported_mg, transformed_mg=transformed_mg),
    )


def dissolved_mass_in_domain_mg(concentration_mg_l_profile, dx_m: float, cross_section_area_m2: float,
                                 porosity: float) -> float:
    """checklist #18: M_water = integral of C(x) * n * A dx over the modeled
    domain (mg), via the trapezoid rule on a supplied concentration profile
    (mg/L at evenly spaced points along the domain, spacing dx_m).

    Each meter of domain length holds (cross_section_area_m2 * porosity)
    m^3 of pore water, i.e. that many *1000 liters; multiplying the
    concentration profile by that constant before integrating turns
    integral(C dx) [mg/L * m] into mg directly.
    """
    import numpy as np
    c = np.asarray(concentration_mg_l_profile, dtype=float)
    liters_of_pore_water_per_meter = cross_section_area_m2 * porosity * 1000.0
    return float(np.trapezoid(c * liters_of_pore_water_per_meter, dx=dx_m))


def sorbed_mass_mg(dissolved_mass_mg: float, retardation_factor: float) -> float:
    """checklist #19: sorbed mass implied by the retardation factor, given
    the dissolved mass already computed for the same domain slice:
    R = 1 + M_sorbed/M_dissolved (equivalent form of S=Kd*C combined with
    M_sorbed = S*M_soil), so M_sorbed = (R-1) * M_dissolved."""
    return (retardation_factor - 1.0) * dissolved_mass_mg

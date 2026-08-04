"""
source_term.py
===============
Checklist sections 1, 15, 16: how much contaminant is in the ash, how much of
it can plausibly become mobile, and at what concentration it enters pore
water. Two release models are implemented, matching this project's own
earlier-established design (see chat_1_ in the source project notes):

  1. `leach_simple`   -- M_leached = M0 * f_L               (checklist #15, "simplest method")
  2. `leach_kinetic`  -- first-order source depletion, M(t) = M0 * exp(-k_L t)  (checklist #15, "more advanced method")

`project_source_release` combines them the way the project's own prior
design specified: at each time step, release is the MINIMUM of the
kinetic-release estimate and a concentration-controlled estimate (you cannot
release more mass in a time step than the leachate concentration times the
water that moved through), and never more than what remains. This prevents
the model from inventing mass out of nowhere when k_L is set too high.
"""
from dataclasses import dataclass
import numpy as np

from units import mg_to_kg, kg_to_mg, m3_to_liters


@dataclass
class AshSourceGeometry:
    pond_area_m2: float
    ash_thickness_m: float
    ash_bulk_density_kg_m3: float
    ash_porosity: float

    @property
    def ash_volume_m3(self) -> float:
        return self.pond_area_m2 * self.ash_thickness_m

    @property
    def ash_mass_kg(self) -> float:
        return self.ash_volume_m3 * self.ash_bulk_density_kg_m3

    @property
    def pore_water_volume_m3(self) -> float:
        """Checklist #21: water stored inside the saturated ash itself."""
        return self.ash_porosity * self.ash_volume_m3


@dataclass
class ContaminantSource:
    name: str
    solid_concentration_mg_kg: float   # C_solid, checklist #1
    leaching_fraction: float           # f_L, dimensionless [0-1]
    data_quality: str = "literature"


def initial_contaminant_mass_mg(mass_ash_kg: float, c_solid_mg_kg: float) -> float:
    """M0 = M_ash * C_solid   (checklist #1)."""
    return mass_ash_kg * c_solid_mg_kg


def leachable_mass_mg(m0_mg: float, leaching_fraction: float) -> float:
    """M_L = M0 * f_L   (checklist #1)."""
    if not 0.0 <= leaching_fraction <= 1.0:
        raise ValueError(f"leaching_fraction must be in [0,1], got {leaching_fraction}")
    return m0_mg * leaching_fraction


def leach_simple(m0_mg: float, leaching_fraction: float) -> float:
    """Checklist #15, simplest method: M_leached = M0 * f_L, released once."""
    return leachable_mass_mg(m0_mg, leaching_fraction)


def leach_kinetic_remaining(m0_mg: float, k_l_per_year: float, t_years: np.ndarray) -> np.ndarray:
    """Checklist #15, advanced method: M(t) = M0 * exp(-k_L t).
    Returns the mass STILL IN THE ASH (not yet leached) at each time in t_years.
    """
    t_years = np.asarray(t_years, dtype=float)
    return m0_mg * np.exp(-k_l_per_year * t_years)


def calibrate_k_l_from_target_fraction(leaching_fraction_target: float, horizon_years: float) -> float:
    """
    There is no calibration dataset in this project's inputs, so k_L cannot be
    derived from field/column data (checklist #15's own warning: 'never use
    an uncalibrated k_L for a regulatory conclusion'). As an explicit,
    documented stand-in, k_L is instead solved so that the KINETIC model
    reproduces the same total leached fraction as the SIMPLE model over one
    chosen horizon: 1 - exp(-k_L * horizon) = leaching_fraction_target.
    This keeps the two release models mutually consistent for comparison
    rather than pretending k_L is a measured rate constant.
    """
    if not 0.0 < leaching_fraction_target < 1.0:
        raise ValueError("leaching_fraction_target must be strictly between 0 and 1")
    return -np.log(1.0 - leaching_fraction_target) / horizon_years


def project_source_release(
    m0_mg: float,
    k_l_per_year: float,
    leachate_concentration_mg_l: float,
    infiltration_m3_per_year: float,
    years: np.ndarray,
) -> dict:
    """
    Year-by-year source depletion combining a kinetic estimate with a
    concentration-controlled ceiling, following this project's own earlier
    design (see chat_1_ / source notes `project_constituent`):

        kinetic_kg   = remaining * (1 - exp(-k_L * dt))
        conc_limited_kg = C_leachate[mg/L] * V_infiltration[m3/yr] * 1000[L/m3] / 1e6[mg/kg]
        released_kg = min(remaining, kinetic_kg, conc_limited_kg)

    Returns arrays (same length as `years`) for remaining mass, released mass
    per step, and cumulative released mass, all in mg for consistency with
    the rest of this package's mg-based mass balance.
    """
    years = np.asarray(years, dtype=float)
    dt = np.diff(years, prepend=0.0)
    remaining = m0_mg
    remaining_hist, released_hist, cumulative_hist = [], [], []
    conc_limited_mg_per_year = leachate_concentration_mg_l * infiltration_m3_per_year * 1000.0
    cumulative = 0.0
    for step_dt in dt:
        if step_dt <= 0:
            remaining_hist.append(remaining)
            released_hist.append(0.0)
            cumulative_hist.append(cumulative)
            continue
        kinetic_mg = remaining * (1.0 - np.exp(-k_l_per_year * step_dt))
        conc_limited_mg = conc_limited_mg_per_year * step_dt
        released = min(remaining, kinetic_mg, conc_limited_mg)
        remaining -= released
        cumulative += released
        remaining_hist.append(remaining)
        released_hist.append(released)
        cumulative_hist.append(cumulative)
    return dict(
        years=years,
        remaining_mass_mg=np.array(remaining_hist),
        released_mass_mg=np.array(released_hist),
        cumulative_released_mg=np.array(cumulative_hist),
    )


def leachate_concentration_mg_l(leached_mass_mg: float, water_volume_m3: float) -> float:
    """C_leachate = M_leached / V_water   (checklist #16), with the explicit
    m3->L conversion the spec calls out as the most common place to
    introduce a unit error."""
    water_volume_l = m3_to_liters(water_volume_m3)
    if water_volume_l <= 0:
        raise ValueError("water_volume_m3 must be positive")
    return leached_mass_mg / water_volume_l

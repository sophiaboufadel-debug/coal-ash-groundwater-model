"""
sorption.py
===========
Checklist sections 10-11, 26: sorption (Kd), retardation, and a directional
pH adjustment. Metal-specific Kd values live in reference_data.py; this
module only implements the arithmetic and the pH-directionality rule.

UNIT NOTE (read this before changing anything): the retardation formula
    R = 1 + (rho_b * Kd) / n
is only dimensionally clean when rho_b and Kd are expressed in a matched
pair -- kg/L bulk density with L/kg Kd (numerically identical to the more
familiar g/cm^3 with mL/g). Mixing kg/m^3 bulk density with L/kg Kd silently
introduces a factor-of-1000 error, which is exactly the kind of mistake the
original checklist (#45) warns will "destroy the model." This module always
converts bulk density to kg/L internally -- see units.bulk_density_kg_m3_to_kg_l.
"""
from dataclasses import dataclass

from units import bulk_density_kg_m3_to_kg_l
from reference_data import (
    METAL_KD_L_KG, OXYANION_FORMING_METALS, CATIONIC_METALS,
    SCREENING_GROUNDWATER_PH, HONEST_GAP_RANK_CONTAMINANTS, DataQuality,
)


def retardation_factor(bulk_density_kg_m3: float, kd_l_kg: float, porosity: float) -> float:
    """R = 1 + (rho_b * Kd) / n   (checklist #11)."""
    if porosity <= 0:
        raise ValueError("porosity must be positive")
    if kd_l_kg < 0:
        raise ValueError("Kd must be non-negative")
    rho_b_kg_l = bulk_density_kg_m3_to_kg_l(bulk_density_kg_m3)
    return 1.0 + (rho_b_kg_l * kd_l_kg) / porosity


def retarded_velocity(pore_velocity: float, retardation_factor_value: float) -> float:
    """v_c = v / R_f   (checklist #11)."""
    if retardation_factor_value <= 0:
        raise ValueError("retardation_factor_value must be positive")
    return pore_velocity / retardation_factor_value


def retarded_travel_time_years(distance_m: float, velocity_m_s: float, retardation_factor_value: float) -> float:
    """t_c = x * R_f / v   (checklist #11)."""
    from units import SECONDS_PER_YEAR
    if velocity_m_s <= 0:
        raise ValueError("velocity_m_s must be positive")
    return (distance_m * retardation_factor_value / velocity_m_s) / SECONDS_PER_YEAR


@dataclass
class KdLookupResult:
    kd_l_kg: float
    quality: str
    ph_adjustment_note: str


def get_metal_kd(metal: str, ph: float = SCREENING_GROUNDWATER_PH,
                  low_high_range: tuple | None = None) -> KdLookupResult:
    """
    Returns a screening Kd for `metal`, with an explicit, DIRECTIONAL
    pH note rather than a fabricated pH-response curve (checklist #26 asks
    for "experimentally measured contaminant-specific relationships rather
    than inventing a multiplier" -- this project's compiled research did not
    turn up a clean per-metal pH-Kd curve for most of these elements, so this
    function is honest about giving direction, not magnitude, when pH is
    off the EPA screening default of 6.8).

    If `low_high_range` is given (used by uncertainty.py for Monte Carlo /
    sensitivity runs), returns a Kd sampled/scaled from that range instead of
    the point default.
    """
    if metal not in METAL_KD_L_KG:
        raise KeyError(f"No Kd entry for '{metal}'. Add it to reference_data.METAL_KD_L_KG.")
    kd_default, quality = METAL_KD_L_KG[metal]
    kd = kd_default if low_high_range is None else low_high_range

    note = ""
    if abs(ph - SCREENING_GROUNDWATER_PH) > 0.3:
        if metal in OXYANION_FORMING_METALS:
            direction = "more mobile (lower effective Kd)" if ph > SCREENING_GROUNDWATER_PH else \
                        "less mobile (higher effective Kd)"
            note = (f"pH={ph:.1f} vs. screening default {SCREENING_GROUNDWATER_PH}: {metal} is an "
                    f"oxyanion-forming element, expected DIRECTION is {direction} at higher pH "
                    "(CaO/alkalinity mechanism, see reference_data.py Section 6). Magnitude not quantified.")
        elif metal in CATIONIC_METALS:
            direction = "less mobile (higher effective Kd)" if ph > SCREENING_GROUNDWATER_PH else \
                        "more mobile (lower effective Kd)"
            note = (f"pH={ph:.1f} vs. screening default {SCREENING_GROUNDWATER_PH}: {metal} behaves as a "
                    f"simple cationic metal, expected DIRECTION is {direction} at higher pH. "
                    "Magnitude not quantified.")
    return KdLookupResult(kd_l_kg=kd, quality=quality, ph_adjustment_note=note)


def is_honest_gap(metal: str) -> bool:
    """True if this project's source research explicitly flagged `metal` as
    lacking a clean rank-comparative / literature leaching dataset (see
    reference_data.HONEST_GAP_RANK_CONTAMINANTS). Callers should propagate
    this flag into any output table rather than silently treating the metal
    like every other one."""
    return metal in HONEST_GAP_RANK_CONTAMINANTS

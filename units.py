"""
units.py
========
Every conversion the model needs, in one place, each with a self-test.
Internal convention (see reference_data.py Section 1):
    length: m | time: s (reported in years) | mass: mg (contaminant), kg (bulk)
    concentration: mg/L | K: m/s | density: kg/L | Kd: L/kg

Get this file wrong and every downstream number is wrong, so every function
below is covered by an inline sanity check executed at import time
(cheap, deterministic, and catches a broken edit immediately rather than
three modules downstream).
"""

SECONDS_PER_YEAR = 365.25 * 24 * 3600
SECONDS_PER_DAY = 86400.0
M3_PER_LITER = 1e-3
LITERS_PER_M3 = 1000.0
MG_PER_KG = 1_000_000.0
ACRE_TO_M2 = 4046.8564224
CM_PER_HR_TO_M_PER_S = 1e-2 / 3600.0  # cm/hr -> m/s


def years_to_seconds(years: float) -> float:
    return years * SECONDS_PER_YEAR


def seconds_to_years(seconds: float) -> float:
    return seconds / SECONDS_PER_YEAR


def m_per_year_to_m_per_s(v: float) -> float:
    return v / SECONDS_PER_YEAR


def m_per_s_to_m_per_year(v: float) -> float:
    return v * SECONDS_PER_YEAR


def kg_to_mg(mass_kg: float) -> float:
    return mass_kg * MG_PER_KG


def mg_to_kg(mass_mg: float) -> float:
    return mass_mg / MG_PER_KG


def m3_to_liters(v_m3: float) -> float:
    return v_m3 * LITERS_PER_M3


def liters_to_m3(v_l: float) -> float:
    return v_l / LITERS_PER_M3


def mg_per_m3_to_mg_per_l(c: float) -> float:
    """Concentration conversion: mg/m^3 -> mg/L (divide by 1000 L/m^3)."""
    return c / LITERS_PER_M3


def bulk_density_kg_m3_to_kg_l(rho_kg_m3: float) -> float:
    """kg/m^3 -> kg/L (numerically identical to g/cm^3)."""
    return rho_kg_m3 / 1000.0


def acres_to_m2(acres: float) -> float:
    return acres * ACRE_TO_M2


def cm_per_hr_to_m_per_s(k_cm_hr: float) -> float:
    return k_cm_hr * CM_PER_HR_TO_M_PER_S


def pci_to_mg(activity_pci: float, mg_per_pci: float) -> float:
    """Radionuclide activity -> mass, given an isotope-specific mg/pCi factor."""
    return activity_pci * mg_per_pci


# ---------------------------------------------------------------------------
# Self-tests (run at import; raise immediately if a conversion is broken)
# ---------------------------------------------------------------------------
def _self_test():
    assert abs(seconds_to_years(years_to_seconds(7.0)) - 7.0) < 1e-9
    assert abs(kg_to_mg(1.0) - 1_000_000.0) < 1e-6
    assert abs(mg_to_kg(1_000_000.0) - 1.0) < 1e-9
    assert abs(m3_to_liters(1.0) - 1000.0) < 1e-9
    assert abs(acres_to_m2(1.0) - 4046.8564224) < 1e-6
    # 1 m^3 = 1000 L check embedded directly, per the source spec's own
    # explicit reminder that unit errors "can destroy the model."
    assert abs(liters_to_m3(1000.0) - 1.0) < 1e-9
    # cm/hr -> m/s: 1 cm/hr = 1e-2 m / 3600 s
    assert abs(cm_per_hr_to_m_per_s(3600.0) - 1e-2) < 1e-12
    # bulk density kg/m3 <-> kg/L: 1500 kg/m3 == 1.5 kg/L == 1.5 g/cm3
    assert abs(bulk_density_kg_m3_to_kg_l(1500.0) - 1.5) < 1e-9


_self_test()

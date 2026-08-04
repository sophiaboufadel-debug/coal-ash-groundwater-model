"""
tests/test_core.py
===================
Not a formality -- this is what stands between "the code runs" and "the
numbers are right." Run with: python3 tests/test_core.py
(kept dependency-free of pytest so it runs in any environment; exits nonzero
on first failure).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import numpy as np
from scipy.special import erfc

import units
import flow
import sorption
import source_term
import transport
import mass_balance
import regulatory


def _check(name, condition):
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {name}")
    if not condition:
        raise AssertionError(name)


def test_transport_formula_derivation():
    print("test_transport_formula_derivation")
    import sympy as sp
    x, t, v, R, D = sp.symbols("x t v R D", positive=True)
    v_eff, D_eff = v / R, D / R
    a_direct = (x - v_eff * t) / (2 * sp.sqrt(D_eff * t))
    a_claimed = (R * x - v * t) / (2 * sp.sqrt(D * R * t))
    b_direct = (x + v_eff * t) / (2 * sp.sqrt(D_eff * t))
    b_claimed = (R * x + v * t) / (2 * sp.sqrt(D * R * t))
    exp_direct = v_eff * x / D_eff
    exp_claimed = v * x / D
    _check("a-argument substitution is algebraically exact", sp.simplify(a_direct - a_claimed) == 0)
    _check("b-argument substitution is algebraically exact", sp.simplify(b_direct - b_claimed) == 0)
    _check("exponent substitution is algebraically exact", sp.simplify(exp_direct - exp_claimed) == 0)


def test_transport_numerical_stability_matches_naive():
    print("test_transport_numerical_stability_matches_naive")
    # Small-enough parameters that the naive (non-erfcx) formula does NOT
    # overflow, so we can compare stabilized vs. naive directly.
    v, D_L, R, C0 = 1e-8, 1e-8, 3.0, 10.0
    x = np.array([1.0, 5.0, 20.0])
    t_years = 5.0
    t_s = t_years * units.SECONDS_PER_YEAR

    stabilized = transport.concentration_distance_profile(x, t_years, v, D_L, R, C0)

    a = (R * x - v * t_s) / (2 * np.sqrt(D_L * R * t_s))
    b = (R * x + v * t_s) / (2 * np.sqrt(D_L * R * t_s))
    naive = 0.5 * C0 * (erfc(a) + np.exp(v * x / D_L) * erfc(b))

    _check("stabilized matches naive erfc/exp computation (max diff < 1e-9)",
           np.max(np.abs(stabilized - naive)) < 1e-9)


def test_transport_boundary_and_initial_conditions():
    print("test_transport_boundary_and_initial_conditions")
    v, D_L, R, C0 = 5e-8, 5e-8, 2.5, 10.0

    # Initial condition: C(x,0) = 0 for all x > 0
    c_t0 = transport.concentration_distance_profile(np.array([1.0, 100.0, 1000.0]), 0.0, v, D_L, R, C0)
    _check("C(x, t=0) = 0 everywhere", np.allclose(c_t0, 0.0))

    # Boundary condition: C(0, t) -> C0 for t > 0 (within a small distance of
    # the source, since x=0 exactly gives erfc(0)=1 term1 and a finite term2)
    c_x0 = transport.concentration_time_series(0.0, np.array([1.0, 10.0, 50.0]), v, D_L, R, C0)
    _check("C(x=0, t) approx C0 for all t>0", np.allclose(c_x0, C0, rtol=1e-6))

    # Long-time limit: at fixed finite x, C(x,t) -> C0 as t -> infinity
    c_late = transport.concentration_time_series(50.0, np.array([1e6]), v, D_L, R, C0)[0]
    _check("C(x, t->inf) -> C0 at fixed finite x", abs(c_late - C0) / C0 < 1e-3)

    # Concentration should never exceed the source concentration or go negative
    x_grid = np.linspace(0.1, 2000, 40)
    t_grid = np.linspace(0.1, 50, 20)
    grid = transport.retarded_transport_concentration(
        x_grid[:, None], t_grid[None, :], v, D_L, R, C0
    )
    _check("0 <= C(x,t) <= C0 everywhere on a grid", np.all(grid >= -1e-9) and np.all(grid <= C0 * (1 + 1e-9)))

    # Monotonic in time at fixed distance (constant boundary source => plume
    # only builds up, never recedes, in this idealized constant-C0 model)
    c_series = transport.concentration_time_series(200.0, np.linspace(0.5, 40, 30), v, D_L, R, C0)
    _check("C(x, t) is non-decreasing in time at fixed x", np.all(np.diff(c_series) >= -1e-9))


def test_retardation_factor():
    print("test_retardation_factor")
    r_no_sorption = sorption.retardation_factor(bulk_density_kg_m3=1500, kd_l_kg=0.0, porosity=0.4)
    _check("R = 1 when Kd = 0 (no sorption -> no retardation)", abs(r_no_sorption - 1.0) < 1e-9)

    r_lead = sorption.retardation_factor(bulk_density_kg_m3=1500, kd_l_kg=900.0, porosity=0.4)
    # rho_b in kg/L = 1.5; R = 1 + (1.5*900)/0.4 = 1 + 3375 = 3376
    _check("R for Lead (Kd=900) matches hand calc (3376)", abs(r_lead - 3376.0) < 1e-6)

    r_high_kd = sorption.retardation_factor(bulk_density_kg_m3=1500, kd_l_kg=1_800_000.0, porosity=0.4)
    _check("R increases with Kd (Cr-III >> Lead)", r_high_kd > r_lead)

    v_retarded = sorption.retarded_velocity(pore_velocity=1e-7, retardation_factor_value=100.0)
    _check("retarded velocity = pore velocity / R", abs(v_retarded - 1e-9) < 1e-15)


def test_darcy_and_travel_time():
    print("test_darcy_and_travel_time")
    q = flow.darcy_flux(k_m_s=5e-7, gradient=0.01)
    _check("Darcy flux q=Ki matches hand calc", abs(q - 5e-9) < 1e-15)

    v = flow.pore_velocity(darcy_flux_m_s=5e-9, effective_porosity=0.541)
    v_year = v * units.SECONDS_PER_YEAR
    _check("pore velocity ~0.29 m/yr matches this project's own worked example",
           abs(v_year - 0.2917) < 0.001)

    t = flow.unretarded_travel_time_years(distance_m=10.0, velocity_m_s=v)
    _check("10 m unretarded travel time ~34.3 yr matches this project's own worked example",
           abs(t - 34.3) < 0.1)


def test_source_term_mass_conservation():
    print("test_source_term_mass_conservation")
    m0 = source_term.initial_contaminant_mass_mg(mass_ash_kg=1e6, c_solid_mg_kg=50.0)
    _check("M0 = M_ash * C_solid", abs(m0 - 5e7) < 1e-3)

    m_l = source_term.leachable_mass_mg(m0, leaching_fraction=0.03)
    _check("leachable mass = 3% of M0", abs(m_l - 0.03 * m0) < 1e-3)

    # project_source_release: cumulative released mass must never exceed M0,
    # and remaining + cumulative must always equal M0 (mass conservation).
    years = np.arange(0, 21)
    k_l = source_term.calibrate_k_l_from_target_fraction(0.20, horizon_years=20)
    result = source_term.project_source_release(
        m0_mg=m0, k_l_per_year=k_l, leachate_concentration_mg_l=1000.0,
        infiltration_m3_per_year=50.0, years=years,
    )
    total = result["remaining_mass_mg"] + result["cumulative_released_mg"]
    _check("remaining + cumulative released = M0 at every step (mass balance)",
           np.allclose(total, m0, rtol=1e-9))
    _check("cumulative released is non-decreasing", np.all(np.diff(result["cumulative_released_mg"]) >= -1e-6))
    _check("cumulative released never exceeds M0", result["cumulative_released_mg"][-1] <= m0 + 1e-6)


def test_mass_balance_module():
    print("test_mass_balance_module")
    check = mass_balance.mass_balance_check(
        initial_mg=1000.0, remaining_in_ash_mg=600.0, dissolved_mg=150.0,
        sorbed_mg=200.0, exported_mg=50.0, transformed_mg=0.0,
    )
    _check("balanced case reports ~0% error", check.percent_error < 1e-6)
    _check("balanced case passes threshold check", check.passes)

    bad = mass_balance.mass_balance_check(
        initial_mg=1000.0, remaining_in_ash_mg=600.0, dissolved_mg=150.0,
        sorbed_mg=200.0, exported_mg=50.0, transformed_mg=100.0,  # +100 extra mass from nowhere
    )
    _check("unbalanced case (extra mass) reports nonzero error", bad.percent_error > 5.0)
    _check("unbalanced case fails threshold check", not bad.passes)


def test_regulatory_exceedance():
    print("test_regulatory_exceedance")
    result = regulatory.exceedance_ratio("Arsenic", predicted_mg_l=0.05)
    _check("exceedance ratio = predicted / standard", abs(result.ratio - 5.0) < 1e-6)
    _check("ratio > 1 flags exceedance", result.exceeds)

    result_ok = regulatory.exceedance_ratio("Arsenic", predicted_mg_l=0.001)
    _check("below-standard case does not flag exceedance", not result_ok.exceeds)


def test_molybdenum_and_cobalt_wired_in():
    print("test_molybdenum_and_cobalt_wired_in")
    import reference_data as rd

    # Kd lookups must resolve without KeyError, and match the documented
    # source values (EPA/600/R-05/074 for Mo; EPA RSL table for Co).
    mo_kd = sorption.get_metal_kd("Molybdenum")
    _check("Molybdenum Kd = 12.6 L/kg (EPA/600/R-05/074 Table 3 median)",
           abs(mo_kd.kd_l_kg - 12.6) < 1e-9)
    co_kd = sorption.get_metal_kd("Cobalt")
    _check("Cobalt Kd = 45.0 L/kg (EPA RSL default)", abs(co_kd.kd_l_kg - 45.0) < 1e-9)

    # Bulk ash concentrations must resolve without KeyError.
    _check("Molybdenum has a bulk ash concentration on file (HONEST_GAP)",
           "Molybdenum" in rd.BULK_ASH_CONCENTRATION_MG_KG)
    _check("Cobalt has a bulk ash concentration on file",
           "Cobalt" in rd.BULK_ASH_CONCENTRATION_MG_KG)

    # The 4 new CCR groundwater protection standards must all be on file.
    for metal, expected_mg_l in [("Cobalt", 0.006), ("Lithium", 0.04),
                                  ("Molybdenum", 0.1), ("Thallium", 0.002)]:
        _check(f"{metal} regulatory standard = {expected_mg_l} mg/L (40 CFR 257.95(h))",
               abs(rd.REGULATORY_STANDARDS[metal].value_mg_l - expected_mg_l) < 1e-9)

    # Full end-to-end scenario run must complete without error for both new
    # active metals and produce physically sane (finite, non-negative)
    # output -- this is the same pipeline run_model.py drives for every
    # region x metal x tier combination.
    import scenarios as sc
    for metal in ["Molybdenum", "Cobalt"]:
        result = sc.run_scenario("South", metal, "Bituminous", "central")
        _check(f"{metal} scenario retardation factor > 1", result.retardation_factor > 1.0)
        _check(f"{metal} scenario c0 is finite and non-negative",
               np.isfinite(result.c0_mg_l) and result.c0_mg_l >= 0.0)
        for x, series in result.concentration_by_receptor.items():
            for year, c in series.items():
                _check(f"{metal} concentration at {x}m/yr{year} is finite and in [0, c0]",
                       np.isfinite(c) and -1e-9 <= c <= result.c0_mg_l * (1 + 1e-6))


def run_all():
    tests = [
        test_transport_formula_derivation,
        test_transport_numerical_stability_matches_naive,
        test_transport_boundary_and_initial_conditions,
        test_retardation_factor,
        test_darcy_and_travel_time,
        test_source_term_mass_conservation,
        test_mass_balance_module,
        test_regulatory_exceedance,
        test_molybdenum_and_cobalt_wired_in,
    ]
    for t in tests:
        t()
    print("\nALL TESTS PASSED")


if __name__ == "__main__":
    run_all()

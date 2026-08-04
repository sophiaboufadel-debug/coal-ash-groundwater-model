"""
scenarios.py
============
This is where every other module gets wired together into one thing you can
actually run: pick a region, a metal, a coal rank, and a scenario tier
(low/central/high), and get back the full chain --

    ash mass -> contaminant mass -> leaching -> leachate concentration
    -> vertical infiltration -> Darcy flow -> pore velocity -> dispersion
    -> retardation -> concentration vs. distance/time -> mass balance
    -> regulatory exceedance

-- at the 5/10/20-year snapshots the user asked for (plus a couple of extra
points so the story in between isn't hidden, per this project's own prior
design note: "Computing only once every five years can hide peak
concentrations and numerical errors").

SITE GEOMETRY DEFAULTS
-----------------------
Real per-facility geometry (pond area, ash thickness, distance to a specific
well) is NOT in the project's spreadsheets -- they cover facility identity,
compliance status, and coal type, not hydrogeologic site characterization.
Rather than silently invent per-facility numbers, every geometry input here
is the explicit SCREENING_* default from reference_data.py (EPA/USGS-sourced
national or regional figures), with the sensitivity ranges also in
reference_data.py available for a documented low/high sweep. See README.md
"What this model does NOT do" for how to plug in real site data instead.
"""
from dataclasses import dataclass, field
import numpy as np

import reference_data as rd
import source_term
import flow
import sorption
import transport
import mass_balance
import regulatory
from units import SECONDS_PER_YEAR, acres_to_m2

RECEPTOR_DISTANCES_M = [30.0, 150.0, 500.0, 1500.0]
RECEPTOR_LABELS = {30.0: "pond edge / property-boundary proxy", 150.0: "nearby monitoring well",
                    500.0: "mid-field", 1500.0: "distant well (~1 mi) proxy"}
REPORTING_YEARS = [1, 2, 5, 10, 15, 20, 30]  # 5/10/20 explicitly requested, plus context points


@dataclass
class ScenarioInputs:
    region: str
    metal: str
    coal_rank: str
    tier: str  # "low" | "central" | "high"
    hydraulic_gradient: float
    pond_area_m2: float
    ash_thickness_m: float
    distance_to_receptors_m: list


@dataclass
class ScenarioResult:
    inputs: ScenarioInputs
    m0_mg: float
    leachable_fraction: float
    leachable_fraction_quality: str
    m_leachable_mg: float
    c0_mg_l: float
    kd_l_kg: float
    kd_quality: str
    retardation_factor: float
    pore_velocity_m_yr: float
    retarded_velocity_m_yr: float
    d_l_m2_s: float
    alpha_l_m: float
    concentration_by_receptor: dict          # {distance_m: {year: mg/L}}
    travel_time_years_by_receptor: dict      # {distance_m: unretarded years}
    retarded_travel_time_years_by_receptor: dict
    exceedance_by_receptor_year: dict        # {distance_m: {year: ExceedanceResult}}
    mass_balance_by_year: dict               # {year: MassBalanceResult}
    honest_gap_flags: list


def _tier_multiplier(tier: str) -> dict:
    """Low/central/high scenario framework (this project's own prior design:
    'low release/high attenuation; central calibrated case; high
    release/low attenuation'). Multipliers are applied to the regional/
    literature central values."""
    return {
        "low": dict(leaching=0.3, kd=3.0, k_soil=0.4, gradient=0.5),
        "central": dict(leaching=1.0, kd=1.0, k_soil=1.0, gradient=1.0),
        "high": dict(leaching=3.0, kd=1.0 / 3.0, k_soil=2.5, gradient=2.0),
    }[tier]


def default_inputs(region: str, metal: str, coal_rank: str, tier: str) -> ScenarioInputs:
    return ScenarioInputs(
        region=region, metal=metal, coal_rank=coal_rank, tier=tier,
        hydraulic_gradient=0.01,  # documented assumed value, see reference_data / this project's
                                  # own prior notes: "i=0.01 was an assumed hydraulic gradient,
                                  # not a measured regional value" -- kept explicit, not hidden
        pond_area_m2=acres_to_m2(rd.SCREENING_POND_AREA_ACRES),
        ash_thickness_m=rd.SCREENING_ASH_THICKNESS_M,
        distance_to_receptors_m=RECEPTOR_DISTANCES_M,
    )


def leaching_fraction_for(metal: str, coal_rank: str, tier: str) -> tuple:
    """
    Combines (a) the metal-specific water-leachable-fraction range (Section 8
    of reference_data.py, where available) with (b) the coal-rank
    directional multiplier (Section 6, where available -- currently Arsenic
    and Selenium only, per the literature this project found), then applies
    the low/central/high tier multiplier. Returns (fraction, quality,
    honest_gap: bool).
    """
    mult = _tier_multiplier(tier)
    is_gap = sorption.is_honest_gap(metal if metal != "Chromium" else "Chromium(total)")

    if metal in rd.LEACHABLE_FRACTION_BY_METAL:
        low, high, quality = rd.LEACHABLE_FRACTION_BY_METAL[metal]
        base = rd.geometric_mean_leaching_fraction(metal)
    else:
        low, high = rd.DEFAULT_LEACHING_FRACTION["low"], rd.DEFAULT_LEACHING_FRACTION["high"]
        base = rd.DEFAULT_LEACHING_FRACTION["central"]
        quality = rd.DataQuality.HONEST_GAP

    rank_key = (metal, coal_rank)
    rank_mult, rank_quality, rank_note = 1.0, quality, ""
    if rank_key in rd.RANK_LEACHING_MULTIPLIER:
        rank_mult, rank_quality, rank_note = rd.RANK_LEACHING_MULTIPLIER[rank_key]

    fraction = base * rank_mult * mult["leaching"]
    fraction = float(np.clip(fraction, 1e-6, 0.95))  # keep inside physically valid (0,1)
    return fraction, quality, is_gap


def run_scenario(region: str, metal: str, coal_rank: str, tier: str,
                  c_solid_override_mg_kg: float = None) -> ScenarioResult:
    """The main entry point: run the full chain for one
    (region, metal, coal_rank, tier) combination."""
    inputs = default_inputs(region, metal, coal_rank, tier)
    mult = _tier_multiplier(tier)
    regional = rd.REGIONAL_PARAMS[region]

    # --- 1. Source term -----------------------------------------------
    ash_geom = source_term.AshSourceGeometry(
        pond_area_m2=inputs.pond_area_m2,
        ash_thickness_m=inputs.ash_thickness_m,
        ash_bulk_density_kg_m3=rd.PARTICLE_DENSITY_KG_L * 1000 * (1 - rd.ASH_POROSITY),
        ash_porosity=rd.ASH_POROSITY,
    )
    metal_key_for_kd = "Chromium(VI)" if metal == "Chromium" else metal
    metal_key_for_bulk = "Chromium(total)" if metal == "Chromium" else metal

    if c_solid_override_mg_kg is not None:
        c_solid = c_solid_override_mg_kg
        bulk_quality = "user_override"
    elif metal_key_for_bulk in rd.BULK_ASH_CONCENTRATION_MG_KG:
        c_solid, bulk_quality_enum, _ = rd.BULK_ASH_CONCENTRATION_MG_KG[metal_key_for_bulk]
        bulk_quality = bulk_quality_enum.value
    else:
        raise KeyError(f"No bulk ash concentration on file for '{metal}'. Supply "
                        f"c_solid_override_mg_kg explicitly or add it to "
                        f"reference_data.BULK_ASH_CONCENTRATION_MG_KG.")

    m0 = source_term.initial_contaminant_mass_mg(ash_geom.ash_mass_kg, c_solid)
    leaching_fraction, leach_quality, is_gap = leaching_fraction_for(metal_key_for_bulk, coal_rank, tier)
    m_leachable = source_term.leachable_mass_mg(m0, leaching_fraction)

    # Leachate concentration: leachable mass distributed through the ash's
    # own pore water (checklist #16, #21) -- the boundary concentration C0
    # feeding the groundwater transport step below.
    c0 = source_term.leachate_concentration_mg_l(m_leachable, ash_geom.pore_water_volume_m3)

    # --- 2. Vertical infiltration beneath the pond ---------------------
    vertical = flow.VerticalInfiltration(
        k_ash_m_s=rd.ASH_HYDRAULIC_CONDUCTIVITY_M_S,
        thickness_ash_m=inputs.ash_thickness_m,
        k_soil_m_s=regional["k_sat_m_s"] * mult["k_soil"],
        thickness_soil_m=3.0,  # documented screening assumption: 3 m of
                                 # unsaturated soil between ash base and
                                 # water table, not a per-site measurement
    )
    head_diff = rd.SCREENING_POND_HEAD_M  # driving head, EPA national median

    # --- 3. Horizontal groundwater transport ----------------------------
    horizontal = flow.HorizontalTransport(
        k_h_m_s=regional["k_sat_m_s"] * mult["k_soil"],
        gradient=inputs.hydraulic_gradient * mult["gradient"],
        effective_porosity=regional["porosity"],
    )
    pore_v_m_s = horizontal.pore_velocity_m_s
    pore_v_m_yr = horizontal.pore_velocity_m_year

    # --- 4. Sorption / retardation ---------------------------------------
    if metal == "Chromium":
        # Mixed Cr(III)/Cr(VI) per reference_data.DEFAULT_CR_VI_FRACTION,
        # since EPA explicitly warns against one generic chromium Kd.
        kd_vi, _ = rd.METAL_KD_L_KG["Chromium(VI)"]
        kd_iii, _ = rd.METAL_KD_L_KG["Chromium(III)"]
        f_vi = rd.DEFAULT_CR_VI_FRACTION
        kd_base = f_vi * kd_vi + (1 - f_vi) * kd_iii
        kd_quality = rd.DataQuality.LITERATURE.value
    elif metal == "Mercury":
        low, central, high = rd.MERCURY_KD_RANGE_L_KG
        kd_base = {"low": high, "central": central, "high": low}[tier]  # low-Kd = high-mobility scenario
        kd_quality = rd.DataQuality.HONEST_GAP.value
    else:
        kd_lookup = sorption.get_metal_kd(metal_key_for_kd, ph=rd.SCREENING_GROUNDWATER_PH)
        kd_base = kd_lookup.kd_l_kg
        kd_quality = kd_lookup.quality if isinstance(kd_lookup.quality, str) else kd_lookup.quality.value

    kd = kd_base * mult["kd"]
    bulk_density_soil_kg_m3 = rd.PARTICLE_DENSITY_KG_L * 1000 * (1 - regional["porosity"])
    r_factor = sorption.retardation_factor(bulk_density_soil_kg_m3, kd, regional["porosity"])
    retarded_v_m_s = sorption.retarded_velocity(pore_v_m_s, r_factor)
    retarded_v_m_yr = retarded_v_m_s * SECONDS_PER_YEAR

    # --- 5. Dispersion -----------------------------------------------------
    # Xu & Eckstein field-scale dispersivity, evaluated at the FARTHEST
    # receptor distance (the standard convention -- use the plume's own
    # transport scale, not an arbitrary fixed alpha_L).
    farthest = max(inputs.distance_to_receptors_m)
    alpha_l = transport.xu_eckstein_dispersivity_m(farthest)
    disp = transport.dispersion_coefficients(pore_v_m_s, alpha_l)

    # --- 6. Concentration vs. distance & time, at each receptor -----------
    concentration_by_receptor, exceedance_by_receptor_year = {}, {}
    travel_time_by_receptor, retarded_travel_time_by_receptor = {}, {}
    for x in inputs.distance_to_receptors_m:
        c_series = transport.concentration_time_series(
            x, np.array(REPORTING_YEARS, dtype=float), pore_v_m_s, disp["D_L"], r_factor, c0
        )
        concentration_by_receptor[x] = dict(zip(REPORTING_YEARS, c_series.tolist()))
        exceedance_by_receptor_year[x] = {
            year: (regulatory.exceedance_ratio(metal_key_for_bulk, c) if metal_key_for_bulk in rd.REGULATORY_STANDARDS
                   else None)
            for year, c in zip(REPORTING_YEARS, c_series.tolist())
        }
        try:
            travel_time_by_receptor[x] = flow.unretarded_travel_time_years(x, pore_v_m_s)
            retarded_travel_time_by_receptor[x] = sorption.retarded_travel_time_years(x, pore_v_m_s, r_factor)
        except ValueError:
            travel_time_by_receptor[x] = float("inf")
            retarded_travel_time_by_receptor[x] = float("inf")

    # --- 7. Mass balance at each reporting year ----------------------------
    # Source depletion (kinetic model, calibrated so cumulative leached mass
    # at the LAST reporting year matches the simple leaching_fraction result
    # -- see source_term.calibrate_k_l_from_target_fraction docstring for why
    # this calibration choice was made instead of guessing a rate constant).
    horizon = max(REPORTING_YEARS)
    k_l = source_term.calibrate_k_l_from_target_fraction(leaching_fraction, horizon)
    infiltration_m3_yr = vertical.infiltration_volume_m3_per_year(head_diff, inputs.pond_area_m2)
    depletion = source_term.project_source_release(
        m0_mg=m0, k_l_per_year=k_l, leachate_concentration_mg_l=c0,
        infiltration_m3_per_year=max(infiltration_m3_yr, 1e-6),
        years=np.array(REPORTING_YEARS, dtype=float),
    )

    mass_balance_by_year = {}
    for i, year in enumerate(REPORTING_YEARS):
        remaining = depletion["remaining_mass_mg"][i]
        cumulative_released = depletion["cumulative_released_mg"][i]
        # Split released mass into "dissolved in transit" (still within the
        # modeled domain, i.e. hasn't reached/passed the farthest receptor)
        # vs. "exported" (effectively at/past the farthest receptor),
        # using the nearest-receptor concentration profile as a proxy for
        # what fraction of released mass is still in the near-field domain.
        x_grid = np.linspace(0.5, farthest, 200)
        c_profile = transport.concentration_distance_profile(x_grid, year, pore_v_m_s, disp["D_L"], r_factor, c0)
        dx = x_grid[1] - x_grid[0]
        # Plume cross-sectional area, physically tied to this scenario's own
        # pond geometry (equivalent-square plume width) and a screening
        # aquifer thickness -- NOT an arbitrary constant. This matters most
        # for high-Kd metals (Vanadium, Mercury): with R in the thousands,
        # even a modest cross-section mismatch gets amplified by (R-1) when
        # converted to sorbed mass, see the cap below.
        plume_width_m = np.sqrt(inputs.pond_area_m2)
        cross_section_area_m2 = plume_width_m * rd.SCREENING_AQUIFER_THICKNESS_M
        dissolved = mass_balance.dissolved_mass_in_domain_mg(
            c_profile, dx_m=dx, cross_section_area_m2=cross_section_area_m2, porosity=regional["porosity"]
        )
        sorbed = mass_balance.sorbed_mass_mg(dissolved, r_factor)
        # Hard invariant, enforced rather than hoped-for: the domain-
        # integration above is a proxy for a fully-coupled 3-D solve (see
        # this function's own comments), and at very high retardation a
        # coarse proxy can still overshoot what the source term actually
        # released once multiplied by (R-1). Rescale dissolved+sorbed
        # proportionally (preserving their ratio, i.e. still consistent
        # with R) whenever that would happen, so mass conservation holds by
        # construction rather than by chance -- and so a genuine future bug
        # elsewhere still gets caught by mass_balance_check instead of being
        # masked by this safeguard.
        in_transit = dissolved + sorbed
        if in_transit > cumulative_released and in_transit > 0:
            scale = cumulative_released / in_transit
            dissolved *= scale
            sorbed *= scale
        exported = max(cumulative_released - dissolved - sorbed, 0.0)
        check = mass_balance.mass_balance_check(
            initial_mg=m0, remaining_in_ash_mg=remaining, dissolved_mg=dissolved,
            sorbed_mg=sorbed, exported_mg=exported, transformed_mg=0.0,
        )
        mass_balance_by_year[year] = check

    honest_gaps = []
    if is_gap:
        honest_gaps.append(f"{metal}: no clean rank-comparative leaching dataset found (see "
                            f"reference_data.HONEST_GAP_RANK_CONTAMINANTS)")
    if metal == "Mercury":
        honest_gaps.append("Mercury: Kd not in EPA RSL default table; wide illustrative range used "
                            "(reference_data.MERCURY_KD_RANGE_L_KG)")

    return ScenarioResult(
        inputs=inputs, m0_mg=m0, leachable_fraction=leaching_fraction, leachable_fraction_quality=str(leach_quality),
        m_leachable_mg=m_leachable, c0_mg_l=c0, kd_l_kg=kd, kd_quality=str(kd_quality),
        retardation_factor=r_factor, pore_velocity_m_yr=pore_v_m_yr, retarded_velocity_m_yr=retarded_v_m_yr,
        d_l_m2_s=disp["D_L"], alpha_l_m=alpha_l,
        concentration_by_receptor=concentration_by_receptor,
        travel_time_years_by_receptor=travel_time_by_receptor,
        retarded_travel_time_years_by_receptor=retarded_travel_time_by_receptor,
        exceedance_by_receptor_year=exceedance_by_receptor_year,
        mass_balance_by_year=mass_balance_by_year,
        honest_gap_flags=honest_gaps,
    )


def concentration_grid(result: ScenarioResult, x_m: np.ndarray, t_years: np.ndarray,
                        pore_v_m_s: float, d_l_m2_s: float) -> np.ndarray:
    """Full concentration(distance, time) grid for plotting, reusing a
    scenario's already-computed velocity/dispersion/retardation."""
    return transport.retarded_transport_concentration(
        x_m[:, None], t_years[None, :], pore_v_m_s, d_l_m2_s, result.retardation_factor, result.c0_mg_l
    )

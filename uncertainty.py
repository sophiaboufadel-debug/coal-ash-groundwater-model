"""
uncertainty.py
==============
Checklist sections 41-43: sensitivity analysis and Monte Carlo. Both reuse
the exact same physics functions as scenarios.py (source_term, flow,
sorption, transport) through one shared, lightweight evaluator --
`_evaluate` -- rather than re-implementing the chain, so a fix to the
physics in one place is automatically reflected everywhere.

Monte Carlo is run against a single (region, metal, coal_rank, distance,
year) target rather than the full regional table, since exploring how
UNCERTAIN a single prediction is, is a different question from exploring the
prediction itself (which is what scenarios.py's low/central/high tiers are
for) -- running full Monte Carlo across every region x metal x distance x
year combination in this project would be thousands of times more
computation for output nobody would read in full; run it for the
combination(s) you actually care about via `monte_carlo_run(...)`.
"""
from dataclasses import dataclass
import numpy as np

import reference_data as rd
import source_term
import flow
import sorption
import transport
from units import SECONDS_PER_YEAR, acres_to_m2


@dataclass
class CoreParams:
    """The uncertain parameters checklist #42-43 ask this project to vary,
    collected in one place so Monte Carlo and sensitivity analysis sample/
    perturb the exact same set."""
    pond_area_m2: float
    ash_thickness_m: float
    ash_porosity: float
    c_solid_mg_kg: float
    leaching_fraction: float
    k_soil_m_s: float
    porosity: float
    hydraulic_gradient: float
    kd_l_kg: float
    dispersivity_m: float


def _evaluate(p: CoreParams, distance_m: float, year: float) -> float:
    """One concentration prediction (mg/L) at (distance_m, year) for a given
    set of core parameters -- the shared inner loop for both Monte Carlo and
    sensitivity analysis."""
    ash_geom = source_term.AshSourceGeometry(
        pond_area_m2=p.pond_area_m2, ash_thickness_m=p.ash_thickness_m,
        ash_bulk_density_kg_m3=rd.PARTICLE_DENSITY_KG_L * 1000 * (1 - p.ash_porosity),
        ash_porosity=p.ash_porosity,
    )
    m0 = source_term.initial_contaminant_mass_mg(ash_geom.ash_mass_kg, p.c_solid_mg_kg)
    # Sampled/perturbed leaching_fraction can exceed the physically valid
    # [0,1] range for metals whose base fraction is already large (e.g.
    # Boron, base ~35%) once a Monte Carlo multiplier is applied on top --
    # clip here rather than let source_term raise mid-simulation, since a
    # clipped-at-1.0 "all of it leached" is still a physically meaningful
    # (if extreme) sample to keep in the distribution.
    leaching_fraction = float(np.clip(p.leaching_fraction, 0.0, 1.0))
    m_leachable = source_term.leachable_mass_mg(m0, leaching_fraction)
    c0 = source_term.leachate_concentration_mg_l(m_leachable, ash_geom.pore_water_volume_m3)

    horizontal = flow.HorizontalTransport(k_h_m_s=p.k_soil_m_s, gradient=p.hydraulic_gradient,
                                           effective_porosity=p.porosity)
    v = horizontal.pore_velocity_m_s

    bulk_density_kg_m3 = rd.PARTICLE_DENSITY_KG_L * 1000 * (1 - p.porosity)
    r_factor = sorption.retardation_factor(bulk_density_kg_m3, p.kd_l_kg, p.porosity)
    disp = transport.dispersion_coefficients(v, p.dispersivity_m)

    return float(transport.concentration_time_series(distance_m, np.array([year]), v, disp["D_L"],
                                                       r_factor, c0)[0])


def default_core_params(region: str, metal: str, tier: str = "central") -> CoreParams:
    """Builds a CoreParams from the same regional/metal defaults
    scenarios.py uses, at the given scenario tier, as a sane starting point
    to perturb for sensitivity analysis or center a Monte Carlo distribution
    on."""
    import scenarios as sc
    regional = rd.REGIONAL_PARAMS[region]
    mult = sc._tier_multiplier(tier)
    metal_key = "Chromium(total)" if metal == "Chromium" else metal
    c_solid, _, _ = rd.BULK_ASH_CONCENTRATION_MG_KG[metal_key]
    leaching_fraction, _, _ = sc.leaching_fraction_for(metal_key, "Bituminous", tier)
    if metal == "Mercury":
        kd = rd.MERCURY_KD_RANGE_L_KG[1]
    else:
        kd_lookup = sorption.get_metal_kd(metal, ph=rd.SCREENING_GROUNDWATER_PH)
        kd = kd_lookup.kd_l_kg
    return CoreParams(
        pond_area_m2=acres_to_m2(rd.SCREENING_POND_AREA_ACRES),
        ash_thickness_m=rd.SCREENING_ASH_THICKNESS_M,
        ash_porosity=rd.ASH_POROSITY,
        c_solid_mg_kg=c_solid,
        leaching_fraction=leaching_fraction,
        k_soil_m_s=regional["k_sat_m_s"] * mult["k_soil"],
        porosity=regional["porosity"],
        hydraulic_gradient=0.01 * mult["gradient"],
        kd_l_kg=kd * mult["kd"],
        dispersivity_m=transport.xu_eckstein_dispersivity_m(500.0),
    )


# ---------------------------------------------------------------------------
# Sensitivity analysis (checklist #41): one-at-a-time perturbation
# ---------------------------------------------------------------------------
SENSITIVITY_PARAMETERS = [
    "k_soil_m_s", "kd_l_kg", "porosity", "hydraulic_gradient",
    "leaching_fraction", "dispersivity_m", "ash_thickness_m", "c_solid_mg_kg",
]


def sensitivity_analysis(base: CoreParams, distance_m: float, year: float,
                          perturbation: float = 0.20) -> list:
    """
    Change one parameter at a time by +/-`perturbation` (default +/-20%),
    holding everything else at `base`, and measure the effect on predicted
    concentration -- checklist #41. Returns a list of dicts sorted by
    impact (largest range first), ready for a tornado chart.
    """
    baseline = _evaluate(base, distance_m, year)
    results = []
    for param in SENSITIVITY_PARAMETERS:
        base_val = getattr(base, param)
        low_params = CoreParams(**{**base.__dict__, param: base_val * (1 - perturbation)})
        high_params = CoreParams(**{**base.__dict__, param: base_val * (1 + perturbation)})
        c_low = _evaluate(low_params, distance_m, year)
        c_high = _evaluate(high_params, distance_m, year)
        results.append(dict(
            parameter=param, baseline_value=base_val,
            concentration_at_low=c_low, concentration_at_high=c_high,
            concentration_baseline=baseline,
            swing=abs(c_high - c_low),
            # signed direction: does increasing the parameter increase or
            # decrease predicted concentration? (useful for the tornado
            # chart's bar direction, not just its magnitude)
            increases_with_param=c_high >= c_low,
        ))
    results.sort(key=lambda r: r["swing"], reverse=True)
    return results


# ---------------------------------------------------------------------------
# Monte Carlo (checklist #43)
# ---------------------------------------------------------------------------
# Distribution shapes are chosen per-parameter based on how the literature
# typically reports uncertainty for that kind of quantity, not applied
# uniformly:
#   - hydraulic conductivity: log-uniform (K is well known to be
#     log-normally distributed in the field -- Freeze & Cherry -- and can
#     span orders of magnitude; sampling log(K) uniformly is the standard
#     screening-level treatment when only a plausible range, not a fitted
#     distribution, is available)
#   - Kd: log-uniform for the same reason (Kd also commonly spans orders of
#     magnitude across sites, as EPA's own Kd guidance stresses)
#   - porosity, leaching_fraction, gradient: triangular (bounded, with a
#     most-likely central value -- appropriate when a low/central/high
#     estimate is what's actually available, which is exactly this
#     project's own data situation, rather than a fitted statistical
#     distribution)
DEFAULT_MC_RANGES = {
    "k_soil_m_s": ("log-uniform", 0.3, 3.0),        # multiplicative low/high vs. central
    "kd_l_kg": ("log-uniform", 0.3, 3.0),
    "porosity": ("triangular", 0.85, 1.15),
    "hydraulic_gradient": ("triangular", 0.5, 2.0),
    "leaching_fraction": ("triangular", 0.3, 3.0),
    "dispersivity_m": ("triangular", 0.5, 2.0),
}


def _sample_one(rng: np.random.Generator, kind: str, central: float, low_mult: float, high_mult: float) -> float:
    if kind == "log-uniform":
        log_low, log_high = np.log(central * low_mult), np.log(central * high_mult)
        return float(np.exp(rng.uniform(log_low, log_high)))
    elif kind == "triangular":
        return float(rng.triangular(central * low_mult, central, central * high_mult))
    raise ValueError(f"Unknown distribution kind: {kind}")


def monte_carlo_run(base: CoreParams, distance_m: float, year: float,
                     n_iterations: int = 2000, threshold_mg_l: float = None,
                     ranges: dict = None, seed: int = 42) -> dict:
    """
    Checklist #43: sample the parameters in `ranges` (default
    DEFAULT_MC_RANGES) around `base`, run `n_iterations` simulations, and
    report the distribution of predicted concentration at (distance_m,
    year) -- median, 5th/95th percentile, and (if `threshold_mg_l` is
    given) the probability of exceeding it.
    """
    ranges = ranges or DEFAULT_MC_RANGES
    rng = np.random.default_rng(seed)
    samples = np.empty(n_iterations)
    param_samples = {k: np.empty(n_iterations) for k in ranges}

    for i in range(n_iterations):
        kwargs = dict(base.__dict__)
        for param, (kind, low_mult, high_mult) in ranges.items():
            central = getattr(base, param)
            sampled = _sample_one(rng, kind, central, low_mult, high_mult)
            kwargs[param] = sampled
            param_samples[param][i] = sampled
        samples[i] = _evaluate(CoreParams(**kwargs), distance_m, year)

    result = dict(
        n_iterations=n_iterations,
        median_mg_l=float(np.median(samples)),
        mean_mg_l=float(np.mean(samples)),
        p5_mg_l=float(np.percentile(samples, 5)),
        p95_mg_l=float(np.percentile(samples, 95)),
        min_mg_l=float(np.min(samples)),
        max_mg_l=float(np.max(samples)),
        samples=samples,
        param_samples=param_samples,
    )
    if threshold_mg_l is not None:
        result["threshold_mg_l"] = threshold_mg_l
        result["probability_exceeds_threshold"] = float(np.mean(samples > threshold_mg_l))
    return result

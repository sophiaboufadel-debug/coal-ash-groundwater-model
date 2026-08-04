"""
visualize.py
============
Checklist section 40. Every function here takes already-computed results
(from scenarios.py / uncertainty.py) and a matplotlib Axes or output path --
no physics happens in this file, on purpose, so a plotting bug can never
silently change a number in a results table.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

import reference_data as rd
import scenarios as sc
import transport

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "font.size": 10, "axes.titlesize": 12, "axes.titleweight": "bold",
    "axes.grid": True, "grid.alpha": 0.25, "axes.spines.top": False, "axes.spines.right": False,
})
REGION_COLORS = {"Northeast": "#4C72B0", "South": "#DD8452", "Midwest": "#55A868", "West": "#C44E52"}


def plot_concentration_vs_distance(result: sc.ScenarioResult, years_to_show: list, out_path: str, title: str = None):
    """checklist #40: concentration vs. distance, one line per year."""
    x = np.linspace(0.5, max(result.inputs.distance_to_receptors_m), 200)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    cmap = plt.get_cmap("viridis")
    pore_v_m_s = result.pore_velocity_m_yr / rd.SECONDS_PER_YEAR
    for i, year in enumerate(years_to_show):
        c = transport.concentration_distance_profile(
            x, year, pore_v_m_s, result.d_l_m2_s, result.retardation_factor, result.c0_mg_l,
        )
        ax.plot(x, c, label=f"{year:g} yr", color=cmap(i / max(len(years_to_show) - 1, 1)), lw=2)
    ax.set_xlabel("Distance from ash pond (m)")
    ax.set_ylabel("Predicted concentration (mg/L)")
    ax.set_title(title or f"{result.inputs.metal} — concentration vs. distance ({result.inputs.region}, "
                           f"{result.inputs.tier} scenario)")
    ax.legend(title="Time", frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_concentration_vs_time(result: sc.ScenarioResult, out_path: str, title: str = None):
    """checklist #40: concentration vs. time, one line per receptor distance."""
    fig, ax = plt.subplots(figsize=(7, 4.5))
    years = sorted(next(iter(result.concentration_by_receptor.values())).keys())
    for x_m, series in result.concentration_by_receptor.items():
        c = [series[y] for y in years]
        ax.plot(years, c, marker="o", ms=4, lw=2,
                label=f"{x_m:g} m ({sc.RECEPTOR_LABELS.get(x_m, '')})")
    std = rd.REGULATORY_STANDARDS.get(result.inputs.metal.replace("Chromium", "Chromium(total)"))
    if std:
        ax.axhline(std.value_mg_l, color="crimson", ls="--", lw=1.5,
                    label=f"{std.kind}: {std.value_mg_l:g} mg/L")
    ax.set_xlabel("Years since t=0")
    ax.set_ylabel("Predicted concentration (mg/L)")
    ax.set_title(title or f"{result.inputs.metal} — concentration vs. time ({result.inputs.region}, "
                           f"{result.inputs.tier} scenario)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_mass_vs_time(result: sc.ScenarioResult, out_path: str, title: str = None):
    """checklist #40: contaminant mass vs. time (remaining in ash vs.
    cumulative released), from the mass-balance components already computed
    for each reporting year."""
    years = sorted(result.mass_balance_by_year.keys())
    remaining = [result.mass_balance_by_year[y].components["remaining_in_ash_mg"] for y in years]
    dissolved = [result.mass_balance_by_year[y].components["dissolved_mg"] for y in years]
    sorbed = [result.mass_balance_by_year[y].components["sorbed_mg"] for y in years]
    exported = [result.mass_balance_by_year[y].components["exported_mg"] for y in years]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.stackplot(years, remaining, dissolved, sorbed, exported,
                 labels=["remaining in ash", "dissolved (in transit)", "sorbed to soil", "exported past domain"],
                 colors=["#8c8c8c", "#4C72B0", "#DD8452", "#55A868"], alpha=0.85)
    ax.set_xlabel("Years since t=0")
    ax.set_ylabel("Contaminant mass (mg)")
    ax.set_title(title or f"{result.inputs.metal} — mass balance over time ({result.inputs.region})")
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_concentration_by_metal(results_by_metal: dict, year: float, distance_m: float, out_path: str,
                                 title: str = None):
    """checklist #40: concentration by metal, one bar per metal at a fixed
    (distance, year), with the regulatory standard overlaid where one
    exists."""
    metals = list(results_by_metal.keys())
    concentrations = [results_by_metal[m].concentration_by_receptor[distance_m][year] for m in metals]
    ratios = []
    for m in metals:
        exc = results_by_metal[m].exceedance_by_receptor_year[distance_m][year]
        ratios.append(exc.ratio if exc else np.nan)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    colors = ["crimson" if (r is not None and not np.isnan(r) and r > 1) else "#4C72B0" for r in ratios]
    ax.bar(metals, concentrations, color=colors)
    ax.set_ylabel("Predicted concentration (mg/L)")
    ax.set_yscale("symlog", linthresh=1e-4)
    ax.set_title(title or f"Concentration by metal at {distance_m:g} m, year {year:g}")
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_concentration_by_region(results_by_region: dict, year: float, distance_m: float, out_path: str,
                                  metal: str = None):
    """checklist #40: concentration by region, one bar per region."""
    regions = list(results_by_region.keys())
    concentrations = [results_by_region[r].concentration_by_receptor[distance_m][year] for r in regions]
    colors = [REGION_COLORS.get(r, "#888") for r in regions]

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.bar(regions, concentrations, color=colors)
    ax.set_ylabel("Predicted concentration (mg/L)")
    ax.set_title(f"{metal or ''} — concentration by region at {distance_m:g} m, year {year:g}".strip(" —"))
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_kd_vs_concentration(kd_values: np.ndarray, result_template: sc.ScenarioResult, distance_m: float,
                              year: float, out_path: str):
    """checklist #40: Kd vs. predicted concentration, holding everything
    else at a scenario's already-computed values (sensitivity-style sweep,
    used to visually show why sorption strength matters so much)."""
    import sorption as sorp
    concentrations = []
    for kd in kd_values:
        r_factor = sorp.retardation_factor(
            rd.PARTICLE_DENSITY_KG_L * 1000 * (1 - rd.REGIONAL_PARAMS[result_template.inputs.region]["porosity"]),
            kd, rd.REGIONAL_PARAMS[result_template.inputs.region]["porosity"],
        )
        v = result_template.pore_velocity_m_yr / (365.25 * 24 * 3600)
        c = transport.concentration_time_series(distance_m, np.array([year]), v, result_template.d_l_m2_s,
                                                  r_factor, result_template.c0_mg_l)[0]
        concentrations.append(c)

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(kd_values, concentrations, lw=2, color="#4C72B0")
    ax.set_xscale("log")
    ax.set_xlabel("Kd (L/kg)")
    ax.set_ylabel(f"Predicted concentration at {distance_m:g} m, year {year:g} (mg/L)")
    ax.set_title(f"{result_template.inputs.metal} — Kd sensitivity ({result_template.inputs.region})")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_hydraulic_conductivity_vs_travel_time(k_values_m_s: np.ndarray, gradient: float, porosity: float,
                                                distance_m: float, out_path: str):
    """checklist #40: hydraulic conductivity vs. travel time."""
    import flow
    travel_times = []
    for k in k_values_m_s:
        v = flow.pore_velocity(flow.darcy_flux(k, gradient), porosity)
        travel_times.append(flow.unretarded_travel_time_years(distance_m, v))

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(k_values_m_s, travel_times, lw=2, color="#DD8452")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Hydraulic conductivity, K (m/s)")
    ax.set_ylabel(f"Unretarded travel time to {distance_m:g} m (years)")
    ax.set_title("Hydraulic conductivity vs. travel time")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_rainfall_vs_leached_mass(infiltration_fractions: np.ndarray, result: sc.ScenarioResult,
                                   pond_area_m2: float, head_diff_m: float, horizon_years: float,
                                   out_path: str):
    """checklist #40: rainfall/infiltration vs. leached mass, via the
    concentration-limited term of the source-release model."""
    import source_term as st
    leached = []
    for frac in infiltration_fractions:
        infiltration_m3_yr = frac * pond_area_m2  # simple proxy: fraction *
        # area treated as an equivalent annual infiltration volume scaling
        k_l = st.calibrate_k_l_from_target_fraction(result.leachable_fraction, horizon_years)
        depletion = st.project_source_release(
            m0_mg=result.m0_mg, k_l_per_year=k_l, leachate_concentration_mg_l=result.c0_mg_l,
            infiltration_m3_per_year=max(infiltration_m3_yr, 1e-6),
            years=np.array([horizon_years]),
        )
        leached.append(depletion["cumulative_released_mg"][0])

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(infiltration_fractions, leached, lw=2, color="#55A868", marker="o", ms=4)
    ax.set_xlabel("Infiltration rate proxy (m water / m pond area / yr)")
    ax.set_ylabel(f"Cumulative mass leached by year {horizon_years:g} (mg)")
    ax.set_title(f"{result.inputs.metal} — rainfall/infiltration vs. leached mass")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_tornado(sensitivity_results: list, out_path: str, title: str = "Sensitivity analysis (tornado)"):
    """checklist #41 visualization: horizontal bar chart of parameter
    swings, sorted largest-impact first (already sorted by
    uncertainty.sensitivity_analysis)."""
    params = [r["parameter"] for r in sensitivity_results]
    lows = [r["concentration_at_low"] for r in sensitivity_results]
    highs = [r["concentration_at_high"] for r in sensitivity_results]
    baseline = sensitivity_results[0]["concentration_baseline"] if sensitivity_results else 0.0

    fig, ax = plt.subplots(figsize=(7, 0.5 * len(params) + 1.5))
    y_pos = np.arange(len(params))
    for i, (lo, hi) in enumerate(zip(lows, highs)):
        left, right = min(lo, hi), max(lo, hi)
        ax.barh(i, right - left, left=left, color="#4C72B0", alpha=0.8, height=0.6)
    ax.axvline(baseline, color="black", lw=1, ls="--", label="baseline")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(params)
    ax.invert_yaxis()
    ax.set_xlabel("Predicted concentration (mg/L)")
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_monte_carlo_histogram(mc_result: dict, out_path: str, metal: str = "", threshold_label: str = None):
    """checklist #43 visualization: Monte Carlo output distribution with
    median/P5/P95 marked, and the exceedance threshold if one was supplied."""
    samples = mc_result["samples"]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(samples, bins=60, color="#4C72B0", alpha=0.75, edgecolor="white", linewidth=0.3)
    ax.axvline(mc_result["median_mg_l"], color="black", lw=1.5, label=f"median = {mc_result['median_mg_l']:.4g}")
    ax.axvline(mc_result["p5_mg_l"], color="gray", lw=1, ls="--", label=f"P5 = {mc_result['p5_mg_l']:.4g}")
    ax.axvline(mc_result["p95_mg_l"], color="gray", lw=1, ls="--", label=f"P95 = {mc_result['p95_mg_l']:.4g}")
    if "threshold_mg_l" in mc_result:
        label = threshold_label or f"standard = {mc_result['threshold_mg_l']:.4g}"
        ax.axvline(mc_result["threshold_mg_l"], color="crimson", lw=1.5, ls=":", label=label)
        ax.set_title(f"{metal} Monte Carlo (n={mc_result['n_iterations']}) — "
                     f"P(exceeds standard) = {mc_result['probability_exceeds_threshold']:.1%}")
    else:
        ax.set_title(f"{metal} Monte Carlo (n={mc_result['n_iterations']})")
    ax.set_xlabel("Predicted concentration (mg/L)")
    ax.set_ylabel("Count")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

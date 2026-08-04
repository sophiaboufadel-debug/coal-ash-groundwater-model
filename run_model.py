#!/usr/bin/env python3
"""
run_model.py
============
The single entry point: runs the full coal-ash leaching / groundwater
transport screening model across every region x metal x scenario-tier
combination, at the 5/10/20-year horizons (plus context years), and writes:

  outputs/tables/full_results_long.csv       every (region, metal, tier,
                                              receptor, year) prediction --
                                              the complete raw output
  outputs/tables/summary_5_10_20yr.csv       the pivoted, at-a-glance table
                                              for the years explicitly asked
                                              for, at the nearest and
                                              farthest receptor
  outputs/tables/mass_balance_check.csv      the mass-balance %error for
                                              every run (checklist #49) --
                                              read this before trusting
                                              anything else in this folder
  outputs/tables/facilities_clean.csv        the real, cleaned 367-facility
                                              table from data_io.py
  outputs/plots/*.png                        a curated set of the
                                              checklist's required graphs

Usage:  python3 run_model.py
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

import numpy as np
import pandas as pd

import reference_data as rd
import data_io
import scenarios as sc
import uncertainty as unc
import visualize as viz
import regulatory
import mass_balance

ROOT = Path(__file__).resolve().parent
TABLES_DIR = ROOT / "outputs" / "tables"
PLOTS_DIR = ROOT / "outputs" / "plots"
TABLES_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

METALS = ["Arsenic", "Selenium", "Boron", "Chromium", "Lead", "Cadmium", "Mercury", "Vanadium", "Aluminum",
          "Molybdenum", "Cobalt"]
# Lithium and Thallium have CCR groundwater protection standards on file
# (reference_data.REGULATORY_STANDARDS) but no defensible Kd/transport
# parameter set yet -- deliberately left out of the active grid rather than
# guessing a retardation behavior. See reference_data.py Section 5 notes.
REGIONS = list(rd.REGIONAL_PARAMS.keys())
TIERS = rd.SCENARIO_NAMES
REPRESENTATIVE_RANK_BY_REGION = {  # dominant rank per region, per this
    # project's own documented USGS COALQUAL regional pattern (see
    # reference_data.COAL_RANK_PRIOR_BY_STATE) -- used only to pick ONE
    # representative rank per region for the headline regional table; the
    # rank-specific multiplier itself is applied per-metal inside scenarios.py
    "Northeast": "Bituminous", "South": "Bituminous", "Midwest": "Subbituminous", "West": "Subbituminous",
}


def run_full_grid() -> pd.DataFrame:
    """Runs every (region, metal, tier) scenario and flattens the results
    into one long-format DataFrame -- one row per (region, metal, tier,
    receptor distance, reporting year)."""
    rows = []
    for region in REGIONS:
        rank = REPRESENTATIVE_RANK_BY_REGION[region]
        for metal in METALS:
            for tier in TIERS:
                result = sc.run_scenario(region, metal, rank, tier)
                for x_m, series in result.concentration_by_receptor.items():
                    for year, c in series.items():
                        exc = result.exceedance_by_receptor_year[x_m][year]
                        mb = result.mass_balance_by_year[year]
                        rows.append(dict(
                            region=region, coal_rank=rank, metal=metal, tier=tier,
                            distance_m=x_m, receptor_label=sc.RECEPTOR_LABELS.get(x_m, ""),
                            year=year, concentration_mg_l=c,
                            regulatory_standard_mg_l=exc.standard_mg_l if exc else np.nan,
                            standard_kind=exc.standard_kind if exc else "",
                            exceedance_ratio=exc.ratio if exc else np.nan,
                            exceeds=exc.exceeds if exc else False,
                            retardation_factor=result.retardation_factor,
                            kd_l_kg=result.kd_l_kg, kd_quality=result.kd_quality,
                            pore_velocity_m_yr=result.pore_velocity_m_yr,
                            retarded_velocity_m_yr=result.retarded_velocity_m_yr,
                            c0_mg_l=result.c0_mg_l, m0_mg=result.m0_mg,
                            leachable_fraction=result.leachable_fraction,
                            unretarded_travel_time_yr=result.travel_time_years_by_receptor[x_m],
                            retarded_travel_time_yr=result.retarded_travel_time_years_by_receptor[x_m],
                            mass_balance_pct_error=mb.percent_error, mass_balance_passes=mb.passes,
                            honest_gap_flags="; ".join(result.honest_gap_flags),
                        ))
    return pd.DataFrame(rows)


def build_5_10_20_summary(full_df: pd.DataFrame) -> pd.DataFrame:
    """checklist's explicit ask: results at 5/10/20-year horizons, one row
    per (region, metal, tier), at the nearest (30 m) and farthest (1500 m)
    receptor so both the near-field and far-field story are visible without
    having to open the full long-format table."""
    keep_years = [5, 10, 20]
    keep_distances = [30.0, 1500.0]
    sub = full_df[full_df["year"].isin(keep_years) & full_df["distance_m"].isin(keep_distances)].copy()
    pivot = sub.pivot_table(
        index=["region", "coal_rank", "metal", "tier"],
        columns=["distance_m", "year"],
        values="concentration_mg_l",
    )
    pivot.columns = [f"C_{int(d)}m_yr{int(y)}_mgL" for d, y in pivot.columns]
    pivot = pivot.reset_index()

    # attach the standard + central-tier exceedance ratio at the 30 m / year
    # 20 point, the single most-referenced cell for a quick screening read
    ref = sub[(sub["distance_m"] == 30.0) & (sub["year"] == 20)][
        ["region", "metal", "tier", "regulatory_standard_mg_l", "exceedance_ratio", "exceeds"]
    ].rename(columns={"regulatory_standard_mg_l": "standard_mg_l",
                       "exceedance_ratio": "exceedance_ratio_30m_yr20",
                       "exceeds": "exceeds_standard_30m_yr20"})
    pivot = pivot.merge(ref, on=["region", "metal", "tier"], how="left")
    return pivot.sort_values(["metal", "region", "tier"]).reset_index(drop=True)


def run_curated_plots():
    """A representative, readable set of plots rather than every possible
    region x metal x tier combination (4x9x3=108 scenarios would mean
    thousands of plots nobody would look at). Picks a few cases that best
    illustrate the checklist's required graph types (#40) and the more
    interesting/illustrative findings from this run."""
    print("Generating plots...")

    # Flagship case: Boron in the South (fastest-migrating metal, most
    # data-rich region match) -- concentration vs distance/time, mass balance
    boron_south = sc.run_scenario("South", "Boron", "Bituminous", "central")
    viz.plot_concentration_vs_distance(boron_south, [5, 10, 20, 30], str(PLOTS_DIR / "boron_south_concentration_vs_distance.png"))
    viz.plot_concentration_vs_time(boron_south, str(PLOTS_DIR / "boron_south_concentration_vs_time.png"))
    viz.plot_mass_vs_time(boron_south, str(PLOTS_DIR / "boron_south_mass_balance_stackplot.png"))

    # Slow-migrating contrast case: Arsenic (heavily retarded) in the same region
    arsenic_south = sc.run_scenario("South", "Arsenic", "Bituminous", "central")
    viz.plot_concentration_vs_time(arsenic_south, str(PLOTS_DIR / "arsenic_south_concentration_vs_time.png"))

    # Concentration by metal, one region, fixed distance/year
    results_by_metal = {m: sc.run_scenario("South", m, "Bituminous", "central") for m in METALS}
    viz.plot_concentration_by_metal(results_by_metal, year=20, distance_m=30.0,
                                     out_path=str(PLOTS_DIR / "concentration_by_metal_south_30m_yr20.png"))

    # Concentration by region, fixed metal/distance/year
    results_by_region = {r: sc.run_scenario(r, "Boron", REPRESENTATIVE_RANK_BY_REGION[r], "central") for r in REGIONS}
    viz.plot_concentration_by_region(results_by_region, year=20, distance_m=30.0,
                                      out_path=str(PLOTS_DIR / "boron_concentration_by_region_30m_yr20.png"),
                                      metal="Boron")

    # Kd vs concentration, hydraulic conductivity vs travel time (checklist #40)
    viz.plot_kd_vs_concentration(np.logspace(-1, 6, 40), boron_south, distance_m=150.0, year=20,
                                  out_path=str(PLOTS_DIR / "kd_vs_concentration_boron_south.png"))
    viz.plot_hydraulic_conductivity_vs_travel_time(
        np.logspace(-8, -3, 40), gradient=0.01, porosity=rd.REGIONAL_PARAMS["South"]["porosity"],
        distance_m=150.0, out_path=str(PLOTS_DIR / "hydraulic_conductivity_vs_travel_time.png"),
    )

    # Rainfall/infiltration vs leached mass
    viz.plot_rainfall_vs_leached_mass(
        np.linspace(0.05, 2.0, 25), boron_south, pond_area_m2=boron_south.inputs.pond_area_m2,
        head_diff_m=rd.SCREENING_POND_HEAD_M, horizon_years=20,
        out_path=str(PLOTS_DIR / "rainfall_vs_leached_mass_boron_south.png"),
    )

    # Sensitivity tornado + Monte Carlo, for the Boron/South case at the
    # nearest receptor (where there's actually visible dynamics within the
    # reporting horizon -- see uncertainty.py module docstring for why the
    # choice of distance/year matters for a legible demo)
    base = unc.default_core_params("South", "Boron", "central")
    sens = unc.sensitivity_analysis(base, distance_m=30.0, year=20.0)
    viz.plot_tornado(sens, str(PLOTS_DIR / "sensitivity_tornado_boron_south.png"),
                      title="Boron (South, 30 m, year 20) — sensitivity analysis")
    mc = unc.monte_carlo_run(base, distance_m=30.0, year=20.0, n_iterations=3000,
                              threshold_mg_l=rd.REGULATORY_STANDARDS["Boron"].value_mg_l, seed=7)
    viz.plot_monte_carlo_histogram(mc, str(PLOTS_DIR / "monte_carlo_boron_south.png"), metal="Boron")

    print(f"  {len(list(PLOTS_DIR.glob('*.png')))} plot(s) written to {PLOTS_DIR}")
    return dict(boron_south=boron_south, arsenic_south=arsenic_south, sensitivity=sens, monte_carlo=mc)


def main():
    print("=" * 78)
    print("COAL ASH POND CONTAMINATION / GROUNDWATER TRANSPORT SCREENING MODEL")
    print("=" * 78)

    print("\n[1/4] Loading and cleaning facility data...")
    master = data_io.build_master_facility_table()
    contaminants = data_io.load_contaminant_reference()
    data_io.save_processed(master, contaminants, ROOT / "data" / "processed")
    master.to_csv(TABLES_DIR / "facilities_clean.csv", index=False)
    print(f"  {len(master)} facilities loaded, region breakdown:")
    print(master["region"].value_counts().to_string().replace("\n", "\n  "))

    print(f"\n[2/4] Running {len(REGIONS)} regions x {len(METALS)} metals x {len(TIERS)} tiers "
          f"= {len(REGIONS)*len(METALS)*len(TIERS)} scenarios...")
    full_df = run_full_grid()
    full_df.to_csv(TABLES_DIR / "full_results_long.csv", index=False)
    print(f"  {len(full_df)} rows -> {TABLES_DIR / 'full_results_long.csv'}")

    mb_check = full_df[["region", "metal", "tier", "year", "mass_balance_pct_error", "mass_balance_passes"]].drop_duplicates()
    mb_check.to_csv(TABLES_DIR / "mass_balance_check.csv", index=False)
    n_fail = (~mb_check["mass_balance_passes"]).sum()
    print(f"  mass-balance check: {len(mb_check)} runs, {n_fail} failed "
          f"(threshold {mass_balance.MASS_BALANCE_ERROR_THRESHOLD_PCT}%) -> {TABLES_DIR / 'mass_balance_check.csv'}")

    print("\n[3/4] Building the 5/10/20-year summary table...")
    summary = build_5_10_20_summary(full_df)
    summary.to_csv(TABLES_DIR / "summary_5_10_20yr.csv", index=False)
    print(f"  {len(summary)} rows -> {TABLES_DIR / 'summary_5_10_20yr.csv'}")

    print("\n[4/4] Generating plots...")
    run_curated_plots()

    print("\n" + "=" * 78)
    print("DONE. See outputs/tables/ and outputs/plots/, and README.md for methodology.")
    print("=" * 78)


if __name__ == "__main__":
    main()

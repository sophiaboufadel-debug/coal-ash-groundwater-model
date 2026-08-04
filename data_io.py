"""
data_io.py
==========
Loads the actual project data (data/raw/*.xlsx) and produces one clean,
analysis-ready facility table. Every field traces back to a named source
file/sheet; nothing here is fabricated. Fields the raw data doesn't cover
(hydrogeology, Kd, etc.) are deliberately NOT filled in here -- that happens
downstream in scenarios.py by joining reference_data.py, so it stays visible
which columns are real facility data and which are regional literature
defaults.
"""
from pathlib import Path
import pandas as pd
import numpy as np

from reference_data import (
    CENSUS_REGION_BY_STATE, STATE_ABBREV_TO_NAME, COAL_RANK_PRIOR_BY_STATE,
    COAL_RANK_PRIOR_NOTE,
)

DATA_RAW = Path(__file__).resolve().parent.parent / "data" / "raw"


def _state_to_region(state_raw: str) -> str:
    if pd.isna(state_raw):
        return "Unknown"
    s = str(state_raw).strip()
    full = STATE_ABBREV_TO_NAME.get(s, s)  # accepts either "GA" or "Georgia"
    return CENSUS_REGION_BY_STATE.get(full, "Unknown")


def _normalize_rank(raw: str):
    """Map an EIA coal-type string (which can list multiple ranks, e.g.
    'Lignite,Subbituminous') to a single dominant rank the leaching model
    understands, keeping the raw string for provenance."""
    if pd.isna(raw):
        return None
    raw_l = str(raw).lower()
    if "bituminous" in raw_l and "sub" not in raw_l:
        return "Bituminous"
    if "subbituminous" in raw_l:
        return "Subbituminous"
    if "lignite" in raw_l:
        return "Lignite"
    return None  # e.g. "Not currently coal-fired", "Refined coal", "Retired..."


def load_ccr_facilities() -> pd.DataFrame:
    """EPA's official list of 367 CCR-regulated facilities (compliance-website
    registry). Source: US_Coal_Ash_CCR_Facilities_1.xlsx, 'All Facilities'."""
    df = pd.read_excel(DATA_RAW / "US_Coal_Ash_CCR_Facilities_1.xlsx", sheet_name="All Facilities")
    df = df.rename(columns={
        "Legacy Surface Impoundment Flag (per EPA)": "legacy_impoundment_flag",
        "Facility CCR Compliance Website (source)": "compliance_website",
        "Plant Name": "plant_name",
        "City": "city",
        "State": "state",
    })
    df["region"] = df["state"].apply(_state_to_region)
    return df


def load_pond_inventory():
    """EPA's 2009-2014 structural-integrity census of ash impoundments.
    Source: US_Coal_Ash_Pond_Inventory.xlsx."""
    ponds = pd.read_excel(DATA_RAW / "US_Coal_Ash_Pond_Inventory.xlsx", sheet_name="Pond-Level Detail")
    ponds = ponds.rename(columns={
        "Facility": "plant_name", "State": "state",
        "Unit Name": "unit_name", "Hazard Rating": "hazard_rating",
        "EPA Condition Assessment": "condition_assessment",
    })
    counts = pd.read_excel(DATA_RAW / "US_Coal_Ash_Pond_Inventory.xlsx", sheet_name="Facilities With Ponds")
    counts = counts.rename(columns={
        "State": "state", "Facility Name": "plant_name",
        "Number of Assessed Ponds/Units": "n_ponds_assessed",
    })
    return ponds, counts


def load_census_power_plants():
    """The 53-plant table with real EIA-860 coal type, impoundment liner
    status, and county demographics.
    Source: census_data_towns_and_power_plants_1.xlsx, 'Power Plants (53)'."""
    df = pd.read_excel(DATA_RAW / "census_data_towns_and_power_plants_1.xlsx", sheet_name="Power Plants (53)")
    df = df.rename(columns={
        "State": "state", "Plant Name (as given)": "plant_name", "County": "county",
        "Coal Ash Impoundment?": "has_impoundment", "Ash Impoundment Lined?": "liner_status",
        "Ash Impoundment Status": "impoundment_status",
        "County Total Population": "county_population",
        "County Median HH Income ($)": "county_median_income",
        "County Poverty Rate %": "county_poverty_rate_pct",
        "2024 Coal Burned (tons)": "coal_burned_tons_2024",
        "Coal Type(s), 2024": "coal_type_raw",
    })
    df["coal_rank"] = df["coal_type_raw"].apply(_normalize_rank)
    df["region"] = df["state"].apply(_state_to_region)
    return df


def load_contaminant_reference() -> pd.DataFrame:
    """The project's own structured leaching/content reference table.
    Source: coal_ash_contaminant_data.xlsx."""
    df = pd.read_excel(DATA_RAW / "coal_ash_contaminant_data.xlsx", sheet_name="coal_ash_contaminant_data")
    df = df.drop(columns=[c for c in df.columns if str(c).startswith("Unnamed")], errors="ignore")
    return df


def build_master_facility_table() -> pd.DataFrame:
    """
    Joins the CCR facility registry (identity/region/legacy status) with the
    census workbook's real coal-type and liner data wherever plant names
    match, and falls back to the documented regional coal-rank prior
    everywhere else. This is the single facility table every scenario run in
    this project starts from.
    """
    facilities = load_ccr_facilities()
    census = load_census_power_plants()
    ponds, pond_counts = load_pond_inventory()

    # Fuzzy-ish join on (state, normalized plant name) -- plant names are
    # written inconsistently across EPA/EIA sources ("Plant Bowen" vs
    # "Bowen"), so match on a simplified key rather than requiring an exact
    # string match.
    def _key(s):
        if pd.isna(s):
            return ""
        s = str(s).lower()
        for junk in ["plant ", "power station", "generating station", "steam plant",
                     "station", "power plant", "electric generating plant"]:
            s = s.replace(junk, "")
        return "".join(ch for ch in s if ch.isalnum())

    facilities["_key"] = facilities["plant_name"].apply(_key)
    facilities["_state_key"] = facilities["state"].apply(
        lambda s: STATE_ABBREV_TO_NAME.get(str(s).strip(), str(s).strip()) if pd.notna(s) else ""
    )
    census["_key"] = census["plant_name"].apply(_key)
    census["_state_key"] = census["state"].apply(
        lambda s: STATE_ABBREV_TO_NAME.get(str(s).strip(), str(s).strip()) if pd.notna(s) else ""
    )

    census_slim = census[["_key", "_state_key", "coal_rank", "coal_type_raw", "liner_status",
                           "impoundment_status", "county_population", "county_median_income",
                           "county_poverty_rate_pct", "coal_burned_tons_2024"]].drop_duplicates(
        subset=["_key", "_state_key"])

    merged = facilities.merge(census_slim, on=["_key", "_state_key"], how="left")

    pond_counts_slim = pond_counts.copy()
    pond_counts_slim["_key"] = pond_counts_slim["plant_name"].apply(_key)
    pond_counts_slim["_state_key"] = pond_counts_slim["state"].apply(
        lambda s: STATE_ABBREV_TO_NAME.get(str(s).strip(), str(s).strip()) if pd.notna(s) else ""
    )
    merged = merged.merge(
        pond_counts_slim[["_key", "_state_key", "n_ponds_assessed"]].drop_duplicates(subset=["_key", "_state_key"]),
        on=["_key", "_state_key"], how="left",
    )

    # Coal-rank fallback: real EIA data where matched, else the documented
    # regional prior, clearly flagged in a companion column so downstream
    # code (and the README) can tell real data from a prior apart.
    merged["coal_rank_source"] = np.where(merged["coal_rank"].notna(), "EIA-860 (matched)", "regional_prior")
    prior_rank = merged["state"].apply(
        lambda s: COAL_RANK_PRIOR_BY_STATE.get(STATE_ABBREV_TO_NAME.get(str(s).strip(), str(s).strip()))
        if pd.notna(s) else None
    )
    merged["coal_rank"] = merged["coal_rank"].fillna(prior_rank)
    merged["coal_rank_note"] = np.where(
        merged["coal_rank_source"] == "regional_prior", COAL_RANK_PRIOR_NOTE, "Matched to census_data_towns_and_power_plants_1.xlsx (EIA-860 2024)."
    )

    # Liner status fallback: real value where known, else the national 94%-
    # unlined directional prior (explicit, not silent).
    merged["liner_status_source"] = np.where(merged["liner_status"].notna(), "matched", "national_prior_94pct_unlined")

    merged = merged.drop(columns=["_key", "_state_key"])
    return merged


def save_processed(master: pd.DataFrame, contaminants: pd.DataFrame, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    master.to_csv(out_dir / "facilities_clean.csv", index=False)
    contaminants.to_csv(out_dir / "contaminant_reference.csv", index=False)


if __name__ == "__main__":
    master = build_master_facility_table()
    contaminants = load_contaminant_reference()
    out_dir = Path(__file__).resolve().parent.parent / "data" / "processed"
    save_processed(master, contaminants, out_dir)
    print(f"Facilities: {len(master)} rows -> {out_dir/'facilities_clean.csv'}")
    print(master["region"].value_counts())
    print()
    print("Coal rank source breakdown:")
    print(master["coal_rank_source"].value_counts())
    print()
    print("Coal rank breakdown:")
    print(master["coal_rank"].value_counts(dropna=False))

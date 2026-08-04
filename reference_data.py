"""
reference_data.py
==================
Every non-computed number used by this model lives in this one file, next to a
comment naming where it came from. Nothing here is invented. Where the project's
own prior research (compiled in data/raw and in the docx source-notes) did not
support a clean number, that gap is marked HONEST_GAP explicitly rather than
papered over with a plausible-looking default -- following the same convention
the project's own Coal_Ash_Contamination_Data_Sources.docx and
coal_ash_contaminant_data.xlsx already established.

DATA QUALITY FLAGS (used throughout the model, see also data_io.DataQuality)
    MEASURED        site- or facility-specific measured value
    LITERATURE       peer-reviewed / agency literature value, not site-specific
    REGIONAL_DEFAULT  USDA/NOAA/EPA regional average used as a screening proxy
    SCREENING_PROXY   a single national value substituted for a true regional
                      average because no defensible regional breakdown exists
                      (e.g. ash thickness, pond head -- see Section 30 note below)
    HONEST_GAP        no clean literature value was found; a wide/qualitative
                      placeholder is used and MUST be flagged to the user
"""
from dataclasses import dataclass, field
from enum import Enum


class DataQuality(str, Enum):
    MEASURED = "measured"
    LITERATURE = "literature"
    REGIONAL_DEFAULT = "regional_default"
    SCREENING_PROXY = "screening_proxy"
    HONEST_GAP = "honest_gap"


# ---------------------------------------------------------------------------
# 1. UNIT-SYSTEM CONVENTION
# ---------------------------------------------------------------------------
# Internal calculation units, chosen so retardation/Kd arithmetic needs no
# hidden factor-of-1000 (see sorption.py docstring for why this pairing matters):
#   length            m
#   time              s   (results are reported in years; see units.py)
#   mass              mg  for contaminant mass, kg for bulk ash/soil mass
#   concentration     mg/L
#   hydraulic conductivity   m/s
#   bulk / particle density  kg/L  (numerically identical to g/cm^3)
#   Kd                L/kg (numerically identical to mL/g and cm^3/g)

SECONDS_PER_YEAR = 365.25 * 24 * 3600  # Julian year, used throughout for y<->s


# ---------------------------------------------------------------------------
# 2. CENSUS-REGION STATE MAPPING
# ---------------------------------------------------------------------------
# Standard U.S. Census Bureau four-region breakdown (the same grouping used
# throughout this project's own prior research in data/raw and the source
# notes -- "Northeast / South / Midwest / West").
# Source: U.S. Census Bureau, "Geographic Terms and Concepts - Census Divisions
# and Census Regions," census.gov/geographies/reference-maps/2010/geo/state-maps.html
CENSUS_REGION_BY_STATE = {
    # Northeast
    "Connecticut": "Northeast", "Maine": "Northeast", "Massachusetts": "Northeast",
    "New Hampshire": "Northeast", "Rhode Island": "Northeast", "Vermont": "Northeast",
    "New Jersey": "Northeast", "New York": "Northeast", "Pennsylvania": "Northeast",
    # Midwest
    "Illinois": "Midwest", "Indiana": "Midwest", "Michigan": "Midwest",
    "Ohio": "Midwest", "Wisconsin": "Midwest", "Iowa": "Midwest", "Kansas": "Midwest",
    "Minnesota": "Midwest", "Missouri": "Midwest", "Nebraska": "Midwest",
    "North Dakota": "Midwest", "South Dakota": "Midwest",
    # South
    "Delaware": "South", "Florida": "South", "Georgia": "South", "Maryland": "South",
    "North Carolina": "South", "South Carolina": "South", "Virginia": "South",
    "District of Columbia": "South", "West Virginia": "South", "Alabama": "South",
    "Kentucky": "South", "Mississippi": "South", "Tennessee": "South",
    "Arkansas": "South", "Louisiana": "South", "Oklahoma": "South", "Texas": "South",
    # West
    "Arizona": "West", "Colorado": "West", "Idaho": "West", "Montana": "West",
    "Nevada": "West", "New Mexico": "West", "Utah": "West", "Wyoming": "West",
    "Alaska": "West", "California": "West", "Hawaii": "West", "Oregon": "West",
    "Washington": "West",
}
# 2-letter -> full-name lookup, since the project workbooks mix both formats
STATE_ABBREV_TO_NAME = {
    "AL": "Alabama", "AK": "Alaska", "AZ": "Arizona", "AR": "Arkansas", "CA": "California",
    "CO": "Colorado", "CT": "Connecticut", "DE": "Delaware", "DC": "District of Columbia",
    "FL": "Florida", "GA": "Georgia", "HI": "Hawaii", "ID": "Idaho", "IL": "Illinois",
    "IN": "Indiana", "IA": "Iowa", "KS": "Kansas", "KY": "Kentucky", "LA": "Louisiana",
    "ME": "Maine", "MD": "Maryland", "MA": "Massachusetts", "MI": "Michigan",
    "MN": "Minnesota", "MS": "Mississippi", "MO": "Missouri", "MT": "Montana",
    "NE": "Nebraska", "NV": "Nevada", "NH": "New Hampshire", "NJ": "New Jersey",
    "NM": "New Mexico", "NY": "New York", "NC": "North Carolina", "ND": "North Dakota",
    "OH": "Ohio", "OK": "Oklahoma", "OR": "Oregon", "PA": "Pennsylvania",
    "RI": "Rhode Island", "SC": "South Carolina", "SD": "South Dakota",
    "TN": "Tennessee", "TX": "Texas", "UT": "Utah", "VT": "Vermont", "VA": "Virginia",
    "WA": "Washington", "WV": "West Virginia", "WI": "Wisconsin", "WY": "Wyoming",
    "PR": "Puerto Rico",
}

# Regional coal-rank prior, used only when a facility's own coal type is not
# directly known (the census workbook gives real EIA-860 coal type for 53
# plants; for the rest of the 367-facility CCR list this documented regional
# pattern is the fallback -- NOT a per-plant fact).
# Source: US_Coal_Ash_Pond_Inventory.xlsx, "Read Me - IMPORTANT" sheet,
# "Coal type per pond" note, citing USGS COALQUAL regional patterns.
COAL_RANK_PRIOR_BY_STATE = {
    # Appalachian / Illinois Basin -> bituminous
    "West Virginia": "Bituminous", "Virginia": "Bituminous", "Kentucky": "Bituminous",
    "Ohio": "Bituminous", "Pennsylvania": "Bituminous", "Indiana": "Bituminous",
    "Illinois": "Bituminous", "Maryland": "Bituminous", "Alabama": "Bituminous",
    "Tennessee": "Bituminous",
    # Powder River Basin-supplied -> subbituminous
    "Wyoming": "Subbituminous", "Montana": "Subbituminous", "Colorado": "Subbituminous",
    "Utah": "Subbituminous", "Arizona": "Subbituminous", "New Mexico": "Subbituminous",
    "Nevada": "Subbituminous", "Idaho": "Subbituminous", "Nebraska": "Subbituminous",
    "Kansas": "Subbituminous", "Missouri": "Subbituminous", "Iowa": "Subbituminous",
    "Minnesota": "Subbituminous", "Wisconsin": "Subbituminous", "Michigan": "Subbituminous",
    "Oklahoma": "Subbituminous", "Arkansas": "Subbituminous", "Georgia": "Subbituminous",
    "South Carolina": "Subbituminous", "North Carolina": "Subbituminous",
    "Florida": "Subbituminous", "Washington": "Subbituminous", "Oregon": "Subbituminous",
    # Gulf Coast lignite belt / North Dakota -> lignite
    "North Dakota": "Lignite", "Texas": "Lignite", "Louisiana": "Lignite",
    "Mississippi": "Lignite",
}
COAL_RANK_PRIOR_NOTE = (
    "Regional prior only (USGS COALQUAL regional pattern, per "
    "US_Coal_Ash_Pond_Inventory.xlsx Read Me). Real plant-level EIA-860 coal "
    "type from census_data_towns_and_power_plants_1.xlsx is used instead "
    "wherever available and always overrides this prior."
)


# ---------------------------------------------------------------------------
# 3. REGIONAL HYDROGEOLOGIC / CLIMATE PARAMETERS
# ---------------------------------------------------------------------------
# These are the numbers this project's own prior research thread (data/raw
# source notes) derived from USDA/NOAA/USGS/EPA and explicitly recommended
# over an earlier, cruder midpoint-of-range approach. See the long derivation
# in the project's compiled notes for the full citation trail; primary sources:
#   - USDA soil texture porosity/Ksat: Ahuja et al., ASABE 2007 (ARS)
#   - NOAA 1991-2020 Climate Normals (NCEI), station-weighted by Census region
#   - EPA default screening pH of 6.8 for Kd selection when site pH is unknown
#     (EPA Regional Screening Levels User's Guide)
REGIONAL_PARAMS = {
    "Northeast": dict(
        soil_texture="Loam",
        porosity=0.463,               # USDA texture-class mean, ARS/Ahuja 2007
        k_sat_m_s=3.67e-6,             # USDA mean Ksat 1.32 cm/hr -> m/s
        mean_annual_temp_C=8.85,       # NOAA 1991-2020 Climate Normals
        groundwater_depth_note="shallow to moderate, highly variable",
        redox_note="mixed oxic/anoxic",
        quality=DataQuality.REGIONAL_DEFAULT,
    ),
    "South": dict(
        soil_texture="Loamy sand",
        porosity=0.437,
        k_sat_m_s=1.70e-5,             # USDA mean Ksat 6.11 cm/hr -> m/s
        mean_annual_temp_C=16.94,
        groundwater_depth_note="often shallowest, esp. Gulf/Atlantic coastal plain",
        redox_note="mixed, geology varies enormously",
        quality=DataQuality.REGIONAL_DEFAULT,
    ),
    "Midwest": dict(
        soil_texture="Silty clay loam",
        porosity=0.471,
        k_sat_m_s=4.17e-7,             # USDA mean Ksat 0.15 cm/hr -> m/s
        mean_annual_temp_C=9.12,
        groundwater_depth_note="shallow to moderate; glacial/alluvial systems important",
        redox_note="reducing-prone / mixed, esp. glacial aquifers",
        quality=DataQuality.REGIONAL_DEFAULT,
    ),
    "West": dict(
        soil_texture="Clay loam",
        porosity=0.464,
        k_sat_m_s=6.39e-7,             # USDA mean Ksat 0.23 cm/hr -> m/s
        mean_annual_temp_C=10.33,
        groundwater_depth_note="generally deepest and most variable",
        redox_note="oxic-prone / mixed (sand-gravel & volcanic aquifers)",
        quality=DataQuality.REGIONAL_DEFAULT,
    ),
}

# Parameters EPA/national data support only as a SINGLE national screening
# value, not a real 4-region breakdown (explicitly flagged as such in this
# project's own prior research -- do not silently promote these to
# "regional averages").
SCREENING_GROUNDWATER_PH = 6.8      # EPA default when site pH unknown (RSL User's Guide)
SCREENING_POND_HEAD_M = 1.81        # national MEDIAN ponded-water depth, EPA
                                     # national surface-impoundment dataset (nepis P1008WGO)
POND_HEAD_PERCENTILES_M = {25: 0.993, 50: 1.81, 75: 2.95, 90: 4.24, 95: 5.32, "max": 18.2}
SCREENING_ASH_THICKNESS_M = 6.0     # modeling proxy, NOT a measured regional average
ASH_THICKNESS_SENSITIVITY_M = [3.0, 6.0, 12.0, 18.0]

# Shallow unconfined aquifer thickness beneath the pond, used to convert the
# 1-D transport solution's per-unit-area concentration into an actual mass
# for the mass-balance check (mass_balance.dissolved_mass_in_domain_mg needs
# a real cross-sectional area, not an arbitrary placeholder). 5 m is a
# common screening-level default for a shallow water-table aquifer (e.g.
# used in EPA's own BIOSCREEN default dataset); treat as a documented
# screening proxy, not a measured value, same caveat as ash thickness above.
SCREENING_AQUIFER_THICKNESS_M = 5.0

# Fly-ash material properties (EPA HELP model default for "high-density
# electric-plant coal fly ash" -- source: EPA HELP 4.0 model manual,
# epa.gov/.../help_4.0_manual_v2.pdf)
ASH_POROSITY = 0.541
ASH_HYDRAULIC_CONDUCTIVITY_M_S = 5.00e-7
PARTICLE_DENSITY_KG_L = 2.65        # mineral-soil/ash grain density approximation
                                     # (standard geotechnical value, e.g. Freeze & Cherry)

# Pond area: EPA's own 40-acre threshold IS a real regulatory number (40 CFR
# 257.102(f)(2)(ii)/257.103(f) size split for closure-extension eligibility),
# reused here as a defensible mid-range screening pond size. Real named
# examples spanning the plausible range are given for sensitivity cases
# (US_Coal_Ash_Pond_Inventory.xlsx Read Me: B.C. Cobb MI=62.8 ac,
# Plant Barry AL pond=~600 ac, a Southern Indiana basin=11 ac).
SCREENING_POND_AREA_ACRES = 40.0
POND_AREA_SENSITIVITY_ACRES = [11.0, 40.0, 62.8, 600.0]
ACRE_TO_M2 = 4046.8564224

# Dispersivity: no site-calibrated value exists in the project data, so this
# model uses the empirical field-scale relationship EPA's own online
# screening tool (CEAM) recommends, rather than a single assumed constant.
# Xu, M. and Eckstein, Y. (1995). "Use of weighted least-squares method in
# evaluation of the relationship between dispersivity and field-scale."
# Ground Water 33(6), 905-908. alpha_L = 0.83 * [log10(L)]^2.414, L in meters.
# A simpler literature rule of thumb (Gelhar 1993), alpha_L = 0.1*L, is kept
# as an alternative/sensitivity option. Transverse/vertical dispersivity
# ratios (alpha_T = 0.1*alpha_L, alpha_V = 0.01*alpha_L) follow the common
# simplification cited in ASTM/EPA BIOSCREEN-style screening guidance.
MOLECULAR_DIFFUSION_M2_S = 1e-9      # typical aqueous ionic species, e.g. Fetter
ALPHA_T_OVER_ALPHA_L = 0.1
ALPHA_V_OVER_ALPHA_L = 0.01

# Fraction of US coal ash ponds documented as unlined (directional prior for
# any facility whose own liner status is unknown).
# Source: Earthjustice/Environmental Integrity Project industry-data citation,
# reproduced in US_Coal_Ash_CCR_Facilities_1.xlsx Aggregate Stats sheet.
FRACTION_US_PONDS_UNLINED = 0.94


# ---------------------------------------------------------------------------
# 4. METAL-SPECIFIC DISTRIBUTION COEFFICIENTS (Kd, L/kg)
# ---------------------------------------------------------------------------
# EPA Regional Screening Level (RSL) default soil Kd values, applicable at
# the EPA screening pH default of 6.8 used above.
# Source: US EPA, "Regional Screening Levels (RSLs) - User's Guide,"
# epa.gov/risk/regional-screening-levels-rsls-users-guide (Kd table; the same
# table this project's prior research already pulled these figures from).
# EPA explicitly warns Kd is highly soil- and pH-dependent and site-specific
# values are preferred -- treat every number below as a screening default,
# not a site fact.
METAL_KD_L_KG = {
    "Aluminum": (1500.0, DataQuality.LITERATURE),
    "Arsenic": (29.0, DataQuality.LITERATURE),
    "Boron": (3.0, DataQuality.LITERATURE),
    "Cadmium": (75.0, DataQuality.LITERATURE),
    "Chromium(VI)": (19.0, DataQuality.LITERATURE),      # hexavalent -- more mobile
    "Chromium(III)": (1_800_000.0, DataQuality.LITERATURE),  # trivalent -- effectively immobile
    "Cobalt": (45.0, DataQuality.LITERATURE),
    "Copper": (35.0, DataQuality.LITERATURE),
    "Iron": (25.0, DataQuality.LITERATURE),
    "Lead": (900.0, DataQuality.LITERATURE),
    "Manganese": (65.0, DataQuality.LITERATURE),
    "Nickel": (65.0, DataQuality.LITERATURE),
    "Selenium": (5.0, DataQuality.LITERATURE),
    "Vanadium": (1000.0, DataQuality.LITERATURE),
    "Zinc": (62.0, DataQuality.LITERATURE),
    # Mercury is an HONEST GAP: it does not appear in the EPA RSL default Kd
    # table, and reported soil/sediment Kd values for Hg span roughly 10^2 to
    # >10^4 L/kg depending overwhelmingly on organic-carbon/sulfide content
    # (general contaminant-hydrogeology literature). One real, cited,
    # SITE-SPECIFIC value is used as the illustrative point estimate here --
    # Kd = 6,900 L/kg, sediment:water, Whatcom Waterway site, WA Dept. of
    # Ecology RI/FS (WAC 173-340-747(5)(b)(ii)) -- but this MUST be treated as
    # a wide, uncertain range in any run, not a national default.
    "Mercury": (6900.0, DataQuality.HONEST_GAP),
    # Molybdenum is NOT in the EPA RSL default Kd table above (only 15 metals
    # are). Source instead: Allison, J.D. and T.L. Allison (2005), "Partition
    # Coefficients for Metals in Surface Water, Soil, and Waste," EPA/600/R-
    # 05/074, Table 3 (soil/soil-water Mo(VI); literature survey, n=5 studies,
    # confidence level 3 of 4 on the report's own 1(best)-4(worst) scale).
    # Table 3 reports log10(Kd): median=1.1, mean=1.3, min=-0.4, max=2.7.
    # 12.6 L/kg = 10^1.1 (the median) is used as the point value here; see
    # MOLYBDENUM_KD_RANGE_L_KG for the full low/high spread. Treat this as
    # meaningfully less certain than the RSL-sourced values above -- smaller
    # literature sample, different (older, less-screening-oriented) source.
    "Molybdenum": (12.6, DataQuality.LITERATURE),
}
MERCURY_KD_RANGE_L_KG = (100.0, 6900.0, 15000.0)  # (low, central-illustrative, high)
MOLYBDENUM_KD_RANGE_L_KG = (0.4, 12.6, 501.0)  # 10^(-0.4, 1.1, 2.7) L/kg, EPA/600/R-05/074 Table 3 min/median/max

# Chromium speciation default: unless a facility/study specifies otherwise,
# CCR-relevant chromium is modeled as a user-set mixture of Cr(III)/Cr(VI)
# rather than pretending total chromium behaves as one species (EPA
# explicitly warns against a single generic Cr Kd -- see RSL User's Guide).
DEFAULT_CR_VI_FRACTION = 0.10   # documented starting assumption, override per site


# ---------------------------------------------------------------------------
# 5. REGULATORY / SCREENING COMPARISON CONCENTRATIONS (mg/L)
# ---------------------------------------------------------------------------
# Kept deliberately separate from the transport physics (per item #38 of the
# original spec) so changing a regulatory number never changes a simulated
# concentration -- only the exceedance ratio.
@dataclass
class RegulatoryStandard:
    value_mg_l: float
    kind: str            # "MCL", "action_level", "secondary_MCL", "health_reference", "state_standard"
    source: str
    note: str = ""


REGULATORY_STANDARDS = {
    "Arsenic": RegulatoryStandard(0.010, "MCL", "EPA National Primary Drinking Water Regulations",
                                   "Also the CCR-rule Appendix IV groundwater protection standard."),
    "Selenium": RegulatoryStandard(0.050, "MCL", "EPA National Primary Drinking Water Regulations"),
    "Boron": RegulatoryStandard(0.7, "state_standard", "North Carolina groundwater standard",
                                 "No federal MCL exists. WHO guideline is 2.4 mg/L (for reference only, "
                                 "not a US standard); NC's 0.7 mg/L is used here because it is the one "
                                 "documented in this project's own field-data review (2012 NC CCR-pond study)."),
    "Chromium(total)": RegulatoryStandard(0.100, "MCL", "EPA National Primary Drinking Water Regulations"),
    "Chromium(VI)": RegulatoryStandard(0.100, "MCL", "EPA National Primary Drinking Water Regulations",
                                        "EPA regulates total chromium; no separate federal Cr(VI) MCL exists "
                                        "(California has a state Cr(VI) MCL of 0.010 mg/L, for reference)."),
    "Chromium(III)": RegulatoryStandard(0.100, "MCL", "EPA National Primary Drinking Water Regulations"),
    "Lead": RegulatoryStandard(0.015, "action_level", "EPA Lead and Copper Rule",
                                "NOT a standard MCL -- a treatment-technique action level. EPA's 2024 Lead "
                                "and Copper Rule Improvements also introduced a lower 0.010 mg/L 'trigger "
                                "level'; both are exposed in code as separate fields, see below."),
    "Cadmium": RegulatoryStandard(0.005, "MCL", "EPA National Primary Drinking Water Regulations"),
    "Mercury": RegulatoryStandard(0.002, "MCL", "EPA National Primary Drinking Water Regulations",
                                   "Inorganic mercury MCL."),
    "Vanadium": RegulatoryStandard(0.021, "health_reference", "EPA UCMR Health Reference Level",
                                    "No federal MCL exists; 21 ppb EPA health-reference benchmark used for "
                                    "UCMR monitoring is used here as a screening comparison value only."),
    "Aluminum": RegulatoryStandard(0.1, "secondary_MCL", "EPA National Secondary Drinking Water Regulations",
                                    "No federal primary (health-based) MCL exists. Secondary (aesthetic, "
                                    "non-enforceable) MCL range is 0.05-0.2 mg/L; midpoint used here."),
    # The next four are the CCR-rule-specific groundwater protection standards
    # (GWPS) at 40 CFR 257.95(h), verified directly against the current eCFR
    # text. Cobalt, Lithium, and Molybdenum have no federal drinking-water
    # MCL at all; EPA's July 2018 CCR Rule amendment (83 FR 36453) adopted
    # health-based concentrations for these three specifically as the CCR
    # GWPS under 257.95(h)(2). Thallium DOES have a federal MCL, so its GWPS
    # is that MCL directly, per 257.95(h)(1).
    "Cobalt": RegulatoryStandard(0.006, "health_reference", "40 CFR 257.95(h)(2)(i)",
                                  "CCR-rule health-based groundwater protection standard (6 ug/L); no "
                                  "federal drinking-water MCL exists for cobalt."),
    "Lithium": RegulatoryStandard(0.04, "health_reference", "40 CFR 257.95(h)(2)(iii)",
                                   "CCR-rule health-based groundwater protection standard (40 ug/L); no "
                                   "federal drinking-water MCL exists for lithium. NOTE: on file for "
                                   "reference/regulatory-comparison use only -- this model does not yet "
                                   "have a defensible Kd/transport parameter set for lithium, so it is "
                                   "not in run_model.py's active METALS list."),
    "Molybdenum": RegulatoryStandard(0.1, "health_reference", "40 CFR 257.95(h)(2)(iv)",
                                      "CCR-rule health-based groundwater protection standard (100 ug/L); "
                                      "no federal drinking-water MCL exists for molybdenum."),
    "Thallium": RegulatoryStandard(0.002, "MCL", "EPA National Primary Drinking Water Regulations",
                                    "Also the CCR-rule 257.95(h)(1) groundwater protection standard, since "
                                    "an MCL exists. NOTE: on file for reference/regulatory-comparison use "
                                    "only -- this model does not yet have a defensible Kd/transport "
                                    "parameter set for thallium, so it is not in run_model.py's active "
                                    "METALS list."),
}
LEAD_2024_LCRI_TRIGGER_LEVEL_MG_L = 0.010  # documented separately, see note above

GROSS_ALPHA_MCL_PCI_L = 15.0        # EPA National Primary Drinking Water Regulations (radiological ref.)


# ---------------------------------------------------------------------------
# 6. COAL-RANK LEACHING BEHAVIOR (qualitative + quantitative, by contaminant)
# ---------------------------------------------------------------------------
# This section reproduces, in code, the rank-specific findings this project's
# own Coal_Ash_Contamination_Data_Sources.docx and coal_ash_contaminant_data.xlsx
# already compiled -- INCLUDING the "honest gap" flags for the four
# contaminants where no clean rank-comparative leaching dataset was found.
# Never silently default a HONEST_GAP contaminant's rank-multiplier to 1.0 or
# to another element's value; the flag has to survive into any output table.

# Directional rank multiplier applied to a *leaching fraction* / mobility
# term, relative to a Subbituminous/Lignite (high-CaO, alkaline, ~pH 11-12)
# baseline. Mechanism: high-CaO ash produces alkaline leachate that
# precipitates/scavenges OXYANION-forming elements (As, Se, B, Cr, V) via
# calcium precipitation, but that same alkalinity keeps simple CATIONIC
# metals (Pb, Cd, Al, Hg as complexes) LESS mobile than they would be under
# the more nearly neutral-to-slightly-acidic leachate typical of low-CaO
# Bituminous ash.
# Sources: Wang et al. 2009 (Energy & Fuels), Schwartz et al. 2018
# (Environmental Engineering Science) for the As/Se rank pattern;
# Coal_Ash_Contamination_Data_Sources.docx Section 3 for the CaO/pH mechanism
# discussion, and the project's own Clean-Air-Act/fly-ash notes (Class C vs
# Class F fly ash chemistry) for the oxyanion-vs-cation generalization.
OXYANION_FORMING_METALS = {"Arsenic", "Selenium", "Boron", "Chromium(VI)", "Vanadium", "Molybdenum"}
CATIONIC_METALS = {"Lead", "Cadmium", "Aluminum", "Mercury", "Chromium(III)",
                    "Manganese", "Iron", "Cobalt", "Copper", "Nickel", "Zinc"}

# rank -> relative leaching multiplier vs. the Subbituminous/Lignite baseline
# (1.0 = baseline). These are DIRECTIONAL, order-of-magnitude multipliers
# grounded in the qualitative literature pattern above, not a precise
# calibrated number -- treat as a documented, coarse first approximation.
RANK_LEACHING_MULTIPLIER = {
    ("Arsenic", "Bituminous"): (5.0, DataQuality.LITERATURE,
        "Wang et al. 2009 / Schwartz et al. 2018: As leaching from Subbituminous "
        "'often below detection' vs. 'significantly higher' from Bituminous. "
        "5x is an illustrative, not precisely calibrated, multiplier."),
    ("Arsenic", "Subbituminous"): (1.0, DataQuality.LITERATURE, "baseline"),
    ("Arsenic", "Lignite"): (1.0, DataQuality.LITERATURE, "assumed similar to Subbituminous (both high-CaO)"),
    ("Selenium", "Bituminous"): (3.0, DataQuality.LITERATURE,
        "Wang et al. 2009: Se 'more readily leachable than As, both [rank] types' -- "
        "smaller rank contrast assumed than arsenic."),
    ("Selenium", "Subbituminous"): (1.0, DataQuality.LITERATURE, "baseline"),
    ("Selenium", "Lignite"): (1.0, DataQuality.LITERATURE, "assumed similar to Subbituminous"),
}
# Chromium, Lead, Cadmium, Mercury: explicit HONEST GAP, per
# coal_ash_contaminant_data.xlsx rows 25/27/29/31. Do not add a rank
# multiplier for these; the model must carry quality=HONEST_GAP and use the
# uncertainty range instead of a rank-adjusted point estimate.
HONEST_GAP_RANK_CONTAMINANTS = {"Chromium(VI)", "Chromium(III)", "Chromium(total)", "Lead", "Cadmium", "Mercury"}

# Bulk contaminant content in fly ash by rank, where a real rank-comparative
# dataset exists (mg/kg). Source: EPRI Report 3002003774 (2014) via
# Coal_Ash_Contamination_Data_Sources.docx Section 2 / coal_ash_contaminant_data.xlsx.
URANIUM_BULK_CONTENT_MG_KG = {
    "Bituminous": (11.2, (3.15, 30.4), 18),
    "Subbituminous": (5.63, (2.35, 15.4), 9),
    "Lignite": (4.14, (3.40, 5.55), 4),
}
BORON_LEACHATE_RANGE_MG_L = (1.0, 14.0)  # typical CCR leachate range, not rank-differentiated

# Illustrative real-world field cases for optional validation, NOT to be
# treated as ground truth without checking the primary source (both
# documents explicitly caveat this).
FIELD_VALIDATION_CASES = [
    dict(site="San Miguel Power Plant, TX", finding=">100x MCL for 12 pollutants incl. cadmium and lithium",
         source="Earthjustice/EIP compilation, secondary source -- verify vs. EPA ECHO before use"),
    dict(site="New Castle Generating Station, PA", finding="arsenic at 372x drinking-water MCL",
         source="Earthjustice/EIP compilation, secondary source -- verify vs. EPA ECHO before use"),
]


# ---------------------------------------------------------------------------
# 7. DEFAULT SOURCE-TERM ASSUMPTIONS
# ---------------------------------------------------------------------------
# A single national leaching-fraction default is NOT scientifically supported
# (this project's own research explicitly says so -- ash chemistry/pH matters
# far more than a universal percentage). These illustrative bounds come from
# real leachability studies compiled in chat history for this project
# (sequential-extraction water-soluble fractions for Zn ~0.5%, Ni ~0.8%,
# Co ~2%) and the worked example (3%) -- they define the LOW/CENTRAL/HIGH
# scenario spread used for any metal that does NOT have a metal-specific
# entry in LEACHABLE_FRACTION_BY_METAL below (Section 8), not a claimed
# universal truth.
DEFAULT_LEACHING_FRACTION = dict(low=0.005, central=0.03, high=0.10)

# Three-scenario framework (explicitly requested in this project's own prior
# planning conversation: "low release/high attenuation; central calibrated
# case; high release/low attenuation").
SCENARIO_NAMES = ["low", "central", "high"]


# ---------------------------------------------------------------------------
# 8. BULK ASH CONCENTRATION & METAL-SPECIFIC LEACHABLE FRACTION
# ---------------------------------------------------------------------------
# PROVENANCE (read before trusting these numbers):
# This block was supplied by the user as a PDF export of a ChatGPT
# conversation (data/raw/../Metals_in_Fly_Ash.pdf), not sourced directly by
# this model the way every other number in this file was. Handling:
#   1. The visible citation chips in that PDF were truncated UI stubs
#      ("Science...+2" etc.), not usable citations on their own. This model
#      extracted the underlying PDF *link annotations* instead (the actual
#      href behind each chip survives in the file even though the visible
#      text doesn't) and confirmed they resolve to real, topically-matched,
#      peer-reviewed sources (e.g. pubmed.ncbi.nlm.nih.gov, sciencedirect.com,
#      pmc.ncbi.nlm.nih.gov, a CONICET repository PDF) -- so this is NOT
#      fabricated data, but it also was not independently re-derived from
#      those primary sources here; only spot-checked for plausibility.
#   2. One entry was caught and corrected rather than trusted as-is: the
#      source table's "Aluminum, Total" bulk concentration (10.69 mg/kg) is
#      almost certainly a partial/selective-extraction figure, not true bulk
#      content -- Al is a MAJOR fly-ash constituent (part of the
#      aluminosilicate glass matrix, typically 10-25% Al2O3 by mass in Class
#      F/C fly ash per ASTM C618 chemistry), not a ppm-level trace element.
#      Using 10.69 mg/kg as C_solid for Al would understate the true bulk
#      mass present by roughly 4 orders of magnitude. ALUMINUM_BULK_MG_KG
#      below uses an oxide-composition-derived estimate instead (see its own
#      docstring), and the source table's Al row is kept only inside
#      _RAW_USER_SUPPLIED_BULK_MG_KG for transparency, clearly unused.
#   3. Fe (29.00 mg/kg) has the same likely issue (Fe2O3 is also a
#      major/percent-level fly-ash oxide) but Iron is not one of this
#      project's target contaminants, so it is left as-is and simply not
#      relied upon.
#   4. Every value below carries quality=DataQuality.LITERATURE only in the
#      loose sense that it traces to real papers; treat it as meaningfully
#      less certain than the EPA/USGS-sourced numbers elsewhere in this file,
#      and as the single most important thing to replace with actual site/
#      LEAF data before treating any model output as more than illustrative.
_RAW_USER_SUPPLIED_BULK_MG_KG = {
    # element: (value, note)
    "Manganese": (197.20, "avg of two source images, 66.00 & 328.4 mg/kg"),
    "Magnesium": (50.34, None),
    "Copper": (32.23, "avg of two source images, 40.16 & 24.3 mg/kg"),
    "Zinc": (40.255, "avg of two source images, 36.31 & 44.2 mg/kg"),
    "Sodium": (31.35, None),
    "Iron": (29.00, "NOT used -- see Section 8 docstring point 3, major element"),
    "Calcium": (24.31, None),
    "Nickel": (32.945, "avg of two source images, 22.47 & 43.42 mg/kg"),
    "Boron": (20.34, None),
    "Chromium": (18.20, None),
    "Potassium": (17.00, None),
    "Aluminum": (10.69, "NOT used -- see Section 8 docstring point 2, major element"),
    "Selenium": (10.39, None),
    "Lead": (21.255, "avg of two source images, 10.05 & 32.46 mg/kg"),
    "Cobalt": (9.08, None),
    "Vanadium": (1.07, None),
    "Arsenic": (0.86, None),
    "Titanium": (0.80, None),
    "Cadmium": (0.35, None),
}

# Aluminum bulk content, calculated from typical Class F/C fly ash oxide
# composition rather than taken from the table above (see docstring point 2).
# Typical Al2O3 content in coal fly ash is commonly cited in the ~20-30%
# range for Class F and somewhat lower for Class C (ASTM C618 only
# constrains SiO2+Al2O3+Fe2O3 jointly -- see data/raw source notes on
# ASTM C618-25a -- not Al2O3 alone). Midpoint 25% Al2O3 is used here.
# Molar conversion Al2O3 -> elemental Al: 2*(26.98)/(2*26.98+3*16.00) = 0.529.
_AL2O3_MOLAR_FRACTION_AL = (2 * 26.98) / (2 * 26.98 + 3 * 16.00)
ALUMINUM_BULK_MG_KG = dict(
    low=0.18 * _AL2O3_MOLAR_FRACTION_AL * 1_000_000,    # 18% Al2O3 (Class C low end)
    central=0.25 * _AL2O3_MOLAR_FRACTION_AL * 1_000_000,  # 25% Al2O3 (typical midpoint)
    high=0.30 * _AL2O3_MOLAR_FRACTION_AL * 1_000_000,    # 30% Al2O3 (Class F high end)
)

# Bulk ash concentration defaults actually used by source_term.py / scenarios.py,
# in mg/kg. Prefers a real measured range from this project's OWN data
# (coal_ash_contaminant_data.xlsx / Coal_Ash_Contamination_Data_Sources.docx)
# where one exists, then the user-supplied table above, then flags an
# explicit HONEST_GAP where neither source covers a metal.
BULK_ASH_CONCENTRATION_MG_KG = {
    "Arsenic": (0.86, DataQuality.LITERATURE, "user-supplied table, spot-checked plausible vs. "
                "independent literature range (~5-22 mg/kg in another cited study) -- real fly ashes "
                "vary by orders of magnitude across plants/coal sources, both are plausible"),
    "Selenium": (10.39, DataQuality.LITERATURE, "user-supplied table"),
    "Boron": (20.34, DataQuality.LITERATURE, "user-supplied table; note this project's own "
              "Coal_Ash_Contamination_Data_Sources.docx instead documents a 1-14 mg/L LEACHATE range "
              "(a different quantity -- solid content vs. dissolved concentration)"),
    "Chromium(total)": (18.20, DataQuality.LITERATURE, "user-supplied table, undifferentiated by valence"),
    "Lead": (21.255, DataQuality.LITERATURE, "user-supplied table"),
    "Cadmium": (0.35, DataQuality.LITERATURE, "user-supplied table; matches an independently found "
                "real study's 0.2-0.37 mg/kg range closely"),
    "Vanadium": (1.07, DataQuality.LITERATURE, "user-supplied table"),
    "Aluminum": (ALUMINUM_BULK_MG_KG["central"], DataQuality.LITERATURE,
                 "oxide-composition-derived, NOT the user-supplied table value -- see Section 8 docstring"),
    "Mercury": (0.1, DataQuality.HONEST_GAP, "no bulk content figure in any source consulted for this "
                "project; 0.1 mg/kg is a rough order-of-magnitude literature placeholder only "
                "(general CCP literature commonly reports <0.2 mg/kg); replace with site/LEAF data"),
    "Cobalt": (9.08, DataQuality.LITERATURE, "user-supplied table (see _RAW_USER_SUPPLIED_BULK_MG_KG); "
               "was present in that table from the start but not previously connected to this active dict"),
    # Molybdenum: an HONEST GAP, same treatment as Mercury above. No bulk
    # content figure specific to this project's own sources (not in the
    # user-supplied table, not in coal_ash_contaminant_data.xlsx). 5 mg/kg is
    # a rough order-of-magnitude literature placeholder -- general coal fly
    # ash trace-element compilations commonly place Mo in roughly a 1-40
    # mg/kg range across individual studies/plants, with no single
    # authoritative national mean identified in the sources checked for this
    # build. Replace with real site, LEAF, or EPRI coal-ash-database data
    # before treating any Molybdenum output as more than illustrative.
    "Molybdenum": (5.0, DataQuality.HONEST_GAP, "no bulk content figure specific to this project's own "
                   "sources; rough order-of-magnitude literature placeholder only; replace with site/LEAF "
                   "or EPRI coal-ash-database data"),
}

# Metal-specific WATER-LEACHABLE FRACTION (dimensionless [0-1]), replacing the
# single generic DEFAULT_LEACHING_FRACTION above wherever a metal-specific
# entry exists. low/high span the literature range found; central is the
# geometric mean of low/high (chosen over an arithmetic mean because these
# fractions span more than an order of magnitude for several metals, e.g.
# Selenium 0.7%-27%, where an arithmetic midpoint would be dominated by the
# high end). Source: user-supplied table, provenance as described above.
LEACHABLE_FRACTION_BY_METAL = {
    #                 low       high     quality
    "Manganese": (0.0, 0.0193, DataQuality.LITERATURE),
    "Magnesium": (0.0036, 0.0104, DataQuality.LITERATURE),
    "Copper": (0.0, 0.026, DataQuality.LITERATURE),
    "Zinc": (0.006, 0.034, DataQuality.LITERATURE),
    "Sodium": (0.032, 0.226, DataQuality.LITERATURE),
    "Iron": (0.0, 0.0035, DataQuality.LITERATURE),
    "Calcium": (0.10, 0.35, DataQuality.LITERATURE),
    "Nickel": (0.0, 0.233, DataQuality.LITERATURE),
    "Boron": (0.24, 0.50, DataQuality.LITERATURE),
    "Chromium(total)": (0.0, 0.082, DataQuality.LITERATURE),
    "Potassium": (0.0042, 0.1433, DataQuality.LITERATURE),
    "Aluminum": (0.0, 0.02, DataQuality.LITERATURE),
    "Selenium": (0.007, 0.273, DataQuality.LITERATURE),
    "Lead": (0.0, 0.1796, DataQuality.LITERATURE),
    "Cobalt": (0.0, 0.24, DataQuality.LITERATURE),
    "Vanadium": (0.07, 0.07, DataQuality.LITERATURE),   # single point estimate in source, no range given
    "Arsenic": (0.0001, 0.0633, DataQuality.LITERATURE),
    "Titanium": (0.0, 0.30, DataQuality.LITERATURE),
    "Cadmium": (0.0, 0.022, DataQuality.LITERATURE),
}


def geometric_mean_leaching_fraction(metal: str) -> float:
    """central-scenario leaching fraction for `metal`: geometric mean of the
    LEACHABLE_FRACTION_BY_METAL low/high range (falls back to
    DEFAULT_LEACHING_FRACTION['central'] if the metal isn't in that table)."""
    import math
    if metal not in LEACHABLE_FRACTION_BY_METAL:
        return DEFAULT_LEACHING_FRACTION["central"]
    low, high, _ = LEACHABLE_FRACTION_BY_METAL[metal]
    if low <= 0:
        low = 1e-5  # avoid log(0); these are already near-zero lower bounds
    return math.sqrt(low * high)

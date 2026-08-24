#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
HIV rapid-ART microsimulation.

Stochastic individual-based (agent-based) dynamic transmission microsimulation
of HIV acquisition, diagnosis, and care in Japan, linked to a cost-effectiveness
analysis from the public-payer perspective. Compares immediate ART (0 month)
with one-month and two-month delays over a 40-year horizon after a 20-year
burn-in. HIV-related mortality is modelled as a CD4-stratum by care-state
standardised mortality ratio applied to the age- and sex-specific
general-population background hazard.

All run conditions are set in the CONFIG section below; no command-line
arguments are used.
"""

import os
os.environ["NUMBA_NUM_THREADS"] = "8"
os.environ["NUMBA_THREADING_LAYER"] = "workqueue"

import os
import numpy as np
try:
    from numba import njit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    raise ImportError(
        "numba is not installed. "
        "Install it with 'pip install numba' to run this model."
    )

import math
import json
import datetime
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import warnings
warnings.filterwarnings('ignore')

# =============================================================================

# =============================================================================
DEFAULT_BASE_OUTPUT_DIR = r"/path/to/20250315_169/144_scenario/core"

BASE_OUTPUT_DIR = os.environ.get("OUTPUT_BASE_DIR", os.environ.get("OUTPUT_DIR", DEFAULT_BASE_OUTPUT_DIR))
OUTPUT_DIR = BASE_OUTPUT_DIR

SAVE_CD4_INDIVIDUAL_CSV = False

# =============================================================================

# =============================================================================
RUN_LABEL = "hiv_rapid_art_model"

SEED = 713
years = 40
months_per_year = 12
simulation_months = years * months_per_year

# =============================================================================

# =============================================================================
USE_DISCOUNTING_QALY = True
DISCOUNT_RATE_QALY = 0.02

# Source: Poku 2025 systematic review (53 studies, 22 countries);
#         MAINTAIN trial EQ-5D-3L; SMART/MMP SF-6D data;
#         adjusted for Japanese EQ-5D-5L value set

UTILITY_CD4_GE500   = 0.87
UTILITY_CD4_350_499 = 0.83
UTILITY_CD4_200_349 = 0.76
UTILITY_CD4_LT200   = 0.65

# Source: Shiroiwa 2021 (depression disutility 0.177);
#         Brandt 2017 (HIV anxiety/depression prevalence 30-50%);
#         Calculation: 0.177 × 0.35 = 0.062 → 0.06

MENTAL_HEALTH_UTIL_LOSS = 0.06

# =============================================================================

# =============================================================================

#

#

#

#

# =============================================================================

FOREIGN_INFLUX_BURNIN_MSM_ANNUAL  = 25
FOREIGN_INFLUX_BURNIN_HETM_ANNUAL = 30
FOREIGN_INFLUX_BURNIN_HETF_ANNUAL = 50

FOREIGN_INFLUX_MAIN_MSM_ANNUAL    = 38
FOREIGN_INFLUX_MAIN_HETM_ANNUAL   = 16
FOREIGN_INFLUX_MAIN_HETF_ANNUAL   = 13

# =============================================================================

# =============================================================================

#

#

#

# =============================================================================
INITIAL_HIV_MSM_RATIO  = 0.75
INITIAL_HIV_HETM_RATIO = 0.125
INITIAL_HIV_HETF_RATIO = 0.125
assert abs(INITIAL_HIV_MSM_RATIO + INITIAL_HIV_HETM_RATIO 
           + INITIAL_HIV_HETF_RATIO - 1.0) < 1e-6, \
       "INITIAL_HIV ratios must sum to 1.0"

# =============================================================================

# =============================================================================

ENABLE_HETEROGENEITY_TRACKING = True

# =============================================================================

# =============================================================================

SCENARIO_0_ENABLED = True
SCENARIO_1_ENABLED = True
SCENARIO_2_ENABLED = True

GRID_SEARCH_ENABLED = False

SCENARIO_CONFIGS = {
    "0 month": {"delay_mean": 0.0, "delay_sd": 0.0, "csv_label": "scenario0"},
    "1 month": {"delay_mean": 1.0, "delay_sd": 0.5, "csv_label": "scenario1"},
    "2 months": {"delay_mean": 2.0, "delay_sd": 0.5, "csv_label": "scenario2"},
}

def get_enabled_scenarios():
    enabled = []
    if SCENARIO_0_ENABLED:
        enabled.append("0 month")
    if SCENARIO_1_ENABLED:
        enabled.append("1 month")
    if SCENARIO_2_ENABLED:
        enabled.append("2 months")
    return enabled

# =============================================================================

# =============================================================================

EVENT_RNG_NAMES = [
    "init",
    "prep",
    "risk_exit",
    "art_delay",
    "art_start",
    "dropout",
    "vl_update",
    "cd4_decline",
    "reconnect",
    "testing",
    "infection",
    "infection_attr",
    "death",
    "newcomer",
]

def create_event_rngs(base_seed: int):
    ss = np.random.SeedSequence(int(base_seed))
    child_ss = ss.spawn(len(EVENT_RNG_NAMES))
    return {name: np.random.default_rng(s) for name, s in zip(EVENT_RNG_NAMES, child_ss)}

# =============================================================================

# =============================================================================
BURN_IN_YEARS = 20
BURN_IN_MONTHS = BURN_IN_YEARS * months_per_year

msm_population_size = 1_000_000
heterosexual_population_size = 1_000_000

NEWCOMER_RATE = 0.02
TARGET_POPULATION = 2_100_000

# =============================================================================

# =============================================================================

#

#

#

#

#

# =============================================================================

# =============================================================================

# =============================================================================

# =============================================================================

#

#

# =============================================================================
initial_total_hiv = 16_500       # hiv_rapid_art_model: keep tune10 burn-in stock; widen post-burn-in scenario separation via main transmission force.
initial_on_art = 7_500

initial_prep_uptake_msm = 0.01
max_prep_uptake_msm = 0.02
prep_growth_years_msm = 10

initial_prep_uptake_het = 0.00
max_prep_uptake_het = 0.005
prep_growth_years_het = 10

prep_effectiveness = 0.95

# =============================================================================

# =============================================================================

prep_dropout_rate_msm = 0.005
prep_dropout_rate_het = 0.005

USE_PRE_ART_CONDITION = True
X = 10
Y = 350.0

# hiv_rapid_art_model: burn-in ART eligibility follows historical CD4 thresholds.
# These rules override the usual CD4/VL criterion during their active years.
BURNIN_START_YEAR = 2004
ART_CD4_WAIT_POLICY_ENABLED = True
ART_CD4_WAIT_POLICY_START_YEAR = 2004
ART_CD4_WAIT_POLICY_END_YEAR = 2013
ART_CD4_WAIT_THRESHOLD = 350.0  # Backward-compatible default; see period table below.
ART_CD4_WAIT_POLICY_PERIODS = (
    (2004, 2008, 350.0),
    (2009, 2013, 500.0),
)

def get_burnin_calendar_year(month):
    """Calendar year corresponding to a burn-in month."""
    return BURNIN_START_YEAR + (int(month) - 1) // months_per_year

def get_art_cd4_wait_threshold(month, is_burnin):
    """Return the historical CD4 ART-start threshold, or None if inactive."""
    if not ART_CD4_WAIT_POLICY_ENABLED or not is_burnin:
        return None
    burnin_year = get_burnin_calendar_year(month)
    for start_year, end_year, threshold in ART_CD4_WAIT_POLICY_PERIODS:
        if start_year <= burnin_year <= end_year:
            return threshold
    return None

def is_art_cd4_wait_policy_active(month, is_burnin):
    """Return True when a historical CD4 ART-start rule is active."""
    return get_art_cd4_wait_threshold(month, is_burnin) is not None

def meets_art_cd4_wait_threshold(cd4_values, month=None, is_burnin=True):
    """CD4 eligibility under the historical ART-start rule."""
    threshold = ART_CD4_WAIT_THRESHOLD
    if month is not None:
        period_threshold = get_art_cd4_wait_threshold(month, is_burnin)
        if period_threshold is not None:
            threshold = period_threshold
    return cd4_values <= threshold

BURNIN_TEST_CAPACITY_POLICY_ON = True
BURNIN_TEST_CAPACITY_SOURCE = "2024_kensa.xlsx administrative tests + calibrated 2011-2014 diagnosis hump"
BURNIN_TEST_CAPACITY_MULTIPLIER = 1.50
BURNIN_TEST_CAPACITY_MONTHLY_BY_YEAR = {
    # Tune7 keeps total testing near the administrative-test scale but moves the
    # burn-in diagnosis peak toward the observed 2013 HIV+AIDS-reporting peak.
    2004: 14_500,
    2005: 15_000,
    2006: 15_500,
    2007: 16_000,
    2008: 16_500,
    2009: 16_500,
    2010: 17_500,
    2011: 18_500,
    2012: 19_500,
    2013: 20_500,
    2014: 20_000,
    2015: 18_500,
    2016: 18_000,
    2017: 17_500,
    2018: 17_000,
    2019: 17_000,
    2020: 14_500,
    2021: 14_000,
    2022: 14_500,
    2023: 15_000,
}
BURNIN_ADMIN_TESTS_ANNUAL = {
    # 2024_kensa.xlsx starts in 2008. The 2004-2007 values are conservative
    # back-casts toward the observed 2008 administrative-testing volume.
    2004: 135_000,
    2005: 145_000,
    2006: 160_000,
    2007: 170_000,
    2008: 177_156,
    2009: 150_252,
    2010: 130_930,
    2011: 131_243,
    2012: 131_235,
    2013: 136_400,
    2014: 145_048,
    2015: 128_241,
    2016: 118_005,
    2017: 123_432,
    2018: 130_759,
    2019: 142_260,
    2020: 68_998,
    2021: 58_172,
    2022: 73_104,
    2023: 106_137,
}

def get_burnin_test_capacity(month):
    """Monthly regular-test capacity during burn-in, anchored to admin tests."""
    burnin_year = get_burnin_calendar_year(month)
    calibrated_capacity = BURNIN_TEST_CAPACITY_MONTHLY_BY_YEAR.get(burnin_year)
    if calibrated_capacity is not None:
        return int(calibrated_capacity)
    admin_tests = BURNIN_ADMIN_TESTS_ANNUAL.get(burnin_year)
    if admin_tests is None:
        return int(TEST_CAPACITY_X)
    return int(round(admin_tests * BURNIN_TEST_CAPACITY_MULTIPLIER / months_per_year))

BURNIN_UNDIAGNOSED_HIV_TEST_BOOST_BY_YEAR = {
    2010: 1.05,
    2011: 1.10,
    2012: 1.20,
    2013: 1.35,
    2014: 1.20,
    2015: 1.10,
    2016: 1.05,
}

def get_undiagnosed_hiv_test_boost(month, is_burnin):
    """Targeted boost for HIV-positive undiagnosed testing during burn-in."""
    if not is_burnin:
        return 1.0
    burnin_year = get_burnin_calendar_year(month)
    return float(BURNIN_UNDIAGNOSED_HIV_TEST_BOOST_BY_YEAR.get(burnin_year, 1.0))

TRANSMISSION_MULTIPLIER_POLICY_ON = True
BURNIN_TRANSMISSION_MULTIPLIER_BY_YEAR = {
    2004: 1.00,
    2005: 1.00,
    2006: 1.00,
    2007: 1.00,
    2008: 1.00,
    2009: 1.02,
    2010: 1.04,
    2011: 1.08,
    2012: 1.12,
    2013: 1.15,
    2014: 1.12,
    2015: 1.09,
    2016: 1.07,
    2017: 1.07,
    2018: 1.07,
    2019: 1.07,
    2020: 1.06,
    2021: 1.06,
    2022: 1.07,
    2023: 1.07,
}
MAIN_TRANSMISSION_MULTIPLIER = 1.18

MORTALITY_RECALIB_FOI_SCALE = 1.05

RECENT_FOI_BOOST = 1.04
RECENT_FOI_BOOST_FROM_YEAR = 2016

def get_transmission_multiplier(month, is_burnin):
    """Calibration multiplier for infection force, keeping initial stock fixed."""
    if not TRANSMISSION_MULTIPLIER_POLICY_ON:
        return 1.0
    if is_burnin:
        burnin_year = get_burnin_calendar_year(month)
        base = BURNIN_TRANSMISSION_MULTIPLIER_BY_YEAR.get(burnin_year, 1.0)
        m = base * MORTALITY_RECALIB_FOI_SCALE
        if burnin_year >= RECENT_FOI_BOOST_FROM_YEAR:
            m *= RECENT_FOI_BOOST
        return float(m)
    return float(MAIN_TRANSMISSION_MULTIPLIER * MORTALITY_RECALIB_FOI_SCALE * RECENT_FOI_BOOST)

USE_CD4_VL_CRITERIA = True

vl_decay_rate_per_day = 0.2 / 365.0
vl_decay_rate_per_month = vl_decay_rate_per_day * 30.0
LOG_DROP = 4.0

dropout_prob_pre_art = 0.01
dropout_prob_on_art = 0.001

# Thresholds for analytical outputs
VL_THRESHOLD_LOG10_5000 = math.log10(5000.0)  # log10 scale viral load threshold for 5000 copies/mL
LONG_TERM_NO_ART_THRESHOLD_MONTHS = 24.0  # definition of "long-term" (2 years) since diagnosis without ART

male_age_upper = np.array(
    [4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 79, 84, 89, 94, 99, 120],
    dtype=np.float64
)
male_mortality_rate = np.array([
    0.5, 0.1, 0.1, 0.3, 0.5, 0.5, 0.6, 0.8, 1.1, 1.8,
    2.9, 4.8, 7.7, 12.8, 20.5, 34.4, 57.4, 104.5, 184.8, 310.6, 478.9
], dtype=np.float64) / (1000.0 * 12.0)

female_age_upper = np.array(
    [4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 79, 84, 89, 94, 99, 120],
    dtype=np.float64
)
female_mortality_rate = np.array([
    0.4, 0.1, 0.1, 0.2, 0.3, 0.3, 0.3, 0.4, 0.7, 1.1,
    1.7, 2.4, 3.4, 5.3, 8.4, 15.3, 28.7, 59.7, 121.5, 235.5, 393.1
], dtype=np.float64) / (1000.0 * 12.0)

# =============================================================================

# =============================================================================

#

#

# =============================================================================
HIV_MORTALITY_MODEL = "multiplier"

CD4_MORT_STRATA_LABELS = [">=500", "350-499", "200-349", "<200"]
HIV_MORT_STATE_LABELS = ["no_art", "on_art_suppressed", "on_art_unsuppressed", "disengaged_dropout"]

_HIV_MORTALITY_MODEL_CODES = {
    "legacy": 0, "multiplier": 1, "background_plus_excess": 2, "combined_max": 3,
    "literature_total_rate": 4,
}
if HIV_MORTALITY_MODEL not in _HIV_MORTALITY_MODEL_CODES:
    raise ValueError("HIV_MORTALITY_MODEL must be one of: "
                     + ", ".join(sorted(_HIV_MORTALITY_MODEL_CODES)))
HIV_MORTALITY_MODE_CODE = _HIV_MORTALITY_MODEL_CODES[HIV_MORTALITY_MODEL]

HIV_MORTALITY_TOTAL_ANNUAL = np.array([
    # >=500, 350-499, 200-349,  <200
    [0.0040, 0.0060, 0.0150, 0.0670],  # no_art
    [0.0024, 0.0038, 0.0069, 0.0200],  # on_art_suppressed
    [0.0030, 0.0048, 0.0100, 0.0350],  # on_art_unsuppressed
    [0.0050, 0.0080, 0.0200, 0.0800],  # disengaged_dropout
], dtype=np.float64)

HIV_MORTALITY_TRUE_EXCESS_ANNUAL = None

#

#   - on ART suppressed / unsuppressed:

HIV_MORTALITY_RR_VERSION = "v2_literature_smr_msm_anchor_with_japan_context_2026_06_18"
HIV_MORTALITY_MULTIPLIERS = np.array([
    # >=500, 350-499, 200-349,  <200
    [1.0, 1.4, 5.0, 10.0],  # no_art: high-CD4 from Lodwick 2010; <350 extrapolated
    [1.0, 1.1, 1.5,  2.6],  # on_art_suppressed: ART-CC 2009 MSM A/B VL<=500, plus COHERE/Konishi support
    [2.0, 2.1, 3.5,  5.5],  # on_art_unsuppressed: ART-CC 2009 MSM A/B VL>500, <200 weighted representative
    [1.5, 2.0, 6.0, 12.0],  # disengaged_dropout: extrapolated from no_art and Japanese interruption/untreated SMR
], dtype=np.float64)

def _build_hiv_mortality_rr_meta():
    meta = {}
    common_reference = "age- and sex-specific general population background mortality hazard"
    common_measure = "standardized mortality ratio or RR used as hazard multiplier"
    common_note = (
        "v2 applies this value to the Japan age- and sex-specific background mortality hazard. "
        "SMR estimates are all-cause mortality ratios and therefore conceptually match multiplier mode."
    )

    def put(state_label, cd4_label, *, lower_ci, upper_ci, calendar_period, age_adjusted,
            adjusted_covariates, source_population, derivation, direct_or_imputed,
            citation, assumption_flag, support_note=""):
        state_i = HIV_MORT_STATE_LABELS.index(state_label)
        cd4_i = CD4_MORT_STRATA_LABELS.index(cd4_label)
        meta[(state_label, cd4_label)] = {
            "point_estimate": float(HIV_MORTALITY_MULTIPLIERS[state_i, cd4_i]),
            "measure": common_measure,
            "reference_group": common_reference,
            "lower_ci": lower_ci,
            "upper_ci": upper_ci,
            "calendar_period": calendar_period,
            "age_adjusted": age_adjusted,
            "adjusted_covariates": adjusted_covariates,
            "source_population": source_population,
            "derivation": derivation,
            "direct_or_imputed": direct_or_imputed,
            "assumption_flag": bool(assumption_flag),
            "citation": citation,
            "support_note": support_note,
            "rr_version": HIV_MORTALITY_RR_VERSION,
            "v2_common_note": common_note,
            "cd4_stratum": cd4_label,
            "hiv_state": state_label,
        }

    # no_art: high CD4 from Lodwick 2010; lower CD4 extrapolated because no direct untreated 4-stratum SMR table was found.
    put(
        "no_art", ">=500",
        lower_ci="not directly estimated after CD4 split",
        upper_ci="not directly estimated after CD4 split",
        calendar_period="1990-2004 pre-ART observations in Europe and North America",
        age_adjusted="yes, SMR standardized by age, sex, country, and calendar year; CD4 split derived from adjusted rate ratio",
        adjusted_covariates="SMR standardized by age, sex, country, and year; CD4 500-699 vs 350-499 adjusted rate ratio used for split",
        source_population="ART-naive patients, MSM, CD4 >350 cells/uL",
        derivation="Lodwick MSM CD4>350 SMR 1.30 and 500-699 vs 350-499 adjusted rate ratio 0.77 imply roughly 1.0 for >=500; rounded to 1.0 to avoid implying excess risk above background when CD4 is high.",
        direct_or_imputed="derived split from direct high-CD4 SMR",
        citation="Lodwick et al. Lancet 2010, PMID 20638118",
        assumption_flag=True,
        support_note="Konishi 2025 reports recent Japan male CD4>=500 SMR around 1.05, supporting a near-background high-CD4 value."
    )
    put(
        "no_art", "350-499",
        lower_ci="not directly estimated for 350-499 alone",
        upper_ci="not directly estimated for 350-499 alone",
        calendar_period="1990-2004 pre-ART observations in Europe and North America",
        age_adjusted="yes, SMR standardized by age, sex, country, and calendar year",
        adjusted_covariates="SMR standardized by age, sex, country, and year; high-CD4 CD4 category effect estimated with Poisson regression",
        source_population="ART-naive patients, MSM, CD4 >350 cells/uL",
        derivation="Anchored to Lodwick MSM CD4>350 SMR 1.30, then assigned a higher value than >=500 using the reported 500-699 vs 350-499 adjusted rate ratio.",
        direct_or_imputed="derived split from direct high-CD4 SMR",
        citation="Lodwick et al. Lancet 2010, PMID 20638118",
        assumption_flag=True,
        support_note="Original all-CD4>350 MSM SMR was 1.30 with 95% CI 1.06-1.58."
    )
    put(
        "no_art", "200-349",
        lower_ci="not available",
        upper_ci="not available",
        calendar_period="not directly observed in the selected untreated high-CD4 paper",
        age_adjusted="not directly estimated for this cell",
        adjusted_covariates="extrapolated, no cell-specific adjustment",
        source_population="untreated or pre-ART PLHIV, CD4 200-349",
        derivation="Extrapolated above ART-CC unsuppressed MSM A/B 200-349 SMR 3.50 and below the <200 no_art value, preserving monotonicity and natural-history severity.",
        direct_or_imputed="extrapolated",
        citation="Lodwick et al. Lancet 2010; ART-CC Int J Epidemiol 2009; modelling assumption",
        assumption_flag=True,
        support_note="No direct age/sex-standardized untreated CD4 200-349 MSM SMR suitable for this 4x4 table was identified."
    )
    put(
        "no_art", "<200",
        lower_ci="not available",
        upper_ci="not available",
        calendar_period="not directly observed in the selected untreated high-CD4 paper",
        age_adjusted="not directly estimated for this cell",
        adjusted_covariates="extrapolated, no cell-specific adjustment",
        source_population="untreated or pre-ART PLHIV, CD4 <200",
        derivation="Extrapolated to represent advanced untreated disease; set above ART-CC unsuppressed weighted <200 representative and below extreme CD4<50 strata.",
        direct_or_imputed="extrapolated",
        citation="ART-CC Int J Epidemiol 2009; Lodwick et al. Lancet 2010; modelling assumption",
        assumption_flag=True,
        support_note="Konishi 2025 reports recent Japan male CD4<200 SMR around 2.59 in all care states, but untreated/disengaged advanced disease is expected to be higher."
    )

    # on_art_suppressed: ART-CC MSM, CDC A/B, VL <=500, plus modern Japan validation.
    put(
        "on_art_suppressed", ">=500",
        lower_ci="not directly estimated after >=350 split",
        upper_ci="not directly estimated after >=350 split",
        calendar_period="ART-CC cohorts assembled in 2007; follow-up from 6 months after ART start",
        age_adjusted="yes, expected deaths matched by age, sex, and country",
        adjusted_covariates="indirect standardization by age, sex, and country; stratified by risk group, CDC stage, 6-month CD4, and 6-month VL",
        source_population="MSM, no AIDS at ART start, on ART, 6-month VL <=500 copies/mL, CD4 >=350 in source table",
        derivation="ART-CC reports SMR 1.05 for CD4 >=350 and VL<=500; split to 1.0 for >=500 based on lower risk at higher CD4 and modern Japan VL<50/CD4>=500 SMR close to 1.",
        direct_or_imputed="derived split from direct ART-CC >=350 cell",
        citation="ART-CC Int J Epidemiol 2009, PMID 19820106; Konishi et al. BMC Infect Dis 2025",
        assumption_flag=True,
        support_note="Konishi 2025 reports 2021-2024 VL<50 SMR around 1.03 and CD4>=500 SMR around 1.05 in men with HIV in Japan."
    )
    put(
        "on_art_suppressed", "350-499",
        lower_ci="0.82 for ART-CC >=350 source cell",
        upper_ci="1.35 for ART-CC >=350 source cell",
        calendar_period="ART-CC cohorts assembled in 2007; follow-up from 6 months after ART start",
        age_adjusted="yes, expected deaths matched by age, sex, and country",
        adjusted_covariates="indirect standardization by age, sex, and country; stratified by risk group, CDC stage, 6-month CD4, and 6-month VL",
        source_population="MSM, no AIDS at ART start, on ART, 6-month VL <=500 copies/mL, CD4 >=350 in source table",
        derivation="ART-CC reports SMR 1.05 for CD4 >=350 and VL<=500; split to 1.1 for 350-499 to preserve a small CD4 gradient.",
        direct_or_imputed="derived split from direct ART-CC >=350 cell",
        citation="ART-CC Int J Epidemiol 2009, PMID 19820106; COHERE PLoS Med 2012, PMID 22448150",
        assumption_flag=True,
        support_note="COHERE supports a residual CD4 gradient during viral suppression."
    )
    put(
        "on_art_suppressed", "200-349",
        lower_ci=1.07,
        upper_ci=2.01,
        calendar_period="ART-CC cohorts assembled in 2007; follow-up from 6 months after ART start",
        age_adjusted="yes, expected deaths matched by age, sex, and country",
        adjusted_covariates="indirect standardization by age, sex, and country; stratified by risk group, CDC stage, 6-month CD4, and 6-month VL",
        source_population="MSM, no AIDS at ART start, on ART, 6-month VL <=500 copies/mL, CD4 200-349",
        derivation="Direct ART-CC SMR 1.47 rounded to 1.5.",
        direct_or_imputed="direct, rounded",
        citation="ART-CC Int J Epidemiol 2009, PMID 19820106; COHERE PLoS Med 2012, PMID 22448150",
        assumption_flag=False,
        support_note="Konishi 2025 reports recent Japan male CD4 200-499 SMR around 0.82, suggesting this older ART-CC anchor may be conservative for modern suppressed care."
    )
    put(
        "on_art_suppressed", "<200",
        lower_ci="not directly estimated for pooled <200 model stratum",
        upper_ci="not directly estimated for pooled <200 model stratum",
        calendar_period="ART-CC cohorts assembled in 2007; follow-up from 6 months after ART start; Japan validation 2021-2024",
        age_adjusted="yes for ART-CC source cells and Konishi SMR",
        adjusted_covariates="ART-CC indirect standardization by age, sex, and country; Konishi SMR based on Japanese male general population",
        source_population="MSM on ART with VL<=500 and CD4 <200 source cells; Japan male in-care CD4<200 validation",
        derivation="Representative <200 value using ART-CC VL<=500 source SMRs for 100-199, 50-99, and 0-49, with most model person-time expected near 100-199; set to 2.6 and supported by recent Japan CD4<200 SMR 2.59.",
        direct_or_imputed="weighted representative from finer CD4 cells",
        citation="ART-CC Int J Epidemiol 2009, PMID 19820106; Konishi et al. BMC Infect Dis 2025; COHERE PLoS Med 2012",
        assumption_flag=True,
        support_note="ART-CC VL<=500 MSM A/B SMRs were 2.27 for 100-199, 3.23 for 50-99, and 3.41 for 0-49."
    )

    # on_art_unsuppressed: ART-CC MSM, CDC A/B, VL >500.
    put(
        "on_art_unsuppressed", ">=500",
        lower_ci="not directly estimated after >=350 split",
        upper_ci="not directly estimated after >=350 split",
        calendar_period="ART-CC cohorts assembled in 2007; follow-up from 6 months after ART start",
        age_adjusted="yes, expected deaths matched by age, sex, and country",
        adjusted_covariates="indirect standardization by age, sex, and country; stratified by risk group, CDC stage, 6-month CD4, and 6-month VL",
        source_population="MSM, no AIDS at ART start, on ART, 6-month VL >500 copies/mL, CD4 >=350 in source table",
        derivation="ART-CC reports SMR 2.06 for CD4 >=350 and VL>500; split to 2.0 for >=500.",
        direct_or_imputed="derived split from direct ART-CC >=350 cell",
        citation="ART-CC Int J Epidemiol 2009, PMID 19820106",
        assumption_flag=True,
        support_note="Konishi 2025 reports 2021-2024 VL>=1000 SMR around 3.37 across CD4 strata, consistent with higher mortality when not suppressed."
    )
    put(
        "on_art_unsuppressed", "350-499",
        lower_ci="1.37 for ART-CC >=350 source cell",
        upper_ci="3.10 for ART-CC >=350 source cell",
        calendar_period="ART-CC cohorts assembled in 2007; follow-up from 6 months after ART start",
        age_adjusted="yes, expected deaths matched by age, sex, and country",
        adjusted_covariates="indirect standardization by age, sex, and country; stratified by risk group, CDC stage, 6-month CD4, and 6-month VL",
        source_population="MSM, no AIDS at ART start, on ART, 6-month VL >500 copies/mL, CD4 >=350 in source table",
        derivation="ART-CC reports SMR 2.06 for CD4 >=350 and VL>500; split to 2.1 for 350-499.",
        direct_or_imputed="derived split from direct ART-CC >=350 cell",
        citation="ART-CC Int J Epidemiol 2009, PMID 19820106",
        assumption_flag=True,
        support_note="Small split preserves monotonic CD4 gradient while staying close to the direct >=350 source cell."
    )
    put(
        "on_art_unsuppressed", "200-349",
        lower_ci=2.33,
        upper_ci=5.27,
        calendar_period="ART-CC cohorts assembled in 2007; follow-up from 6 months after ART start",
        age_adjusted="yes, expected deaths matched by age, sex, and country",
        adjusted_covariates="indirect standardization by age, sex, and country; stratified by risk group, CDC stage, 6-month CD4, and 6-month VL",
        source_population="MSM, no AIDS at ART start, on ART, 6-month VL >500 copies/mL, CD4 200-349",
        derivation="Direct ART-CC SMR 3.50 used without material rounding.",
        direct_or_imputed="direct",
        citation="ART-CC Int J Epidemiol 2009, PMID 19820106",
        assumption_flag=False,
        support_note="This replaces the previous Markov-derived 5.3, which likely overestimated non-suppressed high/mid-CD4 mortality."
    )
    put(
        "on_art_unsuppressed", "<200",
        lower_ci="not directly estimated for pooled <200 model stratum",
        upper_ci="not directly estimated for pooled <200 model stratum",
        calendar_period="ART-CC cohorts assembled in 2007; follow-up from 6 months after ART start; Japan validation 2021-2024",
        age_adjusted="yes for ART-CC source cells and Konishi SMR",
        adjusted_covariates="ART-CC indirect standardization by age, sex, and country; Konishi SMR based on Japanese male general population",
        source_population="MSM on ART with VL>500 and CD4 <200 source cells",
        derivation="Representative <200 value from ART-CC VL>500 source SMRs for 100-199, 50-99, and 0-49, using a moderate value because the model <200 stratum is expected to contain many 100-199 cases and few 0-49 cases.",
        direct_or_imputed="weighted representative from finer CD4 cells",
        citation="ART-CC Int J Epidemiol 2009, PMID 19820106; Konishi et al. BMC Infect Dis 2025",
        assumption_flag=True,
        support_note="ART-CC VL>500 MSM A/B SMRs were 3.52 for 100-199, 6.59 for 50-99, and 23.4 for 0-49."
    )

    # disengaged_dropout: no direct 4-stratum estimate found. Use explicit modelling assumption.
    for cd4_label, deriv, support in [
        (">=500", "Set above no_art high-CD4 value to represent disengagement, but kept below unsuppressed ART because high CD4 implies limited immunodeficiency.", "Konishi 2025 ART interruption/untreated overall SMR was 10.43 over 2007-2024, but not CD4-specific."),
        ("350-499", "Set above no_art 350-499 and near unsuppressed 350-499, preserving monotonicity.", "Direct disengaged CD4-specific SMRs were not identified."),
        ("200-349", "Set above no_art 200-349 and unsuppressed ART 200-349 to represent disengagement and delayed re-entry to care.", "Japanese interruption/untreated overall SMR supports a high-risk disengaged category but cannot directly determine this CD4 cell."),
        ("<200", "Set as upper-risk extrapolation for advanced disease with disengagement, above no_art <200 and above unsuppressed representative <200.", "This cell remains a major sensitivity-analysis target.")
    ]:
        put(
            "disengaged_dropout", cd4_label,
            lower_ci="not available",
            upper_ci="not available",
            calendar_period="no direct 4x4 disengaged CD4-specific SMR estimate identified; Japan validation 2007-2024",
            age_adjusted="not directly estimated for this cell",
            adjusted_covariates="extrapolated, no cell-specific adjustment",
            source_population="disengaged from care, ART interruption, or untreated after diagnosis",
            derivation=deriv,
            direct_or_imputed="extrapolated",
            citation="Konishi et al. BMC Infect Dis 2025; ART-CC Int J Epidemiol 2009; modelling assumption",
            assumption_flag=True,
            support_note=support
        )

    return meta

HIV_MORTALITY_RR_META = _build_hiv_mortality_rr_meta()

HIV_MORTALITY_EXCESS_RATE_MULTIPLIER = 1.0
HIV_MORTALITY_MULTIPLIER_SCALER = 1.0

AGE_MORTALITY_OFFSET = 5.0

def _validate_matrix(name, value, expected_shape=(4, 4), allow_none=False, strictly_positive=False):
    if value is None:
        if allow_none:
            return
        raise ValueError(f"{name} must not be None.")
    arr = np.asarray(value, dtype=np.float64)
    if arr.shape != expected_shape:
        raise ValueError(f"{name} must have shape {expected_shape}, but got {arr.shape}.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains NaN or infinite values.")
    if strictly_positive:
        if not np.all(arr > 0.0):
            raise ValueError(f"{name} must contain strictly positive values.")
    else:
        if not np.all(arr >= 0.0):
            raise ValueError(f"{name} must contain non-negative values.")

    diffs = np.diff(arr, axis=1)
    if np.any(diffs < -1e-12):
        raise ValueError(f"{name} should be non-decreasing as CD4 declines across columns.")

def validate_mortality_config():
    _validate_matrix("HIV_MORTALITY_TOTAL_ANNUAL", HIV_MORTALITY_TOTAL_ANNUAL)
    _validate_matrix("HIV_MORTALITY_MULTIPLIERS", HIV_MORTALITY_MULTIPLIERS, strictly_positive=True)
    _validate_matrix("HIV_MORTALITY_TRUE_EXCESS_ANNUAL", HIV_MORTALITY_TRUE_EXCESS_ANNUAL, allow_none=True)
    if not np.isfinite(float(HIV_MORTALITY_EXCESS_RATE_MULTIPLIER)) or HIV_MORTALITY_EXCESS_RATE_MULTIPLIER < 0.0:
        raise ValueError("HIV_MORTALITY_EXCESS_RATE_MULTIPLIER must be finite and non-negative.")
    if not np.isfinite(float(HIV_MORTALITY_MULTIPLIER_SCALER)) or HIV_MORTALITY_MULTIPLIER_SCALER < 0.0:
        raise ValueError("HIV_MORTALITY_MULTIPLIER_SCALER must be finite and non-negative.")
    if not np.isfinite(float(AGE_MORTALITY_OFFSET)):
        raise ValueError("AGE_MORTALITY_OFFSET must be finite.")

def _format_config_value_for_path(value):
    s = f"{float(value):g}"
    return s.replace("-", "m").replace(".", "p")

def build_run_label():
    return (
        "hiv_rapid_art_model_"
        f"{HIV_MORTALITY_MODEL}_"
        f"rr{_format_config_value_for_path(HIV_MORTALITY_MULTIPLIER_SCALER)}_"
        f"rate{_format_config_value_for_path(HIV_MORTALITY_EXCESS_RATE_MULTIPLIER)}_"
        f"foi{_format_config_value_for_path(MORTALITY_RECALIB_FOI_SCALE)}_"
        f"rfoi{_format_config_value_for_path(RECENT_FOI_BOOST)}_"
        f"ageoff{_format_config_value_for_path(AGE_MORTALITY_OFFSET)}_"
        f"seed{int(SEED)}"
    )

validate_mortality_config()

if HIV_MORTALITY_MODE_CODE in (2, 3):  # background_plus_excess / combined_max
    if HIV_MORTALITY_TRUE_EXCESS_ANNUAL is None:
        raise ValueError(
            "HIV_MORTALITY_MODEL='" + HIV_MORTALITY_MODEL + "' requires a true excess mortality rate."
            "Set HIV_MORTALITY_TRUE_EXCESS_ANNUAL to 'observed all-cause mortality rate - expected background (age-standardised)'."
            "Leaving it unset would double-count background mortality, so execution is halted."
        )
    _ACTIVE_RATE_TABLE = HIV_MORTALITY_TRUE_EXCESS_ANNUAL
else:

    _ACTIVE_RATE_TABLE = HIV_MORTALITY_TOTAL_ANNUAL

RUN_LABEL = build_run_label()
OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, RUN_LABEL)
os.makedirs(OUTPUT_DIR, exist_ok=True)

HIV_MORTALITY_SOURCE_NOTE = (
    "v2 RR/SMR base case: ART-CC Int J Epidemiol 2009 Table 3 MSM A/B CD4-VL SMRs; "
    "Lodwick Lancet 2010 ART-naive high-CD4 MSM SMR and CD4 gradient; "
    "COHERE PLoS Med 2012 suppressed-VL residual CD4 gradient; "
    "Konishi BMC Infect Dis 2025 Japan male HIV SMRs by CD4, viral load, and ART interruption/untreated status"
)

male_age_groups_lower = np.array([20, 30, 40, 50])
male_age_groups_upper = np.array([29, 39, 49, 59])
male_age_groups_prob = np.array([0.3, 0.3, 0.2, 0.2])

female_age_groups_lower = np.array([20, 30, 40, 50])
female_age_groups_upper = np.array([29, 39, 49, 59])
female_age_groups_prob = np.array([0.4, 0.3, 0.2, 0.1])

# =============================================================================

# =============================================================================

#

#

#

#   hiv_rapid_art_model Phase 2:  p_infect_het_bins = [0, 1e-4, 2e-4, 4e-4, 1e-3]    (= base × 2.0)

#   hiv_rapid_art_model Phase 3:  p_infect_het_bins = [0, 1.5e-4, 3e-4, 6e-4, 1.5e-3] (= base × 3.0)

# =============================================================================
VL_BIN_EDGES = np.array([0.0, 2.3, 4.0, 5.0, 6.0, 10.0])
p_infect_msm_bins = np.array([0, 0.001, 0.001, 0.002, 0.005])    # hiv_rapid_art_model: keep baseline per-act infectivity; apply a stronger main-only multiplier.
p_infect_het_bins = np.array([0, 1.5e-4, 3.0e-4, 6.0e-4, 1.5e-3]) # hiv_rapid_art_model: keep baseline per-act infectivity while preserving route balance.

# =============================================================================

# =============================================================================

#

#

#

#   hiv_rapid_art_model/135/136: contact_rate_msm = 11.0, contact_rate_het = 3.0
#   hiv_rapid_art_model Phase 2: contact_rate_msm = 10.0, contact_rate_het = 5.0

# =============================================================================
contact_rate_msm = 10.1   # hiv_rapid_art_model: keep tune5 contact level; incidence is adjusted by transmission multiplier.

contact_rate_het = 7.0    # hiv_rapid_art_model: keep Phase 3 heterosexual contact level; route mix is handled by the bridge layer.

# =============================================================================
# hiv_rapid_art_model: Heterosexual core/high-activity layer
# =============================================================================
# This represents Japan-relevant heterosexual high-activity/bridge
# heterogeneity. It is a latent compartment for partner turnover, concurrency,
# migration/foreign-nationality bridge features, and delayed testing access;
# it is not a directly observed sex-worker population. For this recalibration,
# heterosexual per-contact probabilities remain at the Phase 3 upper calibration.
# The bridge layer is kept within the deep_research_v1 sensitivity range, but
# expressed as a broader and less concentrated bridge group than hiv_rapid_art_model_b.
# The first ACT=9 attempt made heterosexual transmission too weak (Hetero share
# ~5-6%, burn-in incidence 638/year). The subsequent ACT=12 run improved
# reporting calibration (burn-in diagnoses 957/year) but still underfit
# infections (726/year) and heterosexual contribution (40-year Het share
# 8.5-11.2% across ART-delay scenarios). ACT is therefore raised to 14 while
# keeping SHARE=3% and FOREIGN=15%. This preserves the bridge framing while
# moving Het toward the publication target of roughly 17-20%.
HETERO_CORE_SHARE = 0.03
HETERO_CORE_ACTIVITY_MULTIPLIER = 14.0   # bridge concept: ACT=12 still underfit Het share and total incidence
HETERO_CORE_INITIAL_HIV_WEIGHT = 12.0    # keep hiv_rapid_art_model_b conservative initial enrichment; do not re-concentrate HIV in the bridge layer
HETERO_CORE_FOREIGN_SHARE = 0.15         # bridge concept: lower foreign/migrant bridge share within 12-15%

TEST_CAPACITY_POLICY_ON = True
TEST_CAPACITY_X = 14500
TEST_CAPACITY_Y = 18000
TEST_CAPACITY_Z = 5

# =============================================================================
# hiv_rapid_art_model_c publication-calibration anchors from deep_research_v1.md
# =============================================================================
# Primary baseline: 2023-2024. Sensitivity baseline: 2022-2024.
# These are reporting/cascade anchors, not all direct infection-incidence targets.
CALIBRATION_BASELINE = "Primary 2023-2024; sensitivity 2022-2024"
TARGET_HIV_REPORTS_ANNUAL = 654
TARGET_AIDS_REPORTS_ANNUAL = 292
TARGET_TOTAL_REPORTS_ANNUAL = 946
TARGET_AIDS_AT_DIAGNOSIS_SHARE = 0.308
TARGET_ROUTE_MSM_SHARE = 0.80
TARGET_ROUTE_HET_SHARE = 0.20
TARGET_ART_CURRENT = 29_680
TARGET_ART_CURRENT_NARROW_LOW = 28_500
TARGET_ART_CURRENT_NARROW_HIGH = 30_500
TARGET_VL_SUPPRESSED_ON_ART = 0.99
TARGET_UNDIAGNOSED_HIV = 3_500
TARGET_UNDIAGNOSED_HIV_LOW = 3_200
TARGET_UNDIAGNOSED_HIV_HIGH = 3_800
TARGET_PLHIV_SOFT_LOW = 33_000
TARGET_PLHIV_SOFT_HIGH = 36_000

# =============================================================================

# =============================================================================

@njit
def get_mortality_rate_by_age(adj_age, gender, upper_male, rate_male, upper_female, rate_female):
    if gender == 0:
        for i in range(len(upper_male)):
            if adj_age <= upper_male[i]:
                return rate_male[i]
        return rate_male[-1]
    else:
        for i in range(len(upper_female)):
            if adj_age <= upper_female[i]:
                return rate_female[i]
        return rate_female[-1]

@njit
def calc_mortality_risk_single(cd4, age, is_aids, gender, upper_male, rate_male, upper_female, rate_female):
    adj_age = age + AGE_MORTALITY_OFFSET
    mr = get_mortality_rate_by_age(adj_age, gender, upper_male, rate_male, upper_female, rate_female)
    if is_aids:
        mr *= 4
    if cd4 < 100.0:
        mr *= 1.2
    return mr

@njit(parallel=True)
def calc_mortality_risk_batch(cd4_arr, age_arr, is_aids_arr, gender_arr, 
                               upper_male, rate_male, upper_female, rate_female):
    n = len(cd4_arr)
    result = np.empty(n, dtype=np.float64)
    for i in prange(n):
        result[i] = calc_mortality_risk_single(
            cd4_arr[i], age_arr[i], is_aids_arr[i], gender_arr[i],
            upper_male, rate_male, upper_female, rate_female
        )
    return result

@njit(parallel=True)
def calc_mortality_risk_uninfected_batch(age_arr, gender_arr, 
                                          upper_male, rate_male, upper_female, rate_female):
    n = len(age_arr)
    result = np.empty(n, dtype=np.float64)
    for i in prange(n):
        adj_age = age_arr[i] + AGE_MORTALITY_OFFSET
        result[i] = get_mortality_rate_by_age(adj_age, gender_arr[i],
                                               upper_male, rate_male, upper_female, rate_female)
    return result

# =============================================================================

# =============================================================================
@njit
def get_cd4_mortality_stratum_index(cd4):
    if cd4 >= 500.0:
        return 0
    if cd4 >= 350.0:
        return 1
    if cd4 >= 200.0:
        return 2
    return 3

@njit
def get_hiv_mortality_state_index(on_art, vl_suppressed, dropped_out, diagnosed):
    if dropped_out:
        return 3
    if on_art:
        if vl_suppressed:
            return 1
        return 2
    return 0

@njit
def combine_monthly_hazards(background_p, annual_excess_rate):
    if annual_excess_rate < 0.0:
        annual_excess_rate = 0.0
    excess_p = 1.0 - math.exp(-annual_excess_rate / 12.0)
    p = 1.0 - (1.0 - background_p) * (1.0 - excess_p)
    if p < 0.0:
        return 0.0
    if p > 0.95:
        return 0.95
    return p

@njit
def calc_mortality_risk_single_v143(cd4, age, is_aids, gender,
                                    on_art, vl_suppressed, dropped_out, diagnosed,
                                    upper_male, rate_male, upper_female, rate_female,
                                    mode_code, excess_annual, mortality_multipliers,
                                    excess_rate_multiplier, multiplier_scaler):
    adj_age = age + AGE_MORTALITY_OFFSET
    background_p = get_mortality_rate_by_age(adj_age, gender, upper_male, rate_male,
                                             upper_female, rate_female)

    if mode_code == 0:
        mr = background_p
        if is_aids:
            mr *= 4.0
        if cd4 < 100.0:
            mr *= 1.2
        if mr < 0.0:
            return 0.0
        if mr > 0.95:
            return 0.95
        return mr

    cd4_idx = get_cd4_mortality_stratum_index(cd4)
    state_idx = get_hiv_mortality_state_index(on_art, vl_suppressed, dropped_out, diagnosed)

    rr = mortality_multipliers[state_idx, cd4_idx]
    if rr <= 0.0:
        scaled_rr = 0.0
    else:
        scaled_rr = math.exp(multiplier_scaler * math.log(rr))
    bg_p = background_p
    if bg_p >= 1.0:
        bg_p = 0.999999
    bg_hazard = -math.log1p(-bg_p)
    p_multiplier = 1.0 - math.exp(-bg_hazard * scaled_rr)
    if p_multiplier < 0.0:
        p_multiplier = 0.0
    if p_multiplier > 0.95:
        p_multiplier = 0.95

    if mode_code == 1:      # multiplier only
        return p_multiplier

    if mode_code == 4:
        total_rate = excess_annual[state_idx, cd4_idx] * excess_rate_multiplier
        if total_rate < 0.0:
            total_rate = 0.0
        p_total = 1.0 - math.exp(-total_rate / 12.0)
        if p_total > 0.95:
            p_total = 0.95
        return p_total

    annual_excess = excess_annual[state_idx, cd4_idx] * excess_rate_multiplier
    p_excess = combine_monthly_hazards(background_p, annual_excess)
    if mode_code == 2:
        return p_excess

    if p_excess >= p_multiplier:
        return p_excess
    return p_multiplier

@njit(parallel=True)
def calc_mortality_risk_batch_v143(cd4_arr, age_arr, is_aids_arr, gender_arr,
                                   on_art_arr, vl_supp_arr, dropped_out_arr, diagnosed_arr,
                                   upper_male, rate_male, upper_female, rate_female,
                                   mode_code, excess_annual, mortality_multipliers,
                                   excess_rate_multiplier, multiplier_scaler):
    n = len(cd4_arr)
    result = np.empty(n, dtype=np.float64)
    for i in prange(n):
        result[i] = calc_mortality_risk_single_v143(
            cd4_arr[i], age_arr[i], is_aids_arr[i], gender_arr[i],
            on_art_arr[i], vl_supp_arr[i], dropped_out_arr[i], diagnosed_arr[i],
            upper_male, rate_male, upper_female, rate_female,
            mode_code, excess_annual, mortality_multipliers,
            excess_rate_multiplier, multiplier_scaler
        )
    return result

def compute_cd4_state_counts(cd4_arr, on_art_arr, supp_arr, dropped_arr):
    counts = np.zeros((4, 4), dtype=np.int64)
    n = len(cd4_arr)
    if n == 0:
        return counts
    cd4 = np.asarray(cd4_arr)
    strat = np.where(cd4 >= 500.0, 0, np.where(cd4 >= 350.0, 1, np.where(cd4 >= 200.0, 2, 3)))
    on_art = np.asarray(on_art_arr).astype(bool)
    supp = np.asarray(supp_arr).astype(bool)
    dropped = np.asarray(dropped_arr).astype(bool)
    state = np.where(dropped, 3, np.where(on_art & supp, 1, np.where(on_art & ~supp, 2, 0)))
    np.add.at(counts, (state, strat), 1)
    return counts

def compute_cd4_state_weighted(cd4_arr, on_art_arr, supp_arr, dropped_arr, weights):
    out = np.zeros((4, 4), dtype=np.float64)
    n = len(cd4_arr)
    if n == 0:
        return out
    cd4 = np.asarray(cd4_arr)
    strat = np.where(cd4 >= 500.0, 0, np.where(cd4 >= 350.0, 1, np.where(cd4 >= 200.0, 2, 3)))
    on_art = np.asarray(on_art_arr).astype(bool)
    supp = np.asarray(supp_arr).astype(bool)
    dropped = np.asarray(dropped_arr).astype(bool)
    state = np.where(dropped, 3, np.where(on_art & supp, 1, np.where(on_art & ~supp, 2, 0)))
    np.add.at(out, (state, strat), np.asarray(weights, dtype=np.float64))
    return out

@njit
def simulate_viral_load_randomized(vl, time_since_art, time_to_suppression, log_drop, decay_rate,
                                    noise_val, blip_rand):
    if time_since_art <= 0:
        return vl
    t = time_since_art
    tts = time_to_suppression
    if t < tts:
        frac = t / tts
        drop_amount = log_drop * frac
        new_vl = vl - drop_amount
        if new_vl < 0:
            new_vl = 0.0
    else:
        over = t - tts
        init_suppressed = vl - log_drop
        if init_suppressed < 0:
            init_suppressed = 0.0
        new_vl = init_suppressed - (decay_rate * over)
        if new_vl < 0:
            new_vl = 0.0

    new_vl += noise_val
    if new_vl < 0:
        new_vl = 0.0
    p_blip = 0.01
    if blip_rand < p_blip:
        new_vl += 0.5
    return new_vl

@njit(parallel=True)
def simulate_viral_load_batch(vl_arr, time_since_art_arr, time_to_suppression_arr, 
                               log_drop, decay_rate, noise_arr, blip_rand_arr):
    n = len(vl_arr)
    result = np.empty(n, dtype=np.float64)
    for i in prange(n):
        result[i] = simulate_viral_load_randomized(
            vl_arr[i], time_since_art_arr[i], time_to_suppression_arr[i],
            log_drop, decay_rate, noise_arr[i], blip_rand_arr[i]
        )
    return result

@njit
def simulate_cd4_recovery_on_art(cd4_current, cd4_start, cd4_ceiling, suppressed):
    alpha_fast0 = 0.06
    alpha_min = 0.005
    gap0 = max(cd4_ceiling - cd4_start, 1.0)
    gap = max(cd4_ceiling - cd4_current, 0.0)
    if suppressed:
        alpha = alpha_fast0 * (gap / gap0)
        if alpha < alpha_min:
            alpha = alpha_min
    else:
        alpha = alpha_min
    delta = alpha * (cd4_ceiling - cd4_current)
    new_cd4 = cd4_current + (delta if delta > 0.0 else 0.0)
    if new_cd4 > cd4_ceiling:
        new_cd4 = cd4_ceiling
    return new_cd4

@njit(parallel=True)
def simulate_cd4_recovery_batch(cd4_arr, cd4_start_arr, cd4_ceiling_arr, suppressed_arr):
    n = len(cd4_arr)
    result = np.empty(n, dtype=np.float64)
    for i in prange(n):
        result[i] = simulate_cd4_recovery_on_art(
            cd4_arr[i], cd4_start_arr[i], cd4_ceiling_arr[i], suppressed_arr[i]
        )
    return result

@njit(parallel=True)
def update_cd4_tracking(cd4, prev_cd4, 
                        t_first_ge200, t_first_ge350, t_first_ge500,
                        months_at_200, months_at_350, months_at_500,
                        t_sustained_200, t_sustained_350, t_sustained_500, 
                        month):
    n = len(cd4)
    for i in prange(n):
        # First crossing
        if t_first_ge200[i] < 0 and prev_cd4[i] < 200.0 and cd4[i] >= 200.0:
            t_first_ge200[i] = month
        if t_first_ge350[i] < 0 and prev_cd4[i] < 350.0 and cd4[i] >= 350.0:
            t_first_ge350[i] = month
        if t_first_ge500[i] < 0 and prev_cd4[i] < 500.0 and cd4[i] >= 500.0:
            t_first_ge500[i] = month
        
        # Sustained counters
        if cd4[i] >= 200.0:
            months_at_200[i] += 1.0
        else:
            months_at_200[i] = 0.0
        if cd4[i] >= 350.0:
            months_at_350[i] += 1.0
        else:
            months_at_350[i] = 0.0
        if cd4[i] >= 500.0:
            months_at_500[i] += 1.0
        else:
            months_at_500[i] = 0.0
        
        # Mark sustained (6 months)
        if t_sustained_200[i] < 0 and months_at_200[i] >= 6.0:
            t_sustained_200[i] = month
        if t_sustained_350[i] < 0 and months_at_350[i] >= 6.0:
            t_sustained_350[i] = month
        if t_sustained_500[i] < 0 and months_at_500[i] >= 6.0:
            t_sustained_500[i] = month

# =============================================================================

# =============================================================================

def calc_age_exit_probability_vectorized(ages):
    probs = np.zeros_like(ages, dtype=np.float64)
    probs[ages < 50] = 0.0001 / 12.0
    probs[(ages >= 50) & (ages < 60)] = 0.001 / 12.0
    probs[(ages >= 60) & (ages < 70)] = 0.01 / 12.0
    probs[ages >= 70] = 0.03 / 12.0
    return probs

def get_age_distribution_by_year(year_idx, gender="male"):
    if gender == "male":
        base_probs = male_age_groups_prob.copy()
    else:
        base_probs = female_age_groups_prob.copy()
    
    year_factor = min(year_idx / 20, 1.0)
    increase_30s = 0.2 * year_factor
    decrease_others = increase_30s / 3
    
    adjusted_probs = base_probs.copy()
    adjusted_probs[1] += increase_30s
    adjusted_probs[0] -= decrease_others
    adjusted_probs[2] -= decrease_others
    adjusted_probs[3] -= decrease_others
    adjusted_probs = adjusted_probs / adjusted_probs.sum()
    
    return adjusted_probs

def generate_ages_from_distribution(n, age_groups_lower, age_groups_upper, probs, rng):
    bins = np.arange(len(probs))
    chosen_bins = rng.choice(bins, size=n, p=probs)
    ages = rng.uniform(age_groups_lower[chosen_bins], age_groups_upper[chosen_bins])
    return ages

def assign_viral_load_vectorized(cd4_arr, rng):
    result = np.empty_like(cd4_arr)
    mask_high = cd4_arr >= 500
    mask_mid = (cd4_arr >= 200) & (cd4_arr < 500)
    mask_low = cd4_arr < 200
    
    result[mask_high] = rng.uniform(3.0, 4.0, size=mask_high.sum())
    result[mask_mid] = rng.uniform(4.5, 5.5, size=mask_mid.sum())
    result[mask_low] = rng.uniform(5.5, 6.5, size=mask_low.sum())
    return result

def generate_time_to_suppression_skewed(n, rng):
    alpha, beta = 2.0, 5.0
    vals = rng.beta(alpha, beta, size=n)
    return 1.0 + 5.0 * vals

def gamma_delay(mean, sd, size, rng):
    alpha = (mean**2) / (sd**2)
    theta = (sd**2) / mean
    return rng.gamma(shape=alpha, scale=theta, size=size)

def assign_cd4_ceiling_from_nadir_vectorized(nadir_cd4, rng):
    n = len(nadir_cd4)
    result = np.empty(n, dtype=np.float64)
    
    conditions = [
        nadir_cd4 < 100,
        (nadir_cd4 >= 100) & (nadir_cd4 < 200),
        (nadir_cd4 >= 200) & (nadir_cd4 < 350),
        nadir_cd4 >= 350
    ]
    means = [200, 350, 500, 800]
    
    for cond, mean_val in zip(conditions, means):
        count = cond.sum()
        if count > 0:
            vals = rng.normal(loc=mean_val, scale=50, size=count)
            vals = np.clip(vals, 50, 1200)
            result[cond] = vals
    
    return result

def calc_cd4_decline_vectorized(pop, indices, month, rng):
    n = len(indices)
    if n == 0:
        return np.array([])
    
    ti = month - pop['time_of_infection'][indices]
    ti = np.maximum(ti, 0)
    
    cd4_now = pop['cd4'][indices]
    s1 = pop['stage1_end'][indices]
    s2 = pop['stage2_end'][indices]
    s3 = pop['stage3_end'][indices]
    s4 = pop['stage4_end'][indices]
    
    stage1_drop = pop['stage1_drop_ratio'][indices]
    stage2_recover = pop['stage2_recover_ratio'][indices]
    d3 = pop['stage3_decline'][indices]
    d4 = pop['stage4_decline'][indices]
    d5 = pop['stage5_decline'][indices]

    mask_s1 = ti < s1
    mask_s2 = (ti >= s1) & (ti < s2)
    mask_s3 = (ti >= s2) & (ti < s3)
    mask_s4 = (ti >= s3) & (ti < s4)
    mask_s5 = ti >= s4
    
    decline = np.zeros(n, dtype=np.float64)
    decline[mask_s1] = cd4_now[mask_s1] * stage1_drop[mask_s1]
    decline[mask_s2] = -cd4_now[mask_s2] * stage2_recover[mask_s2]
    decline[mask_s3] = d3[mask_s3]
    decline[mask_s4] = d4[mask_s4]
    decline[mask_s5] = d5[mask_s5]

    decline *= pop['cd4_decline_factor'][indices]
    decline *= pop['cd4_decline_base'][indices]

    late_stage = ti >= s3
    decline[late_stage] *= pop['cd4_decline_accelerator'][indices][late_stage]

    decline *= rng.uniform(0.9, 1.1, size=n)
    
    return decline

def generate_new_cd4_for_infection_vectorized(ages, rng):
    n = len(ages)
    out_cd4 = np.empty(n, dtype=np.float64)
    
    base_p = 0.3
    prob_200 = np.where(ages < 30, base_p * 0.8,
                np.where(ages < 50, base_p * 1.0, np.minimum(base_p * 1.5, 1.0)))
    
    is_low_cd4 = rng.random(n) < prob_200

    low_vals = rng.lognormal(mean=np.log(150), sigma=0.5, size=is_low_cd4.sum())
    low_vals = np.minimum(low_vals, 200.0)
    out_cd4[is_low_cd4] = low_vals

    high_count = (~is_low_cd4).sum()
    high_vals = rng.lognormal(mean=np.log(500), sigma=0.5, size=high_count)
    high_vals = np.minimum(high_vals, 1000.0)
    out_cd4[~is_low_cd4] = high_vals
    
    return out_cd4

def get_vl_bin_indices_vectorized(vl_arr):
    return np.clip(np.searchsorted(VL_BIN_EDGES[1:], vl_arr), 0, len(VL_BIN_EDGES) - 2)

def get_age_group_vectorized(ages):
    bins = np.array([0, 30, 40, 50, 60, 70, 200])
    labels = ["20-29", "30-39", "40-49", "50-59", "60-69", "70+"]
    indices = np.digitize(ages, bins[1:])
    return np.array(labels)[np.minimum(indices, len(labels) - 1)]

def get_cd4_category_vectorized(cd4, hiv_positive):
    result = np.empty(len(cd4), dtype='<U10')
    result[~hiv_positive] = "HIV-"
    
    hiv_pos = hiv_positive
    result[hiv_pos & (cd4 < 200)] = "<200"
    result[hiv_pos & (cd4 >= 200) & (cd4 < 350)] = "200-349"
    result[hiv_pos & (cd4 >= 350) & (cd4 < 500)] = "350-499"
    result[hiv_pos & (cd4 >= 500)] = ">=500"
    
    return result

# =============================================================================

# =============================================================================

AGE_GROUP_LABELS = ["20-29", "30-39", "40-49", "50-59", "60-69", "70+"]
CD4_CAT_LABELS = ["HIV-", "<200", "200-349", "350-499", ">=500"]

N_AGE_GROUPS = 6
N_CD4_CATS = 5
N_ART = 2
N_PREP = 2
AGGREGATOR_KEY_MULTIPLIERS = (N_AGE_GROUPS * N_CD4_CATS * N_ART * N_PREP,  # year
                               N_CD4_CATS * N_ART * N_PREP,                 # age
                               N_ART * N_PREP,                              # cd4
                               N_PREP,                                      # art
                               1)                                           # prep

def get_age_group_index_vectorized(ages):
    bins = np.array([0, 30, 40, 50, 60, 70, 200])
    indices = np.digitize(ages, bins[1:])
    return np.minimum(indices, N_AGE_GROUPS - 1).astype(np.int32)

def get_cd4_category_index_vectorized(cd4, hiv_positive):
    result = np.zeros(len(cd4), dtype=np.int32)
    
    hiv_pos = hiv_positive
    result[hiv_pos & (cd4 < 200)] = 1
    result[hiv_pos & (cd4 >= 200) & (cd4 < 350)] = 2
    result[hiv_pos & (cd4 >= 350) & (cd4 < 500)] = 3
    result[hiv_pos & (cd4 >= 500)] = 4
    
    return result

def encode_aggregator_keys(year_idx, age_indices, cd4_indices, on_art, on_prep):
    return (year_idx * AGGREGATOR_KEY_MULTIPLIERS[0] +
            age_indices * AGGREGATOR_KEY_MULTIPLIERS[1] +
            cd4_indices * AGGREGATOR_KEY_MULTIPLIERS[2] +
            on_art.astype(np.int32) * AGGREGATOR_KEY_MULTIPLIERS[3] +
            on_prep.astype(np.int32) * AGGREGATOR_KEY_MULTIPLIERS[4])

def decode_aggregator_key(key):
    year_idx = key // AGGREGATOR_KEY_MULTIPLIERS[0]
    remainder = key % AGGREGATOR_KEY_MULTIPLIERS[0]
    
    age_idx = remainder // AGGREGATOR_KEY_MULTIPLIERS[1]
    remainder = remainder % AGGREGATOR_KEY_MULTIPLIERS[1]
    
    cd4_idx = remainder // AGGREGATOR_KEY_MULTIPLIERS[2]
    remainder = remainder % AGGREGATOR_KEY_MULTIPLIERS[2]
    
    art_idx = remainder // AGGREGATOR_KEY_MULTIPLIERS[3]
    prep_idx = remainder % AGGREGATOR_KEY_MULTIPLIERS[3]
    
    return (year_idx, 
            AGE_GROUP_LABELS[age_idx], 
            CD4_CAT_LABELS[cd4_idx],
            bool(art_idx),
            bool(prep_idx))

def aggregate_counts_vectorized(keys, max_key):
    if len(keys) == 0:
        return np.zeros(max_key + 1, dtype=np.int64)
    return np.bincount(keys, minlength=max_key + 1)

# =============================================================================

# =============================================================================

def create_population_arrays(n):
    pop = {
        'population_type': np.empty(n, dtype=np.int32),
        'gender': np.empty(n, dtype=np.int32),
        'age': np.empty(n, dtype=np.float64),
        'hiv_positive': np.zeros(n, dtype=np.bool_),
        'diagnosed': np.zeros(n, dtype=np.bool_),
        'cd4': np.full(n, 500.0, dtype=np.float64),
        'vl': np.full(n, 5.0, dtype=np.float64),
        'cd4_at_art_start': -np.ones(n, dtype=np.float64),
        'cd4_nadir': np.full(n, np.inf, dtype=np.float64),
        't_first_ge200': -np.ones(n, dtype=np.float64),
        't_first_ge350': -np.ones(n, dtype=np.float64),
        't_first_ge500': -np.ones(n, dtype=np.float64),
        'months_at_or_above_200': np.zeros(n, dtype=np.float64),
        'months_at_or_above_350': np.zeros(n, dtype=np.float64),
        'months_at_or_above_500': np.zeros(n, dtype=np.float64),
        't_sustained_ge200': -np.ones(n, dtype=np.float64),
        't_sustained_ge350': -np.ones(n, dtype=np.float64),
        't_sustained_ge500': -np.ones(n, dtype=np.float64),
        'art_start': -np.ones(n, dtype=np.float64),
        'on_art': np.zeros(n, dtype=np.bool_),
        'time_since_art': np.zeros(n, dtype=np.float64),
        'is_aids': np.zeros(n, dtype=np.bool_),
        'on_prep': np.zeros(n, dtype=np.bool_),
        'transmission_type': np.empty(n, dtype=np.int32),
        'hetero_core': np.zeros(n, dtype=np.bool_),
        'activity_multiplier': np.ones(n, dtype=np.float64),
        'vl_suppressed': np.zeros(n, dtype=np.bool_),
        'dropped_out': np.zeros(n, dtype=np.bool_),
        'time_to_suppression': -np.ones(n, dtype=np.float64),
        'time_of_infection': -np.ones(n, dtype=np.float64),
        'time_of_diagnosis': -np.ones(n, dtype=np.float64),
        'alive': np.ones(n, dtype=np.bool_),
        'time_since_dropout': np.zeros(n, dtype=np.float64),
        'cd4_decline_factor': np.ones(n, dtype=np.float64),
        'cd4_ceiling': np.full(n, 800.0, dtype=np.float64),
        'stage1_end': np.zeros(n, dtype=np.float64),
        'stage2_end': np.zeros(n, dtype=np.float64),
        'stage3_end': np.zeros(n, dtype=np.float64),
        'stage4_end': np.zeros(n, dtype=np.float64),
        'stage1_drop_ratio': np.zeros(n, dtype=np.float64),
        'stage2_recover_ratio': np.zeros(n, dtype=np.float64),
        'stage3_decline': np.zeros(n, dtype=np.float64),
        'stage4_decline': np.zeros(n, dtype=np.float64),
        'stage5_decline': np.zeros(n, dtype=np.float64),
        'cd4_decline_base': np.ones(n, dtype=np.float64),
        'cd4_decline_accelerator': np.ones(n, dtype=np.float64),

        't_art_eligible': -np.ones(n, dtype=np.float64),
        't_art_start_actual': -np.ones(n, dtype=np.float64),

        'infection_year': -np.ones(n, dtype=np.int32),
        'infection_month_index': -np.ones(n, dtype=np.int32),
        'diagnosis_year': -np.ones(n, dtype=np.int32),
        'diagnosis_month_index': -np.ones(n, dtype=np.int32),
        'discounted_qaly': np.zeros(n, dtype=np.float64),
        'discounted_cost': np.zeros(n, dtype=np.float64),
        'discounted_qaly_post_diagnosis': np.zeros(n, dtype=np.float64),
        'discounted_cost_post_diagnosis': np.zeros(n, dtype=np.float64),
    }
    return pop

def initialize_stage_parameters(pop, rng):
    n = len(pop['age'])

    s1 = 3.0 + rng.normal(0, 0.5, size=n)
    s2 = 6.0 + rng.normal(0, 1.0, size=n)
    s3 = 60.0 + rng.normal(0, 15.0, size=n)
    s4 = 90.0 + rng.normal(0, 40.0, size=n)

    stages = np.stack([s1, s2, s3, s4], axis=1)
    stages.sort(axis=1)
    
    pop['stage1_end'] = stages[:, 0]
    pop['stage2_end'] = stages[:, 1]
    pop['stage3_end'] = stages[:, 2]
    pop['stage4_end'] = stages[:, 3]

    pop['stage1_drop_ratio'] = 0.15 * rng.uniform(0.8, 1.2, size=n)
    pop['stage2_recover_ratio'] = 0.05 * rng.uniform(0.8, 1.2, size=n)
    pop['stage3_decline'] = 2.0 * rng.uniform(0.8, 1.2, size=n)
    pop['stage4_decline'] = 5.0 * rng.uniform(0.8, 1.2, size=n)
    pop['stage5_decline'] = 10.0 * rng.uniform(0.8, 1.2, size=n)
    
    pop['cd4_decline_base'] = rng.normal(1.0, 0.3, size=n)
    pop['cd4_decline_accelerator'] = rng.normal(1.0, 0.5, size=n)
    
    return pop

def assign_hetero_core_activity(pop, hetero_indices, rng, core_share=HETERO_CORE_SHARE):
    """Assign a small high-activity heterosexual core group.

    HETERO_CORE_ACTIVITY_MULTIPLIER is treated as the activity ratio of the
    latent high-activity/bridge compartment versus low-activity heterosexuals.
    Because the base core share is only about 1%, the stratum-wide mean contact
    volume increases only modestly while transmission is concentrated among
    persons with higher partner turnover, concurrency, migration/foreign-
    nationality bridge features, or delayed testing access.
    """
    hetero_indices = np.asarray(hetero_indices, dtype=np.int64)
    if len(hetero_indices) == 0:
        return

    pop['hetero_core'][hetero_indices] = False
    pop['activity_multiplier'][hetero_indices] = 1.0

    n_core = int(round(len(hetero_indices) * float(core_share)))
    n_core = min(max(n_core, 0), len(hetero_indices))
    if n_core == 0:
        return

    high_to_low_ratio = max(float(HETERO_CORE_ACTIVITY_MULTIPLIER), 1.0)
    core_indices = rng.choice(hetero_indices, size=n_core, replace=False)
    pop['hetero_core'][core_indices] = True
    pop['activity_multiplier'][core_indices] = high_to_low_ratio

def weighted_choice_without_replacement(indices, size, weights, rng):
    """Choice helper that normalizes weights and handles empty selections."""
    indices = np.asarray(indices, dtype=np.int64)
    if size <= 0:
        return np.array([], dtype=np.int64)

    weights = np.asarray(weights, dtype=np.float64)
    if len(indices) != len(weights):
        raise ValueError("indices and weights must have the same length")

    weight_sum = float(np.sum(weights))
    if weight_sum <= 0.0:
        probs = None
    else:
        probs = weights / weight_sum
    return rng.choice(indices, size=size, replace=False, p=probs)

def choose_initial_hetero_hiv(indices, size, pop, rng):
    """Sample initial heterosexual HIV+ persons with enrichment in core networks."""
    indices = np.asarray(indices, dtype=np.int64)
    weights = np.where(
        pop['hetero_core'][indices],
        HETERO_CORE_INITIAL_HIV_WEIGHT,
        1.0
    )
    return weighted_choice_without_replacement(indices, size, weights, rng)

def initialize_population(rng):
    total_pop = msm_population_size + heterosexual_population_size
    pop = create_population_arrays(total_pop)
    
    # MSM
    pop['population_type'][:msm_population_size] = 0
    pop['gender'][:msm_population_size] = 0
    pop['age'][:msm_population_size] = generate_ages_from_distribution(
        msm_population_size,
        male_age_groups_lower, male_age_groups_upper,
        male_age_groups_prob, rng
    )
    pop['transmission_type'][:msm_population_size] = 2
    
    # Heterosexual
    pop['population_type'][msm_population_size:] = 1
    n_het = heterosexual_population_size
    genders = rng.choice([0, 1], size=n_het, p=[0.5, 0.5])
    pop['gender'][msm_population_size:] = genders
    
    het_indices = np.arange(msm_population_size, total_pop)
    male_het_mask = genders == 0
    female_het_mask = genders == 1
    
    if male_het_mask.sum() > 0:
        pop['age'][het_indices[male_het_mask]] = generate_ages_from_distribution(
            male_het_mask.sum(),
            male_age_groups_lower, male_age_groups_upper,
            male_age_groups_prob, rng
        )
        pop['transmission_type'][het_indices[male_het_mask]] = 1
    
    if female_het_mask.sum() > 0:
        pop['age'][het_indices[female_het_mask]] = generate_ages_from_distribution(
            female_het_mask.sum(),
            female_age_groups_lower, female_age_groups_upper,
            female_age_groups_prob, rng
        )
        pop['transmission_type'][het_indices[female_het_mask]] = 1

    assign_hetero_core_activity(pop, het_indices, rng)
    
    pop = initialize_stage_parameters(pop, rng)
    
    return pop

def initialize_hiv(pop, rng):
    # ----------------------------------------------------------------

    #   MSM:    population_type=0, gender=0

    # ----------------------------------------------------------------
    msm_indices   = np.where((pop['population_type'] == 0) & (pop['gender'] == 0))[0]
    het_m_indices = np.where((pop['population_type'] == 1) & (pop['gender'] == 0))[0]
    het_f_indices = np.where((pop['population_type'] == 1) & (pop['gender'] == 1))[0]

    n_msm_hiv   = int(initial_total_hiv * INITIAL_HIV_MSM_RATIO)
    n_het_m_hiv = int(initial_total_hiv * INITIAL_HIV_HETM_RATIO)
    n_het_f_hiv = initial_total_hiv - n_msm_hiv - n_het_m_hiv

    msm_hiv_indices   = rng.choice(msm_indices,   size=n_msm_hiv,   replace=False)
    het_m_hiv_indices = choose_initial_hetero_hiv(het_m_indices, n_het_m_hiv, pop, rng)
    het_f_hiv_indices = choose_initial_hetero_hiv(het_f_indices, n_het_f_hiv, pop, rng)
    
    hiv_indices = np.concatenate([msm_hiv_indices, het_m_hiv_indices, het_f_hiv_indices])
    pop['hiv_positive'][hiv_indices] = True

    infection_times = -rng.uniform(0, 120, size=len(hiv_indices))
    pop['time_of_infection'][hiv_indices] = infection_times

    #   CD4 <50:     29.62%
    #   CD4 50-200:  28.23%
    #   CD4 200-350: 18.78%
    #   CD4 350-500: 12.38%
    #   CD4 500+:    10.99%
    cd4_percent = [0.2962, 0.2823, 0.1878, 0.1238, 0.1099]
    cd4_ranges = [(0.0, 50.0), (50.0, 200.0), (200.0, 350.0), (350.0, 500.0), (500.0, 1000.0)]

    art_indices = rng.choice(hiv_indices, size=initial_on_art, replace=False)
    pop['on_art'][art_indices] = True
    pop['cd4_at_art_start'][art_indices] = pop['cd4'][art_indices]
    pop['diagnosed'][art_indices] = True
    pop['time_of_diagnosis'][art_indices] = -1.0
    pop['time_since_art'][art_indices] = 6.0
    tts = generate_time_to_suppression_skewed(len(art_indices), rng)
    pop['time_to_suppression'][art_indices] = tts

    untreated_indices = np.setdiff1d(hiv_indices, art_indices)
    untreated_count = len(untreated_indices)
    
    cd4_counts = [int(untreated_count * p) for p in cd4_percent]
    diff = untreated_count - sum(cd4_counts)
    if diff > 0:
        cd4_counts[-1] += diff
    
    rng.shuffle(untreated_indices)
    start_idx = 0
    for cat_i, n_cat in enumerate(cd4_counts):
        if n_cat <= 0:
            continue
        cat_range = cd4_ranges[cat_i]
        cat_indices = untreated_indices[start_idx:start_idx + n_cat]
        start_idx += n_cat
        pop['cd4'][cat_indices] = rng.uniform(cat_range[0], cat_range[1], size=n_cat)
    
    # =======================================================================

    # =======================================================================

    
    undiagnosed_ratio = 0.40
    art_scheduled_ratio = 0.30
    normal_flow_ratio = 0.25
    dropout_ratio = 0.05
    
    n_undiagnosed = int(untreated_count * undiagnosed_ratio)
    n_art_scheduled = int(untreated_count * art_scheduled_ratio)
    n_normal_flow = int(untreated_count * normal_flow_ratio)
    n_dropout = untreated_count - n_undiagnosed - n_art_scheduled - n_normal_flow

    rng.shuffle(untreated_indices)
    idx_start = 0

    undiagnosed_indices = untreated_indices[idx_start:idx_start + n_undiagnosed]
    idx_start += n_undiagnosed

    art_scheduled_indices = untreated_indices[idx_start:idx_start + n_art_scheduled]
    idx_start += n_art_scheduled
    if len(art_scheduled_indices) > 0:
        pop['diagnosed'][art_scheduled_indices] = True
        pop['time_of_diagnosis'][art_scheduled_indices] = -rng.uniform(1, 12, size=len(art_scheduled_indices))

        pop['art_start'][art_scheduled_indices] = rng.uniform(1, 12, size=len(art_scheduled_indices))

    normal_flow_indices = untreated_indices[idx_start:idx_start + n_normal_flow]
    idx_start += n_normal_flow
    if len(normal_flow_indices) > 0:
        pop['diagnosed'][normal_flow_indices] = True
        pop['time_of_diagnosis'][normal_flow_indices] = -rng.uniform(1, 6, size=len(normal_flow_indices))

    dropout_indices = untreated_indices[idx_start:]
    if len(dropout_indices) > 0:
        pop['diagnosed'][dropout_indices] = True
        pop['time_of_diagnosis'][dropout_indices] = -rng.uniform(6, 36, size=len(dropout_indices))
        pop['dropped_out'][dropout_indices] = True
        pop['time_since_dropout'][dropout_indices] = rng.uniform(1, 24, size=len(dropout_indices))
    # =======================================================================

    cd4_vals_art = rng.lognormal(mean=np.log(500), sigma=0.5, size=len(art_indices))
    cd4_vals_art = np.clip(cd4_vals_art, 50, 1000)
    pop['cd4'][art_indices] = cd4_vals_art
    pop['vl'][art_indices] = assign_viral_load_vectorized(cd4_vals_art, rng)

    untreated_cd4s = pop['cd4'][untreated_indices]
    pop['vl'][untreated_indices] = assign_viral_load_vectorized(untreated_cd4s, rng)
    
    pop['is_aids'][hiv_indices] = (pop['cd4'][hiv_indices] < 200.0)
    pop['cd4_decline_factor'][hiv_indices] = rng.uniform(0.8, 1.2, size=len(hiv_indices))

    for idx in art_indices:

        noise_val = rng.uniform(-0.2, 0.2)
        blip_rand = rng.random()
        pop['vl'][idx] = simulate_viral_load_randomized(
            pop['vl'][idx],
            pop['time_since_art'][idx],
            pop['time_to_suppression'][idx],
            LOG_DROP, vl_decay_rate_per_month,
            noise_val, blip_rand
        )
        pop['vl_suppressed'][idx] = (pop['vl'][idx] <= 1.5)
        nadir = pop['cd4'][idx]
        pop['cd4_ceiling'][idx] = assign_cd4_ceiling_from_nadir_vectorized(np.array([nadir]), rng)[0]
    
    return pop

def generate_newcomers(n_new, fraction_msm, year_idx, rng):
    n_msm = int(n_new * fraction_msm)
    n_het = n_new - n_msm
    
    pop = create_population_arrays(n_new)
    
    # MSM
    pop['population_type'][:n_msm] = 0
    pop['gender'][:n_msm] = 0
    msm_age_probs = get_age_distribution_by_year(year_idx, "male")
    pop['age'][:n_msm] = generate_ages_from_distribution(
        n_msm, male_age_groups_lower, male_age_groups_upper, msm_age_probs, rng
    )
    pop['transmission_type'][:n_msm] = 2
    
    # Heterosexual
    if n_het > 0:
        pop['population_type'][n_msm:] = 1
        genders_het = rng.choice([0, 1], size=n_het, p=[0.5, 0.5])
        pop['gender'][n_msm:] = genders_het
        
        male_mask = genders_het == 0
        female_mask = genders_het == 1
        
        if male_mask.sum() > 0:
            male_probs = get_age_distribution_by_year(year_idx, "male")
            pop['age'][n_msm:][male_mask] = generate_ages_from_distribution(
                male_mask.sum(), male_age_groups_lower, male_age_groups_upper, male_probs, rng
            )
            pop['transmission_type'][n_msm:][male_mask] = 1
        
        if female_mask.sum() > 0:
            female_probs = get_age_distribution_by_year(year_idx, "female")
            pop['age'][n_msm:][female_mask] = generate_ages_from_distribution(
                female_mask.sum(), female_age_groups_lower, female_age_groups_upper, female_probs, rng
            )
            pop['transmission_type'][n_msm:][female_mask] = 1

        assign_hetero_core_activity(pop, np.arange(n_msm, n_new), rng)
    
    pop = initialize_stage_parameters(pop, rng)
    
    return pop

def generate_foreigners(n_msm, n_het_m, n_het_f, year_idx, current_month, rng):
    n = n_msm + n_het_m + n_het_f
    pop = create_population_arrays(n)
    
    if n == 0:
        return pop

    pop['hiv_positive'][:] = True
    pop['time_of_infection'][:] = current_month
    
    # ----------------------------------------------------------------

    # ----------------------------------------------------------------
    n_male = n_msm + n_het_m
    pop['gender'][:n_male] = 0
    pop['gender'][n_male:] = 1
    
    pop['population_type'][:n_msm] = 0       # MSM
    pop['population_type'][n_msm:] = 1
    
    pop['transmission_type'][:n_msm] = 0     # MSM
    pop['transmission_type'][n_msm:] = 1

    if n > n_msm:
        assign_hetero_core_activity(
            pop, np.arange(n_msm, n), rng, core_share=HETERO_CORE_FOREIGN_SHARE
        )
    
    # ----------------------------------------------------------------

    # ----------------------------------------------------------------
    if n_male > 0:
        male_probs = get_age_distribution_by_year(year_idx, "male")
        pop['age'][:n_male] = generate_ages_from_distribution(
            n_male, male_age_groups_lower, male_age_groups_upper, male_probs, rng
        )
    if n_het_f > 0:
        female_probs = get_age_distribution_by_year(year_idx, "female")
        pop['age'][n_male:] = generate_ages_from_distribution(
            n_het_f, female_age_groups_lower, female_age_groups_upper, female_probs, rng
        )

    pop['cd4'][:] = rng.uniform(200, 600, size=n)
    pop['vl'][:] = assign_viral_load_vectorized(pop['cd4'], rng)
    pop['is_aids'][:] = pop['cd4'] < 200.0

    pop = initialize_stage_parameters(pop, rng)
    pop['cd4_decline_factor'][:] = rng.uniform(0.8, 1.2, size=n)
    
    return pop

def append_population(pop, new_pop):
    for key in pop:
        pop[key] = np.concatenate([pop[key], new_pop[key]])
    return pop

# =============================================================================

# =============================================================================

# Transmission-source attribution (model-based decomposition of the monthly
# force of infection; the model does NOT retain individual infector-infectee
# pairs). Four mutually exclusive, exhaustive infectious-source states.
TRANSMISSION_SOURCE_CATEGORIES = (
    "undiagnosed",            # source & ~diagnosed
    "diagnosed_pre_art",      # source & diagnosed & ~on_art & ~dropped_out
    "on_art_unsuppressed",    # source & on_art & ~dropped_out
    "disengaged",             # source & dropped_out
)

def _attribute_sources(pop, src_idx, weights, n_inf):
    """Allocate realised infections n_inf across source states in proportion to
    each state's contribution to the (risk-group) force of infection.

    weights[i] = contact_rate * activity_i * p_infect(VL_i) for source i.
    Read-only: consumes no RNG and does not alter the simulation state.
    Returns {category: {"force": float, "alloc": float}}.
    """
    cats = TRANSMISSION_SOURCE_CATEGORIES
    if len(src_idx) == 0:
        return {c: {"force": 0.0, "alloc": 0.0} for c in cats}
    diagnosed = pop['diagnosed'][src_idx]
    on_art    = pop['on_art'][src_idx]
    dropped   = pop['dropped_out'][src_idx]
    masks = {
        "undiagnosed":         ~diagnosed,
        "diagnosed_pre_art":   diagnosed & (~on_art) & (~dropped),
        "on_art_unsuppressed": on_art & (~dropped),
        "disengaged":          dropped,
    }
    F = {c: float(np.sum(weights[masks[c]])) for c in cats}
    Ftot = sum(F.values())
    out = {}
    for c in cats:
        out[c] = {"force": F[c],
                  "alloc": (n_inf * F[c] / Ftot) if Ftot > 0 else 0.0}
    return out

def process_transmissions_in_month(pop, month, label, rng, is_burnin=False):
    alive_mask = pop['alive']

    source_mask = (
        pop['hiv_positive'] & 
        alive_mask & 
        (~pop['vl_suppressed'])
    )

    sus_mask = (~pop['hiv_positive']) & alive_mask

    msm_source_mask = source_mask & (pop['population_type'] == 0)
    het_source_mask = source_mask & (pop['population_type'] == 1)
    msm_sus_mask = sus_mask & (pop['population_type'] == 0)
    het_sus_mask = sus_mask & (pop['population_type'] == 1)
    
    msm_source_indices = np.where(msm_source_mask)[0]
    het_source_indices = np.where(het_source_mask)[0]
    msm_sus_indices = np.where(msm_sus_mask)[0]
    het_sus_indices = np.where(het_sus_mask)[0]
    
    I_msm = len(msm_source_indices)
    I_het = len(het_source_indices)
    S_msm = len(msm_sus_indices)
    S_het = len(het_sus_indices)
    
    new_infected = np.zeros(len(pop['alive']), dtype=np.bool_)
    
    if I_msm == 0 and I_het == 0:
        return new_infected, None

    # =================================================================

    # =================================================================
    if S_msm > 0:
        n_prep_msm = np.sum(pop['on_prep'][msm_sus_indices])
        n_nonprep_msm = S_msm - n_prep_msm
        S_eff_msm = n_nonprep_msm + n_prep_msm * (1.0 - prep_effectiveness)
    else:
        S_eff_msm = 0.0
    
    if S_het > 0:
        het_sus_activity = pop['activity_multiplier'][het_sus_indices]
        het_sus_prep_modifier = np.where(
            pop['on_prep'][het_sus_indices],
            1.0 - prep_effectiveness,
            1.0
        )
        S_activity_het = float(np.sum(het_sus_activity))
        S_eff_het = float(np.sum(het_sus_activity * het_sus_prep_modifier))
    else:
        S_activity_het = 0.0
        S_eff_het = 0.0

    transmission_multiplier = get_transmission_multiplier(month, is_burnin)
    if I_msm > 0 and S_msm > 0:
        vl_bins_msm = get_vl_bin_indices_vectorized(pop['vl'][msm_source_indices])
        p_infect_msm = p_infect_msm_bins[vl_bins_msm]
        total_force_msm = np.sum(contact_rate_msm * p_infect_msm) * transmission_multiplier
        lambda_msm = total_force_msm / max(S_msm, 1)
    else:
        lambda_msm = 0.0

    if I_het > 0 and S_het > 0:
        vl_bins_het = get_vl_bin_indices_vectorized(pop['vl'][het_source_indices])
        p_infect_het = p_infect_het_bins[vl_bins_het]
        het_source_activity = pop['activity_multiplier'][het_source_indices]
        total_force_het = np.sum(contact_rate_het * het_source_activity * p_infect_het) * transmission_multiplier
        lambda_het = total_force_het / max(S_activity_het, 1.0)
    else:
        lambda_het = 0.0

    new_infection_count_msm = rng.poisson(lambda_msm * S_eff_msm) if S_eff_msm > 0 else 0
    new_infection_count_het = rng.poisson(lambda_het * S_eff_het) if S_eff_het > 0 else 0
    
    new_infection_count_msm = min(new_infection_count_msm, S_msm)
    new_infection_count_het = min(new_infection_count_het, S_het)

    if new_infection_count_msm > 0 and S_msm > 0:
        weights = np.where(pop['on_prep'][msm_sus_indices], 1.0 - prep_effectiveness, 1.0)
        weights /= weights.sum()
        chosen = rng.choice(msm_sus_indices, size=new_infection_count_msm, replace=False, p=weights)
        new_infected[chosen] = True

    if new_infection_count_het > 0 and S_het > 0:
        weights = pop['activity_multiplier'][het_sus_indices] * np.where(
            pop['on_prep'][het_sus_indices],
            1.0 - prep_effectiveness,
            1.0
        )
        weights /= weights.sum()
        chosen = rng.choice(het_sus_indices, size=new_infection_count_het, replace=False, p=weights)
        new_infected[chosen] = True

    # -------------------------------------------------------------------------
    # Transmission-source attribution (analytic period only; read-only).
    # Allocate this month's realised infections (MSM and heterosexual separately)
    # across infectious-source states in proportion to each state's contribution
    # to that risk group's force of infection. The transmission_multiplier scales
    # all source states equally and therefore cancels in the within-group shares;
    # it is recorded separately. No RNG is consumed here.
    # -------------------------------------------------------------------------
    source_attr = None
    if not is_burnin:
        if I_msm > 0:
            w_msm = (contact_rate_msm
                     * p_infect_msm_bins[get_vl_bin_indices_vectorized(pop['vl'][msm_source_indices])])
        else:
            w_msm = np.zeros(0, dtype=np.float64)
        if I_het > 0:
            w_het = (contact_rate_het
                     * pop['activity_multiplier'][het_source_indices]
                     * p_infect_het_bins[get_vl_bin_indices_vectorized(pop['vl'][het_source_indices])])
        else:
            w_het = np.zeros(0, dtype=np.float64)
        source_attr = {
            "mult": float(transmission_multiplier),
            "msm": _attribute_sources(pop, msm_source_indices, w_msm, int(new_infection_count_msm)),
            "het": _attribute_sources(pop, het_source_indices, w_het, int(new_infection_count_het)),
        }

    return new_infected, source_attr

# =============================================================================

# =============================================================================

def process_single_month(pop, month, label, rngs, 
                         msm_monthly_growth_prob, hetero_monthly_growth_prob,
                         record_results=True, is_burnin=False, sim_year_idx=1,
                         delay_mean_override=None, delay_sd_override=None):
    prev_cd4 = pop['cd4'].copy()
    year_idx = (month - 1) // 12 + 1
    
    results = {
        'new_infections': 0,
        'new_diagnoses': 0,
        'deaths': 0,
        'deaths_under200': 0,
        'deaths_over200': 0,
        'art_start_under200': 0,
        'art_start_over200': 0,
        'tests': 0,
        'tests_msm': 0,
        'tests_het': 0,
        'tests_hiv_pos': 0,
        'tests_hiv_neg': 0,
        'breakthrough_undiagnosed': 0,
        'newly_infected_indices': np.array([], dtype=np.int64),
        'new_diag_indices': np.array([], dtype=np.int64),
        'dead_hiv_indices': np.array([], dtype=np.int64),
        'art_start_indices': np.array([], dtype=np.int64),
    }

    # -----------------------------------------------------------------

    # -----------------------------------------------------------------
    if isinstance(rngs, np.random.Generator):

        rngs = {name: rngs for name in EVENT_RNG_NAMES}
    elif not isinstance(rngs, dict):
        raise TypeError("rngs must be dict[str, np.random.Generator] or np.random.Generator")

    if "init" not in rngs:
        raise KeyError("rngs must contain key: 'init'")
    rng = rngs["init"]
    
    # =================================================================

    # =================================================================
    if not is_burnin:
        rng = rngs["prep"]
        msm_sus = (pop['population_type'] == 0) & (~pop['hiv_positive']) & pop['alive']
        het_sus = (pop['population_type'] == 1) & (~pop['hiv_positive']) & pop['alive']
        
        msm_sus_indices = np.where(msm_sus)[0]
        het_sus_indices = np.where(het_sus)[0]
        
        # ----- MSM -----
        if len(msm_sus_indices) > 0:

            on_prep_msm_mask = pop['on_prep'][msm_sus_indices]
            on_prep_msm_indices = msm_sus_indices[on_prep_msm_mask]
            if len(on_prep_msm_indices) > 0:
                dropouts = rng.random(len(on_prep_msm_indices)) < prep_dropout_rate_msm
                pop['on_prep'][on_prep_msm_indices[dropouts]] = False

            n_on_prep_msm = np.sum(pop['on_prep'][msm_sus_indices])
            c_current_msm = n_on_prep_msm / len(msm_sus_indices)
            gap_msm = max(max_prep_uptake_msm - c_current_msm, 0.0)
            
            if gap_msm > 0:

                base_rate_msm = 3.0 / (prep_growth_years_msm * 12.0)
                p_month_msm = base_rate_msm * gap_msm

                off_prep_msm_mask = ~pop['on_prep'][msm_sus_indices]
                off_prep_msm_indices = msm_sus_indices[off_prep_msm_mask]
                
                if len(off_prep_msm_indices) > 0:
                    new_prep_msm = rng.random(len(off_prep_msm_indices)) < p_month_msm
                    pop['on_prep'][off_prep_msm_indices[new_prep_msm]] = True
        
        # ----- Hetero -----
        if len(het_sus_indices) > 0:

            on_prep_het_mask = pop['on_prep'][het_sus_indices]
            on_prep_het_indices = het_sus_indices[on_prep_het_mask]
            if len(on_prep_het_indices) > 0:
                dropouts = rng.random(len(on_prep_het_indices)) < prep_dropout_rate_het
                pop['on_prep'][on_prep_het_indices[dropouts]] = False

            n_on_prep_het = np.sum(pop['on_prep'][het_sus_indices])
            c_current_het = n_on_prep_het / len(het_sus_indices)
            gap_het = max(max_prep_uptake_het - c_current_het, 0.0)
            
            if gap_het > 0:
                base_rate_het = 3.0 / (prep_growth_years_het * 12.0)
                p_month_het = base_rate_het * gap_het
                
                off_prep_het_mask = ~pop['on_prep'][het_sus_indices]
                off_prep_het_indices = het_sus_indices[off_prep_het_mask]
                
                if len(off_prep_het_indices) > 0:
                    new_prep_het = rng.random(len(off_prep_het_indices)) < p_month_het
                    pop['on_prep'][off_prep_het_indices[new_prep_het]] = True
    
    # =================================================================

    # =================================================================
    rng = rngs["risk_exit"]
    alive_indices = np.where(pop['alive'])[0]
    exit_probs = calc_age_exit_probability_vectorized(pop['age'][alive_indices])
    exits = rng.random(len(alive_indices)) < exit_probs
    pop['alive'][alive_indices[exits]] = False
    
    # =================================================================

    # =================================================================
    rng = rngs["art_delay"]
    if label == "0 month":
        mask_art = (
            pop['diagnosed'] &
            (~pop['on_art']) &
            (~pop['dropped_out']) &
            (pop['art_start'] < 0) &
            pop['alive']
        )
        if is_art_cd4_wait_policy_active(month, is_burnin):
            mask_art &= meets_art_cd4_wait_threshold(pop['cd4'], month, is_burnin)
        pop['art_start'][mask_art] = month

        first_start_mask = mask_art & (pop['t_art_start_actual'] < 0)
        pop['t_art_eligible'][first_start_mask] = float(month)
    else:

        vl_threshold_log = math.log10(5000.0)
        base_mask = (
            pop['diagnosed'] &
            (~pop['on_art']) &
            (~pop['dropped_out']) &
            (pop['art_start'] < 0) &
            pop['alive']
        )
        
        if is_art_cd4_wait_policy_active(month, is_burnin):
            base_mask &= meets_art_cd4_wait_threshold(pop['cd4'], month, is_burnin)
        elif USE_CD4_VL_CRITERIA:
            base_mask &= (pop['cd4'] < 500.0) | (pop['vl'] >= vl_threshold_log)
        
        eligible_indices = np.where(base_mask)[0]
        if len(eligible_indices) > 0:

            cfg = SCENARIO_CONFIGS.get(label, {"delay_mean": 2.0, "delay_sd": 0.5})
            delay_mean = float(cfg.get("delay_mean", 2.0))
            delay_sd = float(cfg.get("delay_sd", 0.5))
            if delay_mean_override is not None:
                delay_mean = float(delay_mean_override)
            if delay_sd_override is not None:
                delay_sd = float(delay_sd_override)

            first_start_idx = eligible_indices[pop['t_art_start_actual'][eligible_indices] < 0]
            pop['t_art_eligible'][first_start_idx] = float(month)

            delays = gamma_delay(delay_mean, delay_sd, len(eligible_indices), rng)
            pop['art_start'][eligible_indices] = month + delays

    art_ready_mask = (
        pop['diagnosed'] &
        (~pop['dropped_out']) &
        (~pop['on_art']) &
        (pop['art_start'] > 0) &
        (pop['art_start'] <= month)
    )
    art_ready_indices = np.where(art_ready_mask)[0]
    
    if len(art_ready_indices) > 0:
        if is_art_cd4_wait_policy_active(month, is_burnin):
            ready_cd4_ok = meets_art_cd4_wait_threshold(
                pop['cd4'][art_ready_indices], month, is_burnin
            )
            not_yet_eligible = art_ready_indices[~ready_cd4_ok]
            if len(not_yet_eligible) > 0:
                pop['art_start'][not_yet_eligible] = -1.0
                never_started = not_yet_eligible[pop['t_art_start_actual'][not_yet_eligible] < 0]
                pop['t_art_eligible'][never_started] = -1.0
            art_ready_indices = art_ready_indices[ready_cd4_ok]

    if len(art_ready_indices) > 0:
        rng = rngs["art_start"]
        pop['on_art'][art_ready_indices] = True
        pop['cd4_at_art_start'][art_ready_indices] = pop['cd4'][art_ready_indices]

        first_start = art_ready_indices[pop['t_art_start_actual'][art_ready_indices] < 0]
        pop['t_art_start_actual'][first_start] = float(month)
        tts_arr = generate_time_to_suppression_skewed(len(art_ready_indices), rng)
        pop['time_to_suppression'][art_ready_indices] = tts_arr
        pop['cd4_ceiling'][art_ready_indices] = assign_cd4_ceiling_from_nadir_vectorized(
            pop['cd4'][art_ready_indices], rng
        )

        under200_mask = (pop['cd4'][art_ready_indices] < 200.0) | pop['is_aids'][art_ready_indices]
        results['art_start_under200'] = under200_mask.sum()
        results['art_start_over200'] = (~under200_mask).sum()

        results['art_start_indices'] = np.concatenate(
            [results['art_start_indices'], art_ready_indices]
        )
    
    # =================================================================

    # =================================================================
    rng = rngs["dropout"]

    mask_pre_art = pop['diagnosed'] & (~pop['on_art']) & (~pop['dropped_out']) & pop['alive']
    pre_art_indices = np.where(mask_pre_art)[0]
    if len(pre_art_indices) > 0:
        dropouts = rng.random(len(pre_art_indices)) < dropout_prob_pre_art
        drop_idx = pre_art_indices[dropouts]
        pop['dropped_out'][drop_idx] = True
        pop['art_start'][drop_idx] = -1

        never_started = drop_idx[pop['t_art_start_actual'][drop_idx] < 0]
        pop['t_art_eligible'][never_started] = -1.0
        pop['time_since_art'][drop_idx] = 0.0
        pop['vl_suppressed'][drop_idx] = False
        pop['on_art'][drop_idx] = False
        pop['time_to_suppression'][drop_idx] = -1

    mask_on_art = pop['on_art'] & (~pop['dropped_out']) & pop['alive']
    on_art_indices = np.where(mask_on_art)[0]
    if len(on_art_indices) > 0:
        dropouts2 = rng.random(len(on_art_indices)) < dropout_prob_on_art
        drop_idx2 = on_art_indices[dropouts2]
        pop['dropped_out'][drop_idx2] = True
        pop['on_art'][drop_idx2] = False
        pop['time_since_art'][drop_idx2] = 0.0
        pop['vl_suppressed'][drop_idx2] = False
        pop['time_to_suppression'][drop_idx2] = -1
        pop['art_start'][drop_idx2] = -1
    
    # =================================================================

    # =================================================================
    mask_on_art_active = pop['on_art'] & (~pop['dropped_out']) & pop['alive']
    art_active_indices = np.where(mask_on_art_active)[0]
    
    if len(art_active_indices) > 0:
        pop['time_since_art'][art_active_indices] += 1.0

        rng = rngs["vl_update"]
        n_art = len(art_active_indices)
        noise_arr = rng.uniform(-0.2, 0.2, size=n_art)
        blip_rand_arr = rng.random(size=n_art)
        
        new_vl = simulate_viral_load_batch(
            pop['vl'][art_active_indices],
            pop['time_since_art'][art_active_indices],
            pop['time_to_suppression'][art_active_indices],
            LOG_DROP, vl_decay_rate_per_month,
            noise_arr, blip_rand_arr
        )
        pop['vl'][art_active_indices] = new_vl
        pop['vl_suppressed'][art_active_indices] = new_vl <= 1.5

        cd4_start = np.where(
            pop['cd4_at_art_start'][art_active_indices] >= 0,
            pop['cd4_at_art_start'][art_active_indices],
            pop['cd4'][art_active_indices]
        )
        new_cd4 = simulate_cd4_recovery_batch(
            pop['cd4'][art_active_indices],
            cd4_start,
            pop['cd4_ceiling'][art_active_indices],
            pop['vl_suppressed'][art_active_indices]
        )
        pop['cd4'][art_active_indices] = new_cd4
    
    # =================================================================

    # =================================================================
    mask_not_on_art = pop['hiv_positive'] & (~pop['on_art']) & pop['alive']
    not_on_art_indices = np.where(mask_not_on_art)[0]
    
    if len(not_on_art_indices) > 0:
        rng = rngs["cd4_decline"]
        decline_vals = calc_cd4_decline_vectorized(pop, not_on_art_indices, month, rng)
        pop['cd4'][not_on_art_indices] -= decline_vals

    pop['cd4'] = np.maximum(pop['cd4'], 0.0)
    pop['is_aids'] = pop['cd4'] < 200.0
    
    # =================================================================

    # =================================================================
    pop['cd4_nadir'] = np.minimum(pop['cd4_nadir'], pop['cd4'])
    
    update_cd4_tracking(
        pop['cd4'], prev_cd4,
        pop['t_first_ge200'], pop['t_first_ge350'], pop['t_first_ge500'],
        pop['months_at_or_above_200'], pop['months_at_or_above_350'], pop['months_at_or_above_500'],
        pop['t_sustained_ge200'], pop['t_sustained_ge350'], pop['t_sustained_ge500'],
        float(month)
    )
    
    # =================================================================

    # =================================================================
    mask_dropout = pop['dropped_out'] & pop['alive']
    pop['time_since_dropout'][mask_dropout] += 1.0
    
    mask_reconnect = pop['dropped_out'] & pop['hiv_positive'] & (~pop['on_art']) & pop['alive']
    reconnect_indices = np.where(mask_reconnect)[0]
    
    if len(reconnect_indices) > 0:
        rng = rngs["reconnect"]
        cd4_now = pop['cd4'][reconnect_indices]
        base_prob = 0.6 / (1 + np.exp(0.01 * (cd4_now - 250)))
        base_prob = np.maximum(base_prob, 0.01)
        factor = 1.0 - np.exp(-0.05 * pop['time_since_dropout'][reconnect_indices])
        prob_return = base_prob * factor
        
        returns = rng.random(len(reconnect_indices)) < prob_return
        return_indices = reconnect_indices[returns]
        
        pop['dropped_out'][return_indices] = False
        pop['time_since_dropout'][return_indices] = 0.0

        if is_art_cd4_wait_policy_active(month, is_burnin):
            eligible = return_indices[
                meets_art_cd4_wait_threshold(pop['cd4'][return_indices], month, is_burnin)
            ]
        else:
            eligible = return_indices
        
        if len(eligible) > 0:
            delays = rng.normal(2.0, 1.0, size=len(eligible))
            delays = np.clip(delays, 1.0, 4.0)
            pop['art_start'][eligible] = month + delays

            first_start_eligible = eligible[pop['t_art_start_actual'][eligible] < 0]
            pop['t_art_eligible'][first_start_eligible] = float(month)
            pop['time_since_art'][eligible] = 0.0
            pop['vl_suppressed'][eligible] = False
            pop['on_art'][eligible] = False
            pop['time_to_suppression'][eligible] = -1
    
    # =================================================================

    # =================================================================
    rng = rngs["testing"]
    alive_indices = np.where(pop['alive'])[0]
    testing_rate_multiplier = min(1.0 + 0.02 * year_idx, 2.0)
    
    # -----------------------------------------------------------------

    # -----------------------------------------------------------------
    prep_user_mask = pop['on_prep'] & pop['alive'] & (~pop['diagnosed'])
    prep_test_prob = 4.0 / 12.0
    prep_test_mask = (rng.random(len(pop['alive'])) < prep_test_prob) & prep_user_mask
    prep_test_indices = np.where(prep_test_mask)[0]
    
    # -----------------------------------------------------------------

    # -----------------------------------------------------------------

    test_probs = np.zeros(len(pop['alive']), dtype=np.float64)

    undiag_hiv = pop['hiv_positive'] & (~pop['diagnosed']) & pop['alive']
    test_probs[undiag_hiv & (pop['cd4'] < 200)] = 0.3 / 12.0
    test_probs[undiag_hiv & (pop['cd4'] >= 200) & (pop['cd4'] < 350)] = 0.1 / 12.0
    test_probs[undiag_hiv & (pop['cd4'] >= 350)] = 0.05 / 12.0
    undiag_hiv_test_boost = get_undiagnosed_hiv_test_boost(month, is_burnin)
    if undiag_hiv_test_boost != 1.0:
        test_probs[undiag_hiv] *= undiag_hiv_test_boost

    hiv_neg_nonprep = (~pop['hiv_positive']) & pop['alive'] & (~pop['on_prep'])
    test_probs[hiv_neg_nonprep & (pop['population_type'] == 0)] = 0.04 / 12.0
    test_probs[hiv_neg_nonprep & (pop['population_type'] == 1)] = 0.01 / 12.0

    young_mask = pop['age'] < 30
    old_mask = pop['age'] > 50
    test_probs[young_mask | old_mask] *= 0.8

    test_probs[pop['on_prep']] = 0.0
    
    test_probs *= testing_rate_multiplier

    test_probs[pop['diagnosed']] = 0.0

    regular_test_mask = (rng.random(len(test_probs)) < test_probs) & pop['alive'] & (~pop['diagnosed'])
    regular_test_indices = np.where(regular_test_mask)[0]

    if TEST_CAPACITY_POLICY_ON:
        if is_burnin:

            if BURNIN_TEST_CAPACITY_POLICY_ON:
                capacity = get_burnin_test_capacity(month)
            else:
                capacity = TEST_CAPACITY_X
        else:

            if TEST_CAPACITY_Z == 1:
                capacity = TEST_CAPACITY_Y
            elif sim_year_idx >= TEST_CAPACITY_Z:
                capacity = TEST_CAPACITY_Y
            else:
                step = (TEST_CAPACITY_Y - TEST_CAPACITY_X) / float(TEST_CAPACITY_Z - 1)
                capacity = TEST_CAPACITY_X + step * (sim_year_idx - 1)
        capacity = int(capacity)
        
        if len(regular_test_indices) > capacity:
            regular_test_indices = rng.choice(regular_test_indices, capacity, replace=False)
        elif len(regular_test_indices) < capacity:

            pool = np.where(pop['alive'] & (~pop['diagnosed']) & (~pop['on_prep']))[0]
            pool = np.setdiff1d(pool, regular_test_indices)
            needed = min(capacity - len(regular_test_indices), len(pool))
            if needed > 0:
                additional = rng.choice(pool, needed, replace=False)
                regular_test_indices = np.concatenate([regular_test_indices, additional])
    
    # -----------------------------------------------------------------

    # -----------------------------------------------------------------
    test_indices = np.unique(np.concatenate([prep_test_indices, regular_test_indices]))

    results['tests'] = len(test_indices)
    results['tests_msm'] = np.sum(pop['population_type'][test_indices] == 0)
    results['tests_het'] = np.sum(pop['population_type'][test_indices] == 1)
    results['tests_hiv_pos'] = np.sum(pop['hiv_positive'][test_indices])
    results['tests_hiv_neg'] = np.sum(~pop['hiv_positive'][test_indices])

    new_diag_mask = pop['hiv_positive'][test_indices] & (~pop['diagnosed'][test_indices])
    new_diag_indices = test_indices[new_diag_mask]
    pop['diagnosed'][new_diag_indices] = True
    pop['time_of_diagnosis'][new_diag_indices] = month

    # ================================================================

    # ================================================================
    if not is_burnin:
        sim_month_idx = month - BURN_IN_MONTHS
        if sim_month_idx >= 1 and len(new_diag_indices) > 0:
            unset = pop['diagnosis_year'][new_diag_indices] < 0
            if np.any(unset):
                idx_set = new_diag_indices[unset]
                pop['diagnosis_year'][idx_set] = sim_year_idx
                pop['diagnosis_month_index'][idx_set] = sim_month_idx
    results['new_diagnoses'] = len(new_diag_indices)
    results['new_diag_indices'] = new_diag_indices

    pop['on_prep'][new_diag_indices] = False

    # =================================================================

    # =================================================================
    if len(new_diag_indices) > 0:
        if label == "0 month":

            immediate_art_mask = (
                (~pop['on_art'][new_diag_indices]) &
                (~pop['dropped_out'][new_diag_indices]) &
                (pop['art_start'][new_diag_indices] < 0) &
                pop['alive'][new_diag_indices]
            )
            if is_art_cd4_wait_policy_active(month, is_burnin):
                immediate_art_mask &= meets_art_cd4_wait_threshold(
                    pop['cd4'][new_diag_indices], month, is_burnin
                )
            immediate_art_idx = new_diag_indices[immediate_art_mask]
            if len(immediate_art_idx) > 0:
                rng = rngs["art_start"]
                pop['art_start'][immediate_art_idx] = month

                first_start_mask = pop['t_art_start_actual'][immediate_art_idx] < 0
                pop['t_art_eligible'][immediate_art_idx[first_start_mask]] = float(month)

                pop['on_art'][immediate_art_idx] = True
                pop['cd4_at_art_start'][immediate_art_idx] = pop['cd4'][immediate_art_idx]
                first_start_actual = immediate_art_idx[pop['t_art_start_actual'][immediate_art_idx] < 0]
                pop['t_art_start_actual'][first_start_actual] = float(month)
                tts_arr = generate_time_to_suppression_skewed(len(immediate_art_idx), rng)
                pop['time_to_suppression'][immediate_art_idx] = tts_arr
                pop['cd4_ceiling'][immediate_art_idx] = assign_cd4_ceiling_from_nadir_vectorized(
                    pop['cd4'][immediate_art_idx], rng
                )

                under200_mask = (pop['cd4'][immediate_art_idx] < 200.0) | pop['is_aids'][immediate_art_idx]
                results['art_start_under200'] += under200_mask.sum()
                results['art_start_over200'] += (~under200_mask).sum()

                results['art_start_indices'] = np.concatenate(
                    [results['art_start_indices'], immediate_art_idx]
                )
        else:

            vl_threshold_log = math.log10(5000.0)
            base_mask = (
                (~pop['on_art'][new_diag_indices]) &
                (~pop['dropped_out'][new_diag_indices]) &
                (pop['art_start'][new_diag_indices] < 0) &
                pop['alive'][new_diag_indices]
            )
            
            if is_art_cd4_wait_policy_active(month, is_burnin):
                base_mask &= meets_art_cd4_wait_threshold(
                    pop['cd4'][new_diag_indices], month, is_burnin
                )
            elif USE_CD4_VL_CRITERIA:
                base_mask &= (pop['cd4'][new_diag_indices] < 500.0) | (pop['vl'][new_diag_indices] >= vl_threshold_log)
            
            eligible_new_diag = new_diag_indices[base_mask]
            if len(eligible_new_diag) > 0:
                rng = rngs["art_delay"]

                cfg = SCENARIO_CONFIGS.get(label, {"delay_mean": 2.0, "delay_sd": 0.5})
                delay_mean = float(cfg.get("delay_mean", 2.0))
                delay_sd = float(cfg.get("delay_sd", 0.5))
                if delay_mean_override is not None:
                    delay_mean = float(delay_mean_override)
                if delay_sd_override is not None:
                    delay_sd = float(delay_sd_override)

                first_start_idx = eligible_new_diag[pop['t_art_start_actual'][eligible_new_diag] < 0]
                pop['t_art_eligible'][first_start_idx] = float(month)
                
                delays = gamma_delay(delay_mean, delay_sd, len(eligible_new_diag), rng)
                pop['art_start'][eligible_new_diag] = month + delays

    # =================================================================

    # =================================================================
    new_infected, _src_attr = process_transmissions_in_month(
        pop, month, label, rngs["infection"], is_burnin=is_burnin
    )
    results['transmission_source'] = _src_attr
    newly_infected_indices = np.where(new_infected)[0]
    
    if len(newly_infected_indices) > 0:
        rng = rngs["infection_attr"]
        pop['hiv_positive'][newly_infected_indices] = True
        pop['diagnosed'][newly_infected_indices] = False
        pop['time_of_infection'][newly_infected_indices] = month

        # ================================================================

        # ================================================================
        if not is_burnin:
            sim_month_idx = month - BURN_IN_MONTHS
            if sim_month_idx >= 1:
                unset = pop['infection_year'][newly_infected_indices] < 0
                if np.any(unset):
                    idx_set = newly_infected_indices[unset]
                    pop['infection_year'][idx_set] = sim_year_idx
                    pop['infection_month_index'][idx_set] = sim_month_idx
        
        new_cd4 = generate_new_cd4_for_infection_vectorized(
            pop['age'][newly_infected_indices], rng
        )
        pop['cd4'][newly_infected_indices] = new_cd4
        pop['vl'][newly_infected_indices] = assign_viral_load_vectorized(new_cd4, rng)
        pop['on_art'][newly_infected_indices] = False
        pop['time_since_art'][newly_infected_indices] = 0.0
        pop['is_aids'][newly_infected_indices] = new_cd4 < 200.0
        pop['vl_suppressed'][newly_infected_indices] = False
        pop['time_to_suppression'][newly_infected_indices] = -1
        pop['dropped_out'][newly_infected_indices] = False
        pop['time_since_dropout'][newly_infected_indices] = 0.0
        pop['art_start'][newly_infected_indices] = -1

        pop['t_art_eligible'][newly_infected_indices] = -1.0
        pop['t_art_start_actual'][newly_infected_indices] = -1.0
        pop['cd4_decline_factor'][newly_infected_indices] = rng.uniform(0.8, 1.2, size=len(newly_infected_indices))
        pop['cd4_decline_base'][newly_infected_indices] = rng.normal(1.0, 0.3, size=len(newly_infected_indices))
        pop['cd4_decline_accelerator'][newly_infected_indices] = rng.normal(1.0, 0.5, size=len(newly_infected_indices))
    
    results['new_infections'] = len(newly_infected_indices)
    results['newly_infected_indices'] = newly_infected_indices
    
    # =================================================================

    # =================================================================
    rng = rngs["death"]

    hiv_alive_mask = pop['hiv_positive'] & pop['alive']
    hiv_alive_indices = np.where(hiv_alive_mask)[0]
    
    if len(hiv_alive_indices) > 0:

        mr_hiv = calc_mortality_risk_batch_v143(
            pop['cd4'][hiv_alive_indices],
            pop['age'][hiv_alive_indices],
            pop['is_aids'][hiv_alive_indices],
            pop['gender'][hiv_alive_indices],
            pop['on_art'][hiv_alive_indices],
            pop['vl_suppressed'][hiv_alive_indices],
            pop['dropped_out'][hiv_alive_indices],
            pop['diagnosed'][hiv_alive_indices],
            male_age_upper, male_mortality_rate,
            female_age_upper, female_mortality_rate,
            HIV_MORTALITY_MODE_CODE,
            _ACTIVE_RATE_TABLE,
            HIV_MORTALITY_MULTIPLIERS,
            HIV_MORTALITY_EXCESS_RATE_MULTIPLIER,
            HIV_MORTALITY_MULTIPLIER_SCALER,
        )
        deaths_hiv = rng.random(len(hiv_alive_indices)) < mr_hiv
        dead_hiv_indices = hiv_alive_indices[deaths_hiv]
        pop['alive'][dead_hiv_indices] = False
        results['deaths'] += len(dead_hiv_indices)
        
        under200_deaths = (pop['cd4'][dead_hiv_indices] < 200.0) | pop['is_aids'][dead_hiv_indices]
        results['deaths_under200'] = under200_deaths.sum()
        results['deaths_over200'] = (~under200_deaths).sum()
        results['dead_hiv_indices'] = dead_hiv_indices

        results['deaths_by_cd4_state'] = compute_cd4_state_counts(
            pop['cd4'][dead_hiv_indices], pop['on_art'][dead_hiv_indices],
            pop['vl_suppressed'][dead_hiv_indices], pop['dropped_out'][dead_hiv_indices])
        results['pmonths_by_cd4_state'] = compute_cd4_state_counts(
            pop['cd4'][hiv_alive_indices], pop['on_art'][hiv_alive_indices],
            pop['vl_suppressed'][hiv_alive_indices], pop['dropped_out'][hiv_alive_indices])

        results['expected_by_cd4_state'] = compute_cd4_state_weighted(
            pop['cd4'][hiv_alive_indices], pop['on_art'][hiv_alive_indices],
            pop['vl_suppressed'][hiv_alive_indices], pop['dropped_out'][hiv_alive_indices], mr_hiv)
    else:
        results['dead_hiv_indices'] = np.array([], dtype=np.int64)
        results['deaths_by_cd4_state'] = np.zeros((4, 4), dtype=np.int64)
        results['pmonths_by_cd4_state'] = np.zeros((4, 4), dtype=np.int64)
        results['expected_by_cd4_state'] = np.zeros((4, 4), dtype=np.float64)

    uninf_alive_mask = (~pop['hiv_positive']) & pop['alive']
    uninf_alive_indices = np.where(uninf_alive_mask)[0]
    
    if len(uninf_alive_indices) > 0:
        mr_uninf = calc_mortality_risk_uninfected_batch(
            pop['age'][uninf_alive_indices],
            pop['gender'][uninf_alive_indices],
            male_age_upper, male_mortality_rate,
            female_age_upper, female_mortality_rate
        )
        deaths_uninf = rng.random(len(uninf_alive_indices)) < mr_uninf
        dead_uninf_indices = uninf_alive_indices[deaths_uninf]
        pop['alive'][dead_uninf_indices] = False
        results['deaths'] += len(dead_uninf_indices)
    
    # =================================================================

    # =================================================================
    rng = rngs["newcomer"]
    if is_burnin:
        n_msm   = int(rng.poisson(FOREIGN_INFLUX_BURNIN_MSM_ANNUAL  / 12.0))
        n_het_m = int(rng.poisson(FOREIGN_INFLUX_BURNIN_HETM_ANNUAL / 12.0))
        n_het_f = int(rng.poisson(FOREIGN_INFLUX_BURNIN_HETF_ANNUAL / 12.0))
    else:
        n_msm   = int(rng.poisson(FOREIGN_INFLUX_MAIN_MSM_ANNUAL    / 12.0))
        n_het_m = int(rng.poisson(FOREIGN_INFLUX_MAIN_HETM_ANNUAL   / 12.0))
        n_het_f = int(rng.poisson(FOREIGN_INFLUX_MAIN_HETF_ANNUAL   / 12.0))
    
    foreigners_pop = generate_foreigners(
        n_msm, n_het_m, n_het_f, year_idx, month, rng,
    )

    results['foreign_inflow_msm']   = n_msm
    results['foreign_inflow_het_m'] = n_het_m
    results['foreign_inflow_het_f'] = n_het_f

    if n_msm + n_het_m + n_het_f > 0:
        pop = append_population(pop, foreigners_pop)

    # =================================================================

    # =================================================================

    current_alive = np.sum(pop['alive'])

    base_newcomers = int(current_alive * NEWCOMER_RATE / 12)

    if current_alive < TARGET_POPULATION:
        population_deficit = TARGET_POPULATION - current_alive
        catch_up_newcomers = int(population_deficit * 0.02 / 12)
    else:
        catch_up_newcomers = 0
    
    monthly_newcomers = base_newcomers + catch_up_newcomers
    
    if monthly_newcomers > 0:
        newcomers_pop = generate_newcomers(monthly_newcomers, 1/3, year_idx, rng)
        pop = append_population(pop, newcomers_pop)

    # =================================================================

    # =================================================================
    pop['age'] += 1.0 / 12.0

    # =================================================================

    # =================================================================
    if not is_burnin:
        sim_month_idx = month - BURN_IN_MONTHS
        if sim_month_idx >= 1:

            if USE_DISCOUNTING_QALY:
                df = 1.0 / ((1.0 + DISCOUNT_RATE_QALY) ** (sim_month_idx / 12.0))
            else:
                df = 1.0

            alive_hiv = pop['alive'] & pop['hiv_positive']
            if np.any(alive_hiv):
                cd4_vals = pop['cd4'][alive_hiv]
                util = np.where(
                    cd4_vals >= 500.0, UTILITY_CD4_GE500,
                    np.where(
                        cd4_vals >= 350.0, UTILITY_CD4_350_499,
                        np.where(
                            cd4_vals >= 200.0, UTILITY_CD4_200_349,
                            UTILITY_CD4_LT200
                        )
                    )
                )

                diag_not_on_art = pop['diagnosed'][alive_hiv] & (~pop['on_art'][alive_hiv])
                util = np.where(
                    diag_not_on_art,
                    util - MENTAL_HEALTH_UTIL_LOSS,
                    util
                )
                util = np.maximum(util, 0.0)
                pop['discounted_qaly'][alive_hiv] += util * (1.0 / 12.0) * df

            alive_hiv_diag = alive_hiv & pop['diagnosed']
            if np.any(alive_hiv_diag):
                cd4_vals_d = pop['cd4'][alive_hiv_diag]
                util_d = np.where(
                    cd4_vals_d >= 500.0, UTILITY_CD4_GE500,
                    np.where(
                        cd4_vals_d >= 350.0, UTILITY_CD4_350_499,
                        np.where(
                            cd4_vals_d >= 200.0, UTILITY_CD4_200_349,
                            UTILITY_CD4_LT200
                        )
                    )
                )

                not_on_art_d = ~pop['on_art'][alive_hiv_diag]
                util_d = np.where(
                    not_on_art_d,
                    util_d - MENTAL_HEALTH_UTIL_LOSS,
                    util_d
                )
                util_d = np.maximum(util_d, 0.0)
                pop['discounted_qaly_post_diagnosis'][alive_hiv_diag] += util_d * (1.0 / 12.0) * df

    # =================================================================

    # =================================================================
    results['breakthrough_undiagnosed'] = np.sum(
        pop['hiv_positive'] & pop['on_prep'] & (~pop['diagnosed']) & pop['alive']
    )

    return pop, results

# =============================================================================

# =============================================================================

def run_burnin(seed=None):
    if seed is None:
        seed = SEED

    np.random.seed(seed)
    
    rngs = create_event_rngs(seed)
    
    pop = initialize_population(rngs["init"])
    pop = initialize_hiv(pop, rngs["init"])

    msm_monthly_growth_prob = 0.0
    hetero_monthly_growth_prob = 0.0

    BURNIN_POLICY = "2 months"

    burnin_yearly_summary = []
    burnin_new_infections_year = 0
    burnin_new_diagnoses_year = 0
    burnin_deaths_year = 0
    
    for month in range(1, BURN_IN_MONTHS + 1):
        pop, month_results = process_single_month(
            pop, month, BURNIN_POLICY, rngs,
            msm_monthly_growth_prob, hetero_monthly_growth_prob,
            record_results=True,
            is_burnin=True,
            sim_year_idx=0
        )

        burnin_new_infections_year += month_results['new_infections']
        burnin_new_diagnoses_year += month_results['new_diagnoses']
        burnin_deaths_year += month_results['deaths']

        if month % 12 == 0:
            burn_in_year = month // 12

            n_alive = np.sum(pop['alive'])
            n_msm_alive = np.sum((pop['population_type'] == 0) & pop['alive'])
            n_het_alive = np.sum((pop['population_type'] == 1) & pop['alive'])
            n_hiv = np.sum(pop['hiv_positive'] & pop['alive'])
            n_diagnosed = np.sum(pop['diagnosed'] & pop['hiv_positive'] & pop['alive'])
            n_undiagnosed_hiv = np.sum(pop['hiv_positive'] & (~pop['diagnosed']) & pop['alive'])
            n_on_art = np.sum(pop['on_art'] & pop['alive'])
            n_vl_suppressed = np.sum(pop['vl_suppressed'] & pop['on_art'] & pop['alive'])
            n_on_prep = np.sum(pop['on_prep'] & (~pop['hiv_positive']) & pop['alive'])
            n_aids = np.sum(pop['is_aids'] & pop['hiv_positive'] & pop['alive'])
            n_prep_hivpos_undiag = np.sum(pop['on_prep'] & pop['hiv_positive'] & (~pop['diagnosed']) & pop['alive'])

            hiv_alive = pop['hiv_positive'] & pop['alive']
            n_cd4_under200 = np.sum(hiv_alive & (pop['cd4'] < 200))
            n_cd4_200_349 = np.sum(hiv_alive & (pop['cd4'] >= 200) & (pop['cd4'] < 350))
            n_cd4_350_499 = np.sum(hiv_alive & (pop['cd4'] >= 350) & (pop['cd4'] < 500))
            n_cd4_500plus = np.sum(hiv_alive & (pop['cd4'] >= 500))
            
            summary_row = {
                'burn_in_year': burn_in_year,
                'total_alive': n_alive,
                'msm_alive': n_msm_alive,
                'het_alive': n_het_alive,
                'hiv_positive': n_hiv,
                'diagnosed': n_diagnosed,
                'undiagnosed_hiv': n_undiagnosed_hiv,
                'on_art': n_on_art,
                'vl_suppressed': n_vl_suppressed,
                'on_prep': n_on_prep,
                'aids': n_aids,
                'prep_hivpos_undiag': n_prep_hivpos_undiag,
                'cd4_under200': n_cd4_under200,
                'cd4_200_349': n_cd4_200_349,
                'cd4_350_499': n_cd4_350_499,
                'cd4_500plus': n_cd4_500plus,
                'new_infections_year': burnin_new_infections_year,
                'new_diagnoses_year': burnin_new_diagnoses_year,
                'deaths_year': burnin_deaths_year,
            }
            burnin_yearly_summary.append(summary_row)

            burnin_new_infections_year = 0
            burnin_new_diagnoses_year = 0
            burnin_deaths_year = 0

    burnin_df = pd.DataFrame(burnin_yearly_summary)
    burnin_csv_path = os.path.join(OUTPUT_DIR, "burnin_summary.csv")
    burnin_df.to_csv(burnin_csv_path, index=False)
    
    return pop, rngs

def copy_population(pop):
    return {key: arr.copy() for key, arr in pop.items()}

# =============================================================================

# =============================================================================

def run_simulation(label="0 month", initial_pop=None, seed=None, scenario_params=None):
    if seed is None:
        seed = SEED

    np.random.seed(seed)
    
    rngs = create_event_rngs(seed)
    rng = rngs["init"]
    rng_prep = rngs["prep"]
    
    # =========================================================================

    # =========================================================================
    def get_param(name, default_global):
        if scenario_params is not None and name in scenario_params:
            return scenario_params[name]
        return default_global
    
    _initial_prep_uptake_msm = get_param('initial_prep_uptake_msm', initial_prep_uptake_msm)
    _max_prep_uptake_msm = get_param('max_prep_uptake_msm', max_prep_uptake_msm)
    _prep_growth_years_msm = get_param('prep_growth_years_msm', prep_growth_years_msm)
    _initial_prep_uptake_het = get_param('initial_prep_uptake_het', initial_prep_uptake_het)
    _max_prep_uptake_het = get_param('max_prep_uptake_het', max_prep_uptake_het)
    _prep_growth_years_het = get_param('prep_growth_years_het', prep_growth_years_het)
    _prep_effectiveness = get_param('prep_effectiveness', prep_effectiveness)
    _prep_dropout_rate_msm = get_param('prep_dropout_rate_msm', prep_dropout_rate_msm)
    _prep_dropout_rate_het = get_param('prep_dropout_rate_het', prep_dropout_rate_het)
    _TEST_CAPACITY_POLICY_ON = get_param('TEST_CAPACITY_POLICY_ON', TEST_CAPACITY_POLICY_ON)
    _TEST_CAPACITY_X = get_param('TEST_CAPACITY_X', TEST_CAPACITY_X)
    _TEST_CAPACITY_Y = get_param('TEST_CAPACITY_Y', TEST_CAPACITY_Y)
    _TEST_CAPACITY_Z = get_param('TEST_CAPACITY_Z', TEST_CAPACITY_Z)
    
    # =========================================================================

    # =========================================================================
    _globals_backup = None
    if scenario_params is not None:

        _globals_backup = {
            'initial_prep_uptake_msm': initial_prep_uptake_msm,
            'max_prep_uptake_msm': max_prep_uptake_msm,
            'prep_growth_years_msm': prep_growth_years_msm,
            'initial_prep_uptake_het': initial_prep_uptake_het,
            'max_prep_uptake_het': max_prep_uptake_het,
            'prep_growth_years_het': prep_growth_years_het,
            'prep_effectiveness': prep_effectiveness,
            'prep_dropout_rate_msm': prep_dropout_rate_msm,
            'prep_dropout_rate_het': prep_dropout_rate_het,
            'TEST_CAPACITY_POLICY_ON': TEST_CAPACITY_POLICY_ON,
            'TEST_CAPACITY_X': TEST_CAPACITY_X,
            'TEST_CAPACITY_Y': TEST_CAPACITY_Y,
            'TEST_CAPACITY_Z': TEST_CAPACITY_Z,
        }

        import sys
        current_module = sys.modules[__name__]
        for key, value in scenario_params.items():
            if hasattr(current_module, key):
                setattr(current_module, key, value)

    msm_monthly_growth_prob = (_max_prep_uptake_msm - _initial_prep_uptake_msm) / (_prep_growth_years_msm * 12.0)
    hetero_monthly_growth_prob = (_max_prep_uptake_het - _initial_prep_uptake_het) / (_prep_growth_years_het * 12.0)
    
    if initial_pop is None:

        pop = initialize_population(rng)
        pop = initialize_hiv(pop, rng)
        
        BURNIN_POLICY = "2 months"
        burnin_yearly_summary = []
        burnin_new_infections_year = 0
        burnin_new_diagnoses_year = 0
        burnin_deaths_year = 0
        
        for month in range(1, BURN_IN_MONTHS + 1):
            pop, month_results = process_single_month(
                pop, month, BURNIN_POLICY, rngs,
                0.0, 0.0,
                record_results=True,
                is_burnin=True,
                sim_year_idx=0
            )
            burnin_new_infections_year += month_results['new_infections']
            burnin_new_diagnoses_year += month_results['new_diagnoses']
            burnin_deaths_year += month_results['deaths']
            
            if month % 12 == 0:
                burn_in_year = month // 12
                n_alive = np.sum(pop['alive'])
                n_hiv = np.sum(pop['hiv_positive'] & pop['alive'])
                n_on_art = np.sum(pop['on_art'] & pop['alive'])
                
                summary_row = {
                    'burn_in_year': burn_in_year,
                    'total_alive': n_alive,
                    'hiv_positive': n_hiv,
                    'on_art': n_on_art,
                    'new_infections_year': burnin_new_infections_year,
                    'new_diagnoses_year': burnin_new_diagnoses_year,
                    'deaths_year': burnin_deaths_year,
                }
                burnin_yearly_summary.append(summary_row)
                burnin_new_infections_year = 0
                burnin_new_diagnoses_year = 0
                burnin_deaths_year = 0
        
        scenario_label = SCENARIO_CONFIGS[label]["csv_label"]
        burnin_df = pd.DataFrame(burnin_yearly_summary)
        burnin_df.to_csv(os.path.join(OUTPUT_DIR, f"{scenario_label}_burnin_summary.csv"), index=False)
    else:

        pop = copy_population(initial_pop)
    
    # =========================================================================

    # =========================================================================

    pop['t_first_ge200'][:] = -1
    pop['t_first_ge350'][:] = -1
    pop['t_first_ge500'][:] = -1
    pop['months_at_or_above_200'][:] = 0
    pop['months_at_or_above_350'][:] = 0
    pop['months_at_or_above_500'][:] = 0
    pop['t_sustained_ge200'][:] = -1
    pop['t_sustained_ge350'][:] = -1
    pop['t_sustained_ge500'][:] = -1
    
    # =========================================================================

    # =========================================================================

    # MSM
    msm_sus_init = (pop['population_type'] == 0) & (~pop['hiv_positive']) & pop['alive']
    msm_sus_init_indices = np.where(msm_sus_init)[0]
    if len(msm_sus_init_indices) > 0 and initial_prep_uptake_msm > 0:
        n_initial_prep_msm = int(len(msm_sus_init_indices) * initial_prep_uptake_msm)
        if n_initial_prep_msm > 0:
            initial_prep_msm_indices = rng_prep.choice(msm_sus_init_indices, n_initial_prep_msm, replace=False)
            pop['on_prep'][initial_prep_msm_indices] = True
            print(f"  [PrEP init] MSM: {n_initial_prep_msm:,} ({initial_prep_uptake_msm*100:.1f}%) set on PrEP")
    
    # Hetero
    het_sus_init = (pop['population_type'] == 1) & (~pop['hiv_positive']) & pop['alive']
    het_sus_init_indices = np.where(het_sus_init)[0]
    if len(het_sus_init_indices) > 0 and initial_prep_uptake_het > 0:
        n_initial_prep_het = int(len(het_sus_init_indices) * initial_prep_uptake_het)
        if n_initial_prep_het > 0:
            initial_prep_het_indices = rng_prep.choice(het_sus_init_indices, n_initial_prep_het, replace=False)
            pop['on_prep'][initial_prep_het_indices] = True
            print(f"  [PrEP init] Hetero: {n_initial_prep_het:,} ({initial_prep_uptake_het*100:.1f}%) set on PrEP")
    
    # =========================================================================

    # =========================================================================

    yearly_infections = []
    yearly_new_deaths = []
    yearly_on_art_current = []
    yearly_new_diagnoses = []
    yearly_cumulative_infected = []
    yearly_population_total = []
    yearly_population_msm = []
    yearly_population_het = []
    yearly_population_undiagnosed = []
    yearly_new_art_start_under200 = []
    yearly_new_art_start_over200 = []
    yearly_new_deaths_under200 = []
    yearly_new_deaths_over200 = []
    yearly_tests = []
    yearly_tests_msm = []
    yearly_tests_het = []
    yearly_tests_hiv_pos = []
    yearly_tests_hiv_neg = []
    yearly_breakthrough_undiagnosed = []

    monthly_art_start_under200 = np.zeros(simulation_months + 1, dtype=np.int32)
    monthly_art_start_over200 = np.zeros(simulation_months + 1, dtype=np.int32)
    monthly_tests = np.zeros(simulation_months + 1, dtype=np.int32)
    monthly_tests_msm = np.zeros(simulation_months + 1, dtype=np.int32)
    monthly_tests_het = np.zeros(simulation_months + 1, dtype=np.int32)
    monthly_tests_hiv_pos = np.zeros(simulation_months + 1, dtype=np.int32)
    monthly_tests_hiv_neg = np.zeros(simulation_months + 1, dtype=np.int32)
    monthly_breakthrough_undiagnosed = np.zeros(simulation_months + 1, dtype=np.int32)

    # =========================================================================

    # =========================================================================
    monthly_new_infections_msm = np.zeros(simulation_months + 1, dtype=np.int32)
    monthly_new_infections_het = np.zeros(simulation_months + 1, dtype=np.int32)
    monthly_new_deaths_msm = np.zeros(simulation_months + 1, dtype=np.int32)
    monthly_new_deaths_het = np.zeros(simulation_months + 1, dtype=np.int32)
    monthly_art_starts_msm = np.zeros(simulation_months + 1, dtype=np.int32)
    monthly_art_starts_het = np.zeros(simulation_months + 1, dtype=np.int32)

    yearly_new_infections_msm = []
    yearly_new_infections_het = []
    yearly_new_deaths_msm = []
    yearly_new_deaths_het = []
    yearly_art_starts_msm = []
    yearly_art_starts_het = []

    year_new_infections_msm = 0
    year_new_infections_het = 0
    year_new_deaths_msm = 0
    year_new_deaths_het = 0
    year_art_starts_msm = 0
    year_art_starts_het = 0

    # =========================================================================

    # =========================================================================
    monthly_foreign_inflow_msm   = np.zeros(simulation_months + 1, dtype=np.int32)
    monthly_foreign_inflow_het_m = np.zeros(simulation_months + 1, dtype=np.int32)
    monthly_foreign_inflow_het_f = np.zeros(simulation_months + 1, dtype=np.int32)
    
    yearly_foreign_inflow_msm   = []
    yearly_foreign_inflow_het_m = []
    yearly_foreign_inflow_het_f = []
    
    year_foreign_inflow_msm   = 0
    year_foreign_inflow_het_m = 0
    year_foreign_inflow_het_f = 0

    cumulative_infected = np.sum(pop['hiv_positive'] & pop['alive'])

    total_new_infections = 0
    monthly_deaths = 0
    month_new_diagnoses = 0
    monthly_deaths_under200 = 0
    monthly_deaths_over200 = 0
    art_start_under200_this_year = 0
    art_start_over200_this_year = 0
    yearly_test_count = 0
    yearly_test_count_msm = 0
    yearly_test_count_het = 0
    yearly_test_count_hiv_pos = 0
    yearly_test_count_hiv_neg = 0

    max_aggregator_key = (years + 1) * AGGREGATOR_KEY_MULTIPLIERS[0]
    agg_person_months = np.zeros(max_aggregator_key, dtype=np.int64)
    agg_newly_infected = np.zeros(max_aggregator_key, dtype=np.int64)
    agg_newly_diagnosed = np.zeros(max_aggregator_key, dtype=np.int64)
    agg_deaths = np.zeros(max_aggregator_key, dtype=np.int64)

    year_susceptible_msm = np.zeros(years, dtype=np.float64)
    year_susceptible_het = np.zeros(years, dtype=np.float64)
    year_susceptible_onprep = np.zeros(years, dtype=np.float64)
    year_susceptible_offprep = np.zeros(years, dtype=np.float64)
    year_infections_msm = np.zeros(years, dtype=np.float64)
    year_infections_het = np.zeros(years, dtype=np.float64)
    year_infections_onprep = np.zeros(years, dtype=np.float64)
    year_infections_offprep = np.zeros(years, dtype=np.float64)

    # Transmission-source attribution accumulators (analytic period; read-only).
    _src_alloc = {(rg, c): np.zeros(years, dtype=np.float64)
                  for rg in ("msm", "het") for c in TRANSMISSION_SOURCE_CATEGORIES}
    _src_force = {(rg, c): np.zeros(years, dtype=np.float64)
                  for rg in ("msm", "het") for c in TRANSMISSION_SOURCE_CATEGORIES}
    _src_mult = np.zeros(years, dtype=np.float64)
    _year_foreign_inflow = np.zeros(years, dtype=np.float64)

    year_deaths_by_cd4_state = np.zeros((years, 4, 4), dtype=np.float64)
    year_pmonths_by_cd4_state = np.zeros((years, 4, 4), dtype=np.float64)
    year_expected_by_cd4_state = np.zeros((years, 4, 4), dtype=np.float64)

    for sim_month in range(1, simulation_months + 1):

        actual_month = BURN_IN_MONTHS + sim_month
        
        year_idx = (sim_month - 1) // 12 + 1
        yid = year_idx - 1

        if (sim_month - 1) % 12 == 0:
            mask_alive_hivneg = pop['alive'] & (~pop['hiv_positive'])
            year_susceptible_msm[yid] = np.sum(mask_alive_hivneg & (pop['population_type'] == 0))
            year_susceptible_het[yid] = np.sum(mask_alive_hivneg & (pop['population_type'] == 1))
            year_susceptible_onprep[yid] = np.sum(mask_alive_hivneg & pop['on_prep'])
            year_susceptible_offprep[yid] = np.sum(mask_alive_hivneg & (~pop['on_prep']))

        pop, month_results = process_single_month(
            pop, actual_month, label, rngs,
            msm_monthly_growth_prob, hetero_monthly_growth_prob,
            record_results=True,
            is_burnin=False,
            sim_year_idx=year_idx
        )

        monthly_art_start_under200[sim_month] = month_results['art_start_under200']
        monthly_art_start_over200[sim_month] = month_results['art_start_over200']
        art_start_under200_this_year += month_results['art_start_under200']
        art_start_over200_this_year += month_results['art_start_over200']
        
        monthly_tests[sim_month] = month_results['tests']
        monthly_tests_msm[sim_month] = month_results['tests_msm']
        monthly_tests_het[sim_month] = month_results['tests_het']
        monthly_tests_hiv_pos[sim_month] = month_results['tests_hiv_pos']
        monthly_tests_hiv_neg[sim_month] = month_results['tests_hiv_neg']
        monthly_breakthrough_undiagnosed[sim_month] = month_results['breakthrough_undiagnosed']

        yearly_test_count += month_results['tests']
        yearly_test_count_msm += month_results['tests_msm']
        yearly_test_count_het += month_results['tests_het']
        yearly_test_count_hiv_pos += month_results['tests_hiv_pos']
        yearly_test_count_hiv_neg += month_results['tests_hiv_neg']
        
        month_new_diagnoses += month_results['new_diagnoses']
        total_new_infections += month_results['new_infections']
        monthly_deaths += month_results['deaths']
        monthly_deaths_under200 += month_results['deaths_under200']
        monthly_deaths_over200 += month_results['deaths_over200']

        year_deaths_by_cd4_state[yid] += month_results['deaths_by_cd4_state']
        year_pmonths_by_cd4_state[yid] += month_results['pmonths_by_cd4_state']
        year_expected_by_cd4_state[yid] += month_results['expected_by_cd4_state']

        newly_infected_indices = month_results['newly_infected_indices']
        if len(newly_infected_indices) > 0:
            year_infections_msm[yid] += np.sum(pop['population_type'][newly_infected_indices] == 0)
            year_infections_het[yid] += np.sum(pop['population_type'][newly_infected_indices] == 1)
            year_infections_onprep[yid] += np.sum(pop['on_prep'][newly_infected_indices])
            year_infections_offprep[yid] += np.sum(~pop['on_prep'][newly_infected_indices])

        _src = month_results.get('transmission_source')
        if _src is not None:
            for _rg in ("msm", "het"):
                for _c in TRANSMISSION_SOURCE_CATEGORIES:
                    _src_alloc[(_rg, _c)][yid] += _src[_rg][_c]["alloc"]
                    _src_force[(_rg, _c)][yid] += _src[_rg][_c]["force"]
            _src_mult[yid] = _src["mult"]
        _year_foreign_inflow[yid] += (month_results.get('foreign_inflow_msm', 0)
                                      + month_results.get('foreign_inflow_het_m', 0)
                                      + month_results.get('foreign_inflow_het_f', 0))

        # ---------------------------------------------------------------

        # ---------------------------------------------------------------
        if ENABLE_HETEROGENEITY_TRACKING:

            if len(newly_infected_indices) > 0:
                _inf_msm = int(np.sum(pop['population_type'][newly_infected_indices] == 0))
                _inf_het = int(np.sum(pop['population_type'][newly_infected_indices] == 1))
                monthly_new_infections_msm[sim_month] = _inf_msm
                monthly_new_infections_het[sim_month] = _inf_het
                year_new_infections_msm += _inf_msm
                year_new_infections_het += _inf_het

            _dead_hiv_indices = month_results['dead_hiv_indices']
            if len(_dead_hiv_indices) > 0:
                _d_msm = int(np.sum(pop['population_type'][_dead_hiv_indices] == 0))
                _d_het = int(np.sum(pop['population_type'][_dead_hiv_indices] == 1))
                monthly_new_deaths_msm[sim_month] = _d_msm
                monthly_new_deaths_het[sim_month] = _d_het
                year_new_deaths_msm += _d_msm
                year_new_deaths_het += _d_het

            _art_indices = month_results.get('art_start_indices', np.array([], dtype=np.int64))
            if len(_art_indices) > 0:
                _a_msm = int(np.sum(pop['population_type'][_art_indices] == 0))
                _a_het = int(np.sum(pop['population_type'][_art_indices] == 1))
                monthly_art_starts_msm[sim_month] = _a_msm
                monthly_art_starts_het[sim_month] = _a_het
                year_art_starts_msm += _a_msm
                year_art_starts_het += _a_het
        
        # ---------------------------------------------------------------

        # ---------------------------------------------------------------
        _f_msm   = int(month_results.get('foreign_inflow_msm', 0))
        _f_het_m = int(month_results.get('foreign_inflow_het_m', 0))
        _f_het_f = int(month_results.get('foreign_inflow_het_f', 0))
        monthly_foreign_inflow_msm[sim_month]   = _f_msm
        monthly_foreign_inflow_het_m[sim_month] = _f_het_m
        monthly_foreign_inflow_het_f[sim_month] = _f_het_f
        year_foreign_inflow_msm   += _f_msm
        year_foreign_inflow_het_m += _f_het_m
        year_foreign_inflow_het_f += _f_het_f

        alive_indices = np.where(pop['alive'])[0]
        if len(alive_indices) > 0:
            age_idx_alive = get_age_group_index_vectorized(pop['age'][alive_indices])
            cd4_idx_alive = get_cd4_category_index_vectorized(pop['cd4'][alive_indices], 
                                                              pop['hiv_positive'][alive_indices])
            keys_alive = encode_aggregator_keys(
                year_idx, age_idx_alive, cd4_idx_alive,
                pop['on_art'][alive_indices], pop['on_prep'][alive_indices]
            )
            agg_person_months += aggregate_counts_vectorized(keys_alive, max_aggregator_key - 1)
        
        if len(newly_infected_indices) > 0:
            age_idx_inf = get_age_group_index_vectorized(pop['age'][newly_infected_indices])
            cd4_idx_inf = get_cd4_category_index_vectorized(pop['cd4'][newly_infected_indices],
                                                            pop['hiv_positive'][newly_infected_indices])
            keys_inf = encode_aggregator_keys(
                year_idx, age_idx_inf, cd4_idx_inf,
                pop['on_art'][newly_infected_indices], pop['on_prep'][newly_infected_indices]
            )
            agg_newly_infected += aggregate_counts_vectorized(keys_inf, max_aggregator_key - 1)
        
        new_diag_indices = month_results['new_diag_indices']
        if len(new_diag_indices) > 0:
            age_idx_diag = get_age_group_index_vectorized(pop['age'][new_diag_indices])
            cd4_idx_diag = get_cd4_category_index_vectorized(pop['cd4'][new_diag_indices],
                                                             pop['hiv_positive'][new_diag_indices])
            keys_diag = encode_aggregator_keys(
                year_idx, age_idx_diag, cd4_idx_diag,
                pop['on_art'][new_diag_indices], pop['on_prep'][new_diag_indices]
            )
            agg_newly_diagnosed += aggregate_counts_vectorized(keys_diag, max_aggregator_key - 1)
        
        dead_hiv_indices = month_results['dead_hiv_indices']
        if len(dead_hiv_indices) > 0:
            age_idx_dead = get_age_group_index_vectorized(pop['age'][dead_hiv_indices])
            cd4_idx_dead = get_cd4_category_index_vectorized(pop['cd4'][dead_hiv_indices],
                                                             pop['hiv_positive'][dead_hiv_indices])
            keys_dead = encode_aggregator_keys(
                year_idx, age_idx_dead, cd4_idx_dead,
                pop['on_art'][dead_hiv_indices], pop['on_prep'][dead_hiv_indices]
            )
            agg_deaths += aggregate_counts_vectorized(keys_dead, max_aggregator_key - 1)

        if sim_month % 12 == 0:
            yearly_infections.append(total_new_infections)
            yearly_new_deaths.append(monthly_deaths)
            yearly_new_deaths_under200.append(monthly_deaths_under200)
            yearly_new_deaths_over200.append(monthly_deaths_over200)
            yearly_new_art_start_under200.append(art_start_under200_this_year)
            yearly_new_art_start_over200.append(art_start_over200_this_year)
            yearly_on_art_current.append(np.sum(pop['on_art'] & pop['hiv_positive'] & pop['alive']))
            yearly_new_diagnoses.append(month_new_diagnoses)
            yearly_tests.append(yearly_test_count)
            yearly_tests_msm.append(yearly_test_count_msm)
            yearly_tests_het.append(yearly_test_count_het)
            yearly_tests_hiv_pos.append(yearly_test_count_hiv_pos)
            yearly_tests_hiv_neg.append(yearly_test_count_hiv_neg)
            
            cumulative_infected += total_new_infections
            yearly_cumulative_infected.append(cumulative_infected)
            
            yearly_population_total.append(np.sum(pop['alive']))
            yearly_population_msm.append(np.sum((pop['population_type'] == 0) & pop['alive']))
            yearly_population_het.append(np.sum((pop['population_type'] == 1) & pop['alive']))
            yearly_population_undiagnosed.append(np.sum((~pop['hiv_positive'] | ~pop['diagnosed']) & pop['alive']))

            year_start_month = sim_month - 11
            year_end_month = sim_month + 1
            yearly_breakthrough_undiagnosed.append(
                int(np.max(monthly_breakthrough_undiagnosed[year_start_month:year_end_month]))
            )

            total_new_infections = 0
            monthly_deaths = 0
            monthly_deaths_under200 = 0
            monthly_deaths_over200 = 0
            art_start_under200_this_year = 0
            art_start_over200_this_year = 0
            month_new_diagnoses = 0
            yearly_test_count = 0
            yearly_test_count_msm = 0
            yearly_test_count_het = 0
            yearly_test_count_hiv_pos = 0
            yearly_test_count_hiv_neg = 0

            if ENABLE_HETEROGENEITY_TRACKING:
                yearly_new_infections_msm.append(year_new_infections_msm)
                yearly_new_infections_het.append(year_new_infections_het)
                yearly_new_deaths_msm.append(year_new_deaths_msm)
                yearly_new_deaths_het.append(year_new_deaths_het)
                yearly_art_starts_msm.append(year_art_starts_msm)
                yearly_art_starts_het.append(year_art_starts_het)
                year_new_infections_msm = 0
                year_new_infections_het = 0
                year_new_deaths_msm = 0
                year_new_deaths_het = 0
                year_art_starts_msm = 0
                year_art_starts_het = 0

            yearly_foreign_inflow_msm.append(year_foreign_inflow_msm)
            yearly_foreign_inflow_het_m.append(year_foreign_inflow_het_m)
            yearly_foreign_inflow_het_f.append(year_foreign_inflow_het_f)
            year_foreign_inflow_msm = 0
            year_foreign_inflow_het_m = 0
            year_foreign_inflow_het_f = 0

    trans_prob_msm = np.where(year_susceptible_msm > 0,
                              year_infections_msm / year_susceptible_msm, 0.0)
    trans_prob_het = np.where(year_susceptible_het > 0,
                              year_infections_het / year_susceptible_het, 0.0)
    trans_prob_onprep = np.where(year_susceptible_onprep > 0,
                                 year_infections_onprep / year_susceptible_onprep, 0.0)
    trans_prob_offprep = np.where(year_susceptible_offprep > 0,
                                  year_infections_offprep / year_susceptible_offprep, 0.0)

    results = {
        'yearly_new_infections': yearly_infections,
        'yearly_new_deaths': yearly_new_deaths,
        'yearly_on_art_current': yearly_on_art_current,
        'yearly_cumulative_infected': yearly_cumulative_infected,
        'yearly_new_diagnoses': yearly_new_diagnoses,
        'monthly_art_start_under200': monthly_art_start_under200,
        'monthly_art_start_over200': monthly_art_start_over200,
        'yearly_art_start_under200': yearly_new_art_start_under200,
        'yearly_art_start_over200': yearly_new_art_start_over200,
        'yearly_new_deaths_under200': yearly_new_deaths_under200,
        'yearly_new_deaths_over200': yearly_new_deaths_over200,
        'year_deaths_by_cd4_state': year_deaths_by_cd4_state,
        'year_person_months_by_cd4_state': year_pmonths_by_cd4_state,
        'year_expected_by_cd4_state': year_expected_by_cd4_state,
        'yearly_tests': yearly_tests,
        'monthly_tests': monthly_tests,
        'yearly_tests_msm': yearly_tests_msm,
        'yearly_tests_het': yearly_tests_het,
        'yearly_tests_hiv_pos': yearly_tests_hiv_pos,
        'yearly_tests_hiv_neg': yearly_tests_hiv_neg,
        'monthly_tests_msm': monthly_tests_msm,
        'monthly_tests_het': monthly_tests_het,
        'monthly_tests_hiv_pos': monthly_tests_hiv_pos,
        'monthly_tests_hiv_neg': monthly_tests_hiv_neg,
        'yearly_population_total': yearly_population_total,
        'yearly_population_msm': yearly_population_msm,
        'yearly_population_het': yearly_population_het,
        'yearly_population_undiagnosed': yearly_population_undiagnosed,
        'yearly_breakthrough_undiagnosed': yearly_breakthrough_undiagnosed,
        'monthly_breakthrough_undiagnosed': monthly_breakthrough_undiagnosed,
        'year_trans_prob_msm': trans_prob_msm,
        'year_trans_prob_het': trans_prob_het,
        'year_trans_prob_onprep': trans_prob_onprep,
        'year_trans_prob_offprep': trans_prob_offprep,

        'monthly_new_infections_msm': monthly_new_infections_msm,
        'monthly_new_infections_het': monthly_new_infections_het,
        'monthly_new_deaths_msm': monthly_new_deaths_msm,
        'monthly_new_deaths_het': monthly_new_deaths_het,
        'monthly_art_starts_msm': monthly_art_starts_msm,
        'monthly_art_starts_het': monthly_art_starts_het,
        'yearly_new_infections_msm': yearly_new_infections_msm,
        'yearly_new_infections_het': yearly_new_infections_het,
        'yearly_new_deaths_msm': yearly_new_deaths_msm,
        'yearly_new_deaths_het': yearly_new_deaths_het,
        'yearly_art_starts_msm': yearly_art_starts_msm,
        'yearly_art_starts_het': yearly_art_starts_het,

        'monthly_foreign_inflow_msm':   monthly_foreign_inflow_msm,
        'monthly_foreign_inflow_het_m': monthly_foreign_inflow_het_m,
        'monthly_foreign_inflow_het_f': monthly_foreign_inflow_het_f,
        'yearly_foreign_inflow_msm':   yearly_foreign_inflow_msm,
        'yearly_foreign_inflow_het_m': yearly_foreign_inflow_het_m,
        'yearly_foreign_inflow_het_f': yearly_foreign_inflow_het_f,
    }

    scenario_label = SCENARIO_CONFIGS[label]["csv_label"]
    rows_list = []

    nonzero_keys = np.where(
        (agg_person_months > 0) | (agg_newly_infected > 0) | 
        (agg_newly_diagnosed > 0) | (agg_deaths > 0)
    )[0]
    
    for key in nonzero_keys:
        yr, age_grp, cd4_cat, on_art_bool, on_prep_bool = decode_aggregator_key(key)
        rows_list.append({
            "Year": yr,
            "AgeGroup": age_grp,
            "CD4Category": cd4_cat,
            "OnART": on_art_bool,
            "OnPrEP": on_prep_bool,
            "PersonMonthsAlive": int(agg_person_months[key]),
            "NewlyInfected": int(agg_newly_infected[key]),
            "NewlyDiagnosed": int(agg_newly_diagnosed[key]),
            "Deaths": int(agg_deaths[key])
        })
    
    aggregator_df = pd.DataFrame(rows_list)
    aggregator_df.sort_values(by=["Year", "AgeGroup", "CD4Category", "OnART", "OnPrEP"], inplace=True)
    aggregator_df.reset_index(drop=True, inplace=True)
    aggregator_df.to_csv(os.path.join(OUTPUT_DIR, f"{scenario_label}_yearly_summary.csv"), index=False)

    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    _cats = TRANSMISSION_SOURCE_CATEGORIES
    _src_rows = []
    for _y in range(years):
        _rg_force_tot = {rg: sum(_src_force[(rg, c)][_y] for c in _cats) for rg in ("msm", "het")}
        _tot_alloc_cat = {c: _src_alloc[("msm", c)][_y] + _src_alloc[("het", c)][_y] for c in _cats}
        _tot_force_cat = {c: _src_force[("msm", c)][_y] + _src_force[("het", c)][_y] for c in _cats}
        _tot_force = sum(_tot_force_cat.values())
        for rg, label_rg in (("msm", "MSM"), ("het", "Heterosexual")):
            for c in _cats:
                f = _src_force[(rg, c)][_y]
                _src_rows.append({
                    "seed": seed, "scenario": scenario_label, "year": _y + 1,
                    "risk_group": label_rg, "source_category": c,
                    "allocated_domestic_infections": _src_alloc[(rg, c)][_y],
                    "source_force": f,
                    "source_force_share": (f / _rg_force_tot[rg]) if _rg_force_tot[rg] > 0 else float("nan"),
                    "transmission_multiplier": _src_mult[_y],
                })
        for c in _cats:
            _src_rows.append({
                "seed": seed, "scenario": scenario_label, "year": _y + 1,
                "risk_group": "Total", "source_category": c,
                "allocated_domestic_infections": _tot_alloc_cat[c],
                "source_force": _tot_force_cat[c],
                "source_force_share": (_tot_force_cat[c] / _tot_force) if _tot_force > 0 else float("nan"),
                "transmission_multiplier": _src_mult[_y],
            })
        _src_rows.append({
            "seed": seed, "scenario": scenario_label, "year": _y + 1,
            "risk_group": "Total", "source_category": "exogenous_foreign_inflow",
            "allocated_domestic_infections": _year_foreign_inflow[_y],
            "source_force": float("nan"), "source_force_share": float("nan"),
            "transmission_multiplier": _src_mult[_y],
        })
    pd.DataFrame(_src_rows).to_csv(
        os.path.join(OUTPUT_DIR, f"{scenario_label}_transmission_source.csv"), index=False)

    _alloc_total = sum(_src_alloc[(rg, c)].sum() for rg in ("msm", "het") for c in _cats)
    _inf_total = float(year_infections_msm.sum() + year_infections_het.sum())
    if not np.isclose(_alloc_total, _inf_total, rtol=1e-6, atol=1e-6):
        print(f"[WARN] transmission-source allocation total {_alloc_total:.3f} != "
              f"domestic infections {_inf_total:.3f} ({scenario_label})")
    for rg, yr_arr in (("msm", year_infections_msm), ("het", year_infections_het)):
        _a = sum(_src_alloc[(rg, c)].sum() for c in _cats)
        if not np.isclose(_a, float(yr_arr.sum()), rtol=1e-6, atol=1e-6):
            print(f"[WARN] {rg} source allocation {_a:.3f} != {rg} infections "
                  f"{float(yr_arr.sum()):.3f} ({scenario_label})")

    export_incident_cohort_outputs(pop, OUTPUT_DIR, scenario_label)
    
    # =========================================================================

    # =========================================================================
    if _globals_backup is not None:
        import sys
        current_module = sys.modules[__name__]
        for key, value in _globals_backup.items():
            setattr(current_module, key, value)
    
    return results, pop

# =============================================================================

# =============================================================================
def export_incident_cohort_outputs(pop, output_dir, scenario_short):
    from pathlib import Path
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------

    # ------------------------------
    inf_year = pop.get('infection_year', None)
    if inf_year is None:
        return

    mask_inf = inf_year >= 1
    if np.any(mask_inf):
        df_inf = pd.DataFrame({
            "cohort_year": inf_year[mask_inf].astype(int),
            "discounted_qaly": pop['discounted_qaly'][mask_inf],
            "discounted_cost": pop.get('discounted_cost', np.zeros_like(pop['discounted_qaly']))[mask_inf],
        })
        g_inf = df_inf.groupby("cohort_year").agg(
            n_incident=("discounted_qaly", "size"),
            total_discounted_qalys=("discounted_qaly", "sum"),
            mean_discounted_qalys=("discounted_qaly", "mean"),
            median_discounted_qalys=("discounted_qaly", "median"),
            total_discounted_costs=("discounted_cost", "sum"),
            mean_discounted_costs=("discounted_cost", "mean"),
        ).reset_index()
    else:
        g_inf = pd.DataFrame(columns=[
            "cohort_year", "n_incident",
            "total_discounted_qalys", "mean_discounted_qalys", "median_discounted_qalys",
            "total_discounted_costs", "mean_discounted_costs"
        ])

    g_inf.to_csv(out_dir / f"{scenario_short}_incident_infection_cohort_summary.csv", index=False)

    # ------------------------------

    # ------------------------------
    diag_year = pop.get('diagnosis_year', None)
    if diag_year is None:
        return

    mask_diag = diag_year >= 1
    if np.any(mask_diag):
        df_diag = pd.DataFrame({
            "cohort_year": diag_year[mask_diag].astype(int),
            "discounted_qaly_post_diagnosis": pop['discounted_qaly_post_diagnosis'][mask_diag],
            "discounted_cost_post_diagnosis": pop.get('discounted_cost_post_diagnosis', np.zeros_like(pop['discounted_qaly_post_diagnosis']))[mask_diag],
        })
        g_diag = df_diag.groupby("cohort_year").agg(
            n_incident=("discounted_qaly_post_diagnosis", "size"),
            total_discounted_qalys=("discounted_qaly_post_diagnosis", "sum"),
            mean_discounted_qalys=("discounted_qaly_post_diagnosis", "mean"),
            median_discounted_qalys=("discounted_qaly_post_diagnosis", "median"),
            total_discounted_costs=("discounted_cost_post_diagnosis", "sum"),
            mean_discounted_costs=("discounted_cost_post_diagnosis", "mean"),
        ).reset_index()
    else:
        g_diag = pd.DataFrame(columns=[
            "cohort_year", "n_incident",
            "total_discounted_qalys", "mean_discounted_qalys", "median_discounted_qalys",
            "total_discounted_costs", "mean_discounted_costs"
        ])

    g_diag.to_csv(out_dir / f"{scenario_short}_incident_diagnosis_cohort_summary.csv", index=False)

# =============================================================================

# =============================================================================

def export_cd4_threshold_outputs(pop, output_dir, sim_months):
    from pathlib import Path
    
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    n = len(pop['cd4'])
    person_id = np.arange(n, dtype=int)

    df_ind = pd.DataFrame({
        'person_id': person_id,
        'art_start_month': pop['art_start'],
        'on_art_end': pop['on_art'].astype(int),
        'alive_end': pop['alive'].astype(int),
        'cd4_at_art_start': pop['cd4_at_art_start'],
        'cd4_nadir': pop['cd4_nadir'],
        't_first_ge200': pop['t_first_ge200'],
        't_first_ge350': pop['t_first_ge350'],
        't_first_ge500': pop['t_first_ge500'],
        't_sustained_ge200_m6': pop['t_sustained_ge200'],
        't_sustained_ge350_m6': pop['t_sustained_ge350'],
        't_sustained_ge500_m6': pop['t_sustained_ge500'],
    })
    if SAVE_CD4_INDIVIDUAL_CSV:
        df_ind.to_csv(out_dir / "cd4_individual_threshold_times.csv", index=False)

    starters = np.where(pop['art_start'] >= 0)[0]
    if len(starters) > 0:
        max_m = int(sim_months)
        months = np.arange(max_m + 1, dtype=int)
        
        def cum_prop(times):
            reached_by_m = np.zeros(max_m + 1, dtype=float)
            valid = times >= 0
            if not np.any(valid):
                return reached_by_m
            t = times[valid].astype(int)
            counts = np.bincount(np.clip(t, 0, max_m), minlength=max_m + 1)
            cum = np.cumsum(counts)
            return cum / float(len(starters))
        
        df_agg = pd.DataFrame({
            'month': months,
            'first_ge200_prop': cum_prop(pop['t_first_ge200'][starters]),
            'first_ge350_prop': cum_prop(pop['t_first_ge350'][starters]),
            'first_ge500_prop': cum_prop(pop['t_first_ge500'][starters]),
            'sust6_ge200_prop': cum_prop(pop['t_sustained_ge200'][starters]),
            'sust6_ge350_prop': cum_prop(pop['t_sustained_ge350'][starters]),
            'sust6_ge500_prop': cum_prop(pop['t_sustained_ge500'][starters]),
        })
        df_agg.to_csv(out_dir / "cd4_aggregated_threshold_curves.csv", index=False)
        
        def prop_at(col, m):
            m = min(int(m), max_m)
            return float(df_agg.loc[df_agg['month'] == m, col].values[0])
        
        def median_time(times):
            times = times[times >= 0]
            return float(np.median(times)) if len(times) > 0 else -1.0
        
        summary = {
            'prop_first_ge350_12m': prop_at('first_ge350_prop', 12),
            'prop_sust6_ge350_12m': prop_at('sust6_ge350_prop', 12),
            'prop_first_ge500_24m': prop_at('first_ge500_prop', 24),
            'prop_sust6_ge500_24m': prop_at('sust6_ge500_prop', 24),
            'median_first_ge350_months': median_time(pop['t_first_ge350'][starters]),
            'median_sust6_ge350_months': median_time(pop['t_sustained_ge350'][starters]),
            'median_first_ge500_months': median_time(pop['t_first_ge500'][starters]),
            'median_sust6_ge500_months': median_time(pop['t_sustained_ge500'][starters]),
            'art_starters_n': int(len(starters)),
        }
        pd.DataFrame([summary]).to_csv(out_dir / "cd4_threshold_summary.csv", index=False)
    else:
        df_agg = pd.DataFrame(columns=['month', 'first_ge200_prop', 'first_ge350_prop',
                                        'first_ge500_prop', 'sust6_ge200_prop',
                                        'sust6_ge350_prop', 'sust6_ge500_prop'])
        df_agg.to_csv(out_dir / "cd4_aggregated_threshold_curves.csv", index=False)
        summary = {
            'prop_first_ge350_12m': 0.0,
            'prop_sust6_ge350_12m': 0.0,
            'prop_first_ge500_24m': 0.0,
            'prop_sust6_ge500_24m': 0.0,
            'median_first_ge350_months': -1.0,
            'median_sust6_ge350_months': -1.0,
            'median_first_ge500_months': -1.0,
            'median_sust6_ge500_months': -1.0,
            'art_starters_n': 0,
        }
        pd.DataFrame([summary]).to_csv(out_dir / "cd4_threshold_summary.csv", index=False)

# =============================================================================

# =============================================================================

def _run_heterogeneity_consistency_checks(all_results, enabled_scenarios, years):
    print("\n[Heterogeneity Consistency Check]")
    print("=" * 70)
    all_ok = True
    for s in enabled_scenarios:
        r = all_results[s]
        inf_msm = np.array(r["yearly_new_infections_msm"], dtype=np.int64)
        inf_het = np.array(r["yearly_new_infections_het"], dtype=np.int64)
        inf_total = np.array(r["yearly_new_infections"], dtype=np.int64)
        inf_sum_40y = int(inf_msm.sum() + inf_het.sum())
        inf_total_40y = int(inf_total.sum())
        inf_ok = np.array_equal(inf_msm + inf_het, inf_total)
        mark = "OK" if inf_ok else "FAIL"
        print(f"  [{s}] Infections: MSM({int(inf_msm.sum())}) + Het({int(inf_het.sum())}) "
              f"= {inf_sum_40y} vs total {inf_total_40y}  [{mark}]")
        if not inf_ok:
            all_ok = False
        
        d_msm = np.array(r["yearly_new_deaths_msm"], dtype=np.int64)
        d_het = np.array(r["yearly_new_deaths_het"], dtype=np.int64)
        d_hivpos_total = (np.array(r["yearly_new_deaths_under200"], dtype=np.int64) +
                          np.array(r["yearly_new_deaths_over200"], dtype=np.int64))
        d_sum_40y = int(d_msm.sum() + d_het.sum())
        d_total_40y = int(d_hivpos_total.sum())
        d_ok = np.array_equal(d_msm + d_het, d_hivpos_total)
        mark = "OK" if d_ok else "FAIL"
        print(f"  [{s}] Deaths(HIV+): MSM({int(d_msm.sum())}) + Het({int(d_het.sum())}) "
              f"= {d_sum_40y} vs HIV+ deaths {d_total_40y}  [{mark}]")
        if not d_ok:
            all_ok = False
        
        a_msm = np.array(r["yearly_art_starts_msm"], dtype=np.int64)
        a_het = np.array(r["yearly_art_starts_het"], dtype=np.int64)
        a_cd4_total = (np.array(r["yearly_art_start_under200"], dtype=np.int64) +
                       np.array(r["yearly_art_start_over200"], dtype=np.int64))
        a_sum_40y = int(a_msm.sum() + a_het.sum())
        a_total_40y = int(a_cd4_total.sum())
        a_ok = np.array_equal(a_msm + a_het, a_cd4_total)
        mark = "OK" if a_ok else "FAIL"
        print(f"  [{s}] ART starts: MSM({int(a_msm.sum())}) + Het({int(a_het.sum())}) "
              f"= {a_sum_40y} vs CD4 total {a_total_40y}  [{mark}]")
        if not a_ok:
            all_ok = False
    
    if all_ok:
        print("  All consistency checks PASSED.")
    else:
        print("  WARNING: One or more consistency checks FAILED. Execution continues.")
    print("=" * 70)
    print()

def export_mortality_parameter_table(output_dir):

    total_rate_meta = {
        "no_art": dict(
            total_estimate_type="observed all-cause mortality rate, approximate",
            total_source_population="ART-naive or pre-ART",
            total_age_adjusted="no, crude rate",
            total_citation="CASCADE/Lewden summaries; Lodwick et al. Lancet 2010",
            total_assumption_flag=True,
        ),
        "on_art_suppressed": dict(
            total_estimate_type="observed all-cause mortality rate",
            total_source_population="on ART, virally suppressed",
            total_age_adjusted="no, crude rate",
            total_citation="COHERE PLoS Med 2012",
            total_assumption_flag=False,
        ),
        "on_art_unsuppressed": dict(
            total_estimate_type="approximate all-cause mortality rate",
            total_source_population="on ART, not suppressed",
            total_age_adjusted="no, crude rate",
            total_citation="CASCADE/Lewden, interpolated",
            total_assumption_flag=True,
        ),
        "disengaged_dropout": dict(
            total_estimate_type="assumption extrapolated from ART-naive",
            total_source_population="disengaged from care after diagnosis",
            total_age_adjusted="no, crude rate",
            total_citation="extrapolated from no_art, not directly observed",
            total_assumption_flag=True,
        ),
    }
    derivation_note_total = (
        "Used directly as a total all-cause rate in literature_total_rate mode. "
        "Do not add background mortality in that mode. For background_plus_excess, "
        "replace with true excess mortality derived as observed total minus expected age-standardised background."
    )
    derivation_note_rr = (
        "Used in multiplier mode by applying the RR/SMR to the age- and sex-specific background mortality hazard. "
        "v2 values are literature anchored where available; split and extrapolated cells remain assumption flagged and should be retained in sensitivity analyses."
    )
    rows = []
    rr_rows = []
    for state_i, state_label in enumerate(HIV_MORT_STATE_LABELS):
        tmeta = total_rate_meta.get(state_label, {})
        for cd4_i, cd4_label in enumerate(CD4_MORT_STRATA_LABELS):
            total_rate = float(HIV_MORTALITY_TOTAL_ANNUAL[state_i, cd4_i])
            true_excess = ("" if HIV_MORTALITY_TRUE_EXCESS_ANNUAL is None
                           else float(HIV_MORTALITY_TRUE_EXCESS_ANNUAL[state_i, cd4_i]))
            rr = float(HIV_MORTALITY_MULTIPLIERS[state_i, cd4_i])
            rrmeta = HIV_MORTALITY_RR_META.get((state_label, cd4_label), {})
            monthly_prob = 1.0 - math.exp(-total_rate / 12.0)
            if HIV_MORTALITY_MODEL == "multiplier":
                active_type, active_val = "RR applied to background mortality hazard", rr
            elif HIV_MORTALITY_MODEL == "literature_total_rate":
                active_type, active_val = "total all-cause rate per PY", total_rate
            elif HIV_MORTALITY_MODEL == "background_plus_excess":
                active_type, active_val = "true excess rate per PY", true_excess
            elif HIV_MORTALITY_MODEL == "combined_max":
                active_type, active_val = "max of RR-hazard and background plus true_excess", ""
            else:
                active_type, active_val = "legacy AIDS x4 and CD4<100 x1.2", ""
            common = {
                "run_label": RUN_LABEL,
                "mortality_model_active": HIV_MORTALITY_MODEL,
                "hiv_state": state_label,
                "cd4_stratum": cd4_label,
                "active_parameter_type": active_type,
                "active_parameter_value": active_val,
                "total_all_cause_rate": total_rate,
                "total_rate_per_1000_py": total_rate * 1000.0,
                "monthly_probability_if_total_rate": monthly_prob,
                "true_excess_rate": true_excess,
                "background_mortality_RR": rr,
                "excess_rate_multiplier_for_sensitivity": HIV_MORTALITY_EXCESS_RATE_MULTIPLIER,
                "multiplier_scaler_for_sensitivity": HIV_MORTALITY_MULTIPLIER_SCALER,
                "age_mortality_offset": AGE_MORTALITY_OFFSET,
                "source_note": HIV_MORTALITY_SOURCE_NOTE,
            }
            row = dict(common)
            row.update({
                "total_rate_estimate_type": tmeta.get("total_estimate_type", ""),
                "total_rate_source_population": tmeta.get("total_source_population", ""),
                "total_rate_age_adjusted": tmeta.get("total_age_adjusted", ""),
                "total_rate_calendar_period": "TBD",
                "total_rate_lower_ci": "TBD",
                "total_rate_upper_ci": "TBD",
                "total_rate_derivation": derivation_note_total,
                "total_rate_assumption_flag": tmeta.get("total_assumption_flag", ""),
                "total_rate_citation": tmeta.get("total_citation", ""),
                "rr_measure": rrmeta.get("measure", ""),
                "rr_reference_group": rrmeta.get("reference_group", ""),
                "rr_age_adjusted": rrmeta.get("age_adjusted", ""),
                "rr_adjusted_covariates": rrmeta.get("adjusted_covariates", ""),
                "rr_calendar_period": rrmeta.get("calendar_period", ""),
                "rr_lower_ci": rrmeta.get("lower_ci", ""),
                "rr_upper_ci": rrmeta.get("upper_ci", ""),
                "rr_source_population": rrmeta.get("source_population", ""),
                "rr_derivation": rrmeta.get("derivation", derivation_note_rr),
                "rr_direct_or_imputed": rrmeta.get("direct_or_imputed", ""),
                "rr_assumption_flag": rrmeta.get("assumption_flag", ""),
                "rr_citation": rrmeta.get("citation", ""),
                "rr_support_note": rrmeta.get("support_note", ""),
                "rr_version": rrmeta.get("rr_version", ""),
                "rr_v2_common_note": rrmeta.get("v2_common_note", ""),
            })
            rows.append(row)
            rr_rows.append({
                "run_label": RUN_LABEL,
                "hiv_state": state_label,
                "cd4_stratum": cd4_label,
                "rr_point_estimate": rr,
                "rr_measure": rrmeta.get("measure", ""),
                "rr_reference_group": rrmeta.get("reference_group", ""),
                "rr_age_adjusted": rrmeta.get("age_adjusted", ""),
                "rr_adjusted_covariates": rrmeta.get("adjusted_covariates", ""),
                "rr_lower_ci": rrmeta.get("lower_ci", ""),
                "rr_upper_ci": rrmeta.get("upper_ci", ""),
                "rr_source_population": rrmeta.get("source_population", ""),
                "rr_derivation": rrmeta.get("derivation", ""),
                "rr_direct_or_imputed": rrmeta.get("direct_or_imputed", ""),
                "rr_assumption_flag": rrmeta.get("assumption_flag", ""),
                "rr_citation": rrmeta.get("citation", ""),
                "rr_support_note": rrmeta.get("support_note", ""),
                "rr_version": rrmeta.get("rr_version", ""),
                "rr_v2_common_note": rrmeta.get("v2_common_note", ""),
                "coauthor_confirmation_status": "pending",
            })
    pd.DataFrame(rows).to_csv(
        os.path.join(output_dir, RUN_LABEL + "_mortality_parameters_cd4_stratified.csv"),
        index=False,
    )
    pd.DataFrame(rr_rows).to_csv(
        os.path.join(output_dir, RUN_LABEL + "_mortality_rr_parameters_cd4_stratified.csv"),
        index=False,
    )
    _cfg = {
        "run_label": RUN_LABEL,
        "seed": SEED,
        "base_output_dir": BASE_OUTPUT_DIR,
        "output_dir": OUTPUT_DIR,
        "mortality_model": HIV_MORTALITY_MODEL,
        "age_offset": AGE_MORTALITY_OFFSET,
        "excess_rate_multiplier": HIV_MORTALITY_EXCESS_RATE_MULTIPLIER,
        "rr_scaler": HIV_MORTALITY_MULTIPLIER_SCALER,
        "total_rate_table": HIV_MORTALITY_TOTAL_ANNUAL.tolist(),
        "true_excess_table": (None if HIV_MORTALITY_TRUE_EXCESS_ANNUAL is None
                              else HIV_MORTALITY_TRUE_EXCESS_ANNUAL.tolist()),
        "rr_table": HIV_MORTALITY_MULTIPLIERS.tolist(),
        "rr_version": HIV_MORTALITY_RR_VERSION,
        "rr_meta": {f"{k[0]}|{k[1]}": v for k, v in HIV_MORTALITY_RR_META.items()},
        "cd4_strata": CD4_MORT_STRATA_LABELS,
        "care_states": HIV_MORT_STATE_LABELS,
        "source_note": HIV_MORTALITY_SOURCE_NOTE,
        "validation": {
            "matrix_shape_expected": [4, 4],
            "monotonicity_checked_across_cd4_decline": True,
            "true_excess_required_for_background_plus_excess_and_combined_max": True,
        },
    }
    with open(os.path.join(output_dir, RUN_LABEL + "_mortality_model_config.json"), "w", encoding="utf-8") as _f:
        json.dump(_cfg, _f, ensure_ascii=False, indent=2)

def save_all_outputs(all_results, all_pops, enabled_scenarios):

    try:
        export_mortality_parameter_table(OUTPUT_DIR)
    except Exception as _e:
        print(f"[hiv_rapid_art_model] mortality parameter table export skipped: {_e}")

    try:
        def _emit(rows, scen, period, d2, pm2, ex2):
            for _ci, _clab in enumerate(CD4_MORT_STRATA_LABELS):
                for _si, _slab in enumerate(HIV_MORT_STATE_LABELS):
                    deaths = float(d2[_si, _ci])
                    py = float(pm2[_si, _ci]) / 12.0
                    expd = float(ex2[_si, _ci])
                    rate = (deaths / py * 1000.0) if py > 0 else float('nan')
                    oe = (deaths / expd) if expd > 0 else float('nan')
                    rows.append({
                        "scenario": scen, "period": period,
                        "cd4_stratum": _clab, "care_state": _slab,
                        "deaths": deaths, "person_years": py,
                        "mortality_rate_per_1000py": rate,
                        "expected_deaths": expd, "OE_ratio": oe,
                    })
        _mort_rows = []
        for _s in enabled_scenarios:
            _r = all_results.get(_s, {})
            _d = _r.get('year_deaths_by_cd4_state')
            _pm = _r.get('year_person_months_by_cd4_state')
            _ex = _r.get('year_expected_by_cd4_state')
            if _d is None or _pm is None or _ex is None:
                continue
            for _yr in range(years):
                _emit(_mort_rows, _s, f"Y{_yr+1}", _d[_yr], _pm[_yr], _ex[_yr])
            for _a in range(0, years, 5):
                _b = min(_a + 5, years)
                _emit(_mort_rows, _s, f"Y{_a+1}-{_b}",
                      _d[_a:_b].sum(axis=0), _pm[_a:_b].sum(axis=0), _ex[_a:_b].sum(axis=0))
            _emit(_mort_rows, _s, f"Y1-{years}",
                  _d.sum(axis=0), _pm.sum(axis=0), _ex.sum(axis=0))
        if _mort_rows:
            pd.DataFrame(_mort_rows).to_csv(
                os.path.join(OUTPUT_DIR, RUN_LABEL + "_mortality_by_cd4_care_state.csv"), index=False)
    except Exception as _e:
        print(f"[hiv_rapid_art_model] cd4/care-state mortality output skipped: {_e}")
    years_arr = np.arange(1, years + 1)
    months_arr = np.arange(1, simulation_months + 1)

    colors = {
        "0 month": "blue",
        "1 month": "green",
        "2 months": "red"
    }
    linestyles = {
        "0 month": "-",
        "1 month": "--",
        "2 months": ":"
    }
    
    def save_figure(fig, filename):
        fig.savefig(os.path.join(OUTPUT_DIR, filename), dpi=300, format='jpeg')
        plt.close(fig)
    
    # =========================================================================

    # =========================================================================
    # yearly_new_infections.csv
    df_data = {"Year": years_arr}
    for s in enabled_scenarios:
        df_data[s] = all_results[s]["yearly_new_infections"]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "yearly_new_infections.csv"), index=False)
    
    # yearly_new_diagnoses.csv
    df_data = {"Year": years_arr}
    for s in enabled_scenarios:
        df_data[s] = all_results[s]["yearly_new_diagnoses"]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "yearly_new_diagnoses.csv"), index=False)
    
    # yearly_new_deaths.csv
    df_data = {"Year": years_arr}
    for s in enabled_scenarios:
        df_data[s] = all_results[s]["yearly_new_deaths"]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "yearly_new_deaths.csv"), index=False)
    
    # yearly_on_art_current.csv
    df_data = {"Year": years_arr}
    for s in enabled_scenarios:
        df_data[s] = all_results[s]["yearly_on_art_current"]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "yearly_on_art_current.csv"), index=False)
    
    # yearly_cumulative_infected.csv
    df_data = {"Year": years_arr}
    for s in enabled_scenarios:
        df_data[s] = all_results[s]["yearly_cumulative_infected"]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "yearly_cumulative_infected.csv"), index=False)
    
    # yearly_tests.csv
    df_data = {"Year": years_arr}
    for s in enabled_scenarios:
        df_data[s] = all_results[s]["yearly_tests"]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "yearly_tests.csv"), index=False)
    
    # yearly_tests_details.csv
    df_data = {"Year": years_arr}
    for s in enabled_scenarios:
        df_data[f"{s}_MSM"] = all_results[s]["yearly_tests_msm"]
        df_data[f"{s}_HET"] = all_results[s]["yearly_tests_het"]
        df_data[f"{s}_HIV+"] = all_results[s]["yearly_tests_hiv_pos"]
        df_data[f"{s}_HIV-"] = all_results[s]["yearly_tests_hiv_neg"]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "yearly_tests_details.csv"), index=False)
    
    # monthly_tests.csv
    df_data = {"Month": months_arr}
    for s in enabled_scenarios:
        df_data[s] = all_results[s]["monthly_tests"][1:simulation_months+1]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "monthly_tests.csv"), index=False)
    
    # monthly_tests_details.csv
    df_data = {"Month": months_arr}
    for s in enabled_scenarios:
        df_data[f"{s}_MSM"] = all_results[s]["monthly_tests_msm"][1:simulation_months+1]
        df_data[f"{s}_HET"] = all_results[s]["monthly_tests_het"][1:simulation_months+1]
        df_data[f"{s}_HIV+"] = all_results[s]["monthly_tests_hiv_pos"][1:simulation_months+1]
        df_data[f"{s}_HIV-"] = all_results[s]["monthly_tests_hiv_neg"][1:simulation_months+1]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "monthly_tests_details.csv"), index=False)

    df_data = {"Month": months_arr}
    for s in enabled_scenarios:
        df_data[f"{s}_over200"] = all_results[s]["monthly_art_start_over200"][1:simulation_months+1]
        df_data[f"{s}_under200"] = all_results[s]["monthly_art_start_under200"][1:simulation_months+1]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "monthly_art_starts.csv"), index=False)

    df_data = {"Year": years_arr}
    for s in enabled_scenarios:
        df_data[f"{s}_over200"] = all_results[s]["yearly_art_start_over200"]
        df_data[f"{s}_under200"] = all_results[s]["yearly_art_start_under200"]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "yearly_new_art_starts.csv"), index=False)
    
    # yearly_new_deaths_by_cd4.csv
    df_data = {"Year": years_arr}
    for s in enabled_scenarios:
        df_data[f"{s} (CD4<200 or AIDS)"] = all_results[s]["yearly_new_deaths_under200"]
        df_data[f"{s} (CD4>=200)"] = all_results[s]["yearly_new_deaths_over200"]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "yearly_new_deaths_by_cd4.csv"), index=False)
    
    # yearly_deaths_among_hiv_positives.csv
    hiv_pos_deaths = {}
    for s in enabled_scenarios:
        hiv_pos_deaths[s] = [all_results[s]["yearly_new_deaths_under200"][i] + 
                             all_results[s]["yearly_new_deaths_over200"][i] for i in range(years)]
    df_data = {"Year": years_arr}
    for s in enabled_scenarios:
        df_data[s] = hiv_pos_deaths[s]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "yearly_deaths_among_hiv_positives.csv"), index=False)
    
    # yearly_population.csv
    df_data = {"Year": years_arr}
    for s in enabled_scenarios:
        df_data[f"{s}_total"] = all_results[s]["yearly_population_total"]
        df_data[f"{s}_msm"] = all_results[s]["yearly_population_msm"]
        df_data[f"{s}_het"] = all_results[s]["yearly_population_het"]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "yearly_population.csv"), index=False)
    
    # yearly_undiagnosed_population.csv
    df_data = {"Year": years_arr}
    for s in enabled_scenarios:
        df_data[s] = all_results[s]["yearly_population_undiagnosed"]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "yearly_undiagnosed_population.csv"), index=False)

    df_data = {"Year": years_arr}
    for s in enabled_scenarios:
        df_data[s] = all_results[s]["yearly_breakthrough_undiagnosed"]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "yearly_breakthrough_undiagnosed.csv"), index=False)

    df_data = {"Month": months_arr}
    for s in enabled_scenarios:
        df_data[s] = all_results[s]["monthly_breakthrough_undiagnosed"][1:simulation_months+1]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "monthly_breakthrough_undiagnosed.csv"), index=False)

    # yearly_transmission_probabilities.csv
    df_data = {"Year": years_arr}
    for s in enabled_scenarios:
        s_label = s.replace(" ", "_")
        df_data[f"MSM_{s_label}"] = all_results[s]["year_trans_prob_msm"]
        df_data[f"Hetero_{s_label}"] = all_results[s]["year_trans_prob_het"]
        df_data[f"OnPrEP_{s_label}"] = all_results[s]["year_trans_prob_onprep"]
        df_data[f"OffPrEP_{s_label}"] = all_results[s]["year_trans_prob_offprep"]
    pd.DataFrame(df_data).to_csv(os.path.join(OUTPUT_DIR, "yearly_transmission_probabilities.csv"), index=False)
    
    # =========================================================================

    # =========================================================================
    if ENABLE_HETEROGENEITY_TRACKING:

        het_wide = {"Year": years_arr}
        for s in enabled_scenarios:
            r = all_results[s]
            het_wide[f"{s}_new_inf_msm"] = r["yearly_new_infections_msm"]
            het_wide[f"{s}_new_inf_het"] = r["yearly_new_infections_het"]
            het_wide[f"{s}_deaths_msm"] = r["yearly_new_deaths_msm"]
            het_wide[f"{s}_deaths_het"] = r["yearly_new_deaths_het"]
            het_wide[f"{s}_art_starts_msm"] = r["yearly_art_starts_msm"]
            het_wide[f"{s}_art_starts_het"] = r["yearly_art_starts_het"]
        pd.DataFrame(het_wide).to_csv(
            os.path.join(OUTPUT_DIR, "yearly_heterogeneity_wide.csv"), index=False
        )

        het_summary_rows = []
        for s in enabled_scenarios:
            r = all_results[s]
            for grp, key_inf, key_death, key_art in [
                ("MSM", "yearly_new_infections_msm", "yearly_new_deaths_msm", "yearly_art_starts_msm"),
                ("Heterosexual", "yearly_new_infections_het", "yearly_new_deaths_het", "yearly_art_starts_het"),
            ]:
                het_summary_rows.append({
                    "Scenario": s,
                    "Risk_Group": grp,
                    "Cum_New_Infections": int(np.sum(r[key_inf])),
                    "Cum_Deaths_HIV_Pos": int(np.sum(r[key_death])),
                    "Cum_ART_Starts": int(np.sum(r[key_art])),
                })
        pd.DataFrame(het_summary_rows).to_csv(
            os.path.join(OUTPUT_DIR, "heterogeneity_summary.csv"), index=False
        )

        _run_heterogeneity_consistency_checks(all_results, enabled_scenarios, years)
    
    # =========================================================================

    # =========================================================================
    inflow_wide = {"Year": years_arr}
    for s in enabled_scenarios:
        r = all_results[s]
        inflow_wide[f"{s}_inflow_msm"]   = r["yearly_foreign_inflow_msm"]
        inflow_wide[f"{s}_inflow_het_m"] = r["yearly_foreign_inflow_het_m"]
        inflow_wide[f"{s}_inflow_het_f"] = r["yearly_foreign_inflow_het_f"]
        inflow_wide[f"{s}_inflow_total"] = (
            np.array(r["yearly_foreign_inflow_msm"])
            + np.array(r["yearly_foreign_inflow_het_m"])
            + np.array(r["yearly_foreign_inflow_het_f"])
        )
    pd.DataFrame(inflow_wide).to_csv(
        os.path.join(OUTPUT_DIR, "foreign_inflow_breakdown.csv"), index=False
    )

    inflow_summary_rows = []
    for s in enabled_scenarios:
        r = all_results[s]
        cum_msm   = int(np.sum(r["yearly_foreign_inflow_msm"]))
        cum_het_m = int(np.sum(r["yearly_foreign_inflow_het_m"]))
        cum_het_f = int(np.sum(r["yearly_foreign_inflow_het_f"]))
        cum_total = cum_msm + cum_het_m + cum_het_f
        inflow_summary_rows.append({
            "Scenario": s,
            "Cum_Inflow_MSM":   cum_msm,
            "Cum_Inflow_HetM":  cum_het_m,
            "Cum_Inflow_HetF":  cum_het_f,
            "Cum_Inflow_Total": cum_total,
            "Pct_MSM":   round(100.0 * cum_msm   / cum_total, 2) if cum_total > 0 else 0.0,
            "Pct_HetM":  round(100.0 * cum_het_m / cum_total, 2) if cum_total > 0 else 0.0,
            "Pct_HetF":  round(100.0 * cum_het_f / cum_total, 2) if cum_total > 0 else 0.0,
        })
    pd.DataFrame(inflow_summary_rows).to_csv(
        os.path.join(OUTPUT_DIR, "foreign_inflow_summary.csv"), index=False
    )
    
    # =========================================================================

    # =========================================================================

    yearly_new_art_starts_total = {}
    for s in enabled_scenarios:
        yearly_new_art_starts_total[s] = [all_results[s]["yearly_art_start_over200"][i] + 
                                          all_results[s]["yearly_art_start_under200"][i] for i in range(years)]
    
    # -------------------------------------------------------------------------
    # Figure 1: summary_plots.jpg
    # -------------------------------------------------------------------------
    fig_all = plt.figure(figsize=(15, 8))
    
    ax1 = fig_all.add_subplot(2, 3, 1)
    for s in enabled_scenarios:
        ax1.plot(years_arr, all_results[s]['yearly_new_infections'], label=s, color=colors[s])
    ax1.set_title("Yearly New Infections")
    ax1.set_xlabel("Year")
    ax1.grid(True)
    ax1.legend()
    
    ax2 = fig_all.add_subplot(2, 3, 2)
    for s in enabled_scenarios:
        ax2.plot(years_arr, all_results[s]['yearly_new_diagnoses'], label=s, color=colors[s])
    ax2.set_title("Yearly New Diagnoses")
    ax2.set_xlabel("Year")
    ax2.grid(True)
    ax2.legend()
    
    ax3 = fig_all.add_subplot(2, 3, 3)
    for s in enabled_scenarios:
        ax3.plot(years_arr, hiv_pos_deaths[s], label=s, color=colors[s])
    ax3.set_title("Yearly Deaths among HIV positives")
    ax3.set_xlabel("Year")
    ax3.grid(True)
    ax3.legend()
    
    ax4 = fig_all.add_subplot(2, 3, 4)
    for s in enabled_scenarios:
        ax4.plot(years_arr, all_results[s]['yearly_on_art_current'], label=s, color=colors[s])
    ax4.set_title("Currently on ART (Alive)")
    ax4.set_xlabel("Year")
    ax4.grid(True)
    ax4.legend()
    
    ax5 = fig_all.add_subplot(2, 3, 5)
    for s in enabled_scenarios:
        ax5.plot(years_arr, all_results[s]['yearly_cumulative_infected'], label=s, color=colors[s])
    ax5.set_title("Yearly Cumulative HIV Infections")
    ax5.set_xlabel("Year")
    ax5.grid(True)
    ax5.legend()

    ax6 = fig_all.add_subplot(2, 3, 6)
    for s in enabled_scenarios:
        ax6.plot(years_arr, yearly_new_art_starts_total[s], label=s, color=colors[s])
    ax6.set_title("Yearly New ART Starts (Total)")
    ax6.set_xlabel("Year")
    ax6.grid(True)
    ax6.legend()
    
    plt.tight_layout()
    save_figure(fig_all, "summary_plots.jpg")
    
    # -------------------------------------------------------------------------
    # Figure 2: yearly_new_infections.jpg
    # -------------------------------------------------------------------------
    fig1, ax = plt.subplots()
    for s in enabled_scenarios:
        ax.plot(years_arr, all_results[s]["yearly_new_infections"], label=s, color=colors[s])
    ax.set_title("Yearly New Infections")
    ax.set_xlabel("Year")
    ax.grid(True)
    ax.legend()
    save_figure(fig1, "yearly_new_infections.jpg")
    
    # -------------------------------------------------------------------------
    # Figure 3: yearly_new_diagnoses.jpg
    # -------------------------------------------------------------------------
    fig2, ax = plt.subplots()
    for s in enabled_scenarios:
        ax.plot(years_arr, all_results[s]["yearly_new_diagnoses"], label=s, color=colors[s])
    ax.set_title("Yearly New Diagnoses")
    ax.set_xlabel("Year")
    ax.grid(True)
    ax.legend()
    save_figure(fig2, "yearly_new_diagnoses.jpg")
    
    # -------------------------------------------------------------------------
    # Figure 4: yearly_new_deaths.jpg
    # -------------------------------------------------------------------------
    fig3, ax = plt.subplots()
    for s in enabled_scenarios:
        ax.plot(years_arr, all_results[s]["yearly_new_deaths"], label=s, color=colors[s])
    ax.set_title("Yearly New Deaths")
    ax.set_xlabel("Year")
    ax.grid(True)
    ax.legend()
    save_figure(fig3, "yearly_new_deaths.jpg")
    
    # -------------------------------------------------------------------------
    # Figure 5: yearly_on_art_current.jpg
    # -------------------------------------------------------------------------
    fig4, ax = plt.subplots()
    for s in enabled_scenarios:
        ax.plot(years_arr, all_results[s]["yearly_on_art_current"], label=s, color=colors[s])
    ax.set_title("Currently on ART (Alive)")
    ax.set_xlabel("Year")
    ax.grid(True)
    ax.legend()
    save_figure(fig4, "yearly_on_art_current.jpg")
    
    # -------------------------------------------------------------------------
    # Figure 6: yearly_cumulative_infected.jpg
    # -------------------------------------------------------------------------
    fig5, ax = plt.subplots()
    for s in enabled_scenarios:
        ax.plot(years_arr, all_results[s]["yearly_cumulative_infected"], label=s, color=colors[s])
    ax.set_title("Yearly Cumulative HIV Infections")
    ax.set_xlabel("Year")
    ax.grid(True)
    ax.legend()
    save_figure(fig5, "yearly_cumulative_infected.jpg")
    
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    fig6, ax = plt.subplots()
    for s in enabled_scenarios:
        monthly_total = (np.array(all_results[s]["monthly_art_start_over200"][1:simulation_months+1]) +
                        np.array(all_results[s]["monthly_art_start_under200"][1:simulation_months+1]))
        ax.plot(months_arr, monthly_total, label=s, color=colors[s])
    ax.set_title("Monthly New ART Starts (Total)")
    ax.set_xlabel("Month")
    ax.grid(True)
    ax.legend()
    save_figure(fig6, "monthly_art_starts.jpg")
    
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    fig7, ax = plt.subplots()
    for s in enabled_scenarios:
        ax.plot(years_arr, yearly_new_art_starts_total[s], label=s, color=colors[s])
    ax.set_title("Yearly New ART Starts (Total)")
    ax.set_xlabel("Year")
    ax.grid(True)
    ax.legend()
    save_figure(fig7, "yearly_new_art_starts.jpg")
    
    # -------------------------------------------------------------------------
    # Figure 9: yearly_new_deaths_by_cd4.jpg
    # -------------------------------------------------------------------------
    fig8, ax = plt.subplots()
    for s in enabled_scenarios:
        ax.plot(years_arr, all_results[s]["yearly_new_deaths_under200"], 
                label=f"{s} (CD4<200 or AIDS)", color=colors[s], linestyle="--")
        ax.plot(years_arr, all_results[s]["yearly_new_deaths_over200"], 
                label=f"{s} (CD4>=200)", color=colors[s], linestyle="-")
    ax.set_title("Yearly New Deaths by CD4 category")
    ax.set_xlabel("Year")
    ax.grid(True)
    ax.legend()
    save_figure(fig8, "yearly_new_deaths_by_cd4.jpg")
    
    # -------------------------------------------------------------------------
    # Figure 10: yearly_deaths_among_hiv_positives.jpg
    # -------------------------------------------------------------------------
    fig9, ax = plt.subplots()
    for s in enabled_scenarios:
        ax.plot(years_arr, hiv_pos_deaths[s], label=s, color=colors[s])
    ax.set_title("Yearly Deaths among HIV positives")
    ax.set_xlabel("Year")
    ax.grid(True)
    ax.legend()
    save_figure(fig9, "yearly_deaths_among_hiv_positives.jpg")
    
    # -------------------------------------------------------------------------
    # Figure 11: yearly_tests.jpg
    # -------------------------------------------------------------------------
    fig10, ax = plt.subplots()
    for s in enabled_scenarios:
        ax.plot(years_arr, all_results[s]["yearly_tests"], label=s, color=colors[s])
    ax.set_title("Yearly HIV Tests")
    ax.set_xlabel("Year")
    ax.set_ylabel("Number of Tests")
    ax.grid(True)
    ax.legend()
    save_figure(fig10, "yearly_tests.jpg")
    
    # -------------------------------------------------------------------------
    # Figure 12: yearly_tests_by_population.jpg
    # -------------------------------------------------------------------------
    fig11, ax = plt.subplots()
    for s in enabled_scenarios:
        ax.plot(years_arr, all_results[s]["yearly_tests_msm"], 
                label=f"MSM ({s})", color=colors[s], linestyle="-")
        ax.plot(years_arr, all_results[s]["yearly_tests_het"], 
                label=f"Hetero ({s})", color=colors[s], linestyle="--")
    ax.set_title("Yearly HIV Tests by Population Group")
    ax.set_xlabel("Year")
    ax.set_ylabel("Number of Tests")
    ax.grid(True)
    ax.legend()
    save_figure(fig11, "yearly_tests_by_population.jpg")
    
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    first_scenario = enabled_scenarios[0]
    fig12, ax = plt.subplots()
    ax.plot(years_arr, all_results[first_scenario]["yearly_tests_hiv_pos"], 
            label=f"HIV+ ({first_scenario})")
    ax.plot(years_arr, all_results[first_scenario]["yearly_tests_hiv_neg"], 
            label=f"HIV- ({first_scenario})")
    ax.set_title(f"Yearly HIV Tests by HIV Status ({first_scenario})")
    ax.set_xlabel("Year")
    ax.set_ylabel("Number of Tests")
    ax.grid(True)
    ax.legend()
    save_figure(fig12, "yearly_tests_by_hiv_status.jpg")
    
    # -------------------------------------------------------------------------
    # Figure 14: yearly_total_population.jpg
    # -------------------------------------------------------------------------
    fig13, ax = plt.subplots(figsize=(10, 6))
    for s in enabled_scenarios:
        ax.plot(years_arr, all_results[s]["yearly_population_total"], label=s, color=colors[s])
    ax.set_title("Total Population Size (High-Risk)")
    ax.set_xlabel("Year")
    ax.set_ylabel("Number of People")
    ax.grid(True)
    ax.legend()
    save_figure(fig13, "yearly_total_population.jpg")
    
    # -------------------------------------------------------------------------
    # Figure 15: yearly_population_by_group.jpg
    # -------------------------------------------------------------------------
    fig14, ax = plt.subplots(figsize=(10, 6))
    for s in enabled_scenarios:
        ax.plot(years_arr, all_results[s]["yearly_population_msm"], 
                label=f"MSM ({s})", color=colors[s], linestyle="-")
        ax.plot(years_arr, all_results[s]["yearly_population_het"], 
                label=f"Hetero ({s})", color=colors[s], linestyle="--")
    ax.set_title("Population by Risk Group")
    ax.set_xlabel("Year")
    ax.set_ylabel("Number of People")
    ax.grid(True)
    ax.legend()
    save_figure(fig14, "yearly_population_by_group.jpg")
    
    # -------------------------------------------------------------------------
    # Figure 16: yearly_undiagnosed_population.jpg
    # -------------------------------------------------------------------------
    fig15, ax = plt.subplots(figsize=(10, 6))
    for s in enabled_scenarios:
        ax.plot(years_arr, all_results[s]["yearly_population_undiagnosed"], label=s, color=colors[s])
    ax.set_title("Undiagnosed Population (HIV- plus undiagnosed HIV+)")
    ax.set_xlabel("Year")
    ax.set_ylabel("Number of People")
    ax.grid(True)
    ax.legend()
    save_figure(fig15, "yearly_undiagnosed_population.jpg")
    
    # -------------------------------------------------------------------------
    # Figure 17: yearly_testing_rate.jpg
    # -------------------------------------------------------------------------
    fig16, ax = plt.subplots()
    for s in enabled_scenarios:
        test_rate = [all_results[s]["yearly_tests"][i]/max(1, all_results[s]["yearly_population_undiagnosed"][i])*100 
                     for i in range(years)]
        ax.plot(years_arr, test_rate, label=s, color=colors[s])
    ax.set_title("Testing Rate (% of Eligible Population)")
    ax.set_xlabel("Year")
    ax.set_ylabel("Percentage")
    ax.grid(True)
    ax.legend()
    save_figure(fig16, "yearly_testing_rate.jpg")
    
    # -------------------------------------------------------------------------
    # Figure 18: population_tracking.jpg
    # -------------------------------------------------------------------------
    fig_pop = plt.figure(figsize=(15, 10))
    
    ax1_pop = fig_pop.add_subplot(2, 2, 1)
    for s in enabled_scenarios:
        ax1_pop.plot(years_arr, all_results[s]["yearly_population_total"], label=s, color=colors[s])
    ax1_pop.set_title("Total Population Size (High-Risk)")
    ax1_pop.set_xlabel("Year")
    ax1_pop.set_ylabel("Number of People")
    ax1_pop.grid(True)
    ax1_pop.legend()
    
    ax2_pop = fig_pop.add_subplot(2, 2, 2)
    for s in enabled_scenarios:
        ax2_pop.plot(years_arr, all_results[s]["yearly_population_msm"], 
                     label=f"MSM ({s})", color=colors[s], linestyle="-")
        ax2_pop.plot(years_arr, all_results[s]["yearly_population_het"], 
                     label=f"Hetero ({s})", color=colors[s], linestyle="--")
    ax2_pop.set_title("Population by Risk Group")
    ax2_pop.set_xlabel("Year")
    ax2_pop.set_ylabel("Number of People")
    ax2_pop.grid(True)
    ax2_pop.legend()
    
    ax3_pop = fig_pop.add_subplot(2, 2, 3)
    for s in enabled_scenarios:
        ax3_pop.plot(years_arr, all_results[s]["yearly_population_undiagnosed"], label=s, color=colors[s])
    ax3_pop.set_title("Undiagnosed Population (HIV- plus undiagnosed HIV+)")
    ax3_pop.set_xlabel("Year")
    ax3_pop.set_ylabel("Number of People")
    ax3_pop.grid(True)
    ax3_pop.legend()
    
    ax4_pop = fig_pop.add_subplot(2, 2, 4)
    first_s = enabled_scenarios[0]
    test_rate_first = [all_results[first_s]["yearly_tests"][i]/max(1, all_results[first_s]["yearly_population_undiagnosed"][i])*100 
                       for i in range(years)]
    ax4_pop.plot(years_arr, all_results[first_s]["yearly_tests"], label=f"Total Tests ({first_s})")
    ax4_pop.plot(years_arr, all_results[first_s]["yearly_population_undiagnosed"], label=f"Undiagnosed Pop ({first_s})")
    ax4_pop.plot(years_arr, test_rate_first, label=f"Testing Rate % ({first_s})", linestyle="--")
    ax4_pop.set_title("Testing Rate vs. Eligible Population")
    ax4_pop.set_xlabel("Year")
    ax4_pop.set_ylabel("Number / Percentage")
    ax4_pop.grid(True)
    ax4_pop.legend()
    
    plt.tight_layout()
    save_figure(fig_pop, "population_tracking.jpg")
    
    # -------------------------------------------------------------------------
    # Figure 19: yearly_transmission_probabilities.jpg
    # -------------------------------------------------------------------------
    fig_tp, ax = plt.subplots(figsize=(12, 7))
    for s in enabled_scenarios:
        ax.plot(years_arr, all_results[s]["year_trans_prob_msm"], 
                label=f"MSM ({s})", color=colors[s], linestyle="-")
        ax.plot(years_arr, all_results[s]["year_trans_prob_het"], 
                label=f"Hetero ({s})", color=colors[s], linestyle="--")
        ax.plot(years_arr, all_results[s]["year_trans_prob_onprep"], 
                label=f"On PrEP ({s})", color=colors[s], linestyle=":")
        ax.plot(years_arr, all_results[s]["year_trans_prob_offprep"], 
                label=f"Off PrEP ({s})", color=colors[s], linestyle="-.")
    ax.set_title("Yearly Transmission Probability by Group")
    ax.set_xlabel("Year")
    ax.set_ylabel("Probability (annual)")
    ax.grid(True)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    save_figure(fig_tp, "yearly_transmission_probabilities.jpg")
    
    # =========================================================================

    # =========================================================================
    try:
        fig_delay, ax_delay = plt.subplots(figsize=(10, 6))
        delay_summary_rows = []
        has_delay_data = False

        def _scenario_short_label(name: str) -> str:
            if name == "0 month":
                return "scenario0"
            elif name == "1 month":
                return "scenario1"
            else:
                return "scenario2"

        for s in enabled_scenarios:
            pop = all_pops[s]
            diag = pop['time_of_diagnosis']

            if 't_art_start_actual' in pop:
                art = pop['t_art_start_actual']
            else:
                art = np.where(pop['art_start'] >= 0, np.ceil(pop['art_start']), -1.0)
            hiv_pos = pop['hiv_positive']

            mask = (
                hiv_pos &
                (diag >= BURN_IN_MONTHS) &
                (art >= 0) &
                (art >= diag)
            )
            if not np.any(mask):
                continue

            has_delay_data = True
            person_ids = np.arange(len(diag), dtype=int)[mask]
            delays = art[mask] - diag[mask]

            scenario_short = _scenario_short_label(s)
            df_ind = pd.DataFrame({
                "person_id": person_ids,
                "diagnosis_month": diag[mask],
                "art_start_month": art[mask],
                "delay_months": delays,
                "scenario": s,
            })
            df_ind.to_csv(os.path.join(OUTPUT_DIR, f"{scenario_short}_diag_to_art_delay_individual.csv"), index=False)

            delay_int = np.floor(delays + 0.5).astype(int)
            vc = pd.Series(delay_int).value_counts().sort_index()
            hist_df = pd.DataFrame({
                "delay_months": vc.index,
                "count": vc.values,
            })
            hist_df["proportion"] = hist_df["count"] / hist_df["count"].sum()
            hist_df.to_csv(os.path.join(OUTPUT_DIR, f"{scenario_short}_diag_to_art_delay_hist.csv"), index=False)

            hist_df_sorted = hist_df.sort_values("delay_months")
            cdf = hist_df_sorted["proportion"].cumsum()
            ax_delay.step(
                hist_df_sorted["delay_months"],
                cdf,
                where="post",
                label=s,
                color=colors.get(s, None)
            )

            mean_delay = float(np.mean(delays))
            median_delay = float(np.median(delays))
            q10, q25, q50, q75, q90, q95 = np.percentile(delays, [10, 25, 50, 75, 90, 95])
            prop_le_1m = float(np.mean(delays <= 1.0))
            prop_le_3m = float(np.mean(delays <= 3.0))
            prop_le_6m = float(np.mean(delays <= 6.0))

            delay_summary_rows.append({
                "scenario": s,
                "n_art_starters": int(len(delays)),
                "mean_delay_months": mean_delay,
                "median_delay_months": median_delay,
                "p10_delay_months": float(q10),
                "p25_delay_months": float(q25),
                "p50_delay_months": float(q50),
                "p75_delay_months": float(q75),
                "p90_delay_months": float(q90),
                "p95_delay_months": float(q95),
                "prop_delay_le_1m": prop_le_1m,
                "prop_delay_le_3m": prop_le_3m,
                "prop_delay_le_6m": prop_le_6m,
            })

        if has_delay_data:
            ax_delay.set_title("Delay from diagnosis to ART start")
            ax_delay.set_xlabel("Delay (months)")
            ax_delay.set_ylabel("Cumulative proportion")
            ax_delay.grid(True)
            ax_delay.set_ylim(0, 1.0)
            ax_delay.legend()
            save_figure(fig_delay, "diagnosis_to_ART_delay_cdf.jpg")

            df_delay_summary = pd.DataFrame(delay_summary_rows)
            df_delay_summary.to_csv(
                os.path.join(OUTPUT_DIR, "diagnosis_to_ART_delay_summary.csv"),
                index=False
            )
        else:
            plt.close(fig_delay)
    except Exception as e:
        print(f"Diagnosis-to-ART delay analysis failed: {e}")

    # =========================================================================

    # =========================================================================
    try:
        fig_delay2, ax_delay2 = plt.subplots(figsize=(10, 6))
        delay_summary_rows2 = []
        has_delay_data2 = False

        def _scenario_short_label2(name: str) -> str:
            if name == "0 month":
                return "scenario0"
            elif name == "1 month":
                return "scenario1"
            else:
                return "scenario2"

        for s in enabled_scenarios:
            pop = all_pops[s]

            if 't_art_eligible' not in pop:
                continue

            elig = pop['t_art_eligible']

            if 't_art_start_actual' in pop:
                art_actual = pop['t_art_start_actual']
            else:
                art_actual = np.where(pop['art_start'] >= 0, np.ceil(pop['art_start']), -1.0)

            hiv_pos = pop['hiv_positive']

            mask = (
                hiv_pos &
                (elig >= BURN_IN_MONTHS) &
                (art_actual >= 0) &
                (art_actual >= elig)
            )
            if not np.any(mask):
                continue

            has_delay_data2 = True
            person_ids = np.arange(len(elig), dtype=int)[mask]
            delays = art_actual[mask] - elig[mask]

            scenario_short = _scenario_short_label2(s)
            df_ind = pd.DataFrame({
                "person_id": person_ids,
                "eligible_month": elig[mask],
                "art_start_month_actual": art_actual[mask],
                "delay_months": delays,
                "scenario": s,
            })
            df_ind.to_csv(os.path.join(OUTPUT_DIR, f"{scenario_short}_eligible_to_art_delay_individual.csv"), index=False)

            delay_int = np.floor(delays + 0.5).astype(int)
            vc = pd.Series(delay_int).value_counts().sort_index()
            hist_df = pd.DataFrame({
                "delay_months": vc.index,
                "count": vc.values,
            })
            hist_df["proportion"] = hist_df["count"] / hist_df["count"].sum()
            hist_df.to_csv(os.path.join(OUTPUT_DIR, f"{scenario_short}_eligible_to_art_delay_hist.csv"), index=False)

            # CDF
            hist_df_sorted = hist_df.sort_values("delay_months")
            cdf = hist_df_sorted["proportion"].cumsum()
            ax_delay2.step(
                hist_df_sorted["delay_months"],
                cdf,
                where="post",
                label=s,
                color=colors.get(s, None)
            )

            mean_delay = float(np.mean(delays))
            median_delay = float(np.median(delays))
            q10, q25, q50, q75, q90, q95 = np.percentile(delays, [10, 25, 50, 75, 90, 95])
            prop_le_1m = float(np.mean(delays <= 1.0))
            prop_le_3m = float(np.mean(delays <= 3.0))
            prop_le_6m = float(np.mean(delays <= 6.0))

            delay_summary_rows2.append({
                "scenario": s,
                "n_art_starters": int(len(delays)),
                "mean_delay_months": mean_delay,
                "median_delay_months": median_delay,
                "p10_delay_months": float(q10),
                "p25_delay_months": float(q25),
                "p50_delay_months": float(q50),
                "p75_delay_months": float(q75),
                "p90_delay_months": float(q90),
                "p95_delay_months": float(q95),
                "prop_delay_le_1m": prop_le_1m,
                "prop_delay_le_3m": prop_le_3m,
                "prop_delay_le_6m": prop_le_6m,
            })

        if has_delay_data2:
            ax_delay2.set_title("Delay from ART eligibility to ART start")
            ax_delay2.set_xlabel("Delay (months)")
            ax_delay2.set_ylabel("Cumulative proportion")
            ax_delay2.grid(True)
            ax_delay2.set_ylim(0, 1.0)
            ax_delay2.legend()
            save_figure(fig_delay2, "eligible_to_ART_delay_cdf.jpg")

            df_delay_summary2 = pd.DataFrame(delay_summary_rows2)
            df_delay_summary2.to_csv(
                os.path.join(OUTPUT_DIR, "eligible_to_ART_delay_summary.csv"),
                index=False
            )
        else:
            plt.close(fig_delay2)
    except Exception as e:
        print(f"Eligible-to-ART delay analysis failed: {e}")

    # =========================================================================

    # =========================================================================
    try:
        risk_summary_rows = []

        for s in enabled_scenarios:
            pop = all_pops[s]
            current_month = BURN_IN_MONTHS + simulation_months

            diag = pop['time_of_diagnosis']
            hiv_pos = pop['hiv_positive']
            alive = pop['alive']
            on_art = pop['on_art']
            diagnosed = pop['diagnosed']
            cd4 = pop['cd4']
            vl = pop['vl']

            time_since_diag = current_month - diag
            diag_valid = diagnosed & (diag >= 0)

            group_mask = (
                hiv_pos & alive & diag_valid & (~on_art) &
                (cd4 >= 500.0) &
                (vl < VL_THRESHOLD_LOG10_5000) &
                (time_since_diag >= LONG_TERM_NO_ART_THRESHOLD_MONTHS)
            )

            group_indices = np.where(group_mask)[0]
            group_size = int(len(group_indices))
            total_hiv_alive = int(np.sum(hiv_pos & alive))

            if group_size == 0 or total_hiv_alive == 0:
                continue

            frac_hiv = group_size / float(total_hiv_alive)

            msm_sus_mask = (~hiv_pos) & alive & (pop['population_type'] == 0)
            het_sus_mask = (~hiv_pos) & alive & (pop['population_type'] == 1)

            S_msm = int(np.sum(msm_sus_mask))
            S_het = int(np.sum(het_sus_mask))

            if S_msm > 0:
                n_prep_msm = int(np.sum(pop['on_prep'][msm_sus_mask]))
                n_nonprep_msm = S_msm - n_prep_msm
                S_eff_msm = n_nonprep_msm + n_prep_msm * (1.0 - prep_effectiveness)
            else:
                S_eff_msm = 0.0

            if S_het > 0:
                het_sus_indices_snapshot = np.where(het_sus_mask)[0]
                het_sus_activity = pop['activity_multiplier'][het_sus_indices_snapshot]
                het_sus_prep_modifier = np.where(
                    pop['on_prep'][het_sus_indices_snapshot],
                    1.0 - prep_effectiveness,
                    1.0
                )
                S_activity_het = float(np.sum(het_sus_activity))
                S_eff_het = float(np.sum(het_sus_activity * het_sus_prep_modifier))
            else:
                S_activity_het = 0.0
                S_eff_het = 0.0

            source_mask_msm = hiv_pos & alive & (~pop['vl_suppressed']) & (pop['population_type'] == 0)
            source_mask_het = hiv_pos & alive & (~pop['vl_suppressed']) & (pop['population_type'] == 1)

            msm_source_indices = np.where(source_mask_msm)[0]
            het_source_indices = np.where(source_mask_het)[0]

            total_group_infections_month = 0.0

            if len(msm_source_indices) > 0 and S_msm > 0 and S_eff_msm > 0:
                vl_bins_msm = get_vl_bin_indices_vectorized(vl[msm_source_indices])
                p_inf_msm_arr = p_infect_msm_bins[vl_bins_msm]
                total_force_msm = float(np.sum(contact_rate_msm * p_inf_msm_arr))
                if TRANSMISSION_MULTIPLIER_POLICY_ON:
                    total_force_msm *= MAIN_TRANSMISSION_MULTIPLIER

                if total_force_msm > 0:
                    group_on_sources_msm = group_mask[msm_source_indices]
                    group_force_msm = float(np.sum(contact_rate_msm * p_inf_msm_arr[group_on_sources_msm]))
                    if TRANSMISSION_MULTIPLIER_POLICY_ON:
                        group_force_msm *= MAIN_TRANSMISSION_MULTIPLIER
                    lambda_msm = total_force_msm / float(S_msm)
                    expected_infections_msm = lambda_msm * float(S_eff_msm)
                    frac_force_msm = group_force_msm / total_force_msm if total_force_msm > 0 else 0.0
                    total_group_infections_month += expected_infections_msm * frac_force_msm

            if len(het_source_indices) > 0 and S_het > 0 and S_eff_het > 0:
                vl_bins_het = get_vl_bin_indices_vectorized(vl[het_source_indices])
                p_inf_het_arr = p_infect_het_bins[vl_bins_het]
                het_source_activity = pop['activity_multiplier'][het_source_indices]
                total_force_het = float(np.sum(contact_rate_het * het_source_activity * p_inf_het_arr))
                if TRANSMISSION_MULTIPLIER_POLICY_ON:
                    total_force_het *= MAIN_TRANSMISSION_MULTIPLIER

                if total_force_het > 0:
                    group_on_sources_het = group_mask[het_source_indices]
                    group_force_het = float(np.sum(
                        contact_rate_het
                        * het_source_activity[group_on_sources_het]
                        * p_inf_het_arr[group_on_sources_het]
                    ))
                    if TRANSMISSION_MULTIPLIER_POLICY_ON:
                        group_force_het *= MAIN_TRANSMISSION_MULTIPLIER
                    lambda_het = total_force_het / max(S_activity_het, 1.0)
                    expected_infections_het = lambda_het * float(S_eff_het)
                    frac_force_het = group_force_het / total_force_het if total_force_het > 0 else 0.0
                    total_group_infections_month += expected_infections_het * frac_force_het

            group_infections_per_year = total_group_infections_month * 12.0
            per_person_transmission_per_year = group_infections_per_year / float(group_size)

            mr_month_arr = calc_mortality_risk_batch_v143(
                cd4[group_indices],
                pop['age'][group_indices],
                pop['is_aids'][group_indices],
                pop['gender'][group_indices],
                pop['on_art'][group_indices],
                pop['vl_suppressed'][group_indices],
                pop['dropped_out'][group_indices],
                pop['diagnosed'][group_indices],
                male_age_upper, male_mortality_rate,
                female_age_upper, female_mortality_rate,
                HIV_MORTALITY_MODE_CODE,
                _ACTIVE_RATE_TABLE,
                HIV_MORTALITY_MULTIPLIERS,
                HIV_MORTALITY_EXCESS_RATE_MULTIPLIER,
                HIV_MORTALITY_MULTIPLIER_SCALER,
            )

            mean_mr_month = float(np.mean(mr_month_arr))

            approx_yearly_mortality = 1.0 - (1.0 - mean_mr_month) ** 12

            risk_summary_rows.append({
                "scenario": s,
                "group_size": group_size,
                "fraction_of_hiv_positive": frac_hiv,
                "mean_cd4": float(cd4[group_indices].mean()),
                "mean_vl_log10": float(vl[group_indices].mean()),
                "mean_age": float(pop['age'][group_indices].mean()),
                "expected_transmissions_per_year_from_group": float(group_infections_per_year),
                "per_person_transmission_risk_per_year": float(per_person_transmission_per_year),
                "per_100_persons_transmissions_per_year": float(per_person_transmission_per_year * 100.0),
                "mean_monthly_mortality_risk": mean_mr_month,
                "approx_yearly_mortality_risk": approx_yearly_mortality,
                "long_term_threshold_months": LONG_TERM_NO_ART_THRESHOLD_MONTHS,
                "vl_threshold_log10": VL_THRESHOLD_LOG10_5000,
            })

        if risk_summary_rows:
            df_risk = pd.DataFrame(risk_summary_rows)
            df_risk.to_csv(
                os.path.join(OUTPUT_DIR, "high_cd4_low_vl_longterm_noART_risk_summary.csv"),
                index=False
            )

            fig_risk, ax1 = plt.subplots(figsize=(10, 6))
            x_pos = np.arange(len(risk_summary_rows))
            group_sizes = [row["group_size"] for row in risk_summary_rows]
            per100_tx = [row["per_100_persons_transmissions_per_year"] for row in risk_summary_rows]
            labels = [row["scenario"] for row in risk_summary_rows]

            ax1.bar(x_pos - 0.2, group_sizes, width=0.4, label="Group size", zorder=3)
            ax2 = ax1.twinx()
            ax2.bar(x_pos + 0.2, per100_tx, width=0.4, label="Transmissions per 100 persons-year", color="orange", zorder=3)

            ax1.set_xticks(x_pos)
            ax1.set_xticklabels(labels)
            ax1.set_xlabel("Scenario")
            ax1.set_ylabel("Number of people (group size)")
            ax2.set_ylabel("Expected transmissions per 100 persons per year")
            ax1.set_title("High CD4 (>=500) & low VL (<5000) long-term no-ART group: size and risk")

            handles1, labels1 = ax1.get_legend_handles_labels()
            handles2, labels2 = ax2.get_legend_handles_labels()
            ax1.legend(handles1 + handles2, labels1 + labels2, loc="upper left", bbox_to_anchor=(0.01, 0.99))

            ax1.grid(True, axis="y", zorder=0)
            ax1.set_axisbelow(True)
            save_figure(fig_risk, "high_cd4_low_vl_longterm_noART_group.jpg")
    except Exception as e:
        print(f"High CD4 & low VL long-term no-ART group analysis failed: {e}")

    # =========================================================================

    # =========================================================================
    try:
        burnin_csv = os.path.join(OUTPUT_DIR, "burnin_summary.csv")
        if os.path.exists(burnin_csv):
            burn_df = pd.read_csv(burnin_csv)
            if "burn_in_year" in burn_df.columns:
                fig_cal, ax_cal = plt.subplots(figsize=(10, 6))
                x_years = burn_df["burn_in_year"].values
                ax_cal.plot(x_years, burn_df["hiv_positive"], label="HIV positive (model)")
                ax_cal.plot(x_years, burn_df["on_art"], label="On ART (model)")
                if "total_alive" in burn_df.columns:
                    ax_cal.plot(x_years, burn_df["total_alive"], label="Total alive (high-risk)", linestyle=":")

                target_hiv_low = TARGET_PLHIV_SOFT_LOW
                target_hiv_high = TARGET_PLHIV_SOFT_HIGH
                target_on_art = TARGET_ART_CURRENT
                ax_cal.axhspan(
                    target_hiv_low, target_hiv_high,
                    color="gray", alpha=0.12,
                    label="Soft PLHIV target (33k-36k)"
                )
                ax_cal.axhline(
                    target_on_art, color="gray", linestyle="-.",
                    label="Hard On ART target (29,680)"
                )

                ax_cal.set_title("Baseline calibration during burn-in")
                ax_cal.set_xlabel("Burn-in year (before main simulation)")
                ax_cal.set_ylabel("Number of people")
                ax_cal.grid(True)
                ax_cal.legend()
                save_figure(fig_cal, "baseline_calibration.jpg")
    except Exception as e:
        print(f"Baseline calibration plot failed: {e}")

    # =========================================================================

    # =========================================================================
    for s in enabled_scenarios:
        if s == "0 month":
            cd4_dir_name = "cd4_0month"
        elif s == "1 month":
            cd4_dir_name = "cd4_1month"
        else:
            cd4_dir_name = "cd4_2months"

        try:
            export_cd4_threshold_outputs(all_pops[s], os.path.join(OUTPUT_DIR, cd4_dir_name), simulation_months)
        except Exception as e:
            print(f"CD4 export failed for {s}: {e}")

    n_scenarios = len(enabled_scenarios)
    print(f"All outputs saved to: {os.path.abspath(OUTPUT_DIR)}")
    print(f"Total figures: 22, Scenarios: {n_scenarios}")
    print(f"Enabled scenarios: {enabled_scenarios}")

# =============================================================================

# =============================================================================

def _get_art_start_actual_month(pop):
    if 't_art_start_actual' in pop:
        return pop['t_art_start_actual']
    art_sched = pop.get('art_start', None)
    if art_sched is None:
        return None
    return np.where(art_sched >= 0, np.ceil(art_sched), -1.0)

def _delay_stats(delays: np.ndarray):
    delays = np.asarray(delays, dtype=float)
    delays = delays[np.isfinite(delays)]
    if delays.size == 0:
        return {
            'n_art_starters': 0,
            'mean_delay_months': float('nan'),
            'median_delay_months': float('nan'),
            'p10_delay_months': float('nan'),
            'p25_delay_months': float('nan'),
            'p50_delay_months': float('nan'),
            'p75_delay_months': float('nan'),
            'p90_delay_months': float('nan'),
            'p95_delay_months': float('nan'),
            'prop_delay_le_1m': float('nan'),
            'prop_delay_le_3m': float('nan'),
            'prop_delay_le_6m': float('nan'),
        }
    q10, q25, q50, q75, q90, q95 = np.percentile(delays, [10, 25, 50, 75, 90, 95])
    return {
        'n_art_starters': int(delays.size),
        'mean_delay_months': float(np.mean(delays)),
        'median_delay_months': float(np.median(delays)),
        'p10_delay_months': float(q10),
        'p25_delay_months': float(q25),
        'p50_delay_months': float(q50),
        'p75_delay_months': float(q75),
        'p90_delay_months': float(q90),
        'p95_delay_months': float(q95),
        'prop_delay_le_1m': float(np.mean(delays <= 1.0)),
        'prop_delay_le_3m': float(np.mean(delays <= 3.0)),
        'prop_delay_le_6m': float(np.mean(delays <= 6.0)),
    }

def compute_art_delay_metrics(pop, burnin_months=BURN_IN_MONTHS):
    art_actual = _get_art_start_actual_month(pop)
    if art_actual is None:
        return {}

    metrics = {}

    diag = pop.get('time_of_diagnosis', None)
    hiv_pos = pop.get('hiv_positive', None)
    if diag is not None and hiv_pos is not None:
        mask = (
            hiv_pos &
            (diag >= burnin_months) &
            (art_actual >= 0) &
            (art_actual >= diag)
        )
        delays = art_actual[mask] - diag[mask]
        metrics['diagnosis_to_art_start'] = _delay_stats(delays)

    elig = pop.get('t_art_eligible', None)
    if elig is not None and hiv_pos is not None:
        mask = (
            hiv_pos &
            (elig >= burnin_months) &
            (art_actual >= 0) &
            (art_actual >= elig)
        )
        delays = art_actual[mask] - elig[mask]
        metrics['eligible_to_art_start'] = _delay_stats(delays)

    return metrics

def run_calibration_simulation(label, initial_pop, seed=SEED, scenario_params=None,
                              delay_mean_override=None, delay_sd_override=None,
                              calibration_years=5):

    np.random.seed(seed)
    rngs = create_event_rngs(seed)
    rng = rngs["init"]
    rng_prep = rngs["prep"]

    _globals_backup = None
    if scenario_params is not None:
        import sys
        current_module = sys.modules[__name__]
        _globals_backup = {}
        for key, value in scenario_params.items():
            if hasattr(current_module, key):
                _globals_backup[key] = getattr(current_module, key)
                setattr(current_module, key, value)

    try:
        pop = copy_population(initial_pop)

        pop['t_first_ge200'][:] = -1
        pop['t_first_ge350'][:] = -1
        pop['t_first_ge500'][:] = -1
        pop['months_at_or_above_200'][:] = 0
        pop['months_at_or_above_350'][:] = 0
        pop['months_at_or_above_500'][:] = 0
        pop['t_sustained_ge200'][:] = -1
        pop['t_sustained_ge350'][:] = -1
        pop['t_sustained_ge500'][:] = -1

        msm_sus_init = (pop['population_type'] == 0) & (~pop['hiv_positive']) & pop['alive']
        msm_sus_init_indices = np.where(msm_sus_init)[0]
        if len(msm_sus_init_indices) > 0 and initial_prep_uptake_msm > 0:
            n_initial_prep_msm = int(len(msm_sus_init_indices) * initial_prep_uptake_msm)
            if n_initial_prep_msm > 0:
                initial_prep_msm_indices = rng_prep.choice(msm_sus_init_indices, n_initial_prep_msm, replace=False)
                pop['on_prep'][initial_prep_msm_indices] = True

        het_sus_init = (pop['population_type'] == 1) & (~pop['hiv_positive']) & pop['alive']
        het_sus_init_indices = np.where(het_sus_init)[0]
        if len(het_sus_init_indices) > 0 and initial_prep_uptake_het > 0:
            n_initial_prep_het = int(len(het_sus_init_indices) * initial_prep_uptake_het)
            if n_initial_prep_het > 0:
                initial_prep_het_indices = rng_prep.choice(het_sus_init_indices, n_initial_prep_het, replace=False)
                pop['on_prep'][initial_prep_het_indices] = True

        cal_months = int(calibration_years * 12)
        for sim_month in range(1, cal_months + 1):
            actual_month = BURN_IN_MONTHS + sim_month
            year_idx = (sim_month - 1) // 12 + 1

            pop, _ = process_single_month(
                pop, actual_month, label, rngs,
                0.0, 0.0,
                record_results=True,
                is_burnin=False,
                sim_year_idx=year_idx,
                delay_mean_override=delay_mean_override,
                delay_sd_override=delay_sd_override,
            )

        metrics = compute_art_delay_metrics(pop, burnin_months=BURN_IN_MONTHS)
        return metrics, pop

    finally:
        if _globals_backup is not None:
            import sys
            current_module = sys.modules[__name__]
            for key, value in _globals_backup.items():
                setattr(current_module, key, value)

def grid_search_art_delay(label, burnin_pop, delay_mean_grid, target_mean_delay,
                          delay_sd_grid=None, seeds=None,
                          calibration_years=5, scenario_params=None,
                          metric_family='eligible_to_art_start'):
    if seeds is None:
        seeds = [SEED, SEED + 1, SEED + 2]
    if delay_sd_grid is None:
        delay_sd_grid = [0.5]

    rows = []
    best = None

    for dm in delay_mean_grid:
        for ds in delay_sd_grid:
            means = []
            ns = []
            for sd in seeds:
                metrics, _ = run_calibration_simulation(
                    label=label,
                    initial_pop=burnin_pop,
                    seed=sd,
                    scenario_params=scenario_params,
                    delay_mean_override=dm,
                    delay_sd_override=ds,
                    calibration_years=calibration_years,
                )

                fam = metrics.get(metric_family, {})
                means.append(fam.get('mean_delay_months', float('nan')))
                ns.append(fam.get('n_art_starters', 0))

            mean_of_means = float(np.nanmean(means))
            n_total = int(np.sum(ns))
            obj = float((mean_of_means - float(target_mean_delay)) ** 2) if np.isfinite(mean_of_means) else float('inf')

            row = {
                'label': label,
                'delay_mean_override': float(dm),
                'delay_sd_override': float(ds),
                'metric_family': metric_family,
                'mean_delay_months_avg': mean_of_means,
                'target_mean_delay': float(target_mean_delay),
                'objective_sse': obj,
                'n_art_starters_total': n_total,
                'seeds': str(list(seeds)),
                'calibration_years': int(calibration_years),
            }
            rows.append(row)

            if best is None or obj < best['objective_sse']:
                best = row

    return {
        'best': best,
        'all_rows': rows,
    }

def export_model_summary_all_scenarios(all_results, all_pops, enabled_scenarios):
    rows = []

    years_arr = np.arange(1, years + 1)

    for scenario in enabled_scenarios:
        result = all_results[scenario]
        pop = all_pops[scenario]

        for i, y in enumerate(years_arr):
            rows.append({
                "scenario": scenario,
                "table_type": "annual",
                "year": int(y),
                "month": np.nan,
                "metric": "new_infections",
                "value": float(result["yearly_new_infections"][i])
            })
            rows.append({
                "scenario": scenario,
                "table_type": "annual",
                "year": int(y),
                "month": np.nan,
                "metric": "new_diagnoses",
                "value": float(result["yearly_new_diagnoses"][i])
            })
            rows.append({
                "scenario": scenario,
                "table_type": "annual",
                "year": int(y),
                "month": np.nan,
                "metric": "new_deaths",
                "value": float(result["yearly_new_deaths"][i])
            })
            rows.append({
                "scenario": scenario,
                "table_type": "annual",
                "year": int(y),
                "month": np.nan,
                "metric": "on_art_current",
                "value": float(result["yearly_on_art_current"][i])
            })
            rows.append({
                "scenario": scenario,
                "table_type": "annual",
                "year": int(y),
                "month": np.nan,
                "metric": "cumulative_infected",
                "value": float(result["yearly_cumulative_infected"][i])
            })
            rows.append({
                "scenario": scenario,
                "table_type": "annual",
                "year": int(y),
                "month": np.nan,
                "metric": "art_starts_over200",
                "value": float(result["yearly_art_start_over200"][i])
            })
            rows.append({
                "scenario": scenario,
                "table_type": "annual",
                "year": int(y),
                "month": np.nan,
                "metric": "art_starts_under200",
                "value": float(result["yearly_art_start_under200"][i])
            })
            rows.append({
                "scenario": scenario,
                "table_type": "annual",
                "year": int(y),
                "month": np.nan,
                "metric": "tests_total",
                "value": float(result["yearly_tests"][i])
            })
            rows.append({
                "scenario": scenario,
                "table_type": "annual",
                "year": int(y),
                "month": np.nan,
                "metric": "population_total",
                "value": float(result["yearly_population_total"][i])
            })
            rows.append({
                "scenario": scenario,
                "table_type": "annual",
                "year": int(y),
                "month": np.nan,
                "metric": "breakthrough_undiagnosed",
                "value": float(result["yearly_breakthrough_undiagnosed"][i])
            })

        deaths_total = int(np.sum(result["yearly_new_deaths"]))
        deaths_under200 = int(np.sum(result["yearly_new_deaths_under200"]))
        deaths_over200 = int(np.sum(result["yearly_new_deaths_over200"]))
        cumulative_final = float(result["yearly_cumulative_infected"][-1]) if len(result["yearly_cumulative_infected"]) > 0 else np.nan
        art_starts_total = int(np.sum(result["yearly_art_start_over200"]) + np.sum(result["yearly_art_start_under200"]))
        tests_total = int(np.sum(result["yearly_tests"]))
        h_pop = int(np.sum(pop['hiv_positive'] & pop['alive']))
        on_art_pop = int(np.sum(pop['on_art'] & pop['hiv_positive'] & pop['alive']))
        diag_pop = int(np.sum(pop['diagnosed'] & pop['hiv_positive'] & pop['alive']))
        undiagnosed = int(np.sum((~pop['diagnosed']) & pop['hiv_positive'] & pop['alive']))
        incident_infection_total = int(np.sum(pop['infection_year'] >= 1))
        incident_diagnosis_total = int(np.sum(pop['diagnosis_year'] >= 1))

        delay_diag = compute_art_delay_metrics(pop).get('diagnosis_to_art_start', {})
        delay_elig = compute_art_delay_metrics(pop).get('eligible_to_art_start', {})

        summary_items = [
            ("summary_total_new_infections", float(np.sum(result["yearly_new_infections"]))),
            ("summary_total_new_diagnoses", float(np.sum(result["yearly_new_diagnoses"]))),
            ("summary_total_deaths", float(deaths_total)),
            ("summary_deaths_cd4_lt200", float(deaths_under200)),
            ("summary_deaths_cd4_ge200", float(deaths_over200)),
            ("summary_cumulative_infected_final", cumulative_final),
            ("summary_art_starts_total", float(art_starts_total)),
            ("summary_tests_total", float(tests_total)),
            ("summary_final_hiv_positive", float(h_pop)),
            ("summary_final_on_art", float(on_art_pop)),
            ("summary_final_diagnosed", float(diag_pop)),
            ("summary_final_undiagnosed", float(undiagnosed)),
            ("summary_incident_infections_20y_41m", float(incident_infection_total)),
            ("summary_incident_diagnoses_20y_41m", float(incident_diagnosis_total)),
            ("delay_diag_n", float(delay_diag.get("n_art_starters", np.nan))),
            ("delay_diag_mean_months", float(delay_diag.get("mean_delay_months", np.nan))),
            ("delay_diag_median_months", float(delay_diag.get("median_delay_months", np.nan))),
            ("delay_eligible_n", float(delay_elig.get("n_art_starters", np.nan))),
            ("delay_eligible_mean_months", float(delay_elig.get("mean_delay_months", np.nan))),
            ("delay_eligible_median_months", float(delay_elig.get("median_delay_months", np.nan))),
        ]

        for metric_name, value in summary_items:
            rows.append({
                "scenario": scenario,
                "table_type": "summary",
                "year": np.nan,
                "month": np.nan,
                "metric": metric_name,
                "value": float(value) if np.isfinite(value) or np.isscalar(value) else np.nan
            })

        monthly_total_art_starts = np.array(result["monthly_art_start_over200"]) + np.array(result["monthly_art_start_under200"])
        for m in range(1, simulation_months + 1):
            rows.append({
                "scenario": scenario,
                "table_type": "monthly",
                "year": int(((m - 1) // 12) + 1),
                "month": int(m),
                "metric": "monthly_art_starts_total",
                "value": float(monthly_total_art_starts[m])
            })
            rows.append({
                "scenario": scenario,
                "table_type": "monthly",
                "year": int(((m - 1) // 12) + 1),
                "month": int(m),
                "metric": "monthly_tests",
                "value": float(result["monthly_tests"][m])
            })

    summary_df = pd.DataFrame(rows)
    return summary_df

def export_report_ready_summary(all_results, all_pops, enabled_scenarios, output_dir=None):
    if output_dir is None:
        output_dir = OUTPUT_DIR
    import json

    # -------------------------------------------------------------

    # -------------------------------------------------------------
    summary_rows = []
    for scenario in enabled_scenarios:
        result = all_results[scenario]
        pop = all_pops[scenario]

        delay_diag = compute_art_delay_metrics(pop).get('diagnosis_to_art_start', {})
        delay_elig = compute_art_delay_metrics(pop).get('eligible_to_art_start', {})

        summary_rows.append({
            "scenario": scenario,
            "total_new_infections_40y": float(np.sum(result["yearly_new_infections"])),
            "total_new_diagnoses_40y": float(np.sum(result["yearly_new_diagnoses"])),
            "total_new_deaths_40y": float(np.sum(result["yearly_new_deaths"])),
            "deaths_cd4_lt200_40y": float(np.sum(result["yearly_new_deaths_under200"])),
            "deaths_cd4_ge200_40y": float(np.sum(result["yearly_new_deaths_over200"])),
            "cumulative_infected_final": float(result["yearly_cumulative_infected"][-1]),
            "total_art_starts_40y": float(np.sum(result["yearly_art_start_over200"]) + np.sum(result["yearly_art_start_under200"])),
            "total_tests_40y": float(np.sum(result["yearly_tests"])),
            "hiv_positive_end": float(np.sum(pop['hiv_positive'] & pop['alive'])),
            "on_art_end": float(np.sum(pop['on_art'] & pop['hiv_positive'] & pop['alive'])),
            "diagnosed_end": float(np.sum(pop['diagnosed'] & pop['hiv_positive'] & pop['alive'])),
            "undiagnosed_end": float(np.sum((~pop['diagnosed']) & pop['hiv_positive'] & pop['alive'])),
            "incident_infections_after_burnin": float(np.sum(pop['infection_year'] >= 1)),
            "incident_diagnoses_after_burnin": float(np.sum(pop['diagnosis_year'] >= 1)),
            "delay_diag_n": float(delay_diag.get("n_art_starters", np.nan)),
            "delay_diag_mean_months": float(delay_diag.get("mean_delay_months", np.nan)),
            "delay_diag_median_months": float(delay_diag.get("median_delay_months", np.nan)),
            "delay_diag_p50_months": float(delay_diag.get("p50_delay_months", np.nan)),
            "delay_eligible_n": float(delay_elig.get("n_art_starters", np.nan)),
            "delay_eligible_mean_months": float(delay_elig.get("mean_delay_months", np.nan)),
            "delay_eligible_median_months": float(delay_elig.get("median_delay_months", np.nan)),
            "delay_eligible_p50_months": float(delay_elig.get("p50_delay_months", np.nan)),

            "total_new_infections_msm_40y": float(np.sum(result.get("yearly_new_infections_msm", [0]))),
            "total_new_infections_het_40y": float(np.sum(result.get("yearly_new_infections_het", [0]))),
            "total_new_deaths_msm_40y": float(np.sum(result.get("yearly_new_deaths_msm", [0]))),
            "total_new_deaths_het_40y": float(np.sum(result.get("yearly_new_deaths_het", [0]))),
            "total_art_starts_msm_40y": float(np.sum(result.get("yearly_art_starts_msm", [0]))),
            "total_art_starts_het_40y": float(np.sum(result.get("yearly_art_starts_het", [0]))),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_records = summary_df.replace({np.nan: None}).to_dict(orient="records")

    # -------------------------------------------------------------

    # -------------------------------------------------------------
    years_arr = np.arange(1, years + 1)
    annual_rows = {"Year": years_arr}

    metric_map = [
        ("new_infections", "yearly_new_infections"),
        ("new_diagnoses", "yearly_new_diagnoses"),
        ("new_deaths", "yearly_new_deaths"),
        ("on_art_current", "yearly_on_art_current"),
        ("cumulative_infected", "yearly_cumulative_infected"),
        ("tests", "yearly_tests"),
        ("art_start_over200", "yearly_art_start_over200"),
        ("art_start_under200", "yearly_art_start_under200"),
        ("breakthrough_undiagnosed", "yearly_breakthrough_undiagnosed"),
        # Heterogeneity（MSM vs Heterosexual）
        ("new_infections_msm", "yearly_new_infections_msm"),
        ("new_infections_het", "yearly_new_infections_het"),
        ("new_deaths_msm", "yearly_new_deaths_msm"),
        ("new_deaths_het", "yearly_new_deaths_het"),
        ("art_starts_msm", "yearly_art_starts_msm"),
        ("art_starts_het", "yearly_art_starts_het"),
    ]

    for scenario in enabled_scenarios:
        result = all_results[scenario]
        for metric_label, key in metric_map:
            annual_rows[f"{scenario}_{metric_label}"] = result[key]

        annual_rows[f"{scenario}_deaths_cd4_lt200"] = result["yearly_new_deaths_under200"]
        annual_rows[f"{scenario}_deaths_cd4_ge200"] = result["yearly_new_deaths_over200"]

    annual_df = pd.DataFrame(annual_rows)
    annual_records = annual_df.replace({np.nan: None}).to_dict(orient="records")

    # -------------------------------------------------------------

    # -------------------------------------------------------------
    months_arr = np.arange(1, simulation_months + 1)
    monthly_rows = {"Month": months_arr}
    for scenario in enabled_scenarios:
        result = all_results[scenario]
        monthly_rows[f"{scenario}_monthly_art_starts_total"] = (
            np.array(result["monthly_art_start_over200"]) + np.array(result["monthly_art_start_under200"])
        )[1:simulation_months + 1]
        monthly_rows[f"{scenario}_monthly_tests"] = np.array(result["monthly_tests"])[1:simulation_months + 1]

    monthly_df = pd.DataFrame(monthly_rows)
    monthly_records = monthly_df.replace({np.nan: None}).to_dict(orient="records")

    # -------------------------------------------------------------

    # -------------------------------------------------------------
    scenario_long = export_model_summary_all_scenarios(all_results, all_pops, enabled_scenarios)
    scenario_long_records = scenario_long.replace({np.nan: None}).to_dict(orient="records")

    # -------------------------------------------------------------

    # -------------------------------------------------------------
    report_json = {
        "metadata": {
            "script": f"{RUN_LABEL}.py",
            "burn_in_years": BURN_IN_YEARS,
            "simulation_years": years,
            "simulated_months": simulation_months,
            "scenarios": enabled_scenarios,
            "art_cd4_wait_policy": {
                "enabled": ART_CD4_WAIT_POLICY_ENABLED,
                "burnin_start_year": BURNIN_START_YEAR,
                "start_year": ART_CD4_WAIT_POLICY_START_YEAR,
                "end_year": ART_CD4_WAIT_POLICY_END_YEAR,
                "cd4_threshold": ART_CD4_WAIT_THRESHOLD,
                "periods": [
                    {
                        "start_year": int(start_year),
                        "end_year": int(end_year),
                        "cd4_threshold": float(threshold),
                    }
                    for start_year, end_year, threshold in ART_CD4_WAIT_POLICY_PERIODS
                ],
            },
            "burnin_test_capacity_policy": {
                "enabled": BURNIN_TEST_CAPACITY_POLICY_ON,
                "source": BURNIN_TEST_CAPACITY_SOURCE,
                "multiplier": BURNIN_TEST_CAPACITY_MULTIPLIER,
                "monthly_capacity_by_year": BURNIN_TEST_CAPACITY_MONTHLY_BY_YEAR,
                "annual_admin_tests": BURNIN_ADMIN_TESTS_ANNUAL,
                "undiagnosed_hiv_test_boost_by_year": BURNIN_UNDIAGNOSED_HIV_TEST_BOOST_BY_YEAR,
            },
            "transmission_multiplier_policy": {
                "enabled": TRANSMISSION_MULTIPLIER_POLICY_ON,
                "burnin_multiplier_by_year": BURNIN_TRANSMISSION_MULTIPLIER_BY_YEAR,
                "main_multiplier": MAIN_TRANSMISSION_MULTIPLIER,
            },
            "created_at": datetime.datetime.now().isoformat()
        },
        "summary_by_scenario": summary_records,
        "annual": annual_records,
        "monthly": monthly_records,
        "scenario_long": scenario_long_records
    }

    json_path = os.path.join(output_dir, f"{RUN_LABEL}_report_ready.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(report_json, f, ensure_ascii=False, indent=2)

    print(f"Saved report-ready JSON output: {os.path.abspath(json_path)}")

# =============================================================================
# Markdown Summary Report (English)
# =============================================================================
def export_markdown_report(all_results, all_pops, enabled_scenarios, elapsed=None, output_dir=None):
    """
    Generate a human-readable Markdown summary report in English.
    Covers: simulation configuration, per-scenario outcomes, cross-scenario
    comparisons, ART delay metrics, and brief interpretation notes.
    """
    if output_dir is None:
        output_dir = OUTPUT_DIR

    lines = []

    # ------------------------------------------------------------------
    # Title & timestamp
    # ------------------------------------------------------------------
    lines.append(f"# {RUN_LABEL} -- Simulation Summary Report")
    lines.append("")
    lines.append(f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    if elapsed is not None:
        lines.append(f"Wall-clock runtime: {elapsed}")
    lines.append("")

    # ------------------------------------------------------------------
    # 1. Configuration Parameters
    # ------------------------------------------------------------------
    lines.append("## 1. Configuration Parameters")
    lines.append("")
    lines.append("| Parameter | Value |")
    lines.append("|---|---|")
    art_cd4_policy_desc = "; ".join(
        f"{start_year}-{end_year}: CD4 <= {threshold:.0f}"
        for start_year, end_year, threshold in ART_CD4_WAIT_POLICY_PERIODS
    )
    undiag_boost_desc = ", ".join(
        f"{year}:x{boost:.2f}"
        for year, boost in BURNIN_UNDIAGNOSED_HIV_TEST_BOOST_BY_YEAR.items()
    )
    transmission_multiplier_desc = (
        f"Burn-in 2009-2023 x{BURNIN_TRANSMISSION_MULTIPLIER_BY_YEAR[2009]:.2f}-"
        f"{BURNIN_TRANSMISSION_MULTIPLIER_BY_YEAR[2023]:.2f}; "
        f"main x{MAIN_TRANSMISSION_MULTIPLIER:.2f}"
    )
    lines.append(f"| Run label | `{RUN_LABEL}` |")
    lines.append(f"| Random seed | {SEED} |")
    lines.append(f"| Burn-in period | {BURN_IN_YEARS} years ({BURN_IN_MONTHS} months) |")
    lines.append(f"| Burn-in ART CD4 wait policy | {art_cd4_policy_desc} |")
    lines.append(f"| Burn-in test capacity policy | {'Calibrated yearly monthly capacities' if BURNIN_TEST_CAPACITY_POLICY_ON else 'Fixed X'} |")
    lines.append(f"| Burn-in undiagnosed HIV test boost | {undiag_boost_desc} |")
    lines.append(f"| Transmission multiplier | {transmission_multiplier_desc} |")
    lines.append(f"| Main simulation | {years} years ({simulation_months} months) |")
    lines.append(f"| Initial MSM population | {msm_population_size:,} |")
    lines.append(f"| Initial heterosexual population | {heterosexual_population_size:,} |")
    lines.append(f"| Target total population | {TARGET_POPULATION:,} |")
    lines.append(f"| Newcomer rate | {NEWCOMER_RATE:.1%} / year |")
    lines.append(f"| Initial total HIV+ | {initial_total_hiv:,} |")
    lines.append(f"| Initial on ART | {initial_on_art:,} |")
    lines.append(f"| Calibration baseline | {CALIBRATION_BASELINE} |")
    lines.append(f"| Target annual HIV reports | {TARGET_HIV_REPORTS_ANNUAL:,} |")
    lines.append(f"| Target annual AIDS reports | {TARGET_AIDS_REPORTS_ANNUAL:,} |")
    lines.append(f"| Target annual HIV+AIDS reports | {TARGET_TOTAL_REPORTS_ANNUAL:,} |")
    lines.append(f"| Target route mix (MSM/Hetero) | {TARGET_ROUTE_MSM_SHARE:.0%}/{TARGET_ROUTE_HET_SHARE:.0%} |")
    lines.append(f"| Target current ART | {TARGET_ART_CURRENT:,} |")
    lines.append(f"| Target undiagnosed HIV | {TARGET_UNDIAGNOSED_HIV:,} ({TARGET_UNDIAGNOSED_HIV_LOW:,}-{TARGET_UNDIAGNOSED_HIV_HIGH:,}) |")
    lines.append(f"| Soft PLHIV target | {TARGET_PLHIV_SOFT_LOW:,}-{TARGET_PLHIV_SOFT_HIGH:,} |")
    lines.append(f"| Initial HIV+ route mix (MSM/HetM/HetF) | {INITIAL_HIV_MSM_RATIO:.0%}/{INITIAL_HIV_HETM_RATIO:.0%}/{INITIAL_HIV_HETF_RATIO:.0%} |")
    lines.append(f"| Contact rate (MSM/Het) | {contact_rate_msm:.1f}/{contact_rate_het:.1f} |")
    lines.append(f"| Main regular-test capacity X/Y/Z | {TEST_CAPACITY_X:,}/{TEST_CAPACITY_Y:,}/{TEST_CAPACITY_Z} |")
    lines.append(f"| Hetero core share | {HETERO_CORE_SHARE:.1%} |")
    lines.append(f"| Hetero core high-to-low activity ratio | {HETERO_CORE_ACTIVITY_MULTIPLIER:.1f} |")
    lines.append(f"| Hetero core initial HIV+ weight | {HETERO_CORE_INITIAL_HIV_WEIGHT:.1f} |")
    lines.append(f"| Hetero core foreign HIV+ share | {HETERO_CORE_FOREIGN_SHARE:.1%} |")
    lines.append(f"| PrEP effectiveness | {prep_effectiveness:.0%} |")
    lines.append(f"| PrEP max uptake (MSM) | {max_prep_uptake_msm:.1%} |")
    lines.append(f"| PrEP max uptake (Het) | {max_prep_uptake_het:.1%} |")
    lines.append(f"| PrEP dropout rate (monthly) | {prep_dropout_rate_msm:.3f} |")
    lines.append(f"| Discount rate (QALY) | {DISCOUNT_RATE_QALY:.0%} |")
    lines.append(f"| Pre-ART dropout (monthly) | {dropout_prob_pre_art:.3f} |")
    lines.append(f"| On-ART dropout (monthly) | {dropout_prob_on_art:.4f} |")
    lines.append(f"| Grid search enabled | {GRID_SEARCH_ENABLED} |")
    lines.append("")

    # Scenario definitions
    lines.append("### Scenario Definitions")
    lines.append("")
    lines.append("| Scenario | ART delay mean (months) | ART delay SD | Enabled |")
    lines.append("|---|---|---|---|")
    for name, cfg in SCENARIO_CONFIGS.items():
        enabled_flag = name in enabled_scenarios
        lines.append(
            f"| {name} | {cfg['delay_mean']:.1f} | {cfg['delay_sd']:.1f} | {'Yes' if enabled_flag else 'No'} |"
        )
    lines.append("")

    # ------------------------------------------------------------------
    # 2. Per-Scenario Outcome Summary
    # ------------------------------------------------------------------
    lines.append("## 2. Per-Scenario Outcome Summary (40-year main simulation)")
    lines.append("")

    # Build a table with key metrics
    header_cols = ["Metric"] + enabled_scenarios
    lines.append("| " + " | ".join(header_cols) + " |")
    lines.append("|" + "|".join(["---"] * len(header_cols)) + "|")

    # Collect data per scenario
    metric_defs = [
        ("New infections (total)", lambda r, p: float(np.sum(r["yearly_new_infections"]))),
        ("New diagnoses (total)", lambda r, p: float(np.sum(r["yearly_new_diagnoses"]))),
        ("Deaths (total)", lambda r, p: float(np.sum(r["yearly_new_deaths"]))),
        ("Deaths (CD4 < 200)", lambda r, p: float(np.sum(r["yearly_new_deaths_under200"]))),
        ("Deaths (CD4 >= 200)", lambda r, p: float(np.sum(r["yearly_new_deaths_over200"]))),
        ("ART starts (total)", lambda r, p: float(np.sum(r["yearly_art_start_over200"]) + np.sum(r["yearly_art_start_under200"]))),
        ("ART starts (CD4 >= 200)", lambda r, p: float(np.sum(r["yearly_art_start_over200"]))),
        ("ART starts (CD4 < 200)", lambda r, p: float(np.sum(r["yearly_art_start_under200"]))),
        ("Tests (total)", lambda r, p: float(np.sum(r["yearly_tests"]))),
        ("HIV+ alive (end)", lambda r, p: float(np.sum(p['hiv_positive'] & p['alive']))),
        ("On ART (end)", lambda r, p: float(np.sum(p['on_art'] & p['hiv_positive'] & p['alive']))),
        ("Diagnosed (end)", lambda r, p: float(np.sum(p['diagnosed'] & p['hiv_positive'] & p['alive']))),
        ("Undiagnosed (end)", lambda r, p: float(np.sum((~p['diagnosed']) & p['hiv_positive'] & p['alive']))),
        ("Incident infections (post-burnin)", lambda r, p: float(np.sum(p['infection_year'] >= 1))),
        ("Incident diagnoses (post-burnin)", lambda r, p: float(np.sum(p['diagnosis_year'] >= 1))),

        ("New infections (MSM)", lambda r, p: float(np.sum(r.get("yearly_new_infections_msm", [0])))),
        ("New infections (heterosexual)", lambda r, p: float(np.sum(r.get("yearly_new_infections_het", [0])))),
        ("Deaths (MSM, HIV+)", lambda r, p: float(np.sum(r.get("yearly_new_deaths_msm", [0])))),
        ("Deaths (heterosexual, HIV+)", lambda r, p: float(np.sum(r.get("yearly_new_deaths_het", [0])))),
        ("ART starts (MSM)", lambda r, p: float(np.sum(r.get("yearly_art_starts_msm", [0])))),
        ("ART starts (heterosexual)", lambda r, p: float(np.sum(r.get("yearly_art_starts_het", [0])))),
    ]

    scenario_vals = {}  # {metric_label: {scenario: value}}
    for label, fn in metric_defs:
        row_vals = []
        scenario_vals[label] = {}
        for sc in enabled_scenarios:
            v = fn(all_results[sc], all_pops[sc])
            scenario_vals[label][sc] = v
            row_vals.append(f"{v:,.0f}")
        lines.append("| " + label + " | " + " | ".join(row_vals) + " |")
    lines.append("")

    # ------------------------------------------------------------------
    # 3. Cross-Scenario Comparison (if >= 2 scenarios)
    # ------------------------------------------------------------------
    if len(enabled_scenarios) >= 2:
        lines.append("## 3. Cross-Scenario Comparison")
        lines.append("")

        ref = enabled_scenarios[0]
        for comp in enabled_scenarios[1:]:
            lines.append(f"### {comp} vs. {ref}")
            lines.append("")
            lines.append("| Metric | " + ref + " | " + comp + " | Difference | Relative change |")
            lines.append("|---|---|---|---|---|")

            compare_metrics = [
                "New infections (total)",
                "New diagnoses (total)",
                "Deaths (total)",
                "Deaths (CD4 < 200)",
                "ART starts (total)",
                "Tests (total)",
                "Undiagnosed (end)",
            ]
            for m in compare_metrics:
                v_ref = scenario_vals[m][ref]
                v_comp = scenario_vals[m][comp]
                diff = v_comp - v_ref
                rel = (diff / v_ref * 100) if v_ref != 0 else float('nan')
                sign = "+" if diff > 0 else ""
                lines.append(
                    f"| {m} | {v_ref:,.0f} | {v_comp:,.0f} | {sign}{diff:,.0f} | {sign}{rel:.2f}% |"
                )
            lines.append("")

    # ------------------------------------------------------------------
    # 4. ART Delay Metrics
    # ------------------------------------------------------------------
    lines.append("## 4. ART Delay Metrics (Diagnosis to ART Start)")
    lines.append("")
    lines.append("| Scenario | N (ART starters) | Mean delay (months) | Median delay (months) |")
    lines.append("|---|---|---|---|")
    for sc in enabled_scenarios:
        pop = all_pops[sc]
        delay_info = compute_art_delay_metrics(pop).get('diagnosis_to_art_start', {})
        n = delay_info.get("n_art_starters", float('nan'))
        mean_d = delay_info.get("mean_delay_months", float('nan'))
        med_d = delay_info.get("median_delay_months", float('nan'))
        lines.append(f"| {sc} | {n:,.0f} | {mean_d:.2f} | {med_d:.2f} |")
    lines.append("")

    # ------------------------------------------------------------------
    # 5. Annual Trends Snapshot (first 5 and last 5 years)
    # ------------------------------------------------------------------
    lines.append("## 5. Annual New Infections Trend (Years 1-5 and 36-40)")
    lines.append("")
    header = ["Year"] + enabled_scenarios
    lines.append("| " + " | ".join(header) + " |")
    lines.append("|" + "|".join(["---"] * len(header)) + "|")
    snapshot_years = list(range(5)) + list(range(max(0, years - 5), years))
    snapshot_years = sorted(set(snapshot_years))
    for yi in snapshot_years:
        yr_label = str(yi + 1)
        vals = []
        for sc in enabled_scenarios:
            v = float(all_results[sc]["yearly_new_infections"][yi])
            vals.append(f"{v:,.0f}")
        lines.append("| " + yr_label + " | " + " | ".join(vals) + " |")
    lines.append("")

    # ------------------------------------------------------------------
    # 6. Interpretation Notes
    # ------------------------------------------------------------------
    lines.append("## 6. Interpretation Notes")
    lines.append("")
    lines.append("- All counts are from the **post-burn-in main simulation period** "
                 f"({years} years / {simulation_months} months).")
    lines.append("- The burn-in phase is shared across all scenarios to ensure identical "
                 "starting conditions.")
    lines.append("- Common Random Numbers (CRN) approximation is used so that stochastic "
                 "variation across scenarios is minimized; observed differences are primarily "
                 "attributable to ART delay timing.")
    if len(enabled_scenarios) >= 2:
        ref = enabled_scenarios[0]
        lines.append(f"- Cross-scenario comparisons use **{ref}** (immediate ART) as the reference.")
        lines.append("- A positive relative change in infections or deaths indicates worse "
                     "outcomes under delayed ART initiation.")
    lines.append("- QALY discounting uses a " + f"{DISCOUNT_RATE_QALY:.0%}" +
                 " annual rate applied monthly from the point of infection or diagnosis.")
    lines.append(f"- PrEP effectiveness is set at {prep_effectiveness:.0%} with gradual uptake "
                 f"over {prep_growth_years_msm} years (MSM) / {prep_growth_years_het} years (heterosexual).")
    lines.append("")
    lines.append("---")
    lines.append(f"*Report generated by `{RUN_LABEL}.py`*")

    # ------------------------------------------------------------------
    # Write to file
    # ------------------------------------------------------------------
    md_path = os.path.join(output_dir, f"{RUN_LABEL}_summary_report.md")
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Saved Markdown summary report: {os.path.abspath(md_path)}")

def _run_scenario_worker(args):
    scenario, burnin_pop, seed = args
    results, pop = run_simulation(label=scenario, initial_pop=burnin_pop, seed=seed)
    return scenario, results, pop

if __name__ == '__main__':
    start_time = datetime.datetime.now()
    print(f"Simulation started at: {start_time}")
    print(f"Numba available: {NUMBA_AVAILABLE}")
    print(f"Burn-in period: {BURN_IN_YEARS} years")
    print(f"Main simulation: {years} years")

    enabled_scenarios = get_enabled_scenarios()
    print(f"Enabled scenarios: {enabled_scenarios}")
    
    if len(enabled_scenarios) == 0:
        print("ERROR: No scenarios enabled. Please set at least one SCENARIO_X_ENABLED = True")
        exit(1)
    
    # =========================================================================

    # =========================================================================
    print("Running burn-in phase (shared by all scenarios)...")
    burnin_pop, _ = run_burnin(seed=SEED)
    print(f"Burn-in complete. Population state: "
          f"Alive={np.sum(burnin_pop['alive']):,}, "
          f"HIV+={np.sum(burnin_pop['hiv_positive'] & burnin_pop['alive']):,}, "
          f"On ART={np.sum(burnin_pop['on_art'] & burnin_pop['alive']):,}")
    
    # =========================================================================

    # =========================================================================
    all_results = {}
    all_pops = {}

    if len(enabled_scenarios) == 1:

        scenario = enabled_scenarios[0]
        print(f"Running scenario: {scenario}...")
        results, pop = run_simulation(label=scenario, initial_pop=burnin_pop, seed=SEED)
        all_results[scenario] = results
        all_pops[scenario] = pop
        print(f"  Completed: {scenario}")
    else:

        n_workers = len(enabled_scenarios)
        print(f"Running {n_workers} scenarios in parallel...")
        worker_args = [(s, burnin_pop, SEED) for s in enabled_scenarios]

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = list(executor.map(_run_scenario_worker, worker_args))

        for scenario, results, pop in futures:
            all_results[scenario] = results
            all_pops[scenario] = pop
            print(f"  Completed: {scenario}")

    # =========================================================================

    # =========================================================================
    save_all_outputs(all_results, all_pops, enabled_scenarios)
    export_report_ready_summary(all_results, all_pops, enabled_scenarios)

    # =========================================================================

    # =========================================================================
    if GRID_SEARCH_ENABLED:
        print("\n" + "="*70)
        print("Starting Grid Search for ART Delay Parameter Calibration...")
        print("="*70)

        result = grid_search_art_delay(
            label="2 months",
            burnin_pop=burnin_pop,
            delay_mean_grid=[1.0, 1.5, 2.0, 2.5, 3.0],
            target_mean_delay=2.0,
            delay_sd_grid=[0.3, 0.5, 0.7],
            calibration_years=5
        )

        print("\nOptimal parameters:", result['best'])

        grid_search_output_path = os.path.join(OUTPUT_DIR, "calibration_grid_search.csv")
        pd.DataFrame(result['all_rows']).to_csv(grid_search_output_path, index=False)
        print(f"Grid-search result saved: {grid_search_output_path}")
    else:
        print("\n[INFO] Grid search skipped (GRID_SEARCH_ENABLED = False)")

    end_time = datetime.datetime.now()
    elapsed = end_time - start_time
    print(f"\nSimulation ended at: {end_time}")
    print(f"Duration: {elapsed}")

    # =========================================================================

    # =========================================================================
    export_markdown_report(all_results, all_pops, enabled_scenarios, elapsed=elapsed)

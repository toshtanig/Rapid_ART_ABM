#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Cost-effectiveness analysis of the ART-timing scenarios from per-seed model outputs."""

import os
import json
import datetime
import random
import copy
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Any
from enum import Enum

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# --- Data location: set BASE_144 to your local 144_scenario folder ---
BASE_144 = r"/path/to/project/144_scenario"
PROJ_169 = os.path.dirname(BASE_144)
_PROJ_144 = BASE_144

CEA_INPUT_MODE = "mean_input"
if CEA_INPUT_MODE == "mean_input":
    _CEA_INPUT_DIR = os.path.join(_PROJ_144, "cea_param_uncertainty", "mean_input")
elif CEA_INPUT_MODE == "seed713":
    _CEA_INPUT_DIR = os.path.join(_PROJ_144, "core",
                                  "hiv_rapid_art_model_multiplier_mainfoi1p18_rr1_rate1_foi1p05_rfoi1p04_ageoff5_seed713")
else:
    raise ValueError('CEA_INPUT_MODE must be "mean_input" or "seed713"')
# =========================================================================

# =============================================================================

# =============================================================================
RANDOM_SEED = 42

def set_seed(seed: int = RANDOM_SEED) -> None:
    np.random.seed(seed)
    random.seed(seed)

# =============================================================================

# =============================================================================

@dataclass
class ParameterDef:
    name: str
    base_value: float
    se: float
    distribution: str  # "beta", "gamma", "normal", "lognormal"
    low: float = None
    high: float = None
    source: str = ""
    unit: str = ""
    
    def __post_init__(self):
        if self.low is None:
            self.low = self.base_value * 0.8
        if self.high is None:
            self.high = self.base_value * 1.2

@dataclass
class CostParameter(ParameterDef):
    distribution: str = "gamma"
    unit: str = "JPY"

@dataclass
class UtilityParameter(ParameterDef):
    distribution: str = "beta"
    unit: str = ""

class ParameterRegistry:
    
    def __init__(self):
        # ---------------------------------------------------------------------

        # Source: NDB (National Database of Health Insurance Claims), 2023

        # ---------------------------------------------------------------------
        self.costs_on_art = {
            "20-29": CostParameter(
                name="Annual ART cost (20-29)",
                base_value=2375000.0,
                se=475000.0,
                source="NDB FY2013-2020 (published; ref 17)"
            ),
            "30-39": CostParameter(
                name="Annual ART cost (30-39)",
                base_value=2901449.0,
                se=580290.0,
                source="NDB FY2013-2020 (published; ref 17)"
            ),
            "40-49": CostParameter(
                name="Annual ART cost (40-49)",
                base_value=3264767.0,
                se=652953.0,
                source="NDB FY2013-2020 (published; ref 17)"
            ),
            "50-59": CostParameter(
                name="Annual ART cost (50-59)",
                base_value=3433598.0,
                se=686720.0,
                source="NDB FY2013-2020 (published; ref 17)"
            ),
            "60-69": CostParameter(
                name="Annual ART cost (60-69)",
                base_value=3356900.0,
                se=671380.0,
                source="NDB FY2013-2020 (published; ref 17)"
            ),
            "70+": CostParameter(
                name="Annual ART cost (70+)",
                base_value=3164654.0,
                se=632931.0,
                source="NDB FY2013-2020 (published; ref 17)"
            ),
        }
        
        # ---------------------------------------------------------------------

        # Off-ART patients incur no ARV drug costs (~80% of on-ART costs in
        # Japan; Kimura 2002 [PMID 11968794]) but require monitoring visits, lab tests, OI
        # prophylaxis, and hospitalizations as CD4 declines.
        # - High CD4 (>500): monitoring only -> ~0.15-0.20 of on-ART
        # - Low CD4 (<200): OI treatment + hospitalization -> can exceed
        #   on-ART costs (Kimura 2002 [PMID 11968794]: inpatient AIDS CD4 50-199 = 768,000
        #   JPY/mo vs outpatient ART = 216,000 JPY/mo)
        # Time-weighted average across disease trajectory: ~0.40-0.50
        # ART medication is the dominant component of HIV care cost, so off-ART
        # person-time costs substantially less overall (Gebo 2010, HIV Research
        # Network CD4-stratified costs).
        # Base case 0.45 reflects weighted average of monitoring at higher
        # CD4 and increased OI/hospitalization costs at lower CD4.
        # ---------------------------------------------------------------------
        self.off_art_cost_ratio = ParameterDef(
            name="Off-ART cost ratio",
            base_value=0.45,
            se=0.10,
            distribution="beta",
            low=0.20,
            high=0.80,
            source="Kimura 2002 [PMID 11968794] (Japan inpatient/outpatient HIV costs); "
                   "Gebo 2010 (HIVRN CD4-stratified HIV care costs)",
            unit="ratio"
        )
        
        # ---------------------------------------------------------------------

        #

        #

        #   22 countries), MAINTAIN trial (EQ-5D-3L, CD4 375-750: 0.90),
        #   SMART/MMP data (SF-6D: 0.67-0.77, EQ-5D adjustment +0.10-0.13)

        # ---------------------------------------------------------------------

        self.utility_hiv_neg_by_age = {
            "30-39": 0.950,
            "40-49": 0.933,
            "50-59": 0.923,
            "60-69": 0.906,
            "70+":   0.857,
        }

        self.utility_hiv_neg = UtilityParameter(
            name="Utility (HIV-negative, base: age 30-39)",
            base_value=0.950,
            se=0.02,
            low=0.90,
            high=0.98,
            source="Shiroiwa et al. 2021 Value Health; "
                   "Japanese EQ-5D-5L population norms (male)"
        )
        
        self.utilities_by_cd4 = {
            ">=500": UtilityParameter(
                name="Utility (CD4>=500, suppressed on ART)",
                base_value=0.87,
                se=0.03,
                low=0.80,
                high=0.93,
                source="Poku 2025 (systematic review); "
                       "MAINTAIN trial EQ-5D-3L 0.90; "
                       "SMART/MMP SF-6D 0.73-0.77 (+EQ-5D adj); "
                       "adjusted for Japanese EQ-5D-5L value set"
            ),
            "350-499": UtilityParameter(
                name="Utility (CD4 350-499)",
                base_value=0.83,
                se=0.04,
                low=0.74,
                high=0.90,
                source="Poku 2025; SMART/MMP SF-6D 0.70-0.75 "
                       "(+EQ-5D adj); adjusted for Japanese context"
            ),
            "200-349": UtilityParameter(
                name="Utility (CD4 200-349)",
                base_value=0.76,
                se=0.05,
                low=0.65,
                high=0.85,
                source="Poku 2025; OPTIMA trial EQ-5D-3L 0.77 "
                       "(median CD4 108); adjusted for Japanese context"
            ),
            "<200": UtilityParameter(
                name="Utility (CD4<200)",
                base_value=0.65,
                se=0.06,
                low=0.50,
                high=0.78,
                source="Poku 2025; SMART/MMP SF-6D 0.67-0.70; "
                       "pre-ART studies range 0.44-0.85; "
                       "adjusted for Japanese EQ-5D-5L value set"
            ),
        }
        
        # ---------------------------------------------------------------------

        # Source: Shiroiwa et al. 2021 Value Health: Depression disutility
        #         on EQ-5D-5L = 0.177 (Japanese population)
        #         HIV diagnosis anxiety/depression prevalence = 30-50%
        #         (Brandt et al. 2017 Psychosomatics)
        # Calculation: 0.177 * 0.35 (midpoint) = 0.062
        # Rounded to 0.06 for base case
        # ---------------------------------------------------------------------
        self.mental_health_util_loss = UtilityParameter(
            name="Mental health utility loss per month",
            base_value=0.06,
            se=0.02,
            low=0.03,
            high=0.12,
            source="Shiroiwa 2021 (depression disutility 0.177); "
                   "Brandt 2017 (HIV anxiety prevalence 30-50%); "
                   "Calculation: 0.177 x 0.35 = 0.062"
        )
        
        # ---------------------------------------------------------------------

        # ---------------------------------------------------------------------
        self.discount_rate_cost = ParameterDef(
            name="Discount rate (costs)",
            base_value=0.02,
            se=0.005,
            distribution="normal",
            low=0.00,
            high=0.04,
            source="Japanese HTA guidelines (C2H)",
            unit="per year"
        )
        
        self.discount_rate_qaly = ParameterDef(
            name="Discount rate (QALYs)",
            base_value=0.02,
            se=0.005,
            distribution="normal",
            low=0.00,
            high=0.04,
            source="Japanese HTA guidelines (C2H)",
            unit="per year"
        )
    
    def get_all_parameters(self) -> List[ParameterDef]:
        params = [
            self.utility_hiv_neg,
            self.off_art_cost_ratio,
            self.mental_health_util_loss,
            self.discount_rate_cost,
            self.discount_rate_qaly,
        ]
        params.extend(self.costs_on_art.values())
        params.extend(self.utilities_by_cd4.values())
        return params
    
    def generate_parameter_table(self) -> pd.DataFrame:
        records = []
        for param in self.get_all_parameters():
            records.append({
                "Category": "Cost" if isinstance(param, CostParameter) else 
                           "Utility" if isinstance(param, UtilityParameter) else "Other",
                "Parameter": param.name,
                "Base Value": f"{param.base_value:,.4g}" if param.base_value < 1 
                              else f"{param.base_value:,.0f}",
                "SE": f"{param.se:,.4g}" if param.se < 1 else f"{param.se:,.0f}",
                "Distribution": param.distribution.capitalize(),
                "Range (DSA)": f"{param.low:,.4g} - {param.high:,.4g}" if param.low < 1 
                               else f"{param.low:,.0f} - {param.high:,.0f}",
                "Unit": param.unit,
                "Source": param.source,
            })
        return pd.DataFrame(records)

# =============================================================================

# =============================================================================

@dataclass
class ScenarioConfig:
    key: str
    label: str
    agg_file: str
    mean_diag_to_art_delay_months: float

class Config:

    INPUT_DIR = _CEA_INPUT_DIR

    OUTPUT_DIR = os.path.join(_PROJ_144, "bia_results", "cea_analysis")

    N_YEARS = 40

    SCENARIOS: List[ScenarioConfig] = [
        ScenarioConfig(
            key="0m",
            label="Rapid ART (0 month)",
            agg_file="scenario0_yearly_summary.csv",
            mean_diag_to_art_delay_months=0.0,
        ),
        ScenarioConfig(
            key="1m",
            label="Standard care (1 month)",
            agg_file="scenario1_yearly_summary.csv",
            mean_diag_to_art_delay_months=1.0,
        ),
        ScenarioConfig(
            key="2m",
            label="Standard care (2 months)",
            agg_file="scenario2_yearly_summary.csv",
            mean_diag_to_art_delay_months=2.0,
        ),
    ]

    BASELINE_KEY = "2m"

    ADDITIONAL_COMPARISONS = [
        ("0m", "1m"),
    ]

    USE_DISCOUNTING = True
    USE_HALF_CYCLE_CORRECTION = True

    PSA_N_ITERATIONS = 1000

    WTP_THRESHOLDS = [0, 1_000_000, 2_000_000, 3_000_000, 4_000_000, 
                      5_000_000, 6_000_000, 7_000_000, 8_000_000, 
                      9_000_000, 10_000_000, 15_000_000, 20_000_000]

    MAKE_FIGURES = True
    FIGURE_DPI = 300
    
    def __init__(self):
        self.params = ParameterRegistry()

# =============================================================================

# =============================================================================

class PSASampler:
    
    def __init__(self, params: ParameterRegistry, n_iterations: int):
        self.params = params
        self.n_iterations = n_iterations
        self.samples = {}
        
    def _sample_beta(self, mean: float, se: float) -> np.ndarray:
        if se <= 0 or mean <= 0 or mean >= 1:
            return np.full(self.n_iterations, mean)
        
        # Method of moments
        variance = se ** 2
        if variance >= mean * (1 - mean):
            variance = mean * (1 - mean) * 0.9
        
        alpha = mean * (mean * (1 - mean) / variance - 1)
        beta = (1 - mean) * (mean * (1 - mean) / variance - 1)
        
        if alpha <= 0 or beta <= 0:
            return np.full(self.n_iterations, mean)
        
        return stats.beta.rvs(alpha, beta, size=self.n_iterations)
    
    def _sample_gamma(self, mean: float, se: float) -> np.ndarray:
        if se <= 0 or mean <= 0:
            return np.full(self.n_iterations, mean)
        
        # shape = (mean/se)^2, scale = se^2/mean
        shape = (mean / se) ** 2
        scale = (se ** 2) / mean
        
        return stats.gamma.rvs(shape, scale=scale, size=self.n_iterations)
    
    def _sample_normal(self, mean: float, se: float) -> np.ndarray:
        return stats.norm.rvs(loc=mean, scale=se, size=self.n_iterations)
    
    def _sample_lognormal(self, mean: float, se: float) -> np.ndarray:
        if mean <= 0:
            return np.full(self.n_iterations, mean)

        variance = se ** 2
        mu = np.log(mean ** 2 / np.sqrt(variance + mean ** 2))
        sigma = np.sqrt(np.log(1 + variance / mean ** 2))
        
        return stats.lognorm.rvs(s=sigma, scale=np.exp(mu), size=self.n_iterations)
    
    def sample_parameter(self, param: ParameterDef) -> np.ndarray:
        if param.distribution == "beta":
            return self._sample_beta(param.base_value, param.se)
        elif param.distribution == "gamma":
            return self._sample_gamma(param.base_value, param.se)
        elif param.distribution == "normal":
            return self._sample_normal(param.base_value, param.se)
        elif param.distribution == "lognormal":
            return self._sample_lognormal(param.base_value, param.se)
        else:
            return np.full(self.n_iterations, param.base_value)
    
    def generate_all_samples(self) -> Dict[str, np.ndarray]:
        samples = {}

        for age, param in self.params.costs_on_art.items():
            samples[f"cost_on_art_{age}"] = self.sample_parameter(param)

        samples["utility_hiv_neg"] = self.sample_parameter(self.params.utility_hiv_neg)
        
        for cd4, param in self.params.utilities_by_cd4.items():
            samples[f"utility_cd4_{cd4}"] = self.sample_parameter(param)

        samples["off_art_cost_ratio"] = self.sample_parameter(self.params.off_art_cost_ratio)
        samples["mental_health_util_loss"] = self.sample_parameter(self.params.mental_health_util_loss)
        samples["discount_rate_cost"] = self.sample_parameter(self.params.discount_rate_cost)
        samples["discount_rate_qaly"] = self.sample_parameter(self.params.discount_rate_qaly)

        samples["discount_rate_cost"] = np.clip(samples["discount_rate_cost"], 0, 0.1)
        samples["discount_rate_qaly"] = np.clip(samples["discount_rate_qaly"], 0, 0.1)
        
        self.samples = samples
        return samples

# =============================================================================

# =============================================================================

def ensure_output_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def discount_factor(year_index: int, rate: float, use_discounting: bool) -> float:
    if (not use_discounting) or rate <= 0.0:
        return 1.0

    return 1.0 / ((1.0 + rate) ** (year_index - 1))

def load_aggregator(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Aggregator file not found: {path}")
    df = pd.read_csv(path)
    expected_cols = {
        "Year", "AgeGroup", "CD4Category", "OnART", "OnPrEP",
        "PersonMonthsAlive", "NewlyInfected", "NewlyDiagnosed", "Deaths"
    }
    missing = expected_cols.difference(df.columns)
    if missing:
        raise ValueError(f"Aggregator missing columns: {missing}")
    return df

def compute_qalys_and_costs(
    df_agg: pd.DataFrame,
    scen: ScenarioConfig,
    costs_on_art: Dict[str, float],
    off_art_cost_ratio: float,
    utility_hiv_neg: float,
    utilities_by_cd4: Dict[str, float],
    mental_health_util_loss: float,
    discount_rate_cost: float,
    discount_rate_qaly: float,
    use_discounting: bool,
    use_half_cycle_correction: bool,
    n_years: int,
    utility_hiv_neg_by_age: Optional[Dict[str, float]] = None,
) -> Dict:
    df = df_agg.copy()
    
    # =========================================================================

    # =========================================================================

    cd4_utilities = df["CD4Category"].map(utilities_by_cd4)
    cd4_utilities = cd4_utilities.fillna(min(utilities_by_cd4.values())).to_numpy(dtype=float)
    if utility_hiv_neg_by_age is not None:
        hiv_neg_utilities = (
            df["AgeGroup"]
            .map(utility_hiv_neg_by_age)
            .fillna(utility_hiv_neg)
            .to_numpy(dtype=float)
        )
    else:
        hiv_neg_utilities = utility_hiv_neg
    is_hiv_negative = df["CD4Category"].eq("HIV-").to_numpy(dtype=bool)
    df["BaseUtility"] = np.where(is_hiv_negative, hiv_neg_utilities, cd4_utilities)
    df["QALY_base"] = df["PersonMonthsAlive"] * df["BaseUtility"] / 12.0

    df_hiv_pos = df[df["CD4Category"] != "HIV-"]
    diag_by_year = df_hiv_pos.groupby("Year")["NewlyDiagnosed"].sum()
    qaly_mental_loss_by_year = (
        diag_by_year
        * scen.mean_diag_to_art_delay_months
        * mental_health_util_loss
        / 12.0
    )
    
    # =========================================================================

    # =========================================================================
    years = np.sort(df["Year"].unique())
    min_year = int(years.min())
    max_year = int(years.max())

    actual_max_year = min(max_year, n_years)
    year_range = range(min_year, actual_max_year + 1)

    qaly_mental_loss_by_year = qaly_mental_loss_by_year.reindex(year_range, fill_value=0.0)

    qaly_base_by_year = df.groupby("Year")["QALY_base"].sum()
    qaly_total_by_year = qaly_base_by_year.reindex(year_range, fill_value=0.0)
    qaly_total_by_year = qaly_total_by_year - qaly_mental_loss_by_year

    annual_cost_art = df["AgeGroup"].map(costs_on_art).fillna(np.mean(list(costs_on_art.values())))
    annual_cost = np.where(df["OnART"].astype(bool), annual_cost_art, annual_cost_art * off_art_cost_ratio)
    annual_cost = np.where(df["CD4Category"].eq("HIV-"), 0.0, annual_cost)
    df["Cost_hivcare"] = annual_cost * (df["PersonMonthsAlive"] / 12.0)
    cost_by_year = df.groupby("Year")["Cost_hivcare"].sum()
    cost_by_year = cost_by_year.reindex(year_range, fill_value=0.0)

    disc_factors_cost = np.array([
        discount_factor(y, discount_rate_cost, use_discounting)
        for y in year_range
    ])
    disc_factors_qaly = np.array([
        discount_factor(y, discount_rate_qaly, use_discounting)
        for y in year_range
    ])

    qaly_values = qaly_total_by_year.values.copy()
    cost_values = cost_by_year.values.copy()
    
    if use_half_cycle_correction and len(qaly_values) > 1:
        qaly_values[0] *= 0.5
        qaly_values[-1] *= 0.5
        cost_values[0] *= 0.5
        cost_values[-1] *= 0.5

    qaly_disc = qaly_values * disc_factors_qaly
    cost_disc = cost_values * disc_factors_cost
    
    return {
        "scenario_key": scen.key,
        "scenario_label": scen.label,
        "years": np.array(list(year_range)),
        "qaly_per_year": qaly_values,
        "qaly_disc_per_year": qaly_disc,
        "total_qaly": float(qaly_values.sum()),
        "total_qaly_disc": float(qaly_disc.sum()),
        "cost_per_year": cost_values,
        "cost_disc_per_year": cost_disc,
        "total_cost": float(cost_values.sum()),
        "total_cost_disc": float(cost_disc.sum()),
    }

def compute_incremental(
    res_intervention: Dict,
    res_baseline: Dict,
) -> Dict:
    delta_cost = res_intervention["total_cost_disc"] - res_baseline["total_cost_disc"]
    delta_qaly = res_intervention["total_qaly_disc"] - res_baseline["total_qaly_disc"]

    if delta_qaly > 0:
        icer = delta_cost / delta_qaly
    elif delta_qaly < 0:
        icer = np.nan  # Dominated or dominated
    else:
        icer = np.nan

    if delta_cost < 0 and delta_qaly > 0:
        quadrant = "Dominant"
    elif delta_cost > 0 and delta_qaly < 0:
        quadrant = "Dominated"
    elif delta_cost < 0 and delta_qaly < 0:
        quadrant = "SW (less effective, less costly)"
    else:
        quadrant = "NE (more effective, more costly)"
    
    return {
        "delta_cost": delta_cost,
        "delta_qaly": delta_qaly,
        "icer": icer,
        "quadrant": quadrant,
    }

def compute_nmb(delta_cost: float, delta_qaly: float, wtp: float) -> float:
    return delta_qaly * wtp - delta_cost

# =============================================================================

# =============================================================================

def run_psa(
    df_aggs: Dict[str, pd.DataFrame],
    scenarios: List[ScenarioConfig],
    baseline_key: str,
    additional_comparisons: List[Tuple[str, str]],
    config: Config,
) -> pd.DataFrame:
    print(f"\n[PSA] Running {config.PSA_N_ITERATIONS} iterations...")
    
    sampler = PSASampler(config.params, config.PSA_N_ITERATIONS)
    samples = sampler.generate_all_samples()
    
    results = []

    comparison_pairs = []
    for scen in scenarios:
        if scen.key != baseline_key:
            comparison_pairs.append((scen.key, baseline_key))
    for pair in additional_comparisons:
        if pair not in comparison_pairs:
            comparison_pairs.append(pair)
    
    for i in range(config.PSA_N_ITERATIONS):
        if (i + 1) % 200 == 0:
            print(f"  Iteration {i + 1}/{config.PSA_N_ITERATIONS}")

        costs_on_art_i = {
            age: samples[f"cost_on_art_{age}"][i]
            for age in config.params.costs_on_art.keys()
        }
        utilities_by_cd4_i = {
            cd4: samples[f"utility_cd4_{cd4}"][i]
            for cd4 in config.params.utilities_by_cd4.keys()
        }

        psa_ratio = samples["utility_hiv_neg"][i] / config.params.utility_hiv_neg.base_value
        utility_hiv_neg_by_age_i = {
            age: min(val * psa_ratio, 1.0)
            for age, val in config.params.utility_hiv_neg_by_age.items()
        }

        res_by_key = {}
        for scen in scenarios:
            res = compute_qalys_and_costs(
                df_agg=df_aggs[scen.key],
                scen=scen,
                costs_on_art=costs_on_art_i,
                off_art_cost_ratio=samples["off_art_cost_ratio"][i],
                utility_hiv_neg=samples["utility_hiv_neg"][i],
                utilities_by_cd4=utilities_by_cd4_i,
                mental_health_util_loss=samples["mental_health_util_loss"][i],
                discount_rate_cost=samples["discount_rate_cost"][i],
                discount_rate_qaly=samples["discount_rate_qaly"][i],
                use_discounting=config.USE_DISCOUNTING,
                use_half_cycle_correction=config.USE_HALF_CYCLE_CORRECTION,
                n_years=config.N_YEARS,
                utility_hiv_neg_by_age=utility_hiv_neg_by_age_i,
            )
            res_by_key[scen.key] = res

        for intervention_key, comparator_key in comparison_pairs:
            inc = compute_incremental(res_by_key[intervention_key], res_by_key[comparator_key])

            nmb_results = {}
            for wtp in config.WTP_THRESHOLDS:
                nmb = compute_nmb(inc["delta_cost"], inc["delta_qaly"], wtp)
                nmb_results[f"NMB_WTP{wtp // 1_000_000}M"] = nmb
            
            comparison_label = f"{intervention_key}_vs_{comparator_key}"
            
            results.append({
                "iteration": i,
                "comparison": comparison_label,
                "intervention": intervention_key,
                "comparator": comparator_key,
                "total_cost_disc_intervention": res_by_key[intervention_key]["total_cost_disc"],
                "total_qaly_disc_intervention": res_by_key[intervention_key]["total_qaly_disc"],
                "total_cost_disc_comparator": res_by_key[comparator_key]["total_cost_disc"],
                "total_qaly_disc_comparator": res_by_key[comparator_key]["total_qaly_disc"],
                "delta_cost": inc["delta_cost"],
                "delta_qaly": inc["delta_qaly"],
                "icer": inc["icer"],
                "quadrant": inc["quadrant"],
                **nmb_results,
            })
    
    return pd.DataFrame(results)

# =============================================================================

# =============================================================================

def run_dsa(
    df_aggs: Dict[str, pd.DataFrame],
    scenarios: List[ScenarioConfig],
    baseline_key: str,
    comparator_key: str,
    config: Config,
) -> pd.DataFrame:
    print(f"\n[DSA] Running deterministic sensitivity analysis...")
    
    params = config.params

    base_costs = {age: p.base_value for age, p in params.costs_on_art.items()}
    base_utilities = {cd4: p.base_value for cd4, p in params.utilities_by_cd4.items()}
    
    base_params = {
        "costs_on_art": base_costs,
        "off_art_cost_ratio": params.off_art_cost_ratio.base_value,
        "utility_hiv_neg": params.utility_hiv_neg.base_value,
        "utilities_by_cd4": base_utilities,
        "mental_health_util_loss": params.mental_health_util_loss.base_value,
        "discount_rate_cost": params.discount_rate_cost.base_value,
        "discount_rate_qaly": params.discount_rate_qaly.base_value,
        "utility_hiv_neg_by_age": dict(params.utility_hiv_neg_by_age),
    }
    
    def run_with_params(param_dict: Dict) -> Dict:

        dsa_ratio = param_dict["utility_hiv_neg"] / params.utility_hiv_neg.base_value
        scaled_age_utils = {
            age: min(val * dsa_ratio, 1.0)
            for age, val in param_dict.get("utility_hiv_neg_by_age",
                                           params.utility_hiv_neg_by_age).items()
        }
        
        res_by_key = {}
        for scen in scenarios:
            res = compute_qalys_and_costs(
                df_agg=df_aggs[scen.key],
                scen=scen,
                costs_on_art=param_dict["costs_on_art"],
                off_art_cost_ratio=param_dict["off_art_cost_ratio"],
                utility_hiv_neg=param_dict["utility_hiv_neg"],
                utilities_by_cd4=param_dict["utilities_by_cd4"],
                mental_health_util_loss=param_dict["mental_health_util_loss"],
                discount_rate_cost=param_dict["discount_rate_cost"],
                discount_rate_qaly=param_dict["discount_rate_qaly"],
                use_discounting=config.USE_DISCOUNTING,
                use_half_cycle_correction=config.USE_HALF_CYCLE_CORRECTION,
                n_years=config.N_YEARS,
                utility_hiv_neg_by_age=scaled_age_utils,
            )
            res_by_key[scen.key] = res
        
        return compute_incremental(res_by_key[comparator_key], res_by_key[baseline_key])

    base_result = run_with_params(base_params)
    base_icer = base_result["icer"]
    base_nmb = compute_nmb(base_result["delta_cost"], base_result["delta_qaly"], 5_000_000)
    
    results = []

    dsa_params = [
        ("Discount rate (costs)", params.discount_rate_cost, "discount_rate_cost"),
        ("Discount rate (QALYs)", params.discount_rate_qaly, "discount_rate_qaly"),
        ("Off-ART cost ratio", params.off_art_cost_ratio, "off_art_cost_ratio"),
        ("Mental health utility loss", params.mental_health_util_loss, "mental_health_util_loss"),
        ("Utility (HIV-negative)", params.utility_hiv_neg, "utility_hiv_neg"),
    ]

    for cd4, param in params.utilities_by_cd4.items():
        dsa_params.append((f"Utility (CD4 {cd4})", param, f"utility_cd4_{cd4}"))

    for age, param in params.costs_on_art.items():
        dsa_params.append((f"Annual ART cost ({age})", param, f"cost_on_art_{age}"))
    
    for name, param, key in dsa_params:
        for value_type, value in [("low", param.low), ("high", param.high)]:

            test_params = copy.deepcopy(base_params)
            
            if key.startswith("utility_cd4_"):
                cd4 = key.replace("utility_cd4_", "")
                test_params["utilities_by_cd4"][cd4] = value
            elif key.startswith("cost_on_art_"):
                age = key.replace("cost_on_art_", "")
                test_params["costs_on_art"][age] = value
            else:
                test_params[key] = value
            
            result = run_with_params(test_params)
            nmb = compute_nmb(result["delta_cost"], result["delta_qaly"], 5_000_000)
            
            results.append({
                "parameter": name,
                "value_type": value_type,
                "value": value,
                "delta_cost": result["delta_cost"],
                "delta_qaly": result["delta_qaly"],
                "icer": result["icer"],
                "nmb_5m": nmb,
            })
    
    df_dsa = pd.DataFrame(results)

    df_dsa["base_icer"] = base_icer
    df_dsa["base_nmb_5m"] = base_nmb
    
    return df_dsa

# =============================================================================
# CEAC（Cost-Effectiveness Acceptability Curve）
# =============================================================================

def compute_ceac(psa_results: pd.DataFrame, wtp_thresholds: List[float]) -> pd.DataFrame:
    records = []
    
    for comparison in psa_results["comparison"].unique():
        df_comp = psa_results[psa_results["comparison"] == comparison]
        
        for wtp in wtp_thresholds:

            nmb = df_comp["delta_qaly"] * wtp - df_comp["delta_cost"]
            prob_ce = (nmb > 0).mean()
            
            records.append({
                "comparison": comparison,
                "wtp": wtp,
                "probability_cost_effective": prob_ce,
            })
    
    return pd.DataFrame(records)

# =============================================================================

# =============================================================================

def validate_model(res_by_key: Dict[str, Dict], df_aggs: Dict[str, pd.DataFrame]) -> Dict:
    validation = {
        "internal": {},
        "external": {},
        "summary": [],
    }

    for key, df_agg in df_aggs.items():

        hiv_pos = df_agg[df_agg["CD4Category"] != "HIV-"]
        total_person_months = hiv_pos["PersonMonthsAlive"].sum()
        total_deaths = hiv_pos["Deaths"].sum()
        total_new_infections = hiv_pos["NewlyInfected"].sum()
        
        validation["internal"][key] = {
            "total_person_months": total_person_months,
            "total_deaths": total_deaths,
            "total_new_infections": total_new_infections,
            "implied_person_years": total_person_months / 12,
        }

    observed_data = {
        "annual_new_diagnoses_2020": 1095,
        "annual_new_diagnoses_2021": 1057,
        "annual_new_diagnoses_2022": 884,
    }
    
    validation["external"]["observed_data"] = observed_data
    validation["external"]["note"] = "Compare model Year 1 diagnoses with observed data"

    for key, res in res_by_key.items():
        validation["summary"].append({
            "scenario": key,
            "total_qaly_disc": res["total_qaly_disc"],
            "total_cost_disc": res["total_cost_disc"],
            "qaly_per_person_year": res["total_qaly_disc"] / (validation["internal"][key]["implied_person_years"]) 
                                    if validation["internal"][key]["implied_person_years"] > 0 else np.nan,
        })
    
    return validation

# =============================================================================

# =============================================================================

def check_simulation_data_consistency(df_aggs: Dict[str, pd.DataFrame]) -> Dict:
    warnings = []
    stats = {}
    
    for key, df in df_aggs.items():
        hiv_pos = df[df["CD4Category"] != "HIV-"]
        stats[key] = {
            "total_person_months": hiv_pos["PersonMonthsAlive"].sum(),
            "total_newly_infected": hiv_pos["NewlyInfected"].sum(),
            "total_deaths": hiv_pos["Deaths"].sum(),
            "total_newly_diagnosed": hiv_pos["NewlyDiagnosed"].sum(),
        }

    keys = list(stats.keys())
    if len(keys) >= 2:

        pm_values = [stats[k]["total_person_months"] for k in keys]
        max_diff_pct = (max(pm_values) - min(pm_values)) / np.mean(pm_values) * 100
        
        if max_diff_pct > 5:
            warnings.append(
                f"WARNING: Large variation in PersonMonthsAlive across scenarios "
                f"({max_diff_pct:.1f}% difference). This may indicate that scenarios "
                f"were run with different random seeds, introducing stochastic noise "
                f"that could mask the true effect of ART timing."
            )

        ni_values = [stats[k]["total_newly_infected"] for k in keys]
        ni_diff_pct = (max(ni_values) - min(ni_values)) / np.mean(ni_values) * 100
        
        if ni_diff_pct > 10:
            warnings.append(
                f"WARNING: Large variation in NewlyInfected across scenarios "
                f"({ni_diff_pct:.1f}% difference). Consider re-running simulations "
                f"with the same random seed for all scenarios."
            )
    
    return {
        "stats": stats,
        "warnings": warnings,
    }

# =============================================================================

# =============================================================================

def plot_tornado_diagram(dsa_results: pd.DataFrame, output_path: str, dpi: int = 300) -> None:
    comparisons = dsa_results["comparison"].unique()
    n_comparisons = len(comparisons)
    
    fig, axes = plt.subplots(1, n_comparisons, figsize=(12 * n_comparisons // 2, 10))
    
    if n_comparisons == 1:
        axes = [axes]
    
    labels = {
        "0m_vs_2m": "0 vs 2 months",
        "0m_vs_1m": "0 vs 1 month",
    }
    
    for idx, comparison in enumerate(comparisons):
        ax = axes[idx]
        df_comp = dsa_results[dsa_results["comparison"] == comparison]

        df_pivot = df_comp.pivot(index="parameter", columns="value_type", values="nmb_5m")
        df_pivot["base"] = df_comp.groupby("parameter")["base_nmb_5m"].first()
        df_pivot["range"] = abs(df_pivot["high"] - df_pivot["low"])
        df_pivot = df_pivot.sort_values("range", ascending=True)

        y_pos = np.arange(len(df_pivot))
        
        low_bars = (df_pivot["low"] - df_pivot["base"]) / 1e9
        high_bars = (df_pivot["high"] - df_pivot["base"]) / 1e9
        
        ax.barh(y_pos, low_bars, align="center", color="steelblue", alpha=0.7, label="Low")
        ax.barh(y_pos, high_bars, align="center", color="darkorange", alpha=0.7, label="High")
        
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df_pivot.index, fontsize=9)
        ax.set_xlabel("Net monetary benefit vs base case (billion JPY)", fontsize=10)
        ax.set_title(f"Tornado Diagram: {labels.get(comparison, comparison)}", fontsize=12)
        ax.legend(loc="lower right")
        ax.grid(axis="x", alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()

def plot_ce_plane_psa(psa_results: pd.DataFrame, output_path: str, dpi: int = 300) -> None:
    comparisons = psa_results["comparison"].unique()

    fig, ax = plt.subplots(figsize=(10, 8))

    colors = {"0m_vs_2m": "blue", "1m_vs_2m": "green", "0m_vs_1m": "purple"}
    markers = {"0m_vs_2m": "o", "1m_vs_2m": "s", "0m_vs_1m": "^"}
    labels = {"0m_vs_2m": "0 vs 2 months", "1m_vs_2m": "1 vs 2 months", "0m_vs_1m": "0 vs 1 month"}

    for comparison in comparisons:
        df_comp = psa_results[psa_results["comparison"] == comparison]

        color = colors.get(comparison, "blue")
        marker = markers.get(comparison, "o")
        label = labels.get(comparison, comparison)

        ax.scatter(
            df_comp["delta_qaly"],
            df_comp["delta_cost"] / 1e9,
            alpha=0.3,
            s=15,
            color=color,
            marker=marker,
            label=label
        )

    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5)

    x_range = ax.get_xlim()
    x_line = np.linspace(x_range[0], x_range[1], 100)
    y_line = x_line * 5_000_000 / 1e9
    ax.plot(x_line, y_line, "r--", linewidth=1, label="WTP = 5M JPY/QALY")

    ax.set_xlabel("Incremental QALYs", fontsize=11)
    ax.set_ylabel("Incremental Cost (Billion JPY)", fontsize=11)
    ax.set_title("Cost-Effectiveness Plane (PSA)", fontsize=12)
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()

def plot_ceac(ceac_results: pd.DataFrame, output_path: str, dpi: int = 300) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    
    colors = {"0m_vs_2m": "blue", "1m_vs_2m": "green", "0m_vs_1m": "purple"}
    labels = {"0m_vs_2m": "0 vs 2 months", "1m_vs_2m": "1 vs 2 months", "0m_vs_1m": "0 vs 1 month"}
    
    for comparison in ceac_results["comparison"].unique():
        df_comp = ceac_results[ceac_results["comparison"] == comparison]
        color = colors.get(comparison, "blue")
        label = labels.get(comparison, comparison)
        
        ax.plot(
            df_comp["wtp"] / 1e6,
            df_comp["probability_cost_effective"],
            color=color,
            linewidth=2,
            label=label
        )
    
    ax.axvline(5, color="red", linestyle="--", linewidth=1, alpha=0.5, label="WTP = 5M JPY/QALY")
    
    ax.set_xlabel("Willingness-to-Pay Threshold (Million JPY/QALY)", fontsize=11)
    ax.set_ylabel("Probability Cost-Effective", fontsize=11)
    ax.set_title("Cost-Effectiveness Acceptability Curve (CEAC)", fontsize=12)
    ax.set_ylim(0, 1.05)
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()

def plot_psa_histogram(psa_results: pd.DataFrame, output_path: str, dpi: int = 300) -> None:
    comparisons = psa_results["comparison"].unique()
    n_comparisons = len(comparisons)
    
    fig, axes = plt.subplots(2, n_comparisons, figsize=(6 * n_comparisons, 8))
    
    if n_comparisons == 1:
        axes = axes.reshape(-1, 1)
    
    labels = {"0m_vs_2m": "0 vs 2 months", "1m_vs_2m": "1 vs 2 months", "0m_vs_1m": "0 vs 1 month"}
    
    for idx, comparison in enumerate(comparisons):
        df_comp = psa_results[psa_results["comparison"] == comparison]
        label = labels.get(comparison, comparison)

        ax1 = axes[0, idx]
        ax1.hist(df_comp["delta_cost"] / 1e9, bins=50, color="steelblue", alpha=0.7, edgecolor="white")
        ax1.axvline(0, color="red", linestyle="--", linewidth=1)
        ax1.set_xlabel("Incremental cost (billion JPY)", fontsize=11)
        ax1.set_ylabel("Frequency", fontsize=11)
        ax1.set_title(f"{label}: incremental cost distribution", fontsize=12)

        ax2 = axes[1, idx]
        ax2.hist(df_comp["delta_qaly"], bins=50, color="darkorange", alpha=0.7, edgecolor="white")
        ax2.axvline(0, color="red", linestyle="--", linewidth=1)
        ax2.set_xlabel("Incremental QALYs", fontsize=11)
        ax2.set_ylabel("Frequency", fontsize=11)
        ax2.set_title(f"{label}: incremental QALY distribution", fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()

# =============================================================================

# =============================================================================

def generate_summary_tables(
    base_results: Dict[str, Dict],
    psa_results: pd.DataFrame,
    dsa_results: pd.DataFrame,
    baseline_key: str,
    additional_comparisons: List[Tuple[str, str]],
    config: Config,
) -> Dict[str, pd.DataFrame]:
    tables = {}

    base_records = []
    for key, res in base_results.items():
        base_records.append({
            "Scenario": res["scenario_label"],
            "Total QALYs (discounted)": f"{res['total_qaly_disc']:,.2f}",
            "Total Cost (discounted)": f"¥{res['total_cost_disc']:,.0f}",
        })
    tables["base_case_results"] = pd.DataFrame(base_records)

    inc_records = []
    baseline_res = base_results[baseline_key]
    for key, res in base_results.items():
        if key == baseline_key:
            continue
        inc = compute_incremental(res, baseline_res)
        inc_records.append({
            "Comparison": f"{key} vs {baseline_key}",
            "Incremental Cost": f"¥{inc['delta_cost']:,.0f}",
            "Incremental QALYs": f"{inc['delta_qaly']:,.2f}",
            "ICER (JPY/QALY)": f"{inc['icer']:,.0f}" if not np.isnan(inc['icer']) else inc['quadrant'],
            "Quadrant": inc['quadrant'],
        })

    for intervention_key, comparator_key in additional_comparisons:
        inc = compute_incremental(base_results[intervention_key], base_results[comparator_key])
        inc_records.append({
            "Comparison": f"{intervention_key} vs {comparator_key}",
            "Incremental Cost": f"¥{inc['delta_cost']:,.0f}",
            "Incremental QALYs": f"{inc['delta_qaly']:,.2f}",
            "ICER (JPY/QALY)": f"{inc['icer']:,.0f}" if not np.isnan(inc['icer']) else inc['quadrant'],
            "Quadrant": inc['quadrant'],
        })
    tables["incremental_results"] = pd.DataFrame(inc_records)

    psa_summary_records = []
    for comparison in psa_results["comparison"].unique():
        df_comp = psa_results[psa_results["comparison"] == comparison]
        prob_dominant = ((df_comp["delta_cost"] < 0) & (df_comp["delta_qaly"] > 0)).mean()
        prob_ce_5m = (df_comp["delta_qaly"] * 5_000_000 - df_comp["delta_cost"] > 0).mean()
        
        psa_summary_records.append({
            "Comparison": comparison,
            "Mean Incremental Cost": f"¥{df_comp['delta_cost'].mean():,.0f}",
            "95% CI (Cost)": f"(¥{df_comp['delta_cost'].quantile(0.025):,.0f}, ¥{df_comp['delta_cost'].quantile(0.975):,.0f})",
            "Mean Incremental QALYs": f"{df_comp['delta_qaly'].mean():,.2f}",
            "95% CI (QALYs)": f"({df_comp['delta_qaly'].quantile(0.025):,.2f}, {df_comp['delta_qaly'].quantile(0.975):,.2f})",
            "Prob(Dominant)": f"{prob_dominant:.1%}",
            "Prob(CE at WTP=5M)": f"{prob_ce_5m:.1%}",
        })
    tables["psa_summary"] = pd.DataFrame(psa_summary_records)
    
    return tables

# =============================================================================

# =============================================================================

def export_cea_report_json(
    base_results: Dict[str, Dict],
    psa_results: pd.DataFrame,
    dsa_results: pd.DataFrame,
    ceac_results: pd.DataFrame,
    validation: Dict,
    baseline_key: str,
    additional_comparisons: List[Tuple[str, str]],
    config: Config,
    output_path: str,
) -> None:

    def safe_val(v):
        if isinstance(v, (np.integer,)):
            return int(v)
        if isinstance(v, (np.floating,)):
            return None if np.isnan(v) else float(v)
        if isinstance(v, np.ndarray):
            return [safe_val(x) for x in v]
        return v

    report = {}

    # metadata
    report["metadata"] = {
        "script": "cea_analysis_rapid_v3.py",
        "created_at": datetime.datetime.now().isoformat(),
        "n_years": config.N_YEARS,
        "psa_n_iterations": config.PSA_N_ITERATIONS,
        "random_seed": RANDOM_SEED,
        "baseline_key": baseline_key,
        "use_discounting": config.USE_DISCOUNTING,
        "use_half_cycle_correction": config.USE_HALF_CYCLE_CORRECTION,
        "wtp_thresholds": config.WTP_THRESHOLDS,
    }

    # cost_and_utility_parameters
    params = config.params
    param_table = params.generate_parameter_table()
    report["parameters"] = param_table.to_dict(orient="records")

    base_case = {}
    for key, res in base_results.items():
        base_case[key] = {
            "scenario_label": res["scenario_label"],
            "total_qaly": safe_val(res["total_qaly"]),
            "total_qaly_discounted": safe_val(res["total_qaly_disc"]),
            "total_cost": safe_val(res["total_cost"]),
            "total_cost_discounted": safe_val(res["total_cost_disc"]),
            "yearly_qaly": safe_val(res["qaly_per_year"]),
            "yearly_qaly_discounted": safe_val(res["qaly_disc_per_year"]),
            "yearly_cost": safe_val(res["cost_per_year"]),
            "yearly_cost_discounted": safe_val(res["cost_disc_per_year"]),
        }
    report["base_case_results"] = base_case

    incremental = []
    baseline_res = base_results[baseline_key]
    for key, res in base_results.items():
        if key == baseline_key:
            continue
        inc = compute_incremental(res, baseline_res)
        incremental.append({
            "comparison": f"{key}_vs_{baseline_key}",
            "intervention": key,
            "comparator": baseline_key,
            "delta_cost": safe_val(inc["delta_cost"]),
            "delta_qaly": safe_val(inc["delta_qaly"]),
            "icer": safe_val(inc["icer"]),
            "quadrant": inc["quadrant"],
        })

    for intervention_key, comparator_key in additional_comparisons:
        inc = compute_incremental(base_results[intervention_key], base_results[comparator_key])
        incremental.append({
            "comparison": f"{intervention_key}_vs_{comparator_key}",
            "intervention": intervention_key,
            "comparator": comparator_key,
            "delta_cost": safe_val(inc["delta_cost"]),
            "delta_qaly": safe_val(inc["delta_qaly"]),
            "icer": safe_val(inc["icer"]),
            "quadrant": inc["quadrant"],
        })
    report["incremental_results"] = incremental

    # psa_summary
    psa_summary = []
    for comparison in psa_results["comparison"].unique():
        df_comp = psa_results[psa_results["comparison"] == comparison]
        n = len(df_comp)
        n_dominant = int(((df_comp["delta_cost"] < 0) & (df_comp["delta_qaly"] > 0)).sum())
        n_dominated = int(((df_comp["delta_cost"] > 0) & (df_comp["delta_qaly"] < 0)).sum())
        n_ne = int(((df_comp["delta_cost"] > 0) & (df_comp["delta_qaly"] > 0)).sum())
        n_sw = int(((df_comp["delta_cost"] < 0) & (df_comp["delta_qaly"] < 0)).sum())

        nmb_5m = df_comp["delta_qaly"] * 5_000_000 - df_comp["delta_cost"]

        psa_summary.append({
            "comparison": comparison,
            "n_iterations": n,
            "delta_cost_mean": float(df_comp["delta_cost"].mean()),
            "delta_cost_ci_lower": float(df_comp["delta_cost"].quantile(0.025)),
            "delta_cost_ci_upper": float(df_comp["delta_cost"].quantile(0.975)),
            "delta_qaly_mean": float(df_comp["delta_qaly"].mean()),
            "delta_qaly_ci_lower": float(df_comp["delta_qaly"].quantile(0.025)),
            "delta_qaly_ci_upper": float(df_comp["delta_qaly"].quantile(0.975)),
            "prob_dominant": n_dominant / n,
            "prob_dominated": n_dominated / n,
            "prob_ne_quadrant": n_ne / n,
            "prob_sw_quadrant": n_sw / n,
            "nmb_wtp5m_mean": float(nmb_5m.mean()),
            "prob_cost_effective_wtp5m": float((nmb_5m > 0).mean()),
        })
    report["psa_summary"] = psa_summary

    ceac_summary = []
    for comparison in ceac_results["comparison"].unique():
        df_c = ceac_results[ceac_results["comparison"] == comparison]
        for _, row in df_c.iterrows():
            ceac_summary.append({
                "comparison": row["comparison"],
                "wtp": safe_val(row["wtp"]),
                "probability_cost_effective": safe_val(row["probability_cost_effective"]),
            })
    report["ceac"] = ceac_summary

    # dsa_summary
    dsa_list = []
    for _, row in dsa_results.iterrows():
        entry = {}
        for col in dsa_results.columns:
            entry[col] = safe_val(row[col])
        dsa_list.append(entry)
    report["dsa_results"] = dsa_list

    # validation
    val_clean = {
        "internal": {},
        "external": validation.get("external", {}),
        "summary": [],
    }
    for key, v in validation.get("internal", {}).items():
        val_clean["internal"][key] = {k: safe_val(vv) for k, vv in v.items()}
    for item in validation.get("summary", []):
        val_clean["summary"].append({k: safe_val(vv) for k, vv in item.items()})
    report["validation"] = val_clean

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)

    print(f"  Saved: {os.path.basename(output_path)}")

# =============================================================================

# =============================================================================

def main():
    print("=" * 70)
    print("Cost-Effectiveness Analysis: Rapid ART Implementation in Japan")
    print("Version 4.0 (Full population QALY + Year index fix + Seed warning)")
    print("=" * 70)
    
    set_seed(RANDOM_SEED)
    config = Config()
    ensure_output_dir(config.OUTPUT_DIR)
    
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    print("\n[1] Loading aggregator data...")
    
    df_aggs = {}
    available_scenarios = []
    
    for scen in config.SCENARIOS:
        path = os.path.join(config.INPUT_DIR, scen.agg_file)
        if os.path.exists(path):
            df_aggs[scen.key] = load_aggregator(path)
            available_scenarios.append(scen)
            print(f"  Loaded: {scen.agg_file}")
        else:
            print(f"  [SKIP] {scen.agg_file} not found")

    if len(available_scenarios) < 2:
        raise ValueError("At least 2 scenarios are required for comparison. "
                        f"Only {len(available_scenarios)} scenario(s) found.")

    print(f"\n  Available scenarios: {[s.key for s in available_scenarios]}")

    # =========================================================================

    # =========================================================================
    print("\n[1.5] Checking simulation data consistency...")
    consistency_check = check_simulation_data_consistency(df_aggs)
    
    for key, stat in consistency_check["stats"].items():
        print(f"  {key}: PersonMonths={stat['total_person_months']:,}, "
              f"Infections={stat['total_newly_infected']:,}, Deaths={stat['total_deaths']:,}")
    
    if consistency_check["warnings"]:
        print("\n  " + "!" * 60)
        for warning in consistency_check["warnings"]:
            print(f"  {warning}")
        print("  " + "!" * 60)

    available_keys = [s.key for s in available_scenarios]
    if config.BASELINE_KEY in available_keys:
        baseline_key = config.BASELINE_KEY
    elif "1m" in available_keys:
        baseline_key = "1m"
        print(f"  [INFO] Baseline changed to '1m' (original '2m' not available)")
    else:

        baseline_key = [k for k in available_keys if k != "0m"][0]
        print(f"  [INFO] Baseline set to '{baseline_key}'")

    additional_comparisons = [
        (int_key, comp_key)
        for int_key, comp_key in config.ADDITIONAL_COMPARISONS
        if int_key in available_keys and comp_key in available_keys
    ]

    final_additional = []
    for int_key, comp_key in additional_comparisons:

        if comp_key == baseline_key:
            continue
        final_additional.append((int_key, comp_key))
    additional_comparisons = final_additional

    print(f"  Baseline: {baseline_key}")
    if additional_comparisons:
        print(f"  Additional comparisons: {additional_comparisons}")
    
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    print("\n[2] Running base case analysis...")

    base_params = config.params
    base_costs = {age: p.base_value for age, p in base_params.costs_on_art.items()}
    base_utilities = {cd4: p.base_value for cd4, p in base_params.utilities_by_cd4.items()}

    base_results = {}
    for scen in available_scenarios:
        res = compute_qalys_and_costs(
            df_agg=df_aggs[scen.key],
            scen=scen,
            costs_on_art=base_costs,
            off_art_cost_ratio=base_params.off_art_cost_ratio.base_value,
            utility_hiv_neg=base_params.utility_hiv_neg.base_value,
            utilities_by_cd4=base_utilities,
            mental_health_util_loss=base_params.mental_health_util_loss.base_value,
            discount_rate_cost=base_params.discount_rate_cost.base_value,
            discount_rate_qaly=base_params.discount_rate_qaly.base_value,
            use_discounting=config.USE_DISCOUNTING,
            use_half_cycle_correction=config.USE_HALF_CYCLE_CORRECTION,
            n_years=config.N_YEARS,
            utility_hiv_neg_by_age=dict(base_params.utility_hiv_neg_by_age),
        )
        base_results[scen.key] = res
        print(f"  {scen.key}: QALY={res['total_qaly_disc']:,.0f}, Cost={res['total_cost_disc']:,.0f}")

    print("\n  Incremental results (vs baseline):")
    baseline_res = base_results[baseline_key]
    for key, res in base_results.items():
        if key == baseline_key:
            continue
        inc = compute_incremental(res, baseline_res)
        print(f"    {key} vs {baseline_key}: ΔCost={inc['delta_cost']:,.0f}, "
              f"ΔQALY={inc['delta_qaly']:,.2f}, {inc['quadrant']}")

    if additional_comparisons:
        print("\n  Additional comparisons:")
        for intervention_key, comparator_key in additional_comparisons:
            inc = compute_incremental(base_results[intervention_key], base_results[comparator_key])
            print(f"    {intervention_key} vs {comparator_key}: ΔCost={inc['delta_cost']:,.0f}, "
                  f"ΔQALY={inc['delta_qaly']:,.2f}, {inc['quadrant']}")
    
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    print("\n[3] Running Probabilistic Sensitivity Analysis (PSA)...")

    psa_results = run_psa(
        df_aggs, available_scenarios, baseline_key,
        additional_comparisons, config
    )
    psa_path = os.path.join(config.OUTPUT_DIR, "psa_results.csv")
    psa_results.to_csv(psa_path, index=False)
    print(f"  Saved PSA results to: {psa_path}")

    print("\n  PSA Summary:")
    for comparison in psa_results["comparison"].unique():
        df_comp = psa_results[psa_results["comparison"] == comparison]
        prob_dominant = ((df_comp["delta_cost"] < 0) & (df_comp["delta_qaly"] > 0)).mean()
        prob_ce_5m = (df_comp["delta_qaly"] * 5_000_000 - df_comp["delta_cost"] > 0).mean()
        print(f"    {comparison}: Prob(Dominant)={prob_dominant:.1%}, "
              f"Prob(CE at WTP=5M)={prob_ce_5m:.1%}")
    
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    print("\n[4] Running Deterministic Sensitivity Analysis (DSA)...")

    dsa_results_list = []

    dsa_results_main = run_dsa(df_aggs, available_scenarios, baseline_key, "0m", config)
    dsa_results_main["comparison"] = f"0m_vs_{baseline_key}"
    dsa_results_list.append(dsa_results_main)
    print(f"  DSA comparison: 0m vs {baseline_key}")

    for intervention_key, comparator_key in additional_comparisons:
        dsa_results_add = run_dsa(df_aggs, available_scenarios, comparator_key, intervention_key, config)
        dsa_results_add["comparison"] = f"{intervention_key}_vs_{comparator_key}"
        dsa_results_list.append(dsa_results_add)
        print(f"  DSA comparison: {intervention_key} vs {comparator_key}")

    dsa_results = pd.concat(dsa_results_list, ignore_index=True)
    dsa_path = os.path.join(config.OUTPUT_DIR, "dsa_results.csv")
    dsa_results.to_csv(dsa_path, index=False)
    print(f"  Saved DSA results to: {dsa_path}")
    
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    print("\n[5] Computing Cost-Effectiveness Acceptability Curve (CEAC)...")
    
    ceac_results = compute_ceac(psa_results, config.WTP_THRESHOLDS)
    ceac_path = os.path.join(config.OUTPUT_DIR, "ceac_results.csv")
    ceac_results.to_csv(ceac_path, index=False)
    print(f"  Saved CEAC results to: {ceac_path}")
    
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    print("\n[6] Running model validation...")
    
    validation = validate_model(base_results, df_aggs)

    validation_summary = pd.DataFrame(validation["summary"])
    validation_path = os.path.join(config.OUTPUT_DIR, "validation_summary.csv")
    validation_summary.to_csv(validation_path, index=False)
    print(f"  Saved validation summary to: {validation_path}")
    
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    print("\n[7] Generating parameter table (CHEERS 2.0 compliant)...")
    
    param_table = config.params.generate_parameter_table()
    param_path = os.path.join(config.OUTPUT_DIR, "parameter_table.csv")
    param_table.to_csv(param_path, index=False)
    print(f"  Saved parameter table to: {param_path}")
    
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    print("\n[8] Generating summary tables...")

    summary_tables = generate_summary_tables(
        base_results, psa_results, dsa_results, baseline_key,
        additional_comparisons, config
    )
    
    for name, df in summary_tables.items():
        path = os.path.join(config.OUTPUT_DIR, f"{name}.csv")
        df.to_csv(path, index=False)
        print(f"  Saved {name} to: {path}")
    
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    if config.MAKE_FIGURES:
        print("\n[9] Generating figures...")

        tornado_path = os.path.join(config.OUTPUT_DIR, "fig_tornado.png")
        plot_tornado_diagram(dsa_results, tornado_path, config.FIGURE_DPI)
        print(f"  Saved tornado diagram to: {tornado_path}")
        
        # CE plane (PSA)
        ce_plane_path = os.path.join(config.OUTPUT_DIR, "fig_ce_plane_psa.png")
        plot_ce_plane_psa(psa_results, ce_plane_path, config.FIGURE_DPI)
        print(f"  Saved CE plane to: {ce_plane_path}")
        
        # CEAC
        ceac_fig_path = os.path.join(config.OUTPUT_DIR, "fig_ceac.png")
        plot_ceac(ceac_results, ceac_fig_path, config.FIGURE_DPI)
        print(f"  Saved CEAC to: {ceac_fig_path}")

        hist_path = os.path.join(config.OUTPUT_DIR, "fig_psa_histogram.png")
        plot_psa_histogram(psa_results, hist_path, config.FIGURE_DPI)
        print(f"  Saved PSA histogram to: {hist_path}")
    
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    print("\n[10] Exporting report-ready JSON...")
    export_cea_report_json(
        base_results=base_results,
        psa_results=psa_results,
        dsa_results=dsa_results,
        ceac_results=ceac_results,
        validation=validation,
        baseline_key=baseline_key,
        additional_comparisons=additional_comparisons,
        config=config,
        output_path=os.path.join(config.OUTPUT_DIR, "cea_report_ready.json"),
    )

    print("\n" + "=" * 70)
    print("Analysis complete!")
    print(f"All results saved to: {config.OUTPUT_DIR}")
    print("="*70)

if __name__ == "__main__":
    main()

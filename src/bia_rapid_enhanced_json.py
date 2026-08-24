#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Budget-impact analysis of immediate versus delayed ART."""

import os
import json
import datetime
import warnings
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy import stats

import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# --- Data location: set BASE_144 to your local 144_scenario folder ---
BASE_144 = r"/path/to/20250315_169/144_scenario"
PROJ_169 = os.path.dirname(BASE_144)
_PROJ_144 = BASE_144

plt.rcParams['font.family'] = ['DejaVu Sans', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

warnings.filterwarnings('ignore', category=FutureWarning)

# ============================================================

# ============================================================

@dataclass
class Config:

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------

    INPUT_DIR: str = os.path.join(_PROJ_144, "cea_param_uncertainty", "mean_input")
    OUTPUT_DIR: str = os.path.join(_PROJ_144, "bia_results")

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    RANDOM_SEED: int = 12345

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    N_YEARS: int = 40
    N_PSA_ITERATIONS: int = 1000

    SHORT_TERM_YEARS: List[int] = field(default_factory=lambda: [3, 5, 10])

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    SCENARIOS: Dict[str, Dict[str, str]] = field(default_factory=lambda: {
        "0 month": {
            "agg_file": "scenario0_yearly_summary.csv",
            "short_name": "0-month (Rapid ART)",
            "description": "Immediate ART start (no CD4/VL restriction)",
        },
        "1 month": {
            "agg_file": "scenario1_yearly_summary.csv",
            "short_name": "1-month ART",
            "description": "ART start after 1 month (CD4/VL restriction applies)",
        },
        "2 months": {
            "agg_file": "scenario2_yearly_summary.csv",
            "short_name": "2-month ART",
            "description": "ART start after 2 months (CD4/VL restriction applies; current policy)",
        },
    })

    INTERVENTION_LABEL: str = "0 month"

    COMPARATOR_LABELS: List[str] = field(default_factory=lambda: ["1 month", "2 months"])

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    PREP_MONTHLY_COST: float = 0.0
    HIV_TEST_UNIT_COST: float = 3_000.0

    MEAN_ANNUAL_COST_ON_ART: Dict[str, float] = field(default_factory=lambda: {
        "20-29": 2_375_000.0,
        "30-39": 2_901_449.0,
        "40-49": 3_264_767.0,
        "50-59": 3_433_598.0,
        "60-69": 3_356_900.0,
        "70+":   3_164_654.0,
    })

    OFF_ART_COST_RATIO: float = 0.3

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------

    COST_SE_RATIO: float = 0.20

    ART_COVERAGE_ALPHA: float = 90.0
    ART_COVERAGE_BETA: float = 10.0

    OFF_ART_RATIO_ALPHA: float = 30.0
    OFF_ART_RATIO_BETA: float = 70.0

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    HIV_LIFE_EXPECTANCY: Dict[str, float] = field(default_factory=lambda: {
        "20-29": 45.0,
        "30-39": 38.0,
        "40-49": 30.0,
        "50-59": 22.0,
        "60-69": 15.0,
        "70+":   10.0,
    })

    LIFETIME_ART_COVERAGE: float = 0.90

    DISCOUNT_RATE: float = 0.02
    USE_DISCOUNTING: bool = True

    DISCOUNT_RATES_SA: List[float] = field(default_factory=lambda: [0.0, 0.02, 0.04])

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------

    OWSA_RANGE_RATIO: float = 0.25  # ±25%

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------

    PHASED_IMPLEMENTATION: Dict[str, List[float]] = field(default_factory=lambda: {
        "Immediate (100%)": [1.0] * 40,
        "Gradual (3-year)": [0.3, 0.6, 0.9] + [1.0] * 37,
        "Gradual (5-year)": [0.2, 0.4, 0.6, 0.8, 1.0] + [1.0] * 35,
        "Conservative": [0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 1.0] + [1.0] * 33,
    })

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    WTP_THRESHOLDS: List[float] = field(default_factory=lambda: 
        list(np.arange(0, 100_000_001, 5_000_000)))

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    STRICT_MODE: bool = True
    CURRENCY_UNIT: str = "JPY"
    COST_DISPLAY_UNIT: float = 1e6

# ============================================================

# ============================================================

def set_random_seed(seed: int) -> None:
    np.random.seed(seed)

def discount_factor(year: int, rate: float, use_discounting: bool) -> float:
    if (not use_discounting) or rate <= 0.0:
        return 1.0
    return 1.0 / ((1.0 + rate) ** (year - 1))

def gamma_params_from_mean_se(mean: float, se: float) -> Tuple[float, float]:
    if se <= 0 or mean <= 0:
        return (mean, 1.0)
    shape = (mean / se) ** 2
    scale = (se ** 2) / mean
    return (shape, scale)

def format_currency(value: float, unit: float = 1e6, decimals: int = 1) -> str:
    return f"{value / unit:,.{decimals}f}"

def load_aggregator(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Aggregator file not found: {path}")

    df = pd.read_csv(path)

    required_cols = [
        "Year", "AgeGroup", "CD4Category",
        "OnART", "OnPrEP", "PersonMonthsAlive"
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in aggregator: {missing}")

    if df["OnART"].dtype != bool:
        df["OnART"] = df["OnART"].astype(bool)
    if df["OnPrEP"].dtype != bool:
        df["OnPrEP"] = df["OnPrEP"].astype(bool)

    return df

# ============================================================

# ============================================================

def make_cost_dicts_from_means(
    mean_on_art: Dict[str, float],
    off_art_ratio: float
) -> Tuple[Dict[str, float], Dict[str, float]]:
    cost_on_art = dict(mean_on_art)
    cost_off_art = {age: val * off_art_ratio for age, val in mean_on_art.items()}
    return cost_on_art, cost_off_art

def compute_lifetime_cost_by_age(
    annual_cost: Dict[str, float],
    life_expectancy: Dict[str, float],
    art_coverage: float,
    discount_rate: float,
    use_discounting: bool
) -> Dict[str, float]:
    lifetime_cost = {}
    for age_group, annual in annual_cost.items():
        years = life_expectancy.get(age_group, 0.0)
        if years <= 0.0 or annual <= 0.0:
            lifetime_cost[age_group] = 0.0
            continue

        total = 0.0
        for y in range(1, int(years) + 1):
            disc = discount_factor(y, discount_rate, use_discounting)
            total += annual * disc

        lifetime_cost[age_group] = total * art_coverage

    return lifetime_cost

def compute_yearly_costs_from_aggregator(
    df_agg: pd.DataFrame,
    cost_on_art: Dict[str, float],
    cost_off_art: Dict[str, float],
    config: Config,
    discount_rate: Optional[float] = None
) -> Dict[str, np.ndarray]:
    if discount_rate is None:
        discount_rate = config.DISCOUNT_RATE

    years = np.arange(1, config.N_YEARS + 1)
    disc = np.array([
        discount_factor(int(y), discount_rate, config.USE_DISCOUNTING)
        for y in years
    ], dtype=float)

    if "NewlyInfected" in df_agg.columns:
        df_hiv = df_agg[df_agg["CD4Category"] != "HIV-"].copy()
        inf_by_year = df_hiv.groupby("Year")["NewlyInfected"].sum()
        yearly_new_infections = inf_by_year.reindex(years).fillna(0.0).values.astype(float)
    else:
        yearly_new_infections = np.zeros(len(years), dtype=float)

    if "NewlyDiagnosed" in df_agg.columns:
        df_hiv = df_agg[df_agg["CD4Category"] != "HIV-"].copy()
        diag_by_year = df_hiv.groupby("Year")["NewlyDiagnosed"].sum()
        yearly_new_diagnosed = diag_by_year.reindex(years).fillna(0.0).values.astype(float)
    else:
        yearly_new_diagnosed = np.zeros(len(years), dtype=float)

    df_care = df_agg[df_agg["CD4Category"] != "HIV-"].copy()
    if len(df_care) == 0:
        yearly_hivcare_cost = np.zeros(len(years), dtype=float)
    else:
        df_care["CostOnArt"] = df_care["AgeGroup"].map(cost_on_art)
        df_care["CostOffArt"] = df_care["AgeGroup"].map(cost_off_art)
        mean_cost = np.mean(list(cost_on_art.values()))
        df_care["CostOnArt"] = df_care["CostOnArt"].fillna(mean_cost)
        df_care["CostOffArt"] = df_care["CostOffArt"].fillna(
            mean_cost * config.OFF_ART_COST_RATIO
        )

        df_care["PersonYears"] = df_care["PersonMonthsAlive"] / 12.0
        df_care["AnnualCostPerPY"] = np.where(
            df_care["OnART"],
            df_care["CostOnArt"],
            df_care["CostOffArt"],
        )
        df_care["Cost"] = df_care["AnnualCostPerPY"] * df_care["PersonYears"]

        care_by_year = df_care.groupby("Year")["Cost"].sum()
        yearly_hivcare_cost = care_by_year.reindex(years).fillna(0.0).values.astype(float)

    yearly_test_cost = yearly_new_diagnosed * config.HIV_TEST_UNIT_COST * 10

    yearly_total_cost = yearly_hivcare_cost + yearly_test_cost

    yearly_hivcare_cost_disc = yearly_hivcare_cost * disc
    yearly_test_cost_disc = yearly_test_cost * disc
    yearly_total_cost_disc = yearly_total_cost * disc

    return {
        "years": years,
        "new_infections": yearly_new_infections,
        "new_diagnosed": yearly_new_diagnosed,
        "hivcare_cost": yearly_hivcare_cost,
        "test_cost": yearly_test_cost,
        "total_cost": yearly_total_cost,
        "hivcare_cost_disc": yearly_hivcare_cost_disc,
        "test_cost_disc": yearly_test_cost_disc,
        "total_cost_disc": yearly_total_cost_disc,
        "discount_factor": disc,
    }

def get_infections_by_age(df_agg: pd.DataFrame, n_years: int) -> pd.DataFrame:
    if "NewlyInfected" not in df_agg.columns:
        return pd.DataFrame()

    df_hiv = df_agg[df_agg["CD4Category"] != "HIV-"].copy()
    grouped = (
        df_hiv.groupby(["Year", "AgeGroup"])["NewlyInfected"]
        .sum()
        .reset_index()
        .rename(columns={"NewlyInfected": "NewInfections"})
    )
    return grouped

# ============================================================

# ============================================================

def run_deterministic_analysis(
    aggregators: Dict[str, pd.DataFrame],
    config: Config,
    discount_rate: Optional[float] = None
) -> Dict[str, Any]:
    if discount_rate is None:
        discount_rate = config.DISCOUNT_RATE

    cost_on_art, cost_off_art = make_cost_dicts_from_means(
        config.MEAN_ANNUAL_COST_ON_ART,
        config.OFF_ART_COST_RATIO,
    )

    results = {}
    for label, df_agg in aggregators.items():
        res = compute_yearly_costs_from_aggregator(
            df_agg, cost_on_art, cost_off_art, config, discount_rate
        )

        res["cum_infections"] = float(np.sum(res["new_infections"]))
        res["cum_hivcare_cost"] = float(np.sum(res["hivcare_cost"]))
        res["cum_total_cost"] = float(np.sum(res["total_cost"]))
        res["cum_hivcare_cost_disc"] = float(np.sum(res["hivcare_cost_disc"]))
        res["cum_total_cost_disc"] = float(np.sum(res["total_cost_disc"]))

        res["infections_by_age"] = get_infections_by_age(df_agg, config.N_YEARS)

        results[label] = res

    return results

def compute_incremental_results(
    results: Dict[str, Any],
    intervention_label: str,
    comparator_label: str,
    config: Config
) -> Dict[str, float]:
    res_int = results[intervention_label]
    res_comp = results[comparator_label]

    I_int = res_int["cum_infections"]
    I_comp = res_comp["cum_infections"]
    infections_averted = I_comp - I_int

    C_int = res_int["cum_total_cost_disc"]
    C_comp = res_comp["cum_total_cost_disc"]
    delta_cost = C_int - C_comp

    C_care_int = res_int["cum_hivcare_cost_disc"]
    C_care_comp = res_comp["cum_hivcare_cost_disc"]
    delta_care_cost = C_care_int - C_care_comp

    if infections_averted > 0:
        cost_per_infection_averted = delta_cost / infections_averted
    else:
        cost_per_infection_averted = np.nan

    lifetime_cost_by_age = compute_lifetime_cost_by_age(
        config.MEAN_ANNUAL_COST_ON_ART,
        config.HIV_LIFE_EXPECTANCY,
        config.LIFETIME_ART_COVERAGE,
        config.DISCOUNT_RATE,
        config.USE_DISCOUNTING,
    )

    inf_age_int = res_int["infections_by_age"]
    inf_age_comp = res_comp["infections_by_age"]

    if len(inf_age_int) > 0 and len(inf_age_comp) > 0:
        merged = inf_age_comp.merge(
            inf_age_int, on=["Year", "AgeGroup"], suffixes=("_comp", "_int")
        )
        merged["Averted"] = merged["NewInfections_comp"] - merged["NewInfections_int"]
        merged["LifetimeCost"] = merged["AgeGroup"].map(lifetime_cost_by_age)
        mean_lifetime = np.mean(list(lifetime_cost_by_age.values()))
        merged["LifetimeCost"] = merged["LifetimeCost"].fillna(mean_lifetime)
        merged["DiscFactor"] = merged["Year"].apply(
            lambda y: discount_factor(int(y), config.DISCOUNT_RATE, config.USE_DISCOUNTING)
        )
        merged["Savings"] = merged["Averted"] * merged["LifetimeCost"] * merged["DiscFactor"]

        lifetime_savings_disc = float(merged["Savings"].sum())
    else:
        lifetime_savings_disc = infections_averted * np.mean(list(lifetime_cost_by_age.values()))

    # Break-even year
    total_int = res_int["total_cost_disc"]
    total_comp = res_comp["total_cost_disc"]
    diff_yearly = total_int - total_comp
    cum_diff = np.cumsum(diff_yearly)
    idx = np.where(cum_diff <= 0)[0]
    break_even_year = int(idx[0] + 1) if len(idx) > 0 else -1

    return {
        "intervention": intervention_label,
        "comparator": comparator_label,
        "infections_intervention": I_int,
        "infections_comparator": I_comp,
        "infections_averted": infections_averted,
        "cost_intervention_disc": C_int,
        "cost_comparator_disc": C_comp,
        "delta_cost_disc": delta_cost,
        "delta_care_cost_disc": delta_care_cost,
        "cost_per_infection_averted": cost_per_infection_averted,
        "lifetime_savings_disc": lifetime_savings_disc,
        "net_benefit": lifetime_savings_disc + delta_care_cost,
        "break_even_year": break_even_year,
    }

# ============================================================

# ============================================================

def sample_psa_parameters(config: Config) -> Dict[str, Any]:
    params = {}

    sampled_costs = {}
    for age, mean in config.MEAN_ANNUAL_COST_ON_ART.items():
        se = mean * config.COST_SE_RATIO
        shape, scale = gamma_params_from_mean_se(mean, se)
        sampled_costs[age] = float(np.random.gamma(shape, scale))
    params["annual_cost_on_art"] = sampled_costs

    params["art_coverage"] = float(np.random.beta(
        config.ART_COVERAGE_ALPHA,
        config.ART_COVERAGE_BETA
    ))

    params["off_art_ratio"] = float(np.random.beta(
        config.OFF_ART_RATIO_ALPHA,
        config.OFF_ART_RATIO_BETA
    ))

    return params

def run_single_psa_iteration(
    aggregators: Dict[str, pd.DataFrame],
    config: Config,
    psa_params: Dict[str, Any]
) -> Dict[str, Dict[str, float]]:
    cost_on_art, cost_off_art = make_cost_dicts_from_means(
        psa_params["annual_cost_on_art"],
        psa_params["off_art_ratio"],
    )

    scenario_results = {}
    for label, df_agg in aggregators.items():
        res = compute_yearly_costs_from_aggregator(
            df_agg, cost_on_art, cost_off_art, config
        )
        scenario_results[label] = {
            "cum_infections": float(np.sum(res["new_infections"])),
            "cum_total_cost_disc": float(np.sum(res["total_cost_disc"])),
            "cum_hivcare_cost_disc": float(np.sum(res["hivcare_cost_disc"])),
        }

    return scenario_results

def run_psa(
    aggregators: Dict[str, pd.DataFrame],
    config: Config,
    n_iterations: Optional[int] = None
) -> pd.DataFrame:
    if n_iterations is None:
        n_iterations = config.N_PSA_ITERATIONS

    set_random_seed(config.RANDOM_SEED)

    results_list = []
    intervention = config.INTERVENTION_LABEL

    print(f"  Running {n_iterations} PSA iterations...")

    for i in range(n_iterations):
        if (i + 1) % 100 == 0:
            print(f"    Iteration {i + 1}/{n_iterations}")

        psa_params = sample_psa_parameters(config)
        scenario_results = run_single_psa_iteration(aggregators, config, psa_params)

        for comparator in config.COMPARATOR_LABELS:
            res_int = scenario_results[intervention]
            res_comp = scenario_results[comparator]

            infections_averted = res_comp["cum_infections"] - res_int["cum_infections"]
            delta_cost = res_int["cum_total_cost_disc"] - res_comp["cum_total_cost_disc"]
            delta_care = res_int["cum_hivcare_cost_disc"] - res_comp["cum_hivcare_cost_disc"]

            if infections_averted > 0:
                cost_per_averted = delta_cost / infections_averted
            else:
                cost_per_averted = np.nan

            results_list.append({
                "iteration": i + 1,
                "comparator": comparator,
                "infections_averted": infections_averted,
                "delta_cost_disc": delta_cost,
                "delta_care_cost_disc": delta_care,
                "cost_per_infection_averted": cost_per_averted,
            })

    return pd.DataFrame(results_list)

def compute_psa_summary(psa_results: pd.DataFrame) -> pd.DataFrame:
    summary_list = []

    for comparator in psa_results["comparator"].unique():
        df_comp = psa_results[psa_results["comparator"] == comparator]

        for metric in ["infections_averted", "delta_cost_disc", "cost_per_infection_averted"]:
            values = df_comp[metric].dropna()
            if len(values) == 0:
                continue

            summary_list.append({
                "comparator": comparator,
                "metric": metric,
                "mean": values.mean(),
                "median": values.median(),
                "std": values.std(),
                "ci_lower": np.percentile(values, 2.5),
                "ci_upper": np.percentile(values, 97.5),
                "p5": np.percentile(values, 5),
                "p95": np.percentile(values, 95),
            })

    return pd.DataFrame(summary_list)

def compute_ceac(
    psa_results: pd.DataFrame,
    wtp_thresholds: List[float]
) -> pd.DataFrame:
    ceac_list = []

    for comparator in psa_results["comparator"].unique():
        df_comp = psa_results[psa_results["comparator"] == comparator]

        for wtp in wtp_thresholds:
            # Net Monetary Benefit = WTP * ΔEffect - ΔCost

            nmb = wtp * df_comp["infections_averted"] - df_comp["delta_cost_disc"]
            prob_ce = (nmb > 0).mean()

            ceac_list.append({
                "comparator": comparator,
                "wtp": wtp,
                "prob_cost_effective": prob_ce,
            })

    return pd.DataFrame(ceac_list)

# ============================================================

# ============================================================

def run_owsa(
    aggregators: Dict[str, pd.DataFrame],
    config: Config
) -> pd.DataFrame:
    results_list = []
    intervention = config.INTERVENTION_LABEL

    base_results = run_deterministic_analysis(aggregators, config)

    parameters = []

    for age in config.MEAN_ANNUAL_COST_ON_ART.keys():
        parameters.append({
            "name": f"Annual ART cost ({age})",
            "type": "cost",
            "age": age,
            "base": config.MEAN_ANNUAL_COST_ON_ART[age],
        })

    # 2. OFF_ART_COST_RATIO
    parameters.append({
        "name": "Off-ART cost ratio",
        "type": "off_art_ratio",
        "base": config.OFF_ART_COST_RATIO,
    })

    # 3. ART coverage
    parameters.append({
        "name": "ART coverage",
        "type": "art_coverage",
        "base": config.LIFETIME_ART_COVERAGE,
    })

    parameters.append({
        "name": "Discount rate",
        "type": "discount_rate",
        "base": config.DISCOUNT_RATE,
    })

    for param in parameters:
        base_val = param["base"]
        low_val = base_val * (1 - config.OWSA_RANGE_RATIO)
        high_val = base_val * (1 + config.OWSA_RANGE_RATIO)

        for val, scenario in [(low_val, "low"), (high_val, "high")]:

            modified_config = Config()
            modified_config.INPUT_DIR = config.INPUT_DIR
            modified_config.OUTPUT_DIR = config.OUTPUT_DIR
            modified_config.N_YEARS = config.N_YEARS

            if param["type"] == "cost":
                new_costs = dict(config.MEAN_ANNUAL_COST_ON_ART)
                new_costs[param["age"]] = val
                modified_config.MEAN_ANNUAL_COST_ON_ART = new_costs
            elif param["type"] == "off_art_ratio":
                modified_config.OFF_ART_COST_RATIO = val
            elif param["type"] == "art_coverage":
                modified_config.LIFETIME_ART_COVERAGE = val
            elif param["type"] == "discount_rate":
                modified_config.DISCOUNT_RATE = max(0, val)

            try:
                mod_results = run_deterministic_analysis(
                    aggregators, modified_config,
                    discount_rate=modified_config.DISCOUNT_RATE
                )

                for comparator in config.COMPARATOR_LABELS:
                    inc_base = compute_incremental_results(
                        base_results, intervention, comparator, config
                    )
                    inc_mod = compute_incremental_results(
                        mod_results, intervention, comparator, modified_config
                    )

                    results_list.append({
                        "parameter": param["name"],
                        "comparator": comparator,
                        "scenario": scenario,
                        "value": val,
                        "base_value": base_val,
                        "infections_averted": inc_mod["infections_averted"],
                        "delta_cost_disc": inc_mod["delta_cost_disc"],
                        "cost_per_averted": inc_mod["cost_per_infection_averted"],
                        "base_infections_averted": inc_base["infections_averted"],
                        "base_delta_cost": inc_base["delta_cost_disc"],
                        "base_cost_per_averted": inc_base["cost_per_infection_averted"],
                    })
            except Exception as e:
                print(f"Warning: OWSA failed for {param['name']} ({scenario}): {e}")

    return pd.DataFrame(results_list)

# ============================================================

# ============================================================

def run_short_term_bia(
    results: Dict[str, Any],
    config: Config
) -> pd.DataFrame:
    bia_list = []
    intervention = config.INTERVENTION_LABEL

    for horizon in config.SHORT_TERM_YEARS:
        for comparator in config.COMPARATOR_LABELS:
            res_int = results[intervention]
            res_comp = results[comparator]

            idx_end = min(horizon, len(res_int["years"]))

            infections_int = float(np.sum(res_int["new_infections"][:idx_end]))
            infections_comp = float(np.sum(res_comp["new_infections"][:idx_end]))
            infections_averted = infections_comp - infections_int

            cost_int = float(np.sum(res_int["total_cost_disc"][:idx_end]))
            cost_comp = float(np.sum(res_comp["total_cost_disc"][:idx_end]))
            delta_cost = cost_int - cost_comp

            care_int = float(np.sum(res_int["hivcare_cost_disc"][:idx_end]))
            care_comp = float(np.sum(res_comp["hivcare_cost_disc"][:idx_end]))
            delta_care = care_int - care_comp

            bia_list.append({
                "time_horizon_years": horizon,
                "comparator": comparator,
                "cum_infections_intervention": infections_int,
                "cum_infections_comparator": infections_comp,
                "infections_averted": infections_averted,
                "cum_cost_intervention_disc": cost_int,
                "cum_cost_comparator_disc": cost_comp,
                "delta_total_cost_disc": delta_cost,
                "delta_care_cost_disc": delta_care,
                "annual_budget_impact": delta_cost / horizon,
            })

    return pd.DataFrame(bia_list)

# ============================================================

# ============================================================

def run_phased_implementation(
    results: Dict[str, Any],
    config: Config
) -> pd.DataFrame:
    phased_list = []
    intervention = config.INTERVENTION_LABEL

    res_int = results[intervention]

    for comparator in config.COMPARATOR_LABELS:
        res_comp = results[comparator]

        for impl_name, impl_rates in config.PHASED_IMPLEMENTATION.items():

            years = res_int["years"]
            n_years = len(years)

            inf_averted_yearly = res_comp["new_infections"] - res_int["new_infections"]

            cost_diff_yearly = res_int["total_cost_disc"] - res_comp["total_cost_disc"]

            weighted_inf_averted = 0.0
            weighted_cost_diff = 0.0

            for i, year in enumerate(years):
                if i < len(impl_rates):
                    rate = impl_rates[i]
                else:
                    rate = impl_rates[-1]

                weighted_inf_averted += inf_averted_yearly[i] * rate
                weighted_cost_diff += cost_diff_yearly[i] * rate

            phased_list.append({
                "comparator": comparator,
                "implementation_scenario": impl_name,
                "total_infections_averted": weighted_inf_averted,
                "total_delta_cost_disc": weighted_cost_diff,
                "year_1_adoption": impl_rates[0] if impl_rates else 0,
                "year_5_adoption": impl_rates[4] if len(impl_rates) > 4 else impl_rates[-1],
            })

    return pd.DataFrame(phased_list)

# ============================================================

# ============================================================

def plot_tornado_diagram(
    owsa_results: pd.DataFrame,
    config: Config,
    metric: str = "delta_cost_disc",
    output_path: Optional[str] = None
) -> None:
    for comparator in config.COMPARATOR_LABELS:
        df_comp = owsa_results[owsa_results["comparator"] == comparator].copy()

        if len(df_comp) == 0:
            continue

        base_metric = f"base_{metric.replace('delta_', '')}" if "delta" in metric else f"base_{metric}"
        if base_metric not in df_comp.columns:
            base_val = df_comp[metric].mean()
        else:
            base_val = df_comp[base_metric].iloc[0]

        param_ranges = []
        for param in df_comp["parameter"].unique():
            df_param = df_comp[df_comp["parameter"] == param]
            low_val = df_param[df_param["scenario"] == "low"][metric].values
            high_val = df_param[df_param["scenario"] == "high"][metric].values

            if len(low_val) > 0 and len(high_val) > 0:
                low_val = low_val[0]
                high_val = high_val[0]
                range_width = abs(high_val - low_val)
                param_ranges.append({
                    "parameter": param,
                    "low": low_val,
                    "high": high_val,
                    "range": range_width,
                })

        if not param_ranges:
            continue

        df_ranges = pd.DataFrame(param_ranges).sort_values("range", ascending=True)

        fig, ax = plt.subplots(figsize=(10, max(6, len(df_ranges) * 0.4)))

        y_pos = np.arange(len(df_ranges))
        colors = ['#2ecc71', '#e74c3c']  # green for low, red for high

        for i, (_, row) in enumerate(df_ranges.iterrows()):
            low_diff = row["low"] - base_val
            high_diff = row["high"] - base_val

            left = min(low_diff, high_diff, 0)
            right = max(low_diff, high_diff, 0)

            if low_diff < high_diff:
                ax.barh(i, low_diff, color=colors[0], height=0.6, left=0)
                ax.barh(i, high_diff, color=colors[1], height=0.6, left=0)
            else:
                ax.barh(i, high_diff, color=colors[1], height=0.6, left=0)
                ax.barh(i, low_diff, color=colors[0], height=0.6, left=0)

        ax.axvline(x=0, color='black', linewidth=1)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df_ranges["parameter"])

        unit = config.COST_DISPLAY_UNIT
        ax.set_xlabel(f"Change in {metric} (million JPY)")
        ax.set_title(f"Tornado Diagram: 0-month vs {comparator}")

        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=colors[0], label='Low value (-25%)'),
            Patch(facecolor=colors[1], label='High value (+25%)')
        ]
        ax.legend(handles=legend_elements, loc='lower right')

        plt.tight_layout()

        if output_path:
            save_path = output_path.replace(".png", f"_{comparator.replace(' ', '_')}.png")
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"    Saved: {save_path}")
        plt.close()

def plot_ceac(
    ceac_data: pd.DataFrame,
    config: Config,
    output_path: Optional[str] = None
) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = {'1 month': '#3498db', '2 months': '#e74c3c'}
    linestyles = {'1 month': '-', '2 months': '--'}

    for comparator in ceac_data["comparator"].unique():
        df_comp = ceac_data[ceac_data["comparator"] == comparator]
        ax.plot(
            df_comp["wtp"] / 1e6,
            df_comp["prob_cost_effective"],
            label=f"0-month vs {comparator}",
            color=colors.get(comparator, 'gray'),
            linestyle=linestyles.get(comparator, '-'),
            linewidth=2
        )

    ax.set_xlabel("Willingness-to-pay per infection averted (million JPY)")
    ax.set_ylabel("Probability of being cost-effective")
    ax.set_title("Cost-Effectiveness Acceptability Curve (CEAC)")
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, ceac_data["wtp"].max() / 1e6)
    ax.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5)
    ax.axhline(y=0.95, color='gray', linestyle=':', alpha=0.5)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"    Saved: {output_path}")
    plt.close()

def plot_cost_effectiveness_plane(
    psa_results: pd.DataFrame,
    config: Config,
    output_path: Optional[str] = None
) -> None:
    fig, ax = plt.subplots(figsize=(10, 8))

    colors = {'1 month': '#3498db', '2 months': '#e74c3c'}

    for comparator in psa_results["comparator"].unique():
        df_comp = psa_results[psa_results["comparator"] == comparator]
        ax.scatter(
            df_comp["infections_averted"],
            df_comp["delta_cost_disc"] / 1e9,
            alpha=0.3,
            s=10,
            label=f"0-month vs {comparator}",
            color=colors.get(comparator, 'gray')
        )

    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.axvline(x=0, color='black', linewidth=0.5)

    ax.set_xlabel("Infections averted (Rapid ART vs Comparator)")
    ax.set_ylabel("Incremental cost (billion JPY)")
    ax.set_title("Cost-Effectiveness Plane")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax.text(0.95, 0.95, "NE\n(More costly,\nMore effective)",
            transform=ax.transAxes, ha='right', va='top', fontsize=9, alpha=0.7)
    ax.text(0.05, 0.95, "NW\n(More costly,\nLess effective)",
            transform=ax.transAxes, ha='left', va='top', fontsize=9, alpha=0.7)
    ax.text(0.95, 0.05, "SE\n(Less costly,\nMore effective)\nDOMINANT",
            transform=ax.transAxes, ha='right', va='bottom', fontsize=9, alpha=0.7)
    ax.text(0.05, 0.05, "SW\n(Less costly,\nLess effective)",
            transform=ax.transAxes, ha='left', va='bottom', fontsize=9, alpha=0.7)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"    Saved: {output_path}")
    plt.close()

def plot_cumulative_budget_impact(
    results: Dict[str, Any],
    config: Config,
    output_path: Optional[str] = None
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    intervention = config.INTERVENTION_LABEL
    res_int = results[intervention]
    years = res_int["years"]

    colors = {'1 month': '#3498db', '2 months': '#e74c3c'}

    ax = axes[0]
    for comparator in config.COMPARATOR_LABELS:
        res_comp = results[comparator]
        cum_diff = np.cumsum(res_int["total_cost_disc"] - res_comp["total_cost_disc"]) / 1e9
        ax.plot(years, cum_diff, label=f"vs {comparator}",
                color=colors.get(comparator, 'gray'), linewidth=2)

    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_xlabel("Year")
    ax.set_ylabel("Cumulative cost difference (billion JPY)")
    ax.set_title("Cumulative Budget Impact: 0-month ART vs Comparators")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    for comparator in config.COMPARATOR_LABELS:
        res_comp = results[comparator]
        cum_averted = np.cumsum(res_comp["new_infections"] - res_int["new_infections"])
        ax.plot(years, cum_averted, label=f"vs {comparator}",
                color=colors.get(comparator, 'gray'), linewidth=2)

    ax.set_xlabel("Year")
    ax.set_ylabel("Cumulative infections averted")
    ax.set_title("Cumulative Infections Averted: 0-month ART vs Comparators")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"    Saved: {output_path}")
    plt.close()

def plot_yearly_comparison(
    results: Dict[str, Any],
    config: Config,
    output_path: Optional[str] = None
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    colors = {
        '0 month': '#27ae60',
        '1 month': '#3498db',
        '2 months': '#e74c3c'
    }

    ax = axes[0]
    for label, res in results.items():
        ax.plot(res["years"], res["new_infections"],
                label=config.SCENARIOS[label]["short_name"],
                color=colors.get(label, 'gray'), linewidth=2)
    ax.set_xlabel("Year")
    ax.set_ylabel("New HIV infections")
    ax.set_title("Annual New HIV Infections")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    for label, res in results.items():
        ax.plot(res["years"], res["hivcare_cost"] / 1e9,
                label=config.SCENARIOS[label]["short_name"],
                color=colors.get(label, 'gray'), linewidth=2)
    ax.set_xlabel("Year")
    ax.set_ylabel("HIV care cost (billion JPY)")
    ax.set_title("Annual HIV Care Costs")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    for label, res in results.items():
        ax.plot(res["years"], res["total_cost_disc"] / 1e9,
                label=config.SCENARIOS[label]["short_name"],
                color=colors.get(label, 'gray'), linewidth=2)
    ax.set_xlabel("Year")
    ax.set_ylabel("Total cost, discounted (billion JPY)")
    ax.set_title("Annual Total Costs (Discounted)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"    Saved: {output_path}")
    plt.close()

def plot_short_term_bia_bars(
    bia_df: pd.DataFrame,
    config: Config,
    output_path: Optional[str] = None
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    horizons = bia_df["time_horizon_years"].unique()
    comparators = bia_df["comparator"].unique()

    x = np.arange(len(horizons))
    width = 0.35

    colors = {'1 month': '#3498db', '2 months': '#e74c3c'}

    ax = axes[0]
    for i, comp in enumerate(comparators):
        df_comp = bia_df[bia_df["comparator"] == comp]
        values = df_comp.set_index("time_horizon_years").loc[horizons, "infections_averted"].values
        ax.bar(x + i * width, values, width, label=f"vs {comp}", color=colors.get(comp, 'gray'))

    ax.set_xlabel("Time horizon (years)")
    ax.set_ylabel("Infections averted")
    ax.set_title("Short-term BIA: Infections Averted")
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(horizons)
    ax.legend()
    ax.set_axisbelow(True)
    ax.grid(True, alpha=0.3, axis='y')

    # Budget Impact
    ax = axes[1]
    for i, comp in enumerate(comparators):
        df_comp = bia_df[bia_df["comparator"] == comp]
        values = df_comp.set_index("time_horizon_years").loc[horizons, "delta_total_cost_disc"].values / 1e9
        ax.bar(x + i * width, values, width, label=f"vs {comp}", color=colors.get(comp, 'gray'))

    ax.set_xlabel("Time horizon (years)")
    ax.set_ylabel("Budget impact (billion JPY)")
    ax.set_title("Short-term BIA: Budget Impact (Discounted)")
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(horizons)
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.legend()
    ax.set_axisbelow(True)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"    Saved: {output_path}")
    plt.close()

# ============================================================

# ============================================================

# ============================================================

# ============================================================

def export_yearly_costs_to_csv(
    results: Dict[str, Any],
    config: Config,
    output_dir: str,
    format_type: str = "both"
) -> None:

    if format_type in ["long", "both"]:
        rows = []
        for label, res in results.items():
            scenario_name = config.SCENARIOS[label]["short_name"]
            years = res["years"]

            for i, year in enumerate(years):
                rows.append({
                    "Scenario": label,
                    "Scenario_Name": scenario_name,
                    "Year": int(year),
                    "New_Infections": res["new_infections"][i],
                    "New_Diagnosed": res["new_diagnosed"][i],
                    "HIVCare_Cost": res["hivcare_cost"][i],
                    "Test_Cost": res["test_cost"][i],
                    "Total_Cost": res["total_cost"][i],
                    "HIVCare_Cost_Discounted": res["hivcare_cost_disc"][i],
                    "Test_Cost_Discounted": res["test_cost_disc"][i],
                    "Total_Cost_Discounted": res["total_cost_disc"][i],
                    "Discount_Factor": res["discount_factor"][i],
                })

        df_long = pd.DataFrame(rows)
        long_path = os.path.join(output_dir, "yearly_costs_long.csv")
        df_long.to_csv(long_path, index=False)
        print(f"    Saved: yearly_costs_long.csv")

    if format_type in ["wide", "both"]:

        base_label = list(results.keys())[0]
        years = results[base_label]["years"]

        df_wide = pd.DataFrame({"Year": years.astype(int)})

        for label, res in results.items():
            short_name = label.replace(" ", "_")
            df_wide[f"{short_name}_New_Infections"] = res["new_infections"]
            df_wide[f"{short_name}_New_Diagnosed"] = res["new_diagnosed"]
            df_wide[f"{short_name}_HIVCare_Cost"] = res["hivcare_cost"]
            df_wide[f"{short_name}_Test_Cost"] = res["test_cost"]
            df_wide[f"{short_name}_Total_Cost"] = res["total_cost"]
            df_wide[f"{short_name}_HIVCare_Cost_Disc"] = res["hivcare_cost_disc"]
            df_wide[f"{short_name}_Test_Cost_Disc"] = res["test_cost_disc"]
            df_wide[f"{short_name}_Total_Cost_Disc"] = res["total_cost_disc"]

        wide_path = os.path.join(output_dir, "yearly_costs_wide.csv")
        df_wide.to_csv(wide_path, index=False)
        print(f"    Saved: yearly_costs_wide.csv")

    for label, res in results.items():
        scenario_name = config.SCENARIOS[label]["short_name"]
        years = res["years"]

        df_scenario = pd.DataFrame({
            "Year": years.astype(int),
            "New_Infections": res["new_infections"],
            "New_Diagnosed": res["new_diagnosed"],
            "HIVCare_Cost": res["hivcare_cost"],
            "Test_Cost": res["test_cost"],
            "Total_Cost": res["total_cost"],
            "HIVCare_Cost_Discounted": res["hivcare_cost_disc"],
            "Test_Cost_Discounted": res["test_cost_disc"],
            "Total_Cost_Discounted": res["total_cost_disc"],
            "Discount_Factor": res["discount_factor"],

            "Cumulative_Infections": np.cumsum(res["new_infections"]),
            "Cumulative_HIVCare_Cost": np.cumsum(res["hivcare_cost"]),
            "Cumulative_Total_Cost": np.cumsum(res["total_cost"]),
            "Cumulative_HIVCare_Cost_Disc": np.cumsum(res["hivcare_cost_disc"]),
            "Cumulative_Total_Cost_Disc": np.cumsum(res["total_cost_disc"]),
        })

        safe_label = label.replace(" ", "_").replace("-", "")
        scenario_path = os.path.join(output_dir, f"yearly_costs_{safe_label}.csv")
        df_scenario.to_csv(scenario_path, index=False)
        print(f"    Saved: yearly_costs_{safe_label}.csv")

def generate_main_results_table(
    det_results: Dict[str, Any],
    incremental_results: List[Dict[str, float]],
    config: Config
) -> pd.DataFrame:
    rows = []

    for label, res in det_results.items():
        rows.append({
            "Scenario": config.SCENARIOS[label]["short_name"],
            "Cumulative infections (40y)": f"{res['cum_infections']:,.0f}",
            "Cumulative HIV care cost (billion JPY)": f"{res['cum_hivcare_cost_disc'] / 1e9:,.1f}",
            "Cumulative total cost (billion JPY)": f"{res['cum_total_cost_disc'] / 1e9:,.1f}",
        })

    df = pd.DataFrame(rows)
    return df

def generate_incremental_table(
    incremental_results: List[Dict[str, float]],
    config: Config
) -> pd.DataFrame:
    rows = []
    for inc in incremental_results:
        rows.append({
            "Comparison": f"0-month vs {inc['comparator']}",
            "Infections averted": f"{inc['infections_averted']:,.0f}",
            "Incremental cost (billion JPY)": f"{inc['delta_cost_disc'] / 1e9:,.1f}",
            "Cost per infection averted (million JPY)": f"{inc['cost_per_infection_averted'] / 1e6:,.1f}" if not np.isnan(inc['cost_per_infection_averted']) else "Dominant",
            "Lifetime savings (billion JPY)": f"{inc['lifetime_savings_disc'] / 1e9:,.1f}",
            "Net benefit (billion JPY)": f"{inc['net_benefit'] / 1e9:,.1f}",
            "Break-even year": str(inc['break_even_year']) if inc['break_even_year'] > 0 else "N/A",
        })
    return pd.DataFrame(rows)

def generate_psa_table(
    psa_summary: pd.DataFrame,
    config: Config
) -> pd.DataFrame:
    rows = []

    metric_labels = {
        "infections_averted": "Infections averted",
        "delta_cost_disc": "Incremental cost (billion JPY)",
        "cost_per_infection_averted": "Cost per infection averted (million JPY)",
    }

    for _, row in psa_summary.iterrows():
        metric_label = metric_labels.get(row["metric"], row["metric"])

        if "cost" in row["metric"].lower():
            scale = 1e9 if "delta" in row["metric"] else 1e6
            rows.append({
                "Comparison": f"0-month vs {row['comparator']}",
                "Metric": metric_label,
                "Mean": f"{row['mean'] / scale:,.1f}",
                "Median": f"{row['median'] / scale:,.1f}",
                "95% UI": f"({row['ci_lower'] / scale:,.1f}, {row['ci_upper'] / scale:,.1f})",
            })
        else:
            rows.append({
                "Comparison": f"0-month vs {row['comparator']}",
                "Metric": metric_label,
                "Mean": f"{row['mean']:,.0f}",
                "Median": f"{row['median']:,.0f}",
                "95% UI": f"({row['ci_lower']:,.0f}, {row['ci_upper']:,.0f})",
            })

    return pd.DataFrame(rows)

def generate_short_term_bia_table(
    bia_df: pd.DataFrame,
    config: Config
) -> pd.DataFrame:
    rows = []

    for _, row in bia_df.iterrows():
        rows.append({
            "Time horizon": f"{int(row['time_horizon_years'])} years",
            "Comparison": f"0-month vs {row['comparator']}",
            "Infections averted": f"{row['infections_averted']:,.0f}",
            "Budget impact (billion JPY)": f"{row['delta_total_cost_disc'] / 1e9:,.1f}",
            "Annual budget impact (billion JPY/year)": f"{row['annual_budget_impact'] / 1e9:,.2f}",
        })

    return pd.DataFrame(rows)

# ============================================================

# ============================================================

def export_report_json(
    det_results: Dict[str, Any],
    incremental_results: List[Dict[str, float]],
    bia_df: pd.DataFrame,
    phased_df: pd.DataFrame,
    psa_summary: pd.DataFrame,
    ceac_data: pd.DataFrame,
    owsa_results: pd.DataFrame,
    config: Config,
    output_path: str,
) -> None:

    report = {}

    # ------------------------------------------------------------------
    # metadata
    # ------------------------------------------------------------------
    report["metadata"] = {
        "script": "bia_rapid_enhanced.py",
        "created_at": datetime.datetime.now().isoformat(),
        "n_years": config.N_YEARS,
        "discount_rate": config.DISCOUNT_RATE,
        "n_psa_iterations": config.N_PSA_ITERATIONS,
        "random_seed": config.RANDOM_SEED,
        "intervention": config.INTERVENTION_LABEL,
        "comparators": config.COMPARATOR_LABELS,
        "currency": config.CURRENCY_UNIT,
    }

    # ------------------------------------------------------------------
    # cost_parameters
    # ------------------------------------------------------------------
    report["cost_parameters"] = {
        "mean_annual_cost_on_art": config.MEAN_ANNUAL_COST_ON_ART,
        "off_art_cost_ratio": config.OFF_ART_COST_RATIO,
        "hiv_test_unit_cost": config.HIV_TEST_UNIT_COST,
        "hiv_life_expectancy": config.HIV_LIFE_EXPECTANCY,
        "lifetime_art_coverage": config.LIFETIME_ART_COVERAGE,
    }

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    det_summary = {}
    for label, res in det_results.items():
        det_summary[label] = {
            "cumulative_infections": res["cum_infections"],
            "cumulative_hivcare_cost": res["cum_hivcare_cost"],
            "cumulative_total_cost": res["cum_total_cost"],
            "cumulative_hivcare_cost_discounted": res["cum_hivcare_cost_disc"],
            "cumulative_total_cost_discounted": res["cum_total_cost_disc"],
        }
    report["deterministic_summary"] = det_summary

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    yearly_all = []
    for label, res in det_results.items():
        years = res["years"]
        for i, y in enumerate(years):
            yearly_all.append({
                "scenario": label,
                "year": int(y),
                "new_infections": float(res["new_infections"][i]),
                "new_diagnosed": float(res["new_diagnosed"][i]),
                "hivcare_cost": float(res["hivcare_cost"][i]),
                "test_cost": float(res["test_cost"][i]),
                "total_cost": float(res["total_cost"][i]),
                "hivcare_cost_discounted": float(res["hivcare_cost_disc"][i]),
                "total_cost_discounted": float(res["total_cost_disc"][i]),
                "discount_factor": float(res["discount_factor"][i]),
            })
    report["yearly_costs"] = yearly_all

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    inc_list = []
    for inc in incremental_results:
        inc_clean = {}
        for k, v in inc.items():
            if isinstance(v, (np.integer, np.floating)):
                v = float(v)
            if isinstance(v, float) and np.isnan(v):
                v = None
            inc_clean[k] = v
        inc_list.append(inc_clean)
    report["incremental_results"] = inc_list

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    breakeven = {}
    intervention = config.INTERVENTION_LABEL
    res_int = det_results[intervention]
    for comparator in config.COMPARATOR_LABELS:
        res_comp = det_results[comparator]
        diff_yearly = res_int["total_cost_disc"] - res_comp["total_cost_disc"]
        cum_diff = np.cumsum(diff_yearly)
        idx = np.where(cum_diff <= 0)[0]
        be_year = int(idx[0] + 1) if len(idx) > 0 else -1

        yearly_diff_list = []
        for i, y in enumerate(res_int["years"]):
            yearly_diff_list.append({
                "year": int(y),
                "annual_diff_discounted": float(diff_yearly[i]),
                "cumulative_diff_discounted": float(cum_diff[i]),
            })

        breakeven[f"0m_vs_{comparator}"] = {
            "break_even_year": be_year,
            "final_cumulative_diff": float(cum_diff[-1]),
            "yearly_data": yearly_diff_list,
        }
    report["breakeven_analysis"] = breakeven

    # ------------------------------------------------------------------
    # short_term_bia
    # ------------------------------------------------------------------
    report["short_term_bia"] = bia_df.to_dict(orient="records")

    # ------------------------------------------------------------------
    # phased_implementation
    # ------------------------------------------------------------------
    report["phased_implementation"] = phased_df.to_dict(orient="records")

    # ------------------------------------------------------------------
    # psa_summary
    # ------------------------------------------------------------------
    psa_list = []
    for _, row in psa_summary.iterrows():
        entry = {}
        for col in psa_summary.columns:
            val = row[col]
            if isinstance(val, (np.integer, np.floating)):
                val = float(val)
            if isinstance(val, float) and np.isnan(val):
                val = None
            entry[col] = val
        psa_list.append(entry)
    report["psa_summary"] = psa_list

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    ceac_summary = []
    key_wtp = [0, 5e6, 10e6, 20e6, 30e6, 50e6, 100e6]
    for comparator in ceac_data["comparator"].unique():
        df_c = ceac_data[ceac_data["comparator"] == comparator]
        for wtp in key_wtp:
            row_match = df_c[df_c["wtp"] == wtp]
            if len(row_match) > 0:
                ceac_summary.append({
                    "comparator": comparator,
                    "wtp_threshold": float(wtp),
                    "prob_cost_effective": float(row_match["prob_cost_effective"].iloc[0]),
                })
    report["ceac_summary"] = ceac_summary

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    if len(owsa_results) > 0:
        owsa_list = []
        for _, row in owsa_results.iterrows():
            entry = {}
            for col in owsa_results.columns:
                val = row[col]
                if isinstance(val, (np.integer, np.floating)):
                    val = float(val)
                if isinstance(val, float) and np.isnan(val):
                    val = None
                entry[col] = val
            owsa_list.append(entry)
        report["owsa_results"] = owsa_list

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)

    print(f"    Saved: {os.path.basename(output_path)}")

# ============================================================

# ============================================================

def main():
    config = Config()
    os.makedirs(config.OUTPUT_DIR, exist_ok=True)

    print("=" * 70)
    print("Budget Impact Analysis: Rapid ART (0-month) vs Standard Care")
    print("=" * 70)
    print(f"Output directory: {config.OUTPUT_DIR}")
    print(f"PSA iterations: {config.N_PSA_ITERATIONS}")
    print(f"Random seed: {config.RANDOM_SEED}")
    print()

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    print("[1] Loading aggregator data...")
    aggregators = {}
    available_scenarios = []

    for label, info in config.SCENARIOS.items():
        path = os.path.join(config.INPUT_DIR, info["agg_file"])
        if os.path.exists(path):
            aggregators[label] = load_aggregator(path)
            available_scenarios.append(label)
            print(f"    {label}: {len(aggregators[label]):,} rows")
        else:
            print(f"    {label}: FILE NOT FOUND - skipping ({info['agg_file']})")

    if config.INTERVENTION_LABEL not in aggregators:
        raise FileNotFoundError(
            f"Intervention scenario '{config.INTERVENTION_LABEL}' data not found. "
            f"File required: {config.SCENARIOS[config.INTERVENTION_LABEL]['agg_file']}"
        )

    available_comparators = [
        label for label in config.COMPARATOR_LABELS
        if label in aggregators
    ]

    if not available_comparators:
        raise FileNotFoundError(
            "No comparator scenario data found. At least one of the following files is required: "
            + ", ".join([config.SCENARIOS[label]['agg_file'] for label in config.COMPARATOR_LABELS])
        )

    config.COMPARATOR_LABELS = available_comparators

    print()
    print(f"    Available scenarios: {', '.join(available_scenarios)}")
    print(f"    Intervention: {config.INTERVENTION_LABEL}")
    print(f"    Comparators: {', '.join(config.COMPARATOR_LABELS)}")
    print()

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    print("[2] Running deterministic analysis...")
    det_results = run_deterministic_analysis(aggregators, config)

    incremental_results = []
    for comparator in config.COMPARATOR_LABELS:
        inc = compute_incremental_results(
            det_results, config.INTERVENTION_LABEL, comparator, config
        )
        incremental_results.append(inc)
        print(f"    0-month vs {comparator}:")
        print(f"      Infections averted: {inc['infections_averted']:,.0f}")
        print(f"      Incremental cost: {inc['delta_cost_disc'] / 1e9:,.1f} billion JPY")
        if not np.isnan(inc['cost_per_infection_averted']):
            print(f"      Cost per infection averted: {inc['cost_per_infection_averted'] / 1e6:,.1f} million JPY")
        print(f"      Break-even year: {inc['break_even_year'] if inc['break_even_year'] > 0 else 'N/A'}")
    print()

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    print("[2.5] Exporting yearly costs to CSV...")
    export_yearly_costs_to_csv(det_results, config, config.OUTPUT_DIR)
    print()

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    print("[3] Computing short-term BIA...")
    bia_df = run_short_term_bia(det_results, config)
    bia_df.to_csv(
        os.path.join(config.OUTPUT_DIR, "short_term_bia.csv"),
        index=False
    )
    print(f"    Saved: short_term_bia.csv")
    print()

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    print("[4] Computing phased implementation scenarios...")
    phased_df = run_phased_implementation(det_results, config)
    phased_df.to_csv(
        os.path.join(config.OUTPUT_DIR, "phased_implementation.csv"),
        index=False
    )
    print(f"    Saved: phased_implementation.csv")
    print()

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    print("[5] Running Probabilistic Sensitivity Analysis (PSA)...")
    psa_results = run_psa(aggregators, config)
    psa_results.to_csv(
        os.path.join(config.OUTPUT_DIR, "psa_raw_results.csv"),
        index=False
    )
    print(f"    Saved: psa_raw_results.csv")

    psa_summary = compute_psa_summary(psa_results)
    psa_summary.to_csv(
        os.path.join(config.OUTPUT_DIR, "psa_summary.csv"),
        index=False
    )
    print(f"    Saved: psa_summary.csv")

    ceac_data = compute_ceac(psa_results, config.WTP_THRESHOLDS)
    ceac_data.to_csv(
        os.path.join(config.OUTPUT_DIR, "ceac_data.csv"),
        index=False
    )
    print(f"    Saved: ceac_data.csv")
    print()

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    print("[6] Running One-Way Sensitivity Analysis (OWSA)...")
    owsa_results = run_owsa(aggregators, config)
    owsa_results.to_csv(
        os.path.join(config.OUTPUT_DIR, "owsa_results.csv"),
        index=False
    )
    print(f"    Saved: owsa_results.csv")
    print()

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    print("[7] Generating figures...")

    plot_yearly_comparison(
        det_results, config,
        os.path.join(config.OUTPUT_DIR, "fig_yearly_comparison.png")
    )

    plot_cumulative_budget_impact(
        det_results, config,
        os.path.join(config.OUTPUT_DIR, "fig_cumulative_budget_impact.png")
    )

    plot_short_term_bia_bars(
        bia_df, config,
        os.path.join(config.OUTPUT_DIR, "fig_short_term_bia.png")
    )

    # CEAC
    plot_ceac(
        ceac_data, config,
        os.path.join(config.OUTPUT_DIR, "fig_ceac.png")
    )

    # Cost-Effectiveness Plane
    plot_cost_effectiveness_plane(
        psa_results, config,
        os.path.join(config.OUTPUT_DIR, "fig_ce_plane.png")
    )

    # Tornado Diagram
    plot_tornado_diagram(
        owsa_results, config,
        metric="delta_cost_disc",
        output_path=os.path.join(config.OUTPUT_DIR, "fig_tornado.png")
    )
    print()

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    print("[8] Generating paper-ready tables...")

    table1 = generate_main_results_table(det_results, incremental_results, config)
    table1.to_csv(
        os.path.join(config.OUTPUT_DIR, "table1_main_results.csv"),
        index=False
    )
    print(f"    Saved: table1_main_results.csv")

    table2 = generate_incremental_table(incremental_results, config)
    table2.to_csv(
        os.path.join(config.OUTPUT_DIR, "table2_incremental.csv"),
        index=False
    )
    print(f"    Saved: table2_incremental.csv")

    table3 = generate_psa_table(psa_summary, config)
    table3.to_csv(
        os.path.join(config.OUTPUT_DIR, "table3_psa_summary.csv"),
        index=False
    )
    print(f"    Saved: table3_psa_summary.csv")

    table4 = generate_short_term_bia_table(bia_df, config)
    table4.to_csv(
        os.path.join(config.OUTPUT_DIR, "table4_short_term_bia.csv"),
        index=False
    )
    print(f"    Saved: table4_short_term_bia.csv")

    print()
    print("=" * 70)
    print("Analysis complete!")
    print(f"All results saved to: {config.OUTPUT_DIR}")
    print("=" * 70)

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    print("\n[9] Exporting report-ready JSON...")
    export_report_json(
        det_results=det_results,
        incremental_results=incremental_results,
        bia_df=bia_df,
        phased_df=phased_df,
        psa_summary=psa_summary,
        ceac_data=ceac_data,
        owsa_results=owsa_results,
        config=config,
        output_path=os.path.join(config.OUTPUT_DIR, "bia_report_ready.json"),
    )
    print()
    print("=" * 70)
    print("JSON report exported: bia_report_ready.json")
    print("=" * 70)

if __name__ == "__main__":
    main()

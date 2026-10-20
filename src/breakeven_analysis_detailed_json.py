#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Break-even-year analysis of cumulative budget impact."""

import os
import json
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = ['DejaVu Sans', 'sans-serif']

import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# --- Data location: set BASE_144 to your local 144_scenario folder ---
BASE_144 = r"/path/to/project/144_scenario"
PROJ_169 = os.path.dirname(BASE_144)
_PROJ_144 = BASE_144

INPUT_DIR = os.path.join(_PROJ_144, "core", "hiv_rapid_art_model_multiplier_mainfoi1p18_rr1_rate1_foi1p05_rfoi1p04_ageoff5_seed713")
OUTPUT_DIR = os.path.join(_PROJ_144, "bia_results", "Break-even")

MEAN_ANNUAL_COST_ON_ART = {
    "20-29": 2_375_000.0,
    "30-39": 2_901_449.0,
    "40-49": 3_264_767.0,
    "50-59": 3_433_598.0,
    "60-69": 3_356_900.0,
    "70+":   3_164_654.0,
}
OFF_ART_COST_RATIO = 0.3
DISCOUNT_RATE = 0.02
N_YEARS = 40

def load_aggregator(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if df["OnART"].dtype != bool:
        df["OnART"] = df["OnART"].astype(bool)
    return df

def discount_factor(year: int, rate: float) -> float:
    return 1.0 / ((1.0 + rate) ** (year - 1))

def compute_yearly_costs(df_agg: pd.DataFrame) -> np.ndarray:
    years = np.arange(1, N_YEARS + 1)
    
    cost_on_art = dict(MEAN_ANNUAL_COST_ON_ART)
    cost_off_art = {age: val * OFF_ART_COST_RATIO for age, val in MEAN_ANNUAL_COST_ON_ART.items()}
    
    df_care = df_agg[df_agg["CD4Category"] != "HIV-"].copy()
    
    if len(df_care) == 0:
        return np.zeros(len(years))
    
    df_care["CostOnArt"] = df_care["AgeGroup"].map(cost_on_art)
    df_care["CostOffArt"] = df_care["AgeGroup"].map(cost_off_art)
    mean_cost = np.mean(list(cost_on_art.values()))
    df_care["CostOnArt"] = df_care["CostOnArt"].fillna(mean_cost)
    df_care["CostOffArt"] = df_care["CostOffArt"].fillna(mean_cost * OFF_ART_COST_RATIO)
    
    df_care["PersonYears"] = df_care["PersonMonthsAlive"] / 12.0
    df_care["AnnualCostPerPY"] = np.where(
        df_care["OnART"],
        df_care["CostOnArt"],
        df_care["CostOffArt"],
    )
    df_care["Cost"] = df_care["AnnualCostPerPY"] * df_care["PersonYears"]
    
    care_by_year = df_care.groupby("Year")["Cost"].sum()
    yearly_cost = care_by_year.reindex(years).fillna(0.0).values.astype(float)
    
    return yearly_cost

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print("=" * 70)
    print("Break-even Analysis: Detailed Visualization")
    print("=" * 70)

    print("\n[1] Loading data...")
    scenarios = {
        "0 month": "scenario0_yearly_summary.csv",
        "1 month": "scenario1_yearly_summary.csv",
        "2 months": "scenario2_yearly_summary.csv",
    }

    aggregators = {}
    for label, filename in scenarios.items():
        path = os.path.join(INPUT_DIR, filename)
        if os.path.exists(path):
            aggregators[label] = load_aggregator(path)
            print(f"    {label}: loaded")
        else:
            print(f"    {label}: file not found, skipping")

    if "0 month" not in aggregators:
        print("\n    ERROR: 0 month scenario file is required but not found!")
        return

    comparison_targets = [s for s in ["1 month", "2 months"] if s in aggregators]
    if not comparison_targets:
        print("\n    ERROR: At least one of 1 month or 2 months scenario is required!")
        return

    print(f"\n    Available comparisons: 0 month vs {', '.join(comparison_targets)}")

    print("\n[2] Computing yearly costs...")
    years = np.arange(1, N_YEARS + 1)
    disc_factors = np.array([discount_factor(y, DISCOUNT_RATE) for y in years])
    
    yearly_costs = {}
    yearly_costs_disc = {}
    for label, df in aggregators.items():
        costs = compute_yearly_costs(df)
        yearly_costs[label] = costs
        yearly_costs_disc[label] = costs * disc_factors
        print(f"    {label}: Total (undiscounted) = ¥{costs.sum()/1e12:.2f}T")

    print("\n[3] Computing yearly cost differences...")

    comparisons = []
    if "1 month" in aggregators:
        comparisons.append(("0 month", "1 month", "0m vs 1m"))
    if "2 months" in aggregators:
        comparisons.append(("0 month", "2 months", "0m vs 2m"))
    
    results = {}
    for int_label, comp_label, name in comparisons:

        diff_yearly = yearly_costs_disc[int_label] - yearly_costs_disc[comp_label]
        cum_diff = np.cumsum(diff_yearly)
        
        # Break-even year
        idx = np.where(cum_diff <= 0)[0]
        break_even_year = int(idx[0] + 1) if len(idx) > 0 else -1
        
        results[name] = {
            "diff_yearly": diff_yearly,
            "cum_diff": cum_diff,
            "break_even_year": break_even_year,
        }
        
        print(f"\n    {name}:")
        print(f"      Break-even year: {break_even_year}")
        print(f"      Year 1 diff: ¥{diff_yearly[0]/1e9:.2f}B")
        print(f"      Year 5 cum: ¥{cum_diff[4]/1e9:.2f}B")
        print(f"      Year 10 cum: ¥{cum_diff[9]/1e9:.2f}B")
        print(f"      Year 14 cum: ¥{cum_diff[13]/1e9:.2f}B")
        print(f"      Year 15 cum: ¥{cum_diff[14]/1e9:.2f}B")
        print(f"      Year 40 cum: ¥{cum_diff[39]/1e9:.2f}B")

    print("\n[4] Creating figures...")
    
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    n_comparisons = len(results)
    fig, axes = plt.subplots(2, n_comparisons, figsize=(7 * n_comparisons, 10))

    if n_comparisons == 1:
        axes = axes.reshape(2, 1)

    colors = {"0m vs 1m": "#28A745", "0m vs 2m": "#2E86AB"}

    for idx, (name, res) in enumerate(results.items()):

        ax1 = axes[0, idx]
        ax1.bar(years, res["diff_yearly"] / 1e9, color=colors[name], alpha=0.7, width=0.8)
        ax1.axhline(0, color='black', linewidth=1)
        ax1.set_xlabel("Year", fontsize=11)
        ax1.set_ylabel("Annual Cost Difference (Billion JPY)", fontsize=11)
        ax1.set_title(f"{name}: Annual Cost Difference\n(Positive = Rapid ART costs more)", fontsize=12, fontweight='bold')
        ax1.grid(axis='y', alpha=0.3)

        if res["diff_yearly"][0] > 0:
            ax1.annotate('Early: Higher cost\n(more patients on ART)',
                         xy=(5, res["diff_yearly"][4]/1e9),
                         fontsize=9, fontweight='bold', ha='center', color='red')

        ax2 = axes[1, idx]
        ax2.plot(years, res["cum_diff"] / 1e9, color=colors[name], linewidth=2.5)
        ax2.fill_between(years, 0, res["cum_diff"] / 1e9,
                         where=(res["cum_diff"] > 0), alpha=0.3, color='red',
                         label='Investment period')
        ax2.fill_between(years, 0, res["cum_diff"] / 1e9,
                         where=(res["cum_diff"] <= 0), alpha=0.3, color='green',
                         label='Savings period')
        ax2.axhline(0, color='black', linewidth=1.5, linestyle='--')

        # Break-even point
        be_year = res["break_even_year"]
        if be_year > 0:
            ax2.axvline(be_year, color='red', linewidth=2, linestyle=':', label=f'Break-even (Year {be_year})')
            ax2.scatter([be_year], [0], s=150, color='red', zorder=10, marker='o')
            ax2.annotate(f'Break-even\nYear {be_year}',
                         xy=(be_year, 0), xytext=(7, -10),
                         fontsize=11, fontweight='bold', color='red',
                         ha='center',
                         arrowprops=dict(arrowstyle='->', color='red'))

        ax2.set_xlabel("Year", fontsize=11)
        ax2.set_ylabel("Cumulative Cost Difference (Billion JPY)", fontsize=11)
        ax2.set_title(f"{name}: Cumulative Budget Impact\n(Negative = Net Savings)", fontsize=12, fontweight='bold')
        ax2.legend(loc='upper right', fontsize=9)
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim(0, 41)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "fig_breakeven_detailed.png"), dpi=300, bbox_inches='tight')
    plt.close()
    print("    Saved: fig_breakeven_detailed.png")
    
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 7))
    
    for name, res in results.items():
        ax.plot(years, res["cum_diff"] / 1e9, color=colors[name], linewidth=2.5, label=name)
        
        # Break-even point
        be_year = res["break_even_year"]
        if be_year > 0:
            ax.scatter([be_year], [0], s=200, color=colors[name], zorder=10, 
                      marker='o', edgecolor='black', linewidth=2)

    ax.axhspan(0, 50, alpha=0.1, color='red', label='_nolegend_')
    ax.axhspan(-70, 0, alpha=0.1, color='green', label='_nolegend_')

    ax.axhline(0, color='black', linewidth=2, linestyle='--')

    be_years_list = [res["break_even_year"] for res in results.values() if res["break_even_year"] > 0]
    if be_years_list:

        representative_be_year = be_years_list[0]
        ax.axvline(representative_be_year, color='gray', linewidth=1.5, linestyle=':', alpha=0.7)

    ax.text(7, 15, 'Investment Period\n(Short-term cost increase)',
            fontsize=11, ha='center', color='darkred', fontweight='bold')
    ax.text(32, -20, 'Savings Period\n(Long-term cost savings)',
            fontsize=11, ha='center', color='darkgreen', fontweight='bold')

    if be_years_list:
        be_text = f'Break-even\n(Year {representative_be_year})'
        ax.text(representative_be_year + 0.5, 5, be_text,
                fontsize=10, ha='left', color='gray', fontweight='bold')

    for name, res in results.items():
        final_val = res["cum_diff"][39] / 1e9
        ax.annotate(f'{name}\n¥{final_val:.1f}B', 
                   xy=(40, final_val), xytext=(42, final_val),
                   fontsize=10, fontweight='bold', color=colors[name],
                   va='center')
    
    ax.set_xlabel("Time Horizon (Years)", fontsize=12)
    ax.set_ylabel("Cumulative Budget Impact (Billion JPY)", fontsize=12)
    ax.set_title("Break-even Analysis: Rapid ART Implementation\n(Negative values = cost savings)", 
                fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 45)
    ax.set_ylim(-70, 30)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "fig_breakeven_summary.png"), dpi=300, bbox_inches='tight')
    plt.close()
    print("    Saved: fig_breakeven_summary.png")
    
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    print("\n[5] Saving data...")

    yearly_data = {"Year": years}
    for name, res in results.items():
        yearly_data[f"{name}_annual_diff_B"] = res["diff_yearly"] / 1e9
        yearly_data[f"{name}_cumulative_B"] = res["cum_diff"] / 1e9
    
    df_yearly = pd.DataFrame(yearly_data)
    df_yearly.to_csv(os.path.join(OUTPUT_DIR, "yearly_cost_differences.csv"), index=False)
    print("    Saved: yearly_cost_differences.csv")

    summary_data = []
    for name, res in results.items():
        summary_data.append({
            "Comparison": name,
            "Break-even year": res["break_even_year"],
            "Year 3 cumulative discounted (B JPY)": f"{res['cum_diff'][2]/1e9:.2f}",
            "Year 5 cumulative discounted (B JPY)": f"{res['cum_diff'][4]/1e9:.2f}",
            "Year 10 cumulative discounted (B JPY)": f"{res['cum_diff'][9]/1e9:.2f}",
            "Year 14 cumulative discounted (B JPY)": f"{res['cum_diff'][13]/1e9:.2f}",
            "Year 15 cumulative discounted (B JPY)": f"{res['cum_diff'][14]/1e9:.2f}",
            "Year 40 cumulative discounted (B JPY)": f"{res['cum_diff'][39]/1e9:.2f}",
        })
    
    df_summary = pd.DataFrame(summary_data)
    df_summary.to_csv(os.path.join(OUTPUT_DIR, "breakeven_summary.csv"), index=False)
    print("    Saved: breakeven_summary.csv")
    
    print("\n" + "=" * 70)
    print("Break-even Analysis Summary:")
    print("=" * 70)
    print(df_summary.to_string(index=False))
    print("\n" + "=" * 70)
    print(f"All outputs saved to: {OUTPUT_DIR}")
    print("=" * 70)

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    print("\n[6] Exporting report-ready JSON...")

    json_report = {
        "metadata": {
            "script": "breakeven_analysis_detailed.py",
            "created_at": datetime.datetime.now().isoformat(),
            "n_years": N_YEARS,
            "discount_rate": DISCOUNT_RATE,
            "cost_parameters": {
                "mean_annual_cost_on_art": MEAN_ANNUAL_COST_ON_ART,
                "off_art_cost_ratio": OFF_ART_COST_RATIO,
            },
        },
        "yearly_costs_by_scenario": {},
        "breakeven_comparisons": {},
    }

    for label in yearly_costs:
        json_report["yearly_costs_by_scenario"][label] = {
            "yearly_undiscounted": yearly_costs[label].tolist(),
            "yearly_discounted": yearly_costs_disc[label].tolist(),
            "total_undiscounted": float(yearly_costs[label].sum()),
            "total_discounted": float(yearly_costs_disc[label].sum()),
        }

    for name, res in results.items():
        yearly_data = []
        for i, y in enumerate(years):
            yearly_data.append({
                "year": int(y),
                "annual_diff_billion_jpy": float(res["diff_yearly"][i] / 1e9),
                "cumulative_diff_billion_jpy": float(res["cum_diff"][i] / 1e9),
            })

        json_report["breakeven_comparisons"][name] = {
            "break_even_year": res["break_even_year"],
            "final_cumulative_billion_jpy": float(res["cum_diff"][-1] / 1e9),
            "milestones": {
                "year_3": float(res["cum_diff"][2] / 1e9),
                "year_5": float(res["cum_diff"][4] / 1e9),
                "year_10": float(res["cum_diff"][9] / 1e9),
                "year_14": float(res["cum_diff"][13] / 1e9),
                "year_15": float(res["cum_diff"][14] / 1e9),
                "year_20": float(res["cum_diff"][19] / 1e9),
                "year_40": float(res["cum_diff"][39] / 1e9),
            },
            "yearly_data": yearly_data,
        }

    json_path = os.path.join(OUTPUT_DIR, "breakeven_report_ready.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(json_report, f, ensure_ascii=False, indent=2)
    print(f"    Saved: breakeven_report_ready.json")
    print("=" * 70)

if __name__ == "__main__":
    main()

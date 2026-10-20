#!/usr/bin/env python3
"""Visualise yearly cost trajectories by scenario."""

import os
import json
import datetime
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

plt.rcParams['font.family'] = ['DejaVu Sans', 'sans-serif']

import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# --- Data location: set BASE_144 to your local 144_scenario folder ---
BASE_144 = r"/path/to/project/144_scenario"
PROJ_169 = os.path.dirname(BASE_144)
_PROJ_144 = BASE_144

csv_path = os.path.join(_PROJ_144, "bia_results", "enhanced", "yearly_costs_wide.csv")
df = pd.read_csv(csv_path)

output_dir = os.path.dirname(csv_path)
output_path = os.path.join(output_dir, 'healthcare_cost_trends.png')

SCALE = 1e8

colors = {
    '0_month': '#2E86AB',
    '1_month': '#E94F37',
    '2_months': '#F39C12'
}

labels = {
    '0_month': 'Immediate ART (0 months)',
    '1_month': '1 month delay',
    '2_months': '2 months delay'
}

all_scenarios = ['0_month', '1_month', '2_months']
available_scenarios = []
for scenario in all_scenarios:
    test_col = f'{scenario}_Total_Cost_Disc'
    if test_col in df.columns:
        available_scenarios.append(scenario)

comparison_scenarios = [s for s in available_scenarios if s != '0_month']

print(f"Detected scenarios: {available_scenarios}")
if len(available_scenarios) < 3:
    missing = set(all_scenarios) - set(available_scenarios)
    print(f"Note: the following scenarios have no data and are skipped: {missing}")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('HIV Healthcare Cost Trends by ART Initiation Timing', fontsize=14, fontweight='bold')

ax1 = axes[0, 0]
for scenario in available_scenarios:
    col = f'{scenario}_Total_Cost_Disc'
    ax1.plot(df['Year'], df[col] / SCALE,
             color=colors[scenario], linewidth=2, label=labels[scenario])

ax1.set_xlabel('Year')
ax1.set_ylabel('Total Cost (100 million JPY)')
ax1.set_title('A. Discounted Total Healthcare Cost')
ax1.legend(loc='upper right', fontsize=9)
ax1.grid(True, alpha=0.3)
ax1.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f'{x:,.0f}'))

ax2 = axes[0, 1]

cumulative_costs = {}
for scenario in available_scenarios:
    col = f'{scenario}_Total_Cost_Disc'
    cumulative_costs[scenario] = df[col].cumsum()

for scenario in comparison_scenarios:
    diff = (cumulative_costs[scenario] - cumulative_costs['0_month']) / SCALE
    ax2.plot(df['Year'], diff, color=colors[scenario], linewidth=2, label=labels[scenario])

ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax2.set_xlabel('Year')
ax2.set_ylabel('Cumulative Cost Difference (100 million JPY)')
ax2.set_title('B. Cumulative Cost Difference vs Immediate ART')
ax2.legend(loc='upper left', fontsize=9)
ax2.grid(True, alpha=0.3)
ax2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f'{x:,.0f}'))

ax3 = axes[1, 0]
for scenario in available_scenarios:
    col = f'{scenario}_New_Infections'
    ax3.plot(df['Year'], df[col],
             color=colors[scenario], linewidth=2, label=labels[scenario])

ax3.set_xlabel('Year')
ax3.set_ylabel('New Infections (persons)')
ax3.set_title('C. Annual New HIV Infections')
ax3.legend(loc='upper right', fontsize=9)
ax3.grid(True, alpha=0.3)

ax4 = axes[1, 1]

selected_years = [1, 20, 40]
n_scenarios = len(available_scenarios)
bar_width = 0.8 / n_scenarios
x = np.arange(len(selected_years))

for i, scenario in enumerate(available_scenarios):
    care_costs = []
    for year in selected_years:
        row = df[df['Year'] == year].iloc[0]
        care_costs.append(row[f'{scenario}_HIVCare_Cost_Disc'] / SCALE)

    offset = (i - (n_scenarios - 1) / 2) * bar_width
    ax4.bar(x + offset, care_costs, bar_width,
            color=colors[scenario], alpha=0.85)

ax4.set_xlabel('Year')
ax4.set_ylabel('HIV Care Cost (100 million JPY)')
ax4.set_title('D. HIV Care Cost at Selected Years')
ax4.set_xticks(x)
ax4.set_xticklabels([f'Year {y}' for y in selected_years])
ax4.legend([plt.Rectangle((0,0),1,1, color=colors[s]) for s in available_scenarios],
           [labels[s] for s in available_scenarios], loc='upper right', fontsize=9)
ax4.grid(True, alpha=0.3, axis='y')
ax4.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f'{x:,.0f}'))

plt.tight_layout()
plt.savefig(output_path, dpi=150, bbox_inches='tight')
plt.close()

print("=" * 60)
print("HIV Healthcare Cost Analysis Summary")
print("=" * 60)

print("\n1. 40-Year Cumulative Total Cost (Discounted, billion JPY):")
for scenario in available_scenarios:
    col = f'{scenario}_Total_Cost_Disc'
    total = df[col].sum() / 1e9
    print(f"   {labels[scenario]}: {total:,.1f}")

print("\n2. Total New Infections over 40 Years:")
for scenario in available_scenarios:
    col = f'{scenario}_New_Infections'
    total = df[col].sum()
    print(f"   {labels[scenario]}: {total:,.0f}")

print("\n3. Cost Savings with Immediate ART (vs delayed initiation):")
base = df['0_month_Total_Cost_Disc'].sum()
for scenario in comparison_scenarios:
    col = f'{scenario}_Total_Cost_Disc'
    diff = df[col].sum() - base
    print(f"   vs {labels[scenario]}: {diff/1e9:,.1f} billion JPY saved")

print("\n4. Infections Prevented with Immediate ART:")
base_inf = df['0_month_New_Infections'].sum()
for scenario in comparison_scenarios:
    col = f'{scenario}_New_Infections'
    diff = df[col].sum() - base_inf
    print(f"   vs {labels[scenario]}: {diff:,.0f} infections prevented")

print("\n" + "=" * 60)
print(f"Figure saved: {output_path}")
print("=" * 60)

print("\nExporting report-ready JSON...")

json_report = {
    "metadata": {
        "script": "visual_yearly_cost.py",
        "created_at": datetime.datetime.now().isoformat(),
        "data_source": csv_path,
        "scale_unit": "JPY",
        "available_scenarios": available_scenarios,
    },
    "forty_year_summary": {},
    "yearly_trends": [],
    "cost_savings_vs_immediate_art": [],
    "infections_prevented_vs_immediate_art": [],
    "selected_year_costs": {},
}

# 40-year cumulative totals
for scenario in available_scenarios:
    col_cost = f'{scenario}_Total_Cost_Disc'
    col_inf = f'{scenario}_New_Infections'
    col_care = f'{scenario}_HIVCare_Cost_Disc'
    json_report["forty_year_summary"][scenario] = {
        "label": labels[scenario],
        "total_cost_discounted": float(df[col_cost].sum()),
        "total_cost_discounted_billion_jpy": float(df[col_cost].sum() / 1e9),
        "total_hivcare_cost_discounted": float(df[col_care].sum()),
        "total_new_infections": int(df[col_inf].sum()),
    }

# Yearly trends
for _, row in df.iterrows():
    entry = {"year": int(row["Year"])}
    for scenario in available_scenarios:
        entry[f"{scenario}_total_cost_disc"] = float(row[f'{scenario}_Total_Cost_Disc'])
        entry[f"{scenario}_hivcare_cost_disc"] = float(row[f'{scenario}_HIVCare_Cost_Disc'])
        entry[f"{scenario}_new_infections"] = int(row[f'{scenario}_New_Infections'])
    json_report["yearly_trends"].append(entry)

# Cost savings vs immediate ART
base_cost = df['0_month_Total_Cost_Disc'].sum()
base_inf = df['0_month_New_Infections'].sum()
for scenario in comparison_scenarios:
    col_cost = f'{scenario}_Total_Cost_Disc'
    col_inf = f'{scenario}_New_Infections'
    cost_diff = df[col_cost].sum() - base_cost
    inf_diff = df[col_inf].sum() - base_inf
    json_report["cost_savings_vs_immediate_art"].append({
        "comparison": f"0_month_vs_{scenario}",
        "label": f"Immediate ART vs {labels[scenario]}",
        "cost_difference": float(cost_diff),
        "cost_difference_billion_jpy": float(cost_diff / 1e9),
    })
    json_report["infections_prevented_vs_immediate_art"].append({
        "comparison": f"0_month_vs_{scenario}",
        "label": f"Immediate ART vs {labels[scenario]}",
        "infections_prevented": int(inf_diff),
    })

# Selected year costs
for year in [1, 10, 20, 30, 40]:
    row = df[df['Year'] == year]
    if len(row) > 0:
        row = row.iloc[0]
        year_data = {}
        for scenario in available_scenarios:
            year_data[scenario] = {
                "hivcare_cost_disc": float(row[f'{scenario}_HIVCare_Cost_Disc']),
                "total_cost_disc": float(row[f'{scenario}_Total_Cost_Disc']),
                "new_infections": int(row[f'{scenario}_New_Infections']),
            }
        json_report["selected_year_costs"][f"year_{year}"] = year_data

json_path = os.path.join(output_dir, "visual_yearly_cost_report_ready.json")
with open(json_path, "w", encoding="utf-8") as f:
    json.dump(json_report, f, ensure_ascii=False, indent=2)
print(f"JSON report saved: {json_path}")
print("=" * 60)
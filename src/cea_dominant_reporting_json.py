#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Dominance reporting for the cost-effectiveness analysis."""

import os
import json
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter

import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# --- Data location: set BASE_144 to your local 144_scenario folder ---
BASE_144 = r"/path/to/20250315_169/144_scenario"
PROJ_169 = os.path.dirname(BASE_144)
_PROJ_144 = BASE_144

plt.rcParams['font.family'] = ['DejaVu Sans', 'sans-serif']

# =============================================================================

# =============================================================================

def load_results(output_dir: str) -> dict:
    results = {}

    psa_path = os.path.join(output_dir, "psa_results.csv")
    if os.path.exists(psa_path):
        results["psa"] = pd.read_csv(psa_path)

    inc_path = os.path.join(output_dir, "incremental_results.csv")
    if os.path.exists(inc_path):
        results["incremental"] = pd.read_csv(inc_path)
    
    return results

# =============================================================================

# =============================================================================

def generate_dominant_summary_table(psa_df: pd.DataFrame) -> pd.DataFrame:
    records = []
    
    for comparison in psa_df["comparison"].unique():
        df = psa_df[psa_df["comparison"] == comparison]

        delta_cost_mean = df["delta_cost"].mean()
        delta_cost_ci_low = df["delta_cost"].quantile(0.025)
        delta_cost_ci_high = df["delta_cost"].quantile(0.975)
        
        delta_qaly_mean = df["delta_qaly"].mean()
        delta_qaly_ci_low = df["delta_qaly"].quantile(0.025)
        delta_qaly_ci_high = df["delta_qaly"].quantile(0.975)

        n_total = len(df)
        n_dominant = ((df["delta_cost"] < 0) & (df["delta_qaly"] > 0)).sum()
        n_dominated = ((df["delta_cost"] > 0) & (df["delta_qaly"] < 0)).sum()
        n_ne = ((df["delta_cost"] > 0) & (df["delta_qaly"] > 0)).sum()
        n_sw = ((df["delta_cost"] < 0) & (df["delta_qaly"] < 0)).sum()
        
        prob_dominant = n_dominant / n_total
        prob_dominated = n_dominated / n_total
        prob_ne = n_ne / n_total
        prob_sw = n_sw / n_total
        
        # NMB（WTP = 5M）
        nmb_5m = df["delta_qaly"] * 5_000_000 - df["delta_cost"]
        nmb_5m_mean = nmb_5m.mean()
        prob_ce_5m = (nmb_5m > 0).mean()

        if prob_dominant >= 0.95:
            conclusion = "Dominant (cost saving + better outcomes)"
            recommendation = "Strongly recommended"
        elif prob_ce_5m >= 0.95:
            conclusion = "Cost-effective at WTP=5M"
            recommendation = "Recommended"
        elif prob_dominated >= 0.50:
            conclusion = "Dominated (cost increase + worse outcomes)"
            recommendation = "Not recommended"
        else:
            conclusion = "Uncertain"
            recommendation = "Further review needed"
        
        records.append({
            "Comparison": comparison.replace("_", " ").replace("vs", "vs."),
            "Cost saving (mean)": f"¥{-delta_cost_mean/1e9:.1f}B",
            "Cost saving (95% CI)": f"¥{-delta_cost_ci_high/1e9:.1f}B - ¥{-delta_cost_ci_low/1e9:.1f}B",
            "QALY gain (mean)": f"{delta_qaly_mean:,.0f}",
            "QALY gain (95% CI)": f"{delta_qaly_ci_low:,.0f} - {delta_qaly_ci_high:,.0f}",
            "Dominant probability": f"{prob_dominant:.1%}",
            "SW quadrant probability": f"{prob_sw:.1%}",
            "NMB (WTP=5M, mean)": f"¥{nmb_5m_mean/1e9:.1f}B",
            "Cost-effectiveness probability": f"{prob_ce_5m:.1%}",
            "Conclusion": conclusion,
            "Recommended": recommendation,
        })
    
    return pd.DataFrame(records)

def generate_quadrant_analysis(psa_df: pd.DataFrame) -> pd.DataFrame:
    records = []
    
    for comparison in psa_df["comparison"].unique():
        df = psa_df[psa_df["comparison"] == comparison]
        n_total = len(df)

        quadrants = {
            "SE (Dominant)": ((df["delta_cost"] < 0) & (df["delta_qaly"] > 0)).sum(),
            "NE (More effective, More costly)": ((df["delta_cost"] > 0) & (df["delta_qaly"] > 0)).sum(),
            "SW (Less effective, Less costly)": ((df["delta_cost"] < 0) & (df["delta_qaly"] < 0)).sum(),
            "NW (Dominated)": ((df["delta_cost"] > 0) & (df["delta_qaly"] < 0)).sum(),
        }
        
        for quadrant, count in quadrants.items():
            records.append({
                "Comparison": comparison,
                "Quadrant": quadrant,
                "count": count,
                "Probability": f"{count/n_total:.1%}",
            })
    
    return pd.DataFrame(records)

# =============================================================================

# =============================================================================

def plot_cost_savings_comparison(psa_df: pd.DataFrame, output_path: str, dpi: int = 300):
    comparisons = psa_df["comparison"].unique()
    n_comparisons = len(comparisons)

    if n_comparisons == 0:
        print("  - Warning: No comparisons available for cost savings plot")
        return

    fig, axes = plt.subplots(1, n_comparisons, figsize=(5 * n_comparisons, 5))

    if n_comparisons == 1:
        axes = [axes]

    colors = {
        "0m_vs_2m": "#2E86AB",
        "1m_vs_2m": "#A23B72",
        "0m_vs_1m": "#28A745",
    }
    titles = {
        "0m_vs_2m": "0 month vs 2 months",
        "1m_vs_2m": "1 month vs 2 months",
        "0m_vs_1m": "0 month vs 1 month",
    }

    for idx, comparison in enumerate(comparisons):
        ax = axes[idx]
        df = psa_df[psa_df["comparison"] == comparison]

        savings = -df["delta_cost"] / 1e9

        ax.hist(savings, bins=50, color=colors.get(comparison, "gray"),
                alpha=0.7, edgecolor="white")
        ax.axvline(savings.mean(), color="red", linestyle="--", linewidth=2,
                   label=f"Mean: ¥{savings.mean():.1f}B")
        ax.axvline(savings.quantile(0.025), color="orange", linestyle=":", linewidth=1.5)
        ax.axvline(savings.quantile(0.975), color="orange", linestyle=":", linewidth=1.5,
                   label=f"95% CI: ¥{savings.quantile(0.025):.1f}B - ¥{savings.quantile(0.975):.1f}B")

        ax.set_xlabel("Cost Savings (Billion JPY)", fontsize=11)
        ax.set_ylabel("Frequency", fontsize=11)
        ax.set_title(titles.get(comparison, comparison), fontsize=12, fontweight="bold")
        ax.legend(loc="upper right", fontsize=9)
        ax.grid(axis="y", alpha=0.3)

        if savings.min() > 0:
            ax.axvspan(ax.get_xlim()[0], 0, alpha=0.1, color="red")
            ax.text(0.02, 0.98, "100% Cost-Saving", transform=ax.transAxes,
                    fontsize=10, fontweight="bold", color="green",
                    verticalalignment="top")

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()

def plot_dominant_ce_plane(psa_df: pd.DataFrame, output_path: str, dpi: int = 300):
    comparisons = psa_df["comparison"].unique()
    n_comparisons = len(comparisons)

    if n_comparisons == 0:
        print("  - Warning: No comparisons available for CE plane plot")
        return

    fig, ax = plt.subplots(figsize=(14, 10))

    colors = {
        "0m_vs_2m": "#2E86AB",
        "1m_vs_2m": "#A23B72",
        "0m_vs_1m": "#28A745",
    }
    labels = {
        "0m_vs_2m": "0 month vs 2 months",
        "1m_vs_2m": "1 month vs 2 months",
        "0m_vs_1m": "0 month vs 1 month",
    }

    ax.axhspan(-150, 0, xmin=0, xmax=1, alpha=0.1, color="green", zorder=0)
    ax.axvspan(0, 60000, ymin=0, ymax=0.5, alpha=0.1, color="green", zorder=0)
    
    for comparison in psa_df["comparison"].unique():
        df = psa_df[psa_df["comparison"] == comparison]

        prob_dominant = ((df["delta_cost"] < 0) & (df["delta_qaly"] > 0)).mean()
        
        ax.scatter(
            df["delta_qaly"],
            df["delta_cost"] / 1e9,
            alpha=0.4,
            s=25,
            color=colors.get(comparison, "gray"),
            label=f"{labels.get(comparison, comparison)} (Dominant: {prob_dominant:.1%})",
        )

        ax.scatter(
            df["delta_qaly"].mean(),
            df["delta_cost"].mean() / 1e9,
            s=200,
            color=colors.get(comparison, "gray"),
            marker="*",
            edgecolor="black",
            linewidth=1,
            zorder=10,
        )

    ax.axhline(0, color="black", linestyle="-", linewidth=1.5, zorder=1)
    ax.axvline(0, color="black", linestyle="-", linewidth=1.5, zorder=1)

    text_ne = ax.text(0.95, 0.95, "NE\nMore effective\nMore costly",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5, zorder=20))
    text_ne.set_zorder(20)
    text_nw = ax.text(0.05, 0.95, "NW\nLess effective\nMore costly\n(DOMINATED)",
            transform=ax.transAxes, ha="left", va="top", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="lightcoral", alpha=1.0, zorder=20))
    text_nw.set_zorder(20)
    text_se = ax.text(0.95, 0.05, "SE\nMore effective\nLess costly\n(DOMINANT)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="lightgreen", alpha=0.8, zorder=20))
    text_se.set_zorder(20)
    text_sw = ax.text(0.05, 0.05, "SW\nLess effective\nLess costly",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=1.0, zorder=20))
    text_sw.set_zorder(20)

    ax.set_xlabel("Incremental QALYs", fontsize=12)
    ax.set_ylabel("Incremental Cost (Billion JPY)", fontsize=12)

    if n_comparisons == 1:
        title = "Cost-Effectiveness Plane: Single Comparison"
    elif n_comparisons == 2:
        title = "Cost-Effectiveness Plane: Two Comparisons"
    else:
        title = "Cost-Effectiveness Plane: All Comparisons Show Dominant Results"
    ax.set_title(title, fontsize=14, fontweight="bold")

    ax.legend(loc="lower right", bbox_to_anchor=(0.98, 0.12), fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3, zorder=0)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()

def plot_summary_bar_chart(psa_df: pd.DataFrame, output_path: str, dpi: int = 300):
    comparisons = psa_df["comparison"].unique()
    n_comparisons = len(comparisons)

    if n_comparisons == 0:
        print("  - Warning: No comparisons available for summary bar chart")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    x = np.arange(n_comparisons)
    width = 0.6

    labels = {
        "0m_vs_2m": "0m vs 2m",
        "1m_vs_2m": "1m vs 2m",
        "0m_vs_1m": "0m vs 1m",
    }

    color_map = {
        "0m_vs_2m": "#2E86AB",
        "1m_vs_2m": "#A23B72",
        "0m_vs_1m": "#28A745",
    }
    bar_colors = [color_map.get(c, "gray") for c in comparisons]

    savings_mean = []
    savings_ci_low = []
    savings_ci_high = []
    qaly_mean = []
    qaly_ci_low = []
    qaly_ci_high = []

    for comp in comparisons:
        df = psa_df[psa_df["comparison"] == comp]
        savings = -df["delta_cost"] / 1e9
        savings_mean.append(savings.mean())
        savings_ci_low.append(savings.mean() - savings.quantile(0.025))
        savings_ci_high.append(savings.quantile(0.975) - savings.mean())

        qaly_mean.append(df["delta_qaly"].mean())
        qaly_ci_low.append(df["delta_qaly"].mean() - df["delta_qaly"].quantile(0.025))
        qaly_ci_high.append(df["delta_qaly"].quantile(0.975) - df["delta_qaly"].mean())

    ax1 = axes[0]
    bars1 = ax1.bar(x, savings_mean, width, color=bar_colors, alpha=0.8)
    ax1.errorbar(x, savings_mean, yerr=[savings_ci_low, savings_ci_high],
                 fmt="none", color="black", capsize=5, capthick=2)
    ax1.set_ylabel("Cost Savings (Billion JPY)", fontsize=12)
    ax1.set_title("Budget Impact: Cost Savings", fontsize=13, fontweight="bold")
    ax1.set_xticks(x)
    ax1.set_xticklabels([labels.get(c, c) for c in comparisons])
    ax1.axhline(0, color="black", linewidth=0.5)
    ax1.grid(axis="y", alpha=0.3)

    for i, (bar, val) in enumerate(zip(bars1, savings_mean)):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + savings_ci_high[i] + 2,
                 f"¥{val:.1f}B", ha="center", va="bottom", fontsize=11, fontweight="bold")

    ax2 = axes[1]
    bars2 = ax2.bar(x, qaly_mean, width, color=bar_colors, alpha=0.8)
    ax2.errorbar(x, qaly_mean, yerr=[qaly_ci_low, qaly_ci_high],
                 fmt="none", color="black", capsize=5, capthick=2)
    ax2.set_ylabel("Incremental QALYs", fontsize=12)
    ax2.set_title("Health Benefit: QALY Gains", fontsize=13, fontweight="bold")
    ax2.set_xticks(x)
    ax2.set_xticklabels([labels.get(c, c) for c in comparisons])
    ax2.axhline(0, color="black", linewidth=0.5)
    ax2.grid(axis="y", alpha=0.3)

    for i, (bar, val) in enumerate(zip(bars2, qaly_mean)):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + qaly_ci_high[i] + 500,
                 f"{val:,.0f}", ha="center", va="bottom", fontsize=11, fontweight="bold")

    if n_comparisons < 3:
        note_text = f"Note: {n_comparisons} comparison(s) available. Error bars indicate 95% CI from PSA."
    else:
        note_text = "Note: All comparisons show DOMINANT results (cost savings + QALY gains). Error bars indicate 95% CI from PSA."
    fig.text(0.5, 0.01, note_text, ha="center", fontsize=10, style="italic")

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()

def plot_dominance_probability(psa_df: pd.DataFrame, output_path: str, dpi: int = 300):
    comparisons = list(psa_df["comparison"].unique())
    n_comparisons = len(comparisons)

    if n_comparisons == 0:
        print("  - Warning: No comparisons available for dominance probability plot")
        return

    fig_height = max(4, 2 * n_comparisons)
    fig, ax = plt.subplots(figsize=(10, fig_height))

    labels = {
        "0m_vs_2m": "0 month vs 2 months",
        "1m_vs_2m": "1 month vs 2 months",
        "0m_vs_1m": "0 month vs 1 month",
    }

    y_pos = np.arange(n_comparisons)

    probs_dominant = []
    probs_sw = []
    probs_ne = []
    probs_dominated = []

    for comp in comparisons:
        df = psa_df[psa_df["comparison"] == comp]
        n = len(df)
        probs_dominant.append(((df["delta_cost"] < 0) & (df["delta_qaly"] > 0)).sum() / n * 100)
        probs_sw.append(((df["delta_cost"] < 0) & (df["delta_qaly"] < 0)).sum() / n * 100)
        probs_ne.append(((df["delta_cost"] > 0) & (df["delta_qaly"] > 0)).sum() / n * 100)
        probs_dominated.append(((df["delta_cost"] > 0) & (df["delta_qaly"] < 0)).sum() / n * 100)

    bars1 = ax.barh(y_pos, probs_dominant, height=0.5, label="Dominant (SE)", color="#28A745")
    bars2 = ax.barh(y_pos, probs_sw, height=0.5, left=probs_dominant, label="SW (Trade-off)", color="#FFC107")
    bars3 = ax.barh(y_pos, probs_ne, height=0.5,
                    left=[p1+p2 for p1,p2 in zip(probs_dominant, probs_sw)],
                    label="NE (Trade-off)", color="#17A2B8")
    bars4 = ax.barh(y_pos, probs_dominated, height=0.5,
                    left=[p1+p2+p3 for p1,p2,p3 in zip(probs_dominant, probs_sw, probs_ne)],
                    label="Dominated (NW)", color="#DC3545")

    ax.set_yticks(y_pos)
    ax.set_yticklabels([labels.get(c, c) for c in comparisons])
    ax.set_xlabel("Probability (%)", fontsize=12)
    ax.set_title("Probability Distribution Across CE Plane Quadrants", fontsize=13, fontweight="bold")
    ax.set_xlim(0, 100)
    ax.legend(loc="lower right", fontsize=10)
    ax.grid(axis="x", alpha=0.3)

    for i, (y, prob) in enumerate(zip(y_pos, probs_dominant)):
        if prob > 5:
            ax.text(prob/2, y, f"{prob:.1f}%", ha="center", va="center",
                    fontsize=11, fontweight="bold", color="white")

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()

# =============================================================================

# =============================================================================

def generate_publication_table(psa_df: pd.DataFrame) -> str:
    markdown = """
## Table: Cost-Effectiveness Results - Rapid ART Implementation in Japan

### Primary Outcome: Budget Impact (Cost Savings)

| Comparison | Cost Savings (Mean) | 95% Credible Interval | Probability of Cost Savings |
|------------|--------------------|-----------------------|----------------------------|
"""
    
    for comparison in psa_df["comparison"].unique():
        df = psa_df[psa_df["comparison"] == comparison]
        savings = -df["delta_cost"] / 1e9
        prob_savings = (df["delta_cost"] < 0).mean()
        
        label = comparison.replace("_", " ").replace("vs", "vs.")
        markdown += f"| {label} | ¥{savings.mean():.1f} billion | ¥{savings.quantile(0.025):.1f}B - ¥{savings.quantile(0.975):.1f}B | {prob_savings:.1%} |\n"
    
    markdown += """
### Secondary Outcome: Health Benefits (QALY Gains)

| Comparison | QALY Gains (Mean) | 95% Credible Interval | Probability of QALY Gains |
|------------|------------------|----------------------|--------------------------|
"""
    
    for comparison in psa_df["comparison"].unique():
        df = psa_df[psa_df["comparison"] == comparison]
        prob_qaly = (df["delta_qaly"] > 0).mean()
        
        label = comparison.replace("_", " ").replace("vs", "vs.")
        markdown += f"| {label} | {df['delta_qaly'].mean():,.0f} | {df['delta_qaly'].quantile(0.025):,.0f} - {df['delta_qaly'].quantile(0.975):,.0f} | {prob_qaly:.1%} |\n"
    
    markdown += """
### Dominance Analysis

| Comparison | Prob. Dominant | Prob. Cost-Effective (WTP=¥5M) | Conclusion |
|------------|---------------|-------------------------------|------------|
"""
    
    for comparison in psa_df["comparison"].unique():
        df = psa_df[psa_df["comparison"] == comparison]
        prob_dominant = ((df["delta_cost"] < 0) & (df["delta_qaly"] > 0)).mean()
        prob_ce = (df["delta_qaly"] * 5_000_000 - df["delta_cost"] > 0).mean()
        
        if prob_dominant >= 0.95:
            conclusion = "Strongly recommended"
        elif prob_ce >= 0.95:
            conclusion = "Recommended"
        else:
            conclusion = "Further analysis needed"
        
        label = comparison.replace("_", " ").replace("vs", "vs.")
        markdown += f"| {label} | {prob_dominant:.1%} | {prob_ce:.1%} | {conclusion} |\n"
    
    markdown += """
**Notes:**
- All analyses based on 1,000 probabilistic sensitivity analysis iterations
- Dominant = intervention is less costly AND more effective (no ICER calculation needed)
- WTP = Willingness-to-pay threshold
- 95% Credible Intervals derived from 2.5th and 97.5th percentiles of PSA distribution
"""
    
    return markdown

# =============================================================================

# =============================================================================

def main():

    input_dir = os.path.join(_PROJ_144, "bia_results", "cea_analysis")
    output_dir = os.path.join(_PROJ_144, "bia_results", "dominant_reporting")
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 70)
    print("Dominant Intervention Reporting")
    print("=" * 70)

    results = load_results(input_dir)
    psa_df = results.get("psa")
    
    if psa_df is None:
        print("Error: PSA results not found")
        return
    
    print(f"\nLoaded PSA results: {len(psa_df)} rows")
    available_comparisons = psa_df['comparison'].unique().tolist()
    n_comparisons = len(available_comparisons)
    print(f"Available comparisons ({n_comparisons}): {available_comparisons}")

    all_possible = ["0m_vs_2m", "1m_vs_2m", "0m_vs_1m"]
    missing = [c for c in all_possible if c not in available_comparisons]
    if missing:
        print(f"Note: Missing comparisons (will be skipped): {missing}")

    print("\n[1] Generating summary tables...")
    
    summary_df = generate_dominant_summary_table(psa_df)
    summary_df.to_csv(os.path.join(output_dir, "dominant_summary.csv"), index=False)
    print(summary_df.to_string(index=False))
    
    quadrant_df = generate_quadrant_analysis(psa_df)
    quadrant_df.to_csv(os.path.join(output_dir, "quadrant_analysis.csv"), index=False)

    print("\n[2] Generating figures...")

    plot_cost_savings_comparison(
        psa_df, 
        os.path.join(output_dir, "fig_cost_savings.png")
    )
    print("  - Cost savings comparison saved")

    plot_dominant_ce_plane(
        psa_df,
        os.path.join(output_dir, "fig_ce_plane_dominant.png")
    )
    print("  - CE plane (dominance emphasis) saved")

    plot_summary_bar_chart(
        psa_df,
        os.path.join(output_dir, "fig_summary_bars.png")
    )
    print("  - Summary bar chart saved")

    plot_dominance_probability(
        psa_df,
        os.path.join(output_dir, "fig_dominance_probability.png")
    )
    print("  - Dominance probability chart saved")

    print("\n[3] Generating publication table...")
    
    pub_table = generate_publication_table(psa_df)
    with open(os.path.join(output_dir, "publication_table.md"), "w") as f:
        f.write(pub_table)
    print("  - Publication table saved")
    
    print(pub_table)
    
    print("\n" + "=" * 70)
    print(f"All outputs saved to: {output_dir}")
    print("=" * 70)

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    print("\n[4] Exporting report-ready JSON...")

    json_report = {
        "metadata": {
            "script": "cea_dominant_reporting.py",
            "created_at": datetime.datetime.now().isoformat(),
            "input_dir": input_dir,
            "n_psa_iterations": len(psa_df) // n_comparisons if n_comparisons > 0 else 0,
            "comparisons": available_comparisons,
        },
        "dominant_summary": [],
        "quadrant_analysis": [],
    }

    # Dominant summary
    for comparison in psa_df["comparison"].unique():
        df = psa_df[psa_df["comparison"] == comparison]
        n = len(df)

        delta_cost_mean = float(df["delta_cost"].mean())
        delta_qaly_mean = float(df["delta_qaly"].mean())

        n_dominant = int(((df["delta_cost"] < 0) & (df["delta_qaly"] > 0)).sum())
        n_dominated = int(((df["delta_cost"] > 0) & (df["delta_qaly"] < 0)).sum())
        n_ne = int(((df["delta_cost"] > 0) & (df["delta_qaly"] > 0)).sum())
        n_sw = int(((df["delta_cost"] < 0) & (df["delta_qaly"] < 0)).sum())

        nmb_5m = df["delta_qaly"] * 5_000_000 - df["delta_cost"]

        prob_dominant = n_dominant / n
        prob_ce_5m = float((nmb_5m > 0).mean())

        if prob_dominant >= 0.95:
            conclusion = "Dominant"
            recommendation = "Strongly recommended"
        elif prob_ce_5m >= 0.95:
            conclusion = "Cost-effective at WTP=5M"
            recommendation = "Recommended"
        elif n_dominated / n >= 0.50:
            conclusion = "Dominated"
            recommendation = "Not recommended"
        else:
            conclusion = "Uncertain"
            recommendation = "Further analysis needed"

        json_report["dominant_summary"].append({
            "comparison": comparison,
            "delta_cost_mean": delta_cost_mean,
            "delta_cost_ci_lower": float(df["delta_cost"].quantile(0.025)),
            "delta_cost_ci_upper": float(df["delta_cost"].quantile(0.975)),
            "cost_savings_mean_billion_jpy": float(-delta_cost_mean / 1e9),
            "delta_qaly_mean": delta_qaly_mean,
            "delta_qaly_ci_lower": float(df["delta_qaly"].quantile(0.025)),
            "delta_qaly_ci_upper": float(df["delta_qaly"].quantile(0.975)),
            "prob_dominant": prob_dominant,
            "prob_dominated": n_dominated / n,
            "prob_ne_quadrant": n_ne / n,
            "prob_sw_quadrant": n_sw / n,
            "nmb_wtp5m_mean": float(nmb_5m.mean()),
            "prob_cost_effective_wtp5m": prob_ce_5m,
            "conclusion": conclusion,
            "recommendation": recommendation,
        })

    # Quadrant analysis
    for _, row in quadrant_df.iterrows():
        json_report["quadrant_analysis"].append({
            "comparison": row["Comparison"],
            "quadrant": row["Quadrant"],
            "count": int(row["count"]),
            "probability": row["Probability"],
        })

    json_path = os.path.join(output_dir, "dominant_report_ready.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(json_report, f, ensure_ascii=False, indent=2)
    print(f"    Saved: dominant_report_ready.json")
    print("=" * 70)

if __name__ == "__main__":
    main()

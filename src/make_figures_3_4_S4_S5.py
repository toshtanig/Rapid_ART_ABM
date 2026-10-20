#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Redraw Figures 3 and 4 (main text) and Supplementary Figures S4 (PSA) and S5 (CEAC) without titles inside the
artwork and without a box around the graphs, following the journal's figure guidelines.

The data are the published outputs of the 200-seed analysis and of the probabilistic and one-way analyses; nothing
is recomputed:
  - prevention_dominance_per_seed.csv   per-seed incremental QALYs (prevention decomposition) and incremental
                                         cost under the budget-impact convention (off-ART ratio 0.30)   -> Figure 3
  - seed_*/breakeven/breakeven_report_ready.json   per-seed cumulative discounted cost difference      -> Figure 4b
  - dsa_results.csv                      one-way sensitivity analysis                                   -> Figure 4a
  - psa_results.csv, ceac_results.csv    probabilistic analysis and acceptability curve                -> Figures S4, S5

Run as-is in VS Code (no arguments); paths in CONFIG. Panel letters are drawn in the figure (lowercase); the panel
titles are given in the figure legends of the manuscript and Supplement.
"""
import os
import json
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================================
# CONFIG
# =============================================================================
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                                       # repository root
DATA = os.path.join(ROOT, "data")
PER_SEED = os.path.join(DATA, "prevention_dominance_per_seed.csv")
SEED_ROOT = os.path.join(ROOT, "outputs", "core_multiseed")
CEA_DIR = os.path.join(ROOT, "outputs", "cea_param_uncertainty")
OUT_MAIN = os.path.join(ROOT, "outputs", "figures")
OUT_SUPP = os.path.join(ROOT, "outputs", "figures")
DPI = 300
W_MULTI = 7.09      # 180 mm
W_SINGLE = 4.21     # 107 mm
WTP = 5_000_000
COLORS = {"0m_vs_1m": "#1f77b4", "0m_vs_2m": "#d62728", "1m_vs_2m": "#2ca02c"}
NAMES = {"0m_vs_1m": "0 vs 1 month", "0m_vs_2m": "0 vs 2 months", "1m_vs_2m": "1 vs 2 months"}
BE_KEY = {"0m_vs_1m": "0m vs 1m", "0m_vs_2m": "0m vs 2m"}
FIRST_YEAR = 2024
# =============================================================================

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
                     "legend.fontsize": 7.5, "axes.linewidth": 0.6})


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.3, linewidth=0.5)


def letter(ax, s, x=-0.14, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")


def ci(x):
    x = np.asarray(x, float)
    m = x.mean()
    h = 1.96 * x.std(ddof=1) / np.sqrt(len(x))
    return m, h


def save(fig, path_base):
    fig.savefig(path_base + ".png", dpi=DPI)
    try:
        fig.savefig(path_base + ".tiff", dpi=DPI, pil_kwargs={"compression": "tiff_lzw"})
    except Exception as e:  # noqa: BLE001
        print("TIFF export skipped for", path_base, e)
    plt.close(fig)
    print("saved", path_base + ".png")


# ----------------------------------------------------------------------------- Figure 3
def figure3():
    df = pd.read_csv(PER_SEED)
    fig = plt.figure(figsize=(W_MULTI, 8.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.6, 1.0], hspace=0.35, wspace=0.35)
    ax = fig.add_subplot(gs[0, :])
    for comp in ["0m_vs_1m", "0m_vs_2m"]:
        g = df[df.comparison == comp]
        ax.scatter(g.dQALY, g.dCost_JPY / 1e9, s=14, alpha=0.55, color=COLORS[comp], label=f"{NAMES[comp]} (per seed)", lw=0)
        ax.scatter([g.dQALY.mean()], [g.dCost_JPY.mean() / 1e9], s=180, marker="*", color=COLORS[comp],
                   edgecolor="black", lw=0.8, label=f"{NAMES[comp]} mean", zorder=5)
    xs = np.array([-1000, df.dQALY.max() * 1.05])
    ax.plot(xs, xs * WTP / 1e9, ls="--", color="green", lw=1, label="Willingness to pay ¥5 million per QALY")
    ax.axhline(0, color="black", lw=0.8)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Incremental QALYs (immediate ART minus comparator)")
    ax.set_ylabel("Incremental cost (¥ billion)")
    ax.legend(loc="lower left", frameon=False)
    clean(ax)
    letter(ax, "a", x=-0.08)
    for k, (col, ylab, scale, lab) in enumerate([("dQALY", "Incremental QALYs", 1.0, "b"),
                                                 ("dCost_JPY", "Cumulative cost saved (¥ billion)", -1e-9, "c")]):
        ax = fig.add_subplot(gs[1, k])
        xs, ms, hs, cs = [], [], [], []
        for comp in ["0m_vs_1m", "0m_vs_2m"]:
            g = df[df.comparison == comp]
            m, h = ci(g[col] * scale)
            xs.append(NAMES[comp]); ms.append(m); hs.append(h); cs.append(COLORS[comp])
        ax.bar(xs, ms, yerr=hs, color=cs, alpha=0.85, capsize=4, error_kw={"lw": 1})
        ax.set_ylabel(ylab)
        clean(ax)
        ax.grid(False, axis="x")
        ax.set_axisbelow(True)                 # gridlines behind the bars
        letter(ax, lab, x=-0.22)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.96, bottom=0.07)
    save(fig, os.path.join(OUT_MAIN, "Figure3_R1"))


# ----------------------------------------------------------------------------- Figure 4
def breakeven_matrix(comp):
    rows = []
    for d in sorted(glob.glob(os.path.join(SEED_ROOT, "seed_*"))):
        path = os.path.join(d, "breakeven", "breakeven_report_ready.json")
        if not os.path.exists(path):
            continue
        rep = json.load(open(path, encoding="utf-8"))
        c = rep.get("breakeven_comparisons", {}).get(BE_KEY[comp])
        if not c or not c.get("yearly_data"):
            continue
        rows.append([x["cumulative_diff_billion_jpy"] for x in c["yearly_data"]][:40])
    return np.array(rows, float)


def figure4():
    dsa = pd.read_csv(os.path.join(CEA_DIR, "dsa_results.csv"))
    fig = plt.figure(figsize=(W_MULTI, 9.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.0], hspace=0.32, wspace=1.1)
    for k, comp in enumerate(["0m_vs_2m", "0m_vs_1m"]):
        ax = fig.add_subplot(gs[0, k])
        g = dsa[dsa.comparison == comp].copy()
        g["d"] = (g.nmb_5m - g.base_nmb_5m) / 1e9
        piv = g.pivot_table(index="parameter", columns="value_type", values="d")
        piv["span"] = (piv["high"] - piv["low"]).abs()
        piv = piv.sort_values("span")
        y = np.arange(len(piv))
        ax.barh(y, piv["low"], color="#7fa7cc", label="Low value", height=0.7)
        ax.barh(y, piv["high"], color="#f4a742", label="High value", height=0.7)
        ax.set_yticks(y)
        ax.set_yticklabels(piv.index, fontsize=6.5)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_xlabel("Change in NMB (¥ billion)", fontsize=7.5)
        ax.text(0.98, 0.02, NAMES[comp], transform=ax.transAxes, fontsize=7.5, ha="right", va="bottom")
        if k == 0:
            ax.legend(loc="lower right", frameon=False, bbox_to_anchor=(1.0, 0.08))
        clean(ax)
        ax.grid(False, axis="y")
        ax.set_axisbelow(True)                 # gridlines behind the bars
    fig.text(0.02, 0.975, "a", fontsize=11, fontweight="bold", va="top", ha="left")
    fig.text(0.02, 0.43, "b", fontsize=11, fontweight="bold", va="top", ha="left")
    ax = fig.add_subplot(gs[1, :])
    years = np.arange(FIRST_YEAR, FIRST_YEAR + 40)
    for comp in ["0m_vs_1m", "0m_vs_2m"]:
        mat = breakeven_matrix(comp)
        m = mat.mean(axis=0)
        h = 1.96 * mat.std(axis=0, ddof=1) / np.sqrt(mat.shape[0])
        ax.fill_between(years, m - h, m + h, color=COLORS[comp], alpha=0.2, lw=0)
        ax.plot(years, m, color=COLORS[comp], lw=1.6, label=f"{NAMES[comp]} (n = {mat.shape[0]} seeds)")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("Calendar year")
    ax.set_ylabel("Cumulative discounted cost difference,\nimmediate ART minus comparator (¥ billion)")
    ax.legend(loc="lower left", frameon=False)
    clean(ax)
    fig.subplots_adjust(left=0.25, right=0.97, top=0.97, bottom=0.06)
    save(fig, os.path.join(OUT_MAIN, "Figure4_R1"))


# ----------------------------------------------------------------------------- Figure S4 (PSA)
def figure_s4():
    psa = pd.read_csv(os.path.join(CEA_DIR, "psa_results.csv"))
    comps = [c for c in ["0m_vs_1m", "0m_vs_2m"] if (psa.comparison == c).any()]
    fig = plt.figure(figsize=(W_MULTI, 7.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.5, 1.0], hspace=0.35, wspace=0.3)
    ax = fig.add_subplot(gs[0, :])
    for comp in comps:
        g = psa[psa.comparison == comp]
        ax.scatter(g.delta_qaly, g.delta_cost / 1e9, s=8, alpha=0.35, color=COLORS[comp], lw=0, label=f"{NAMES[comp]} (per iteration)")
        ax.scatter([g.delta_qaly.mean()], [g.delta_cost.mean() / 1e9], s=160, marker="*", color=COLORS[comp],
                   edgecolor="black", lw=0.8, zorder=5, label=f"{NAMES[comp]} mean")
    xs = np.array([0, psa.delta_qaly.max() * 1.05])
    ax.plot(xs, xs * WTP / 1e9, ls="--", color="green", lw=1, label="Willingness to pay ¥5 million per QALY")
    ax.axhline(0, color="black", lw=0.8); ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Incremental QALYs (immediate ART minus comparator)")
    ax.set_ylabel("Incremental cost (¥ billion)")
    ax.legend(loc="upper left", frameon=False)
    clean(ax); letter(ax, "a", x=-0.08)
    for k, (col, xlab, scale) in enumerate([("delta_cost", "Incremental cost (¥ billion)", 1e-9), ("delta_qaly", "Incremental QALYs", 1.0)]):
        ax = fig.add_subplot(gs[1, k])
        for comp in comps:
            g = psa[psa.comparison == comp]
            ax.hist(g[col] * scale, bins=40, alpha=0.55, color=COLORS[comp], label=NAMES[comp])
        ax.set_xlabel(xlab); ax.set_ylabel("Iterations")
        clean(ax)
        ax.set_axisbelow(True)                 # gridlines behind the histogram bars
        if k == 0:
            ax.legend(frameon=False)
            letter(ax, "b", x=-0.22)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.96, bottom=0.07)
    save(fig, os.path.join(OUT_SUPP, "FigureS4_PSA_R1"))


# ----------------------------------------------------------------------------- Figure S5 (CEAC)
def figure_s5():
    ceac = pd.read_csv(os.path.join(CEA_DIR, "ceac_results.csv"))
    fig, ax = plt.subplots(figsize=(W_SINGLE, 3.1))
    for comp in ["0m_vs_1m", "0m_vs_2m", "1m_vs_2m"]:
        g = ceac[ceac.comparison == comp].sort_values("wtp")
        if g.empty:
            continue
        ax.plot(g.wtp / 1e6, g.probability_cost_effective, marker="o", ms=3, lw=1.4, color=COLORS[comp], label=NAMES[comp])
    ax.set_xlabel("Willingness to pay (¥ million per QALY)")
    ax.set_ylabel("Probability cost-effective")
    ax.set_ylim(-0.02, 1.05)
    ax.legend(frameon=False, loc="lower right")
    clean(ax)
    fig.tight_layout()
    save(fig, os.path.join(OUT_SUPP, "FigureS5_CEAC_R1"))


if __name__ == "__main__":
    figure3()
    figure4()
    figure_s4()
    figure_s5()

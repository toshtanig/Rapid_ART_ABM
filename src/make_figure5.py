#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Figure 5 for the R1 revision: observed-practice comparators and decomposition of the reform.

Panel a: infections averted (top) and net monetary benefit (bottom) of immediate ART against four
         comparators/scenarios (mean with 95% CI across seeds).
Panel b: factorial decomposition under the base-case and the observed ineligible wait.
Panel c: cumulative discounted cost difference (immediate minus comparator), mean with 95% CI across
         seeds, for the statutory-minimum, observed-practice, and conservative comparators.
Outputs Figure5.png and Figure5.tiff (300 dpi, LZW, 180 mm wide). Run as-is; paths in CONFIG."""

import os
import sys
import glob
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

# =============================================================================
# CONFIG
# =============================================================================
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _paths import BASE_MULTISEED, FIGURES_R1, REVISION_RUNS as RUNS
SUMMARY = os.path.join(RUNS, "summary", "revision_scenarios_summary.csv")
A9II = os.path.join(RUNS, "A9_observed_ineligible_wait", "empirical_pathways_plus_observed_tail")
A10 = os.path.join(RUNS, "A10_conservative_combined", "empirical_tail_behaviour_0.7")
OUT_PNG = os.path.join(FIGURES_R1, "Figure5.png")
OUT_TIFF = os.path.join(FIGURES_R1, "Figure5.tiff")
WIDTH_MM = 180
DPI = 300
FONT = "Arial"
MAX_SEEDS_BASE = 200
# =============================================================================

plt.rcParams.update({"font.family": FONT, "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.spines.top": False, "axes.spines.right": False})
C_BASE, C_OBSWAIT, C_OBS, C_CONS = "#4C72B0", "#8172B2", "#55A868", "#C44E52"


def row(df, analysis_prefix, setting):
    r = df[(df.analysis.str.startswith(analysis_prefix)) & (df.setting == setting)]
    assert len(r) == 1, (analysis_prefix, setting, len(r))
    return r.iloc[0]


def cum_curve(dirs, key="0m vs 1m", n=None):
    mats = []
    for d in sorted(dirs)[:n]:
        p = os.path.join(d, "breakeven", "breakeven_report_ready.json")
        if not os.path.exists(p):
            continue
        try:
            comp = json.load(open(p, encoding="utf-8"))["breakeven_comparisons"][key]["yearly_data"]
        except Exception:
            continue
        mats.append([x["cumulative_diff_billion_jpy"] for x in comp])
    m = np.array(mats)
    mean = m.mean(axis=0); se = m.std(axis=0, ddof=1) / np.sqrt(len(m))
    return mean, mean - 1.96 * se, mean + 1.96 * se, len(m)


def main():
    df = pd.read_csv(SUMMARY)
    base = row(df, "A9_", "base (published model, same seeds)")       # 100-seed same-seed base (for panel a use published 200-seed values)
    a9i = row(df, "A9_", "one_month_plus_observed_tail")
    a9ii = row(df, "A9_", "empirical_pathways_plus_observed_tail")
    a10 = row(df, "A10", "empirical_tail_behaviour_0.7")
    # published base case (200 seeds) for panel a
    pub = {"inf": (2991, 2871, 3112), "nmb": (99.7, 94.8, 104.6)}   # 200-seed base case; NMB CI from the per-seed NMB of the 200 published seeds (CEA convention; Table S10)

    fig = plt.figure(figsize=(WIDTH_MM / 25.4, 165 / 25.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.85], hspace=0.55, wspace=0.35)
    gsa = gs[0, 0].subgridspec(2, 1, hspace=0.35)
    ax_a1 = fig.add_subplot(gsa[0]); ax_a2 = fig.add_subplot(gsa[1])
    ax_b = fig.add_subplot(gs[0, 1]); ax_c = fig.add_subplot(gs[1, :])

    # ---------------- panel a
    labels = ["Statutory\nminimum\n(base case)", "Statutory\nminimum +\nobserved wait", "Observed\npractice", "Conservative\nscenario"]
    cols = [C_BASE, C_OBSWAIT, C_OBS, C_CONS]
    inf = [(pub["inf"][0], pub["inf"][1], pub["inf"][2]),
           (a9i.infections_averted_mean, a9i.infections_averted_lo, a9i.infections_averted_hi),
           (a9ii.infections_averted_mean, a9ii.infections_averted_lo, a9ii.infections_averted_hi),
           (a10.infections_averted_mean, a10.infections_averted_lo, a10.infections_averted_hi)]
    nmb = [pub["nmb"],
           (a9i.NMB_mean / 1e9, a9i.NMB_lo / 1e9, a9i.NMB_hi / 1e9),
           (a9ii.NMB_mean / 1e9, a9ii.NMB_lo / 1e9, a9ii.NMB_hi / 1e9),
           (a10.NMB_mean / 1e9, a10.NMB_lo / 1e9, a10.NMB_hi / 1e9)]
    x = np.arange(4)
    for ax, vals, ylab in ((ax_a1, inf, "Infections averted\nover 40 years"), (ax_a2, nmb, "Net monetary benefit\n(¥ billion)")):
        means = [v[0] for v in vals]; lo = [v[0] - v[1] for v in vals]; hi = [v[2] - v[0] for v in vals]
        ax.bar(x, means, color=cols, width=0.62, edgecolor="none")
        ax.errorbar(x, means, yerr=[lo, hi], fmt="none", ecolor="black", elinewidth=0.8, capsize=2)
        for xi, m, h in zip(x, means, hi):
            ax.text(xi, m + h + (max(means) * 0.03), f"{m:,.0f}" if m > 150 else f"{m:.1f}", ha="center", va="bottom", fontsize=6.5)
        ax.set_ylabel(ylab); ax.set_xticks(x); ax.set_ylim(0, max(means) * 1.32)
        ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)      # gridlines behind the bars
    ax_a1.set_xticklabels([]); ax_a2.set_xticklabels(labels, fontsize=6.5)
    ax_a1.text(-0.12, 1.08, "a", transform=ax_a1.transAxes, fontsize=10, fontweight="bold")

    # ---------------- panel b
    contr = {
        "base": {"Waiting interval removed": (760, 593, 926), "Thresholds removed": (2476, 2310, 2641), "Full reform": (3072, 2893, 3250)},
        "obs": {"Waiting interval removed": (539, 305, 772), "Thresholds removed": (660, 448, 872), "Full reform": (1199, 969, 1429)},
    }
    names = ["Waiting interval removed", "Thresholds removed", "Full reform"]
    xb = np.arange(3); w = 0.36
    for k, (key, col, lab) in enumerate((("base", "#7f7f7f", "Base case: ineligible people wait\nfor laboratory progression"), ("obs", C_OBS, "Observed: ineligible people wait\na median of about 11 months"))):
        vals = [contr[key][n] for n in names]
        means = [v[0] for v in vals]; lo = [v[0] - v[1] for v in vals]; hi = [v[2] - v[0] for v in vals]
        ax_b.bar(xb + (k - 0.5) * w, means, width=w, color=col, label=lab, edgecolor="none")
        ax_b.errorbar(xb + (k - 0.5) * w, means, yerr=[lo, hi], fmt="none", ecolor="black", elinewidth=0.8, capsize=2)
        for xi, m, h in zip(xb + (k - 0.5) * w, means, hi):
            ax_b.text(xi, m + h + 60, f"{m:,.0f}", ha="center", va="bottom", fontsize=6.3)
    ax_b.set_xticks(xb); ax_b.set_xticklabels(["Waiting interval\nremoved", "Thresholds\nremoved", "Full reform"], fontsize=6.5)
    ax_b.set_ylabel("Infections averted over 40 years"); ax_b.set_ylim(0, 4300); ax_b.grid(axis="y", alpha=0.3)
    ax_b.set_axisbelow(True)                                         # gridlines behind the bars
    ax_b.legend(loc="upper left", frameon=False, fontsize=6.3)
    ax_b.text(-0.12, 1.08, "b", transform=ax_b.transAxes, fontsize=10, fontweight="bold")

    # ---------------- panel c
    years = np.arange(1, 41)
    series = [
        (glob.glob(os.path.join(BASE_MULTISEED, "seed_*")), "0m vs 1m", C_BASE, "Statutory-minimum comparator", MAX_SEEDS_BASE),
        (glob.glob(os.path.join(A9II, "seed_*")), "0m vs 1m", C_OBS, "Observed-practice comparator", None),
        (glob.glob(os.path.join(A10, "seed_*")), "0m vs 1m", C_CONS, "Conservative scenario", None),
    ]
    for dirs, key, col, lab, n in series:
        mean, lo, hi, nn = cum_curve(dirs, key, n)
        yrs = years[: len(mean)]
        ax_c.fill_between(yrs, lo, hi, color=col, alpha=0.18, linewidth=0)
        idx = np.where(mean <= 0)[0]
        by = yrs[idx[0]] if len(idx) else None
        ax_c.plot(yrs, mean, color=col, lw=1.6, label=f"{lab} (n = {nn} seeds; break-even year {by})")
        if by is not None:
            ax_c.plot([by], [0], marker="o", color=col, ms=4)
    ax_c.axhline(0, color="black", lw=0.8)
    ax_c.set_xlabel("Year of the analytic horizon (2024 = year 1)"); ax_c.set_ylabel("Cumulative discounted cost difference,\nimmediate ART minus comparator (¥ billion)")
    ax_c.set_xlim(1, 40); ax_c.grid(alpha=0.3); ax_c.legend(loc="lower left", frameon=False, fontsize=6.5)
    ax_c.text(-0.05, 1.05, "c", transform=ax_c.transAxes, fontsize=10, fontweight="bold")

    os.makedirs(FIGURES_R1, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight", pad_inches=0.05)
    im = Image.open(OUT_PNG).convert("RGB")
    im.save(OUT_TIFF, "TIFF", compression="tiff_lzw", dpi=(DPI, DPI))
    print("saved", OUT_PNG, OUT_TIFF, "size px", im.size, "width mm %.0f" % (im.size[0] / DPI * 25.4))


if __name__ == "__main__":
    main()

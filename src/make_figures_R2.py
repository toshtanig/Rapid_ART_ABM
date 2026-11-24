#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Figures 2 to 5 of the main text and Supplementary Figures S2 to S8, from the post-processing outputs.

Nothing is read from the seed folders and nothing is recomputed here. Inputs (all under
outputs/revision_R2_runs/post/; the same files are provided in data/revision_R2/ of the repository, so the
figures can be drawn without running the model by setting POST to that folder):
  tables/figure_data_analytic_period.csv     yearly series by strategy (median, 95% prediction interval)  -> Figure 2
  summary/per_seed/<analysis>__<setting>__0m_vs_1m.csv (and 0m_vs_2m)   per-seed incremental outcomes     -> Figures 3, S7
  economics/<set>/dsa_results.csv            one-way sensitivity analysis                                 -> Figure 4a
  tables/breakeven_curve.csv                 cumulative discounted cost difference by year                -> Figure 4b
  tables/figure_data_burnin.csv              yearly series of the calibration period                      -> Figure S2
  tables/delay_histograms.csv                diagnosis-to-ART and eligibility-to-ART intervals            -> Figure S3
  economics/<set>/psa_results.csv, ceac_results.csv   probabilistic analysis, acceptability curve        -> Figures S4, S5
  tables/survival_weight.csv                 counterfactual survival weight                               -> Figure S6
  tables/subgroup_high_cd4.csv               untreated subgroup with high CD4 count and low viral load    -> Figure S8
  summary/per_seed/<key>.csv                 per-seed contrasts named in FIG5_PANEL_A and FIG5_PANEL_B    -> Figure 5a, 5b
  tables/breakeven_curve.csv                 curves of the settings named in FIG5_PANEL_C                 -> Figure 5c
Run summarize_R2.py, core_tables_R2.py and economics_R2.py first.

There are no titles inside the artwork and no box around the graphs; panel letters are lowercase. PNG and TIFF
(300 dpi, LZW) are written. A figure whose input is not yet available is skipped with a message.

TARGET = "published" draws the same figures from the published outputs (post/validation_published/), to check
the figures against the published ones; the output then goes to post/validation_published/figures/.

Run as-is in VS Code; all settings are in the CONFIG block, no command-line arguments.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MultipleLocator

# =============================================================================
# CONFIG
# =============================================================================
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _paths import FIGURES_R2, REVISION_RUNS
POST = os.path.join(REVISION_RUNS, "post")      # or os.path.join(os.path.dirname(HERE), "data", "revision_R2")

TARGET = "R2"                                    # "R2" or "published"
FIGURES = ["2", "3", "4", "5", "S2", "S3", "S4", "S5", "S6", "S7", "S8"]  # which figures to draw

# The comparison shown: setting folder of the runs (analysis, setting) and its name in post/economics/
ANALYSIS, SETTING, ECONOMICS_SET = "B0_core", "default", "statutory_minimum"
COMPARISONS = ["0m_vs_1m", "0m_vs_2m"]           # slots compared with the 0 month slot (one or two)
OUTPUT_TAG = "_R2"                               # appended to the file names (for example Figure3_R2.png)

OUT_MAIN = FIGURES_R2                                              # Figures 2 to 5
OUT_SUPP = FIGURES_R2                                              # Figures S2 to S8
FIG3_COST_COLUMN = "dCost_BIA"                   # cost of the cost-effectiveness plane: "dCost_BIA" (as published) or "dCost_CEA"
SURVIVAL_WEIGHT_BASE = 0.648
SURVIVAL_WEIGHT_MIN_SHOWN = 0.55
S3_MAX_MONTHS = 12
DPI = 300
W_MULTI = 7.09                                   # 180 mm
W_SINGLE = 4.21                                  # 107 mm
WTP = 5_000_000
FIRST_YEAR = 2024
FONT = ["Arial", "Liberation Sans", "DejaVu Sans"]

# Figure 5 (does not depend on ANALYSIS, SETTING and COMPARISONS above). Per-seed contrasts are the files
# summary/per_seed/<key>.csv written by summarize_R2.py; bars show the mean with its 95% CI across seeds.
FIG5_PANEL_A = [                                 # x label, per-seed contrast, colour
    ("Statutory\nminimum\n(base case)", "B0_core__default__0m_vs_1m", "#4C72B0"),
    ("Statutory\nminimum +\nobserved wait", "A9tail__full", "#8172B2"),
    ("Observed\npractice", "A9_observed_ineligible_wait__empirical_pathways_plus_observed_tail__0m_vs_1m", "#55A868"),
    ("Conservative\nscenario", "A10_conservative_combined__empirical_tail_behaviour_0.3__0m_vs_1m", "#C44E52"),
]
FIG5_PANEL_B = [                                 # legend label ({n} = number of seeds), colour, per-seed contrasts (waiting interval removed, thresholds removed, full reform)
    ("Base case: ineligible people wait\nfor laboratory progression (n = {n} seeds)", "#7f7f7f",
     ("A1__wait_given_criteria", "A1__criteria_given_wait", "A1__full")),
    ("Observed wait of ineligible people\n(median about 11 months; n = {n} seeds)", "#8172B2",
     ("A9tail__wait", "A9tail__criteria", "A9tail__full")),
    ("Observed-practice comparator (n = {n} seeds)", "#55A868",
     ("A9obs__wait", "A9obs__criteria", "A9obs__full")),
]
FIG5_PANEL_C = [                                 # legend label, (analysis, setting) in tables/breakeven_curve.csv, colour
    ("Statutory-minimum comparator", ("B0_core", "default"), "#4C72B0"),
    ("Observed-practice comparator", ("A9_observed_ineligible_wait", "empirical_pathways_plus_observed_tail"), "#55A868"),
    ("Conservative scenario", ("A10_conservative_combined", "empirical_tail_behaviour_0.3"), "#C44E52"),
]
FIG5_HEIGHT = 6.5                                # inches (165 mm)
# Figure S2 (calibration, 2004 to 2023): observed values that are not in tables/figure_data_burnin.csv
S2_OBSERVED_NOTIFICATIONS = {2024: 994}          # not a calibration target; shown for comparison
S2_OBSERVED_ON_ART = {2020: 26452}               # people retained in care at designated core hospitals
S2_UNDIAGNOSED_BAND = (10, 20)                   # estimated national range of the undiagnosed fraction (%)
# =============================================================================

SLOTS = ["0 month", "1 month", "2 months"]
SLOT_LABEL = {"0 month": "Immediate (0 month)", "1 month": "One-month delay", "2 months": "Two-month delay"}
SLOT_COLOR = {"0 month": "#4C72B0", "1 month": "#55A868", "2 months": "#C44E52"}
SLOT_COLOR_S = {"0 month": "#1f77b4", "1 month": "#2ca02c", "2 months": "#d62728"}
COMP_COLOR = {"0m_vs_1m": "#1f77b4", "0m_vs_2m": "#d62728", "1m_vs_2m": "#2ca02c"}
COMP_NAME = {"0m_vs_1m": "0 vs 1 month", "0m_vs_2m": "0 vs 2 months", "1m_vs_2m": "1 vs 2 months"}
COMP_KEY = {"0m_vs_1m": "0m vs 1m", "0m_vs_2m": "0m vs 2m"}

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": FONT, "font.size": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5, "axes.linewidth": 0.6})


def in_dir(*parts):
    root = POST if TARGET == "R2" else os.path.join(POST, "validation_published")
    return os.path.join(root, *parts)


def out_base(main, name):
    if TARGET == "R2":
        folder = OUT_MAIN if main else OUT_SUPP
    else:
        folder = os.path.join(POST, "validation_published", "figures")
    os.makedirs(folder, exist_ok=True)
    return os.path.join(folder, name + (OUTPUT_TAG if TARGET == "R2" else "_published"))


def read(*parts):
    p = in_dir(*parts)
    if not os.path.exists(p):
        raise FileNotFoundError(p)
    return pd.read_csv(p)


def of_setting(df):
    return df[(df["analysis"] == ANALYSIS) & (df["setting"] == SETTING)]


def per_seed(comp):
    return read("summary", "per_seed", f"{ANALYSIS}__{SETTING}__{comp}.csv").sort_values("seed")


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.3, linewidth=0.5)
    ax.set_axisbelow(True)


def letter(ax, s, x=-0.14, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")


def mean_half(x):
    x = np.asarray(x, float)
    return x.mean(), 1.96 * x.std(ddof=1) / np.sqrt(len(x))


def save(fig, path_base):
    fig.savefig(path_base + ".png", dpi=DPI)
    try:
        fig.savefig(path_base + ".tiff", dpi=DPI, pil_kwargs={"compression": "tiff_lzw"})
    except Exception as e:  # noqa: BLE001
        print("  TIFF export skipped for", os.path.basename(path_base), e)
    plt.close(fig)
    print("  saved", path_base + ".png")


# ----------------------------------------------------------------------------- Figure 2
def figure2():
    d = of_setting(read("tables", "figure_data_analytic_period.csv"))
    panels = [("new infections", "Infections per year"), ("new diagnoses", "Diagnoses per year"),
              ("on ART", "People on ART"), ("cumulative new infections since 2024", "Cumulative new infections"),
              ("deaths among people with HIV", "Deaths per year"), ("ART initiations", "ART initiations per year")]
    fs = 7
    rc = {"font.size": fs, "axes.labelsize": fs, "xtick.labelsize": fs - 0.5, "ytick.labelsize": fs - 0.5,
          "legend.fontsize": fs, "xtick.major.width": 0.6, "ytick.major.width": 0.6}
    with plt.rc_context(rc):
        fig, axes = plt.subplots(2, 3, figsize=(W_MULTI, 4.8))
        for k, (series, ylab) in enumerate(panels):
            ax = axes.ravel()[k]
            for slot in SLOTS:
                g = d[(d["series"] == series) & (d["slot"] == slot)].sort_values("year")
                if g.empty:
                    continue
                ax.fill_between(g["year"], g["pi95_lo"], g["pi95_hi"], color=SLOT_COLOR[slot], alpha=0.18, linewidth=0)
                ax.plot(g["year"], g["median"], color=SLOT_COLOR[slot], lw=1.2, label=SLOT_LABEL[slot])
            ax.set_xlabel("Calendar year")
            ax.set_ylabel(ylab)
            clean(ax)
            ax.xaxis.set_major_locator(MultipleLocator(10))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:,.0f}"))
            ax.text(-0.26, 1.03, "abcdef"[k], transform=ax.transAxes, fontsize=fs + 3, fontweight="bold", va="bottom", ha="left")
        handles, labels = axes.ravel()[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=len(labels), frameon=False, bbox_to_anchor=(0.5, 0.0))
        fig.tight_layout(rect=[0, 0.05, 1, 1])
        save(fig, out_base(True, "Figure2"))
    print("  Figure 2: median and 95% prediction interval across", int(d["n_seeds"].max()), "seeds")


# ----------------------------------------------------------------------------- Figure 3
def figure3():
    data = {c: per_seed(c) for c in COMPARISONS}
    fig = plt.figure(figsize=(W_MULTI, 8.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.6, 1.0], hspace=0.35, wspace=0.35)
    ax = fig.add_subplot(gs[0, :])
    xmax = 0.0
    for comp, g in data.items():
        ax.scatter(g["dQALY_decomposition"], g[FIG3_COST_COLUMN] / 1e9, s=14, alpha=0.55, color=COMP_COLOR[comp],
                   label=f"{COMP_NAME[comp]} (per seed)", lw=0)
        ax.scatter([g["dQALY_decomposition"].mean()], [g[FIG3_COST_COLUMN].mean() / 1e9], s=180, marker="*",
                   color=COMP_COLOR[comp], edgecolor="black", lw=0.8, label=f"{COMP_NAME[comp]} mean", zorder=5)
        xmax = max(xmax, g["dQALY_decomposition"].max())
    xs = np.array([-1000, xmax * 1.05])
    ax.plot(xs, xs * WTP / 1e9, ls="--", color="green", lw=1, label="Willingness to pay ¥5 million per QALY")
    ax.axhline(0, color="black", lw=0.8)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Incremental QALYs (immediate ART minus comparator)")
    ax.set_ylabel("Incremental cost (¥ billion)")
    ax.legend(loc="lower left", frameon=False)
    clean(ax)
    letter(ax, "a", x=-0.08)
    for k, (col, ylab, scale, lab) in enumerate([("dQALY_decomposition", "Incremental QALYs", 1.0, "b"),
                                                 (FIG3_COST_COLUMN, "Cumulative cost saved (¥ billion)", -1e-9, "c")]):
        ax = fig.add_subplot(gs[1, k])
        names, ms, hs, cs = [], [], [], []
        for comp, g in data.items():
            m, h = mean_half(g[col] * scale)
            names.append(COMP_NAME[comp]); ms.append(m); hs.append(h); cs.append(COMP_COLOR[comp])
        ax.bar(names, ms, yerr=hs, color=cs, alpha=0.85, capsize=4, error_kw={"lw": 1})
        ax.set_ylabel(ylab)
        clean(ax)
        ax.grid(False, axis="x")
        letter(ax, lab, x=-0.22)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.96, bottom=0.07)
    save(fig, out_base(True, "Figure3"))
    print("  Figure 3:", {c: len(g) for c, g in data.items()}, "seeds; cost column", FIG3_COST_COLUMN)


# ----------------------------------------------------------------------------- Figure 4
def figure4():
    dsa = read("economics", ECONOMICS_SET, "dsa_results.csv")
    curve = of_setting(read("tables", "breakeven_curve.csv"))
    comps = [c for c in ["0m_vs_2m", "0m_vs_1m"] if c in COMPARISONS and (dsa["comparison"] == c).any()]
    fig = plt.figure(figsize=(W_MULTI, 9.2))
    gs = fig.add_gridspec(2, len(comps), height_ratios=[1.35, 1.0], hspace=0.32, wspace=1.1)
    for k, comp in enumerate(comps):
        ax = fig.add_subplot(gs[0, k])
        g = dsa[dsa["comparison"] == comp].copy()
        g["d"] = (g["nmb_5m"] - g["base_nmb_5m"]) / 1e9
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
        ax.text(0.98, 0.02, COMP_NAME[comp], transform=ax.transAxes, fontsize=7.5, ha="right", va="bottom")
        if k == 0:
            ax.legend(loc="lower right", frameon=False, bbox_to_anchor=(1.0, 0.08))
        clean(ax)
        ax.grid(False, axis="y")
    fig.text(0.02, 0.975, "a", fontsize=11, fontweight="bold", va="top", ha="left")
    fig.text(0.02, 0.43, "b", fontsize=11, fontweight="bold", va="top", ha="left")
    ax = fig.add_subplot(gs[1, :])
    for comp in COMPARISONS:
        g = curve[curve["comparison"] == COMP_KEY[comp]].sort_values("year")
        if g.empty:
            continue
        ax.fill_between(g["calendar_year"], g["ci95_lo"], g["ci95_hi"], color=COMP_COLOR[comp], alpha=0.2, lw=0)
        ax.plot(g["calendar_year"], g["cumulative_bn_mean"], color=COMP_COLOR[comp], lw=1.6,
                label=f"{COMP_NAME[comp]} (n = {int(g['n_seeds'].iloc[0])} seeds)")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("Calendar year")
    ax.set_ylabel("Cumulative discounted cost difference,\nimmediate ART minus comparator (¥ billion)")
    ax.legend(loc="lower left", frameon=False)
    clean(ax)
    fig.subplots_adjust(left=0.25, right=0.97, top=0.97, bottom=0.06)
    save(fig, out_base(True, "Figure4"))


# ----------------------------------------------------------------------------- Figure 5
def contrast(key, column, scale=1.0):
    """(mean, half-width of the 95% CI, number of seeds) of one column of summary/per_seed/<key>.csv."""
    x = read("summary", "per_seed", key + ".csv")[column].to_numpy(float) / scale
    m, h = mean_half(x)
    return m, h, len(x)


def bar_values(ax, x, vals, colors, width, nd=0, fs=6.5):
    means = [v[0] for v in vals]
    halves = [v[1] for v in vals]
    ax.bar(x, means, color=colors, width=width, edgecolor="none")
    ax.errorbar(x, means, yerr=halves, fmt="none", ecolor="black", elinewidth=0.8, capsize=2)
    top = ax.get_ylim()[1]
    for xi, m, h in zip(x, means, halves):
        ax.text(xi, m + h + 0.015 * top, f"{m:,.{nd}f}", ha="center", va="bottom", fontsize=fs)


def figure5():
    fig = plt.figure(figsize=(W_MULTI, FIG5_HEIGHT))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.85], hspace=0.36, wspace=0.36)
    gsa = gs[0, 0].subgridspec(2, 1, hspace=0.3)
    ax_a1, ax_a2 = fig.add_subplot(gsa[0]), fig.add_subplot(gsa[1])
    ax_b, ax_c = fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])

    # panel a: infections averted and net monetary benefit against each comparator
    xa = np.arange(len(FIG5_PANEL_A))
    cols = [c for _, _, c in FIG5_PANEL_A]
    for ax, column, scale, ylab, nd in ((ax_a1, "infections_averted", 1.0, "Infections averted\nover 40 years", 0),
                                        (ax_a2, "NMB", 1e9, "Net monetary benefit\n(¥ billion)", 1)):
        vals = [contrast(key, column, scale) for _, key, _ in FIG5_PANEL_A]
        ax.set_ylim(0, max(m + h for m, h, _ in vals) * 1.22)
        bar_values(ax, xa, vals, cols, 0.62, nd)
        ax.set_ylabel(ylab)
        ax.set_xticks(xa)
        clean(ax)
        ax.grid(False, axis="x")
    ax_a1.set_xticklabels([])
    ax_a2.set_xticklabels([lab for lab, _, _ in FIG5_PANEL_A], fontsize=6.5)
    letter(ax_a1, "a", x=-0.2, y=1.04)

    # panel b: decomposition of the reform, by assumption about the comparator
    names = ["Waiting interval\nremoved", "Thresholds\nremoved", "Full reform"]
    xb = np.arange(len(names))
    n_groups = len(FIG5_PANEL_B)
    w = 0.8 / n_groups
    all_vals = [[contrast(key, "infections_averted") for key in keys] for _, _, keys in FIG5_PANEL_B]
    ax_b.set_ylim(0, max(m + h for vals in all_vals for m, h, _ in vals) * 1.55)
    for k, ((lab, col, _), vals) in enumerate(zip(FIG5_PANEL_B, all_vals)):
        pos = xb + (k - (n_groups - 1) / 2) * w
        bar_values(ax_b, pos, vals, col, w, 0, fs=5.6)
        ax_b.bar([0], [0], color=col, label=lab.format(n=vals[-1][2]))                # legend entry
    ax_b.set_xticks(xb)
    ax_b.set_xticklabels(names, fontsize=6.5)
    ax_b.set_ylabel("Infections averted over 40 years")
    ax_b.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax_b.legend(loc="upper left", frameon=False, fontsize=6.0, handlelength=1.2, borderaxespad=0.2)
    clean(ax_b)
    ax_b.grid(False, axis="x")
    letter(ax_b, "b", x=-0.2, y=1.02)
    ax_a1.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))

    # panel c: cumulative discounted cost difference (budget-impact cost convention)
    curve = read("tables", "breakeven_curve.csv")
    for lab, (an, st), col in FIG5_PANEL_C:
        g = curve[(curve["analysis"] == an) & (curve["setting"] == st) & (curve["comparison"] == "0m vs 1m")].sort_values("year")
        if g.empty:
            print(f"  Figure 5c: no curve for {an} / {st}")
            continue
        ax_c.fill_between(g["year"], g["ci95_lo"], g["ci95_hi"], color=col, alpha=0.18, linewidth=0)
        be = g["break_even_year_of_the_mean_curve"].iloc[0]
        be_txt = "none" if pd.isna(be) else f"{int(be)}"
        ax_c.plot(g["year"], g["cumulative_bn_mean"], color=col, lw=1.6,
                  label=f"{lab} (n = {int(g['n_seeds'].iloc[0])} seeds; break-even year {be_txt})")
        if not pd.isna(be):
            ax_c.plot([int(be)], [0], marker="o", color=col, ms=4)
    ax_c.axhline(0, color="black", lw=0.8)
    ax_c.set_xlim(1, 40)
    ax_c.set_xlabel(f"Year of the analytic horizon ({FIRST_YEAR} = year 1)")
    ax_c.set_ylabel("Cumulative discounted cost difference,\nimmediate ART minus comparator (¥ billion)")
    ax_c.legend(loc="lower left", frameon=False, fontsize=6.5)
    clean(ax_c)
    letter(ax_c, "c", x=-0.085, y=1.02)
    fig.subplots_adjust(left=0.105, right=0.985, top=0.96, bottom=0.075)
    save(fig, out_base(True, "Figure5"))


# ----------------------------------------------------------------------------- Figure S3
def figure_s3():
    d = of_setting(read("tables", "delay_histograms.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(W_MULTI, 2.6))
    for ax, kind, lab in zip(axes, ("diagnosis to ART", "eligibility to ART"), "ab"):
        for slot in SLOTS:
            g = d[(d["interval"] == kind) & (d["slot"] == slot)].sort_values("delay_months")
            if g.empty:
                continue
            x = np.concatenate([[0], g["delay_months"].to_numpy()])
            y = np.concatenate([[0], np.cumsum(g["count"].to_numpy()) / g["count"].sum()])
            ax.step(x, y, where="post", color=SLOT_COLOR_S[slot], lw=1.6, label=SLOT_LABEL[slot])
        ax.set_xlim(0, S3_MAX_MONTHS)
        ax.set_ylim(0, 1.02)
        ax.set_xlabel("Diagnosis-to-ART interval (months)" if kind.startswith("diagnosis") else "Eligibility-to-ART interval (months)")
        ax.set_ylabel("Cumulative proportion")
        clean(ax)
        letter(ax, lab, x=-0.2)
        if kind.startswith("diagnosis"):
            ax.legend(loc="lower right", frameon=False)
    fig.tight_layout(w_pad=2.0)
    save(fig, out_base(False, "FigureS3_delays"))
    print("  Figure S3: pooled over", int(d["n_seeds"].max()), "seeds")


# ----------------------------------------------------------------------------- Figure S2 (calibration)
def figure_s2():
    d = read("tables", "figure_data_burnin.csv")
    panels = [("notifications", "Notifications per year"), ("on ART", "People on ART"),
              ("new infections", "Infections per year"), ("undiagnosed fraction (%)", "Undiagnosed (%)")]
    line, obs = "#1f77b4", "#d62728"
    fig, axes = plt.subplots(2, 2, figsize=(W_MULTI, 5.6))
    n = 0
    for ax, (series, ylab), lab in zip(axes.ravel(), panels, "abcd"):
        g = d[d["series"] == series].sort_values("year")
        if g.empty:
            raise FileNotFoundError(f"series '{series}' in tables/figure_data_burnin.csv (run core_tables_R2.py first)")
        n = int(g["n_seeds"].max())
        if series.startswith("undiagnosed"):
            ax.axhspan(S2_UNDIAGNOSED_BAND[0], S2_UNDIAGNOSED_BAND[1], color="green", alpha=0.10, lw=0)
        ax.fill_between(g["year"], g["pi95_lo"], g["pi95_hi"], color=line, alpha=0.25, lw=0)
        ax.plot(g["year"], g["median"], color=line, lw=1.6)
        if series == "notifications":
            o = g.dropna(subset=["observed"])
            ax.plot(o["year"], o["observed"], "o", color=obs, ms=3.5)
            ax.plot(list(S2_OBSERVED_NOTIFICATIONS), list(S2_OBSERVED_NOTIFICATIONS.values()), "o", color=obs, ms=3.5)
        if series == "on ART":
            ax.plot(list(S2_OBSERVED_ON_ART), list(S2_OBSERVED_ON_ART.values()), "o", color=obs, ms=3.5)
        ax.set_xlabel("Calendar year")
        ax.set_ylabel(ylab)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
        clean(ax)
        letter(ax, lab, x=-0.2)
    handles = [Line2D([0], [0], color=line, lw=1.6, label="Model (median)"),
               Patch(facecolor=line, alpha=0.25, label=f"95% prediction interval ({n} seeds)"),
               Line2D([0], [0], marker="o", color=obs, lw=0, ms=3.5, label="Observed (Japan)"),
               Patch(facecolor="green", alpha=0.10, label=f"National estimate ({S2_UNDIAGNOSED_BAND[0]} to {S2_UNDIAGNOSED_BAND[1]}%)")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0.05, 1, 1), w_pad=2.0, h_pad=1.5)
    save(fig, out_base(False, "FigureS2_calibration"))


# ----------------------------------------------------------------------------- Figure S4 (PSA)
def figure_s4():
    psa = read("economics", ECONOMICS_SET, "psa_results.csv")
    comps = [c for c in COMPARISONS if (psa["comparison"] == c).any()]
    fig = plt.figure(figsize=(W_MULTI, 7.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.5, 1.0], hspace=0.35, wspace=0.3)
    ax = fig.add_subplot(gs[0, :])
    for comp in comps:
        g = psa[psa["comparison"] == comp]
        ax.scatter(g["delta_qaly"], g["delta_cost"] / 1e9, s=8, alpha=0.35, color=COMP_COLOR[comp], lw=0,
                   label=f"{COMP_NAME[comp]} (per iteration)")
        ax.scatter([g["delta_qaly"].mean()], [g["delta_cost"].mean() / 1e9], s=160, marker="*", color=COMP_COLOR[comp],
                   edgecolor="black", lw=0.8, zorder=5, label=f"{COMP_NAME[comp]} mean")
    xs = np.array([0, psa[psa["comparison"].isin(comps)]["delta_qaly"].max() * 1.05])
    ax.plot(xs, xs * WTP / 1e9, ls="--", color="green", lw=1, label="Willingness to pay ¥5 million per QALY")
    ax.axhline(0, color="black", lw=0.8)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Incremental QALYs (immediate ART minus comparator)")
    ax.set_ylabel("Incremental cost (¥ billion)")
    ax.legend(loc="upper left", frameon=False)
    clean(ax)
    letter(ax, "a", x=-0.08)
    for k, (col, xlab, scale) in enumerate([("delta_cost", "Incremental cost (¥ billion)", 1e-9), ("delta_qaly", "Incremental QALYs", 1.0)]):
        ax = fig.add_subplot(gs[1, k])
        for comp in comps:
            g = psa[psa["comparison"] == comp]
            ax.hist(g[col] * scale, bins=40, alpha=0.55, color=COMP_COLOR[comp], label=COMP_NAME[comp])
        ax.set_xlabel(xlab)
        ax.set_ylabel("Iterations")
        clean(ax)
        if k == 0:
            ax.legend(frameon=False)
            letter(ax, "b", x=-0.22)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.96, bottom=0.07)
    save(fig, out_base(False, "FigureS4_PSA"))


# ----------------------------------------------------------------------------- Figure S5 (CEAC)
def figure_s5():
    ceac = read("economics", ECONOMICS_SET, "ceac_results.csv")
    fig, ax = plt.subplots(figsize=(W_SINGLE, 3.1))
    for comp in COMPARISONS + (["1m_vs_2m"] if len(COMPARISONS) > 1 else []):
        g = ceac[ceac["comparison"] == comp].sort_values("wtp")
        if g.empty:
            continue
        ax.plot(g["wtp"] / 1e6, g["probability_cost_effective"], marker="o", ms=3, lw=1.4, color=COMP_COLOR[comp], label=COMP_NAME[comp])
    ax.set_xlabel("Willingness to pay (¥ million per QALY)")
    ax.set_ylabel("Probability cost-effective")
    ax.set_ylim(-0.02, 1.05)
    ax.legend(frameon=False, loc="lower right")
    clean(ax)
    fig.tight_layout()
    save(fig, out_base(False, "FigureS5_CEAC"))


# ----------------------------------------------------------------------------- Figure S6 (survival weight)
def figure_s6():
    d = read("tables", "survival_weight.csv")
    d = d[d["survival_weight"] >= SURVIVAL_WEIGHT_MIN_SHOWN]
    with plt.rc_context({"font.size": 7.5, "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7}):
        fig, axes = plt.subplots(1, 3, figsize=(W_MULTI, 2.45))
        panels = [("qaly_loss_per_averted_infection", None, None, "QALY loss per infection"),
                  ("dQALY_mean", "dQALY_lo", "dQALY_hi", "Incremental QALYs"),
                  ("NMB_mean_bn", "NMB_lo_bn", "NMB_hi_bn", "Net monetary benefit (¥ billion)")]
        for ax, (ycol, lo, hi, ylab), lab in zip(axes, panels, "abc"):
            for comp in COMPARISONS:
                g = d[d["comparison"] == COMP_KEY[comp]].sort_values("survival_weight")
                if g.empty:
                    continue
                ax.plot(g["survival_weight"], g[ycol], marker="o", ms=3.5, lw=1.4, color=COMP_COLOR[comp], label=COMP_NAME[comp])
                if lo is not None:
                    ax.fill_between(g["survival_weight"], g[lo], g[hi], color=COMP_COLOR[comp], alpha=0.18, lw=0)
            ax.axvline(SURVIVAL_WEIGHT_BASE, color="black", ls="--", lw=0.8)
            ax.set_xlabel("Counterfactual survival weight")
            ax.set_ylabel(ylab)
            clean(ax)
            ax.text(-0.3, 1.04, lab, transform=ax.transAxes, fontsize=10, fontweight="bold")
            if lab == "a":
                ax.legend(frameon=False, fontsize=7, loc="lower right")
        fig.tight_layout(w_pad=1.5)
        save(fig, out_base(False, "FigureS6_survival_weight"))


# ----------------------------------------------------------------------------- Figure S7 (convergence)
def running(x):
    x = np.asarray(x, float)
    k = np.arange(1, len(x) + 1)
    m = np.cumsum(x) / k
    var = np.maximum((np.cumsum(x ** 2) - k * m ** 2) / np.maximum(k - 1, 1), 0.0)
    se = np.sqrt(var / k)
    z = np.where(k > 30, 1.96, 2.1)        # approximate t quantile for small k
    lo, hi = m - z * se, m + z * se
    lo[:2] = np.nan
    hi[:2] = np.nan
    return k, m, lo, hi


def figure_s7():
    fig, axes = plt.subplots(1, 3, figsize=(W_MULTI, 2.4))
    lowest = 1.0
    for comp in COMPARISONS:
        g = per_seed(comp)
        col, name = COMP_COLOR[comp], COMP_NAME[comp]
        k, m, lo, hi = running(g["dQALY_decomposition"])
        axes[0].fill_between(k, lo, hi, color=col, alpha=0.18, lw=0)
        axes[0].plot(k, m, color=col, lw=1.4, label=name)
        k, m, lo, hi = running(g["infections_averted"])
        axes[1].fill_between(k, lo, hi, color=col, alpha=0.18, lw=0)
        axes[1].plot(k, m, color=col, lw=1.4, label=name)
        frac = np.cumsum(g["dominant"].astype(float).to_numpy()) / k
        axes[2].plot(k, frac, color=col, lw=1.4, label=name)
        lowest = min(lowest, frac.min())
    axes[0].axhline(0, color="black", lw=0.6)
    axes[0].set_ylabel("Incremental QALYs")
    axes[1].set_ylabel("Infections averted")
    axes[2].set_ylabel("Fraction dominant")
    axes[2].set_ylim(min(0.8, np.floor(lowest * 10) / 10), 1.01)
    for ax, lab in zip(axes, "abc"):
        ax.set_xlabel("Number of seeds")
        clean(ax)
        letter(ax, lab, x=-0.36, y=1.04)
    axes[0].legend(loc="lower right", frameon=False)
    fig.tight_layout(w_pad=2.0)
    save(fig, out_base(False, "FigureS7_convergence"))


# ----------------------------------------------------------------------------- Figure S8 (subgroup)
def figure_s8():
    d = read("tables", "subgroup_high_cd4.csv")
    fig, ax = plt.subplots(figsize=(W_SINGLE, 3.0))
    used = []
    for i, slot in enumerate(SLOTS):
        r = d[d["slot"] == slot]
        if r.empty:
            continue
        r = r.iloc[0]
        ax.bar(i, r["group_size_mean"], color=SLOT_COLOR_S[slot], alpha=0.85,
               yerr=[[r["group_size_mean"] - r["group_size_lo"]], [r["group_size_hi"] - r["group_size_mean"]]],
               capsize=4, error_kw={"lw": 1})
        used.append((i, slot))
    ax.set_xticks([i for i, _ in used])
    ax.set_xticklabels([SLOT_LABEL[s] for _, s in used])
    ax.set_ylabel("Group size (persons)")
    clean(ax)
    ax.grid(False, axis="x")
    fig.tight_layout()
    save(fig, out_base(False, "FigureS8_highcd4"))


DRAW = {"2": figure2, "3": figure3, "4": figure4, "5": figure5, "S2": figure_s2, "S3": figure_s3, "S4": figure_s4, "S5": figure_s5,
        "S6": figure_s6, "S7": figure_s7, "S8": figure_s8}


def main():
    if TARGET not in ("R2", "published"):
        raise SystemExit('TARGET must be "R2" or "published"')
    print(f"figures of {ANALYSIS} / {SETTING} ({TARGET})")
    for name in FIGURES:
        try:
            DRAW[name]()
        except FileNotFoundError as e:
            print(f"  Figure {name} skipped: input not found ({e})")


if __name__ == "__main__":
    main()

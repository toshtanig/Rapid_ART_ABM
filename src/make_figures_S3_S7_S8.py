#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Redraw Supplementary Figures S3 (delay distributions), S7 (Monte Carlo convergence) and S8 (high-CD4 subgroup
size) without titles inside the artwork and without a box around the graphs, following the journal's figure
guidelines. Nothing is recomputed: the figures are drawn from the published per-seed outputs
of the 200-seed analysis.

  S3  per-seed delay histograms scenario{0,1,2}_{diag,eligible}_to_art_delay_hist.csv (delays are whole months,
      so pooling the per-seed counts reproduces the pooled individual-level distribution exactly)
  S7  prevention_dominance_per_seed.csv (running mean, 95% CI and dominance rate, seeds in ascending order)
  S8  high_cd4_low_vl_longterm_noART_risk_summary.csv per seed (mean with 95% CI across seeds)

Run as-is in VS Code (no arguments); paths in CONFIG. The per-seed delay counts are cached in OUT_DIR so that a run
interrupted while reading the 1,200 histogram files can be resumed by running the script again.
"""
import os
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
OUT_DIR = os.path.join(ROOT, "outputs", "figures")
OUT_SUPP = os.path.join(ROOT, "outputs", "figures")
DPI = 300
W_MULTI = 7.09      # 180 mm
W_SINGLE = 4.21     # 107 mm
SCEN = {"scenario0": "Immediate (0 month)", "scenario1": "One-month delay", "scenario2": "Two-month delay"}
SCEN_COLOR = {"scenario0": "#1f77b4", "scenario1": "#2ca02c", "scenario2": "#d62728"}
COMP_COLOR = {"0m_vs_1m": "#1f77b4", "0m_vs_2m": "#d62728"}
COMP_NAME = {"0m_vs_1m": "0 vs 1 month", "0m_vs_2m": "0 vs 2 months"}
MAX_MONTHS_S3 = 12
# =============================================================================

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
                     "legend.fontsize": 7.5, "axes.linewidth": 0.6})


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.3, linewidth=0.5)


def letter(ax, s, x=-0.14, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")


def save(fig, path_base):
    fig.savefig(path_base + ".png", dpi=DPI)
    try:
        fig.savefig(path_base + ".tiff", dpi=DPI, pil_kwargs={"compression": "tiff_lzw"})
    except Exception as e:  # noqa: BLE001
        print("TIFF export skipped for", path_base, e)
    plt.close(fig)
    print("saved", path_base + ".png")


# ----------------------------------------------------------------------------- S3
def pooled_counts():
    """Per-seed delay counts (cached) -> pooled counts by kind, scenario, delay month."""
    os.makedirs(OUT_DIR, exist_ok=True)
    cache = os.path.join(OUT_DIR, "figS3_delay_counts_per_seed.csv")
    done = pd.read_csv(cache) if os.path.exists(cache) else pd.DataFrame(columns=["seed", "kind", "scenario", "delay_months", "count"])
    have = set(done.seed.astype(int)) if len(done) else set()
    rows = []
    for d in sorted(glob.glob(os.path.join(SEED_ROOT, "seed_*"))):
        sid = int(d.split("_")[-1])
        if sid in have:
            continue
        for kind in ("diag", "eligible"):
            for sc in SCEN:
                p = os.path.join(d, f"{sc}_{kind}_to_art_delay_hist.csv")
                if not os.path.exists(p):
                    continue
                h = pd.read_csv(p)
                for m, c in zip(h.delay_months, h["count"]):
                    rows.append({"seed": sid, "kind": kind, "scenario": sc, "delay_months": int(m), "count": int(c)})
        if rows and len({r["seed"] for r in rows}) % 20 == 0:        # checkpoint every 20 seeds
            pd.concat([done, pd.DataFrame(rows)]).to_csv(cache, index=False)
    allrows = pd.concat([done, pd.DataFrame(rows)]) if rows else done
    allrows.to_csv(cache, index=False)
    n_seeds = allrows.seed.nunique()
    pooled = allrows.groupby(["kind", "scenario", "delay_months"], as_index=False)["count"].sum()
    return pooled, n_seeds


def figure_s3():
    pooled, n = pooled_counts()
    fig, axes = plt.subplots(1, 2, figsize=(W_MULTI, 2.6))
    for ax, kind, lab in zip(axes, ("diag", "eligible"), "ab"):
        for sc in SCEN:
            g = pooled[(pooled.kind == kind) & (pooled.scenario == sc)].sort_values("delay_months")
            if g.empty:
                continue
            total = g["count"].sum()
            x = np.concatenate([[0], g.delay_months.to_numpy()])
            y = np.concatenate([[0], np.cumsum(g["count"].to_numpy()) / total])
            ax.step(x, y, where="post", color=SCEN_COLOR[sc], lw=1.6, label=SCEN[sc])
        ax.set_xlim(0, MAX_MONTHS_S3)
        ax.set_ylim(0, 1.02)
        ax.set_xlabel("Diagnosis-to-ART interval (months)" if kind == "diag" else "Eligibility-to-ART interval (months)")
        ax.set_ylabel("Cumulative proportion")
        clean(ax)
        letter(ax, lab, x=-0.2)
        if kind == "diag":
            ax.legend(loc="lower right", frameon=False)
    fig.tight_layout(w_pad=2.0)
    save(fig, os.path.join(OUT_SUPP, "FigureS3_delays_R1"))
    print(f"S3 pooled over {n} seeds")


# ----------------------------------------------------------------------------- S7
def running(x):
    x = np.asarray(x, float)
    k = np.arange(1, len(x) + 1)
    m = np.cumsum(x) / k
    sq = np.cumsum(x ** 2)
    var = np.maximum((sq - k * m ** 2) / np.maximum(k - 1, 1), 0.0)
    se = np.sqrt(var / k)
    z = np.where(k > 30, 1.96, 2.1)        # approximate t for small k
    lo, hi = m - z * se, m + z * se
    lo[:2] = np.nan; hi[:2] = np.nan
    return k, m, lo, hi


def figure_s7():
    df = pd.read_csv(PER_SEED)
    fig, axes = plt.subplots(1, 3, figsize=(W_MULTI, 2.4))
    for comp in ("0m_vs_1m", "0m_vs_2m"):
        g = df[df.comparison == comp].sort_values("seed")
        if g.empty:
            continue
        col, name = COMP_COLOR[comp], COMP_NAME[comp]
        k, m, lo, hi = running(g.dQALY)
        axes[0].fill_between(k, lo, hi, color=col, alpha=0.18, lw=0); axes[0].plot(k, m, color=col, lw=1.4, label=name)
        k, m, lo, hi = running(g.infections_averted)
        axes[1].fill_between(k, lo, hi, color=col, alpha=0.18, lw=0); axes[1].plot(k, m, color=col, lw=1.4, label=name)
        dom = g.dominant_SE.astype(float).to_numpy()
        axes[2].plot(k, np.cumsum(dom) / k, color=col, lw=1.4, label=name)
    axes[0].axhline(0, color="black", lw=0.6)
    axes[0].set_ylabel("Incremental QALYs")
    axes[1].set_ylabel("Infections averted")
    axes[2].set_ylabel("Fraction dominant")
    axes[2].set_ylim(0.8, 1.01)
    for ax, lab in zip(axes, "abc"):
        ax.set_xlabel("Number of seeds")
        clean(ax)
        letter(ax, lab, x=-0.36, y=1.04)
    axes[0].legend(loc="lower right", frameon=False)
    fig.tight_layout(w_pad=2.0)
    save(fig, os.path.join(OUT_SUPP, "FigureS7_convergence_R1"))


# ----------------------------------------------------------------------------- S8
def figure_s8():
    sizes = {s: [] for s in ("0 month", "1 month", "2 months")}
    for d in sorted(glob.glob(os.path.join(SEED_ROOT, "seed_*"))):
        p = os.path.join(d, "high_cd4_low_vl_longterm_noART_risk_summary.csv")
        if not os.path.exists(p):
            continue
        h = pd.read_csv(p)
        for s in sizes:
            r = h[h.scenario == s]
            if len(r):
                sizes[s].append(float(r.group_size.iloc[0]))
    fig, ax = plt.subplots(figsize=(W_SINGLE, 3.0))
    labels = {"0 month": "Immediate (0 month)", "1 month": "One-month delay", "2 months": "Two-month delay"}
    cols = {"0 month": "#1f77b4", "1 month": "#2ca02c", "2 months": "#d62728"}
    for i, s in enumerate(sizes):
        x = np.asarray(sizes[s], float)
        if x.size == 0:
            continue
        m = x.mean(); h = 1.96 * x.std(ddof=1) / np.sqrt(x.size)
        ax.bar(i, m, color=cols[s], alpha=0.85, yerr=h, capsize=4, error_kw={"lw": 1})
    ax.set_xticks(range(3)); ax.set_xticklabels([labels[s] for s in sizes])
    ax.set_ylabel("Group size (persons)")
    clean(ax); ax.grid(False, axis="x")
    ax.set_axisbelow(True)                     # gridlines behind the bars
    fig.tight_layout()
    save(fig, os.path.join(OUT_SUPP, "FigureS8_highcd4_R1"))
    print("S8 seeds:", {s: len(v) for s, v in sizes.items()})


if __name__ == "__main__":
    figure_s3()
    figure_s7()
    figure_s8()

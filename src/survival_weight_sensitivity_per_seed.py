#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""One-way sensitivity to the counterfactual survival weight (SURVIVAL_ADJ) with per-seed incremental costs
(Supplementary Table S4 and Figure S6, revision R1).

The dominance indicator and the net monetary benefit are computed seed by seed from the per-seed CEA
incremental cost (off-ART cost ratio 0.45, per_seed_cea_incremental_cost.csv) and the per-seed
prevention-decomposition QALY at each weight, as in the main analysis. survival_adj_sensitivity.py
(version 1.0.0) applies the mean incremental cost to every seed and is kept for reproducing the
original submission.

Run as-is in VS Code (no arguments). Settings below. Reuses compute_seed_comparison() from
prevention_qaly_dominance.py (SRC_DIR), setting its module-level SURVIVAL_ADJ for each grid value.
"""
import os
import sys
import glob
import numpy as np
import pandas as pd

# =============================================================================
# CONFIG
# =============================================================================
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                                   # repository root
SRC_DIR = HERE                                                 # prevention_qaly_dominance.py lives here
SEED_ROOT = os.path.join(ROOT, "outputs", "core_multiseed")     # the 200 published seeds (run_multiseed.py)
CEA_PER_SEED = os.path.join(ROOT, "data", "per_seed_cea_incremental_cost.csv")
OUT_DIR = os.path.join(ROOT, "outputs", "survival_weight")
FIG_OUT = os.path.join(OUT_DIR, "FigureS6_survival_weight.png")
GRID = [0.45, 0.55, 0.60, 0.648, 0.70, 0.75, 0.80, 0.85]
BASE_SA = 0.648
WTP = 5_000_000
COMPARISONS = [("0m", "1m", "0m_vs_1m"), ("0m", "2m", "0m_vs_2m")]
FIG_DPI = 300
FIG_WIDTH_IN = 7.09          # 180 mm; shown at 160 mm in the Supplement (about 340 dpi)
Z = 1.96                     # normal approximation, as in the other revision tables (200 seeds)
REUSE_CACHE = True           # reuse results/survival_weight_sensitivity_per_seed_values.csv if present (figure tweaks)
# =============================================================================

sys.path.insert(0, SRC_DIR)
import prevention_qaly_dominance as PV  # noqa: E402


def ci(x):
    x = np.asarray(x, float)
    m = x.mean()
    h = Z * x.std(ddof=1) / np.sqrt(len(x))
    return m, m - h, m + h


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    cea = pd.read_csv(CEA_PER_SEED)
    seed_dirs = sorted(glob.glob(os.path.join(SEED_ROOT, "seed_*")))
    rows, per_seed = [], []
    cache = os.path.join(OUT_DIR, "survival_weight_sensitivity_per_seed_values.csv")
    cached = pd.read_csv(cache) if (REUSE_CACHE and os.path.exists(cache)) else None
    for inter, comp, label in COMPARISONS:
        if cached is not None and (cached.comparison == label).any():
            per_seed.append(cached[cached.comparison == label].set_index("seed").sort_index())
            df = per_seed[-1]
            for sa in GRID:
                dq = df[f"dQALY_{sa}"]
                nmb = (WTP * dq - df.dCost_CEA_JPY) / 1e9
                m, lo, hi = ci(dq); nm, nlo, nhi = ci(nmb)
                rows.append({"comparison": label, "SURVIVAL_ADJ": sa, "is_base": sa == BASE_SA, "n_seeds": len(df),
                             "per_infection_QALY_mean": df[f"perinf_{sa}"].mean(),
                             "dQALY_mean": m, "dQALY_lo": lo, "dQALY_hi": hi,
                             "dCost_CEA_mean_JPY": df.dCost_CEA_JPY.mean(),
                             "NMB_bn_mean": nm, "NMB_bn_lo": nlo, "NMB_bn_hi": nhi,
                             "dominant_seeds": int(((dq > 0) & (df.dCost_CEA_JPY < 0)).sum()),
                             "dQALY_positive_seeds": int((dq > 0).sum()), "NMB_positive_seeds": int((nmb > 0).sum())})
            continue
        dcost = cea[cea.comparison == label].set_index("seed")["delta_cost_JPY_offART045"]
        recs = {}
        for d in seed_dirs:
            sid = int(d.split("_")[-1])
            if sid not in dcost.index:
                continue
            rec = {"seed": sid, "comparison": label, "dCost_CEA_JPY": float(dcost.loc[sid])}
            # dQALY and the per-infection QALY loss are exactly linear in SURVIVAL_ADJ for a given seed
            # (prevention component = sum over cohort years of averted x (U_NEG x residual years x SA - observed QALYs);
            # treatment component independent of SA), so two evaluations determine every grid value.
            PV.SURVIVAL_ADJ = 0.0
            r0 = PV.compute_seed_comparison(d, inter, comp)
            PV.SURVIVAL_ADJ = 1.0
            r1 = PV.compute_seed_comparison(d, inter, comp)
            if r0 is None or r1 is None:
                continue
            for sa in GRID:
                rec[f"dQALY_{sa}"] = r0["dQALY"] + sa * (r1["dQALY"] - r0["dQALY"])
                rec[f"perinf_{sa}"] = r0["per_infection_QALY_loss"] + sa * (r1["per_infection_QALY_loss"] - r0["per_infection_QALY_loss"])
            recs[sid] = rec
        df = pd.DataFrame(list(recs.values())).set_index("seed").sort_index()
        per_seed.append(df)
        for sa in GRID:
            dq = df[f"dQALY_{sa}"]
            nmb = (WTP * dq - df.dCost_CEA_JPY) / 1e9
            dominant = int(((dq > 0) & (df.dCost_CEA_JPY < 0)).sum())
            m, lo, hi = ci(dq)
            nm, nlo, nhi = ci(nmb)
            rows.append({"comparison": label, "SURVIVAL_ADJ": sa, "is_base": sa == BASE_SA, "n_seeds": len(df),
                         "per_infection_QALY_mean": df[f"perinf_{sa}"].mean(),
                         "dQALY_mean": m, "dQALY_lo": lo, "dQALY_hi": hi,
                         "dCost_CEA_mean_JPY": df.dCost_CEA_JPY.mean(),
                         "NMB_bn_mean": nm, "NMB_bn_lo": nlo, "NMB_bn_hi": nhi,
                         "dominant_seeds": dominant, "dQALY_positive_seeds": int((dq > 0).sum()),
                         "NMB_positive_seeds": int((nmb > 0).sum())})
    summary = pd.DataFrame(rows)
    summary.to_csv(os.path.join(OUT_DIR, "survival_weight_sensitivity_per_seed_costs.csv"), index=False)
    pd.concat(per_seed).to_csv(os.path.join(OUT_DIR, "survival_weight_sensitivity_per_seed_values.csv"))
    pd.set_option("display.width", 200)
    print(summary.to_string(index=False))

    # ---- Figure S6: three panels, no in-figure title, no box (top/right spines removed) ------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plot = summary[summary.SURVIVAL_ADJ >= 0.55]
    plt.rcParams.update({"font.size": 7.5, "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7})
    fig, axes = plt.subplots(1, 3, figsize=(FIG_WIDTH_IN, 2.45))
    colors = {"0m_vs_1m": "tab:blue", "0m_vs_2m": "tab:red"}
    names = {"0m_vs_1m": "0 vs 1 month", "0m_vs_2m": "0 vs 2 months"}
    panels = [("per_infection_QALY_mean", None, None, "QALY loss per infection"),
              ("dQALY_mean", "dQALY_lo", "dQALY_hi", "Incremental QALYs"),
              ("NMB_bn_mean", "NMB_bn_lo", "NMB_bn_hi", "Net monetary benefit (¥ billion)")]
    for ax, (ycol, lo, hi, ylab), letter in zip(axes, panels, "abc"):
        for label in colors:
            g = plot[plot.comparison == label]
            ax.plot(g.SURVIVAL_ADJ, g[ycol], marker="o", ms=3.5, lw=1.4, color=colors[label], label=names[label])
            if lo is not None:
                ax.fill_between(g.SURVIVAL_ADJ, g[lo], g[hi], color=colors[label], alpha=0.18, lw=0)
        ax.axvline(BASE_SA, color="black", ls="--", lw=0.8)
        ax.set_xlabel("Counterfactual survival weight")
        ax.set_ylabel(ylab)
        ax.grid(True, alpha=0.3)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.text(-0.3, 1.04, letter, transform=ax.transAxes, fontsize=10, fontweight="bold")
        if letter == "a":
            ax.legend(frameon=False, fontsize=7, loc="lower right")
    fig.tight_layout(w_pad=2.0)
    fig.savefig(FIG_OUT, dpi=FIG_DPI)
    try:
        fig.savefig(FIG_OUT[:-4] + ".tiff", dpi=FIG_DPI, pil_kwargs={"compression": "tiff_lzw"})
    except Exception as e:  # noqa: BLE001
        print("TIFF export skipped:", e)
    plt.close(fig)
    print("saved", os.path.join(OUT_DIR, "survival_weight_sensitivity_per_seed_costs.csv"), "and", FIG_OUT)


if __name__ == "__main__":
    main()

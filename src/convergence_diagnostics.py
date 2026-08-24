#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Monte Carlo convergence diagnostics across seeds."""

import os
import glob
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Lancet TIFF export hook (300 dpi, LZW; multipanel 180 mm / single 107 mm) ---
EXPORT_TIFF = True
if EXPORT_TIFF:
    try:
        import os as _os, sys as _sys
        _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
        import lancet_export
        lancet_export.enable_tiff_export()
    except Exception as _e:
        print("[lancet_export] disabled:", _e)

import prevention_qaly_dominance as P

# =============================================================================
# CONFIG
# =============================================================================
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# --- Data location: set BASE_144 to your local 144_scenario folder ---
BASE_144 = r"/path/to/20250315_169/144_scenario"
PROJ_169 = os.path.dirname(BASE_144)
BASE_169 = PROJ_169

RUN_MODE = "validation"
_MS_DIR = "core_multiseed_validation" if RUN_MODE == "validation" else "core_multiseed"
SCAN_DIRS = [
    os.path.join(BASE_169, "144_scenario", _MS_DIR),
]
COMPARISONS = [("0m", "1m"), ("0m", "2m")]
OUTPUT_DIR = os.path.join(BASE_169, "144_scenario", "convergence_diagnostics")
TARGET_HALFWIDTH_PCT = [5.0, 10.0]
DPI = 200
# =============================================================================

def _tval(n):
    tt = {2:12.71,3:4.303,4:3.182,5:2.776,6:2.571,7:2.447,8:2.365,9:2.306,10:2.262,
          11:2.228,12:2.201,13:2.179,14:2.160,15:2.145,16:2.131,18:2.110,20:2.093,
          25:2.064,30:2.045,40:2.021,60:2.000}
    if n in tt:
        return tt[n]
    if n > 60:
        return 1.96
    return tt[min(tt, key=lambda k: abs(k - n))]

def _seed_key(d):
    s = os.path.basename(d).replace("seed_", "")
    try:
        return (0, int(s))
    except ValueError:
        return (1, s)

def collect_per_seed():
    dirs = []
    for base in SCAN_DIRS:
        dirs += glob.glob(os.path.join(base, "seed_*"))
    dirs = sorted(set(dirs), key=_seed_key)
    data = {comp: [] for comp in COMPARISONS}
    for d in dirs:
        for comp in COMPARISONS:
            res = P.compute_seed_comparison(d, comp[0], comp[1])
            if res is None or res["dCost"] is None:
                continue
            nmb = res["dQALY"] * P.WTP - res["dCost"]
            data[comp].append({
                "seed": os.path.basename(d).replace("seed_", ""),
                "dQALY": res["dQALY"],
                "infections_averted": res["infections_averted"],
                "dCost": res["dCost"],
                "NMB": nmb,
                "dominant": (res["dQALY"] > 0 and res["dCost"] < 0),
            })
    return data

def running_stats(values):
    arr = np.asarray(values, float)
    N = len(arr)
    means, los, his, mcses = [], [], [], []
    for k in range(1, N + 1):
        x = arr[:k]
        m = x.mean()
        sd = x.std(ddof=1) if k > 1 else 0.0
        se = sd / math.sqrt(k) if k > 1 else 0.0
        t = _tval(k)
        means.append(m); los.append(m - t * se); his.append(m + t * se); mcses.append(se)
    return np.array(means), np.array(los), np.array(his), np.array(mcses)

def n_needed(sd, mean, pct):
    if mean == 0 or sd == 0:
        return float("nan")
    target = abs(mean) * pct / 100.0
    return (1.96 * sd / target) ** 2

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    data = collect_per_seed()
    colors = {("0m", "1m"): "#1f77b4", ("0m", "2m"): "#d62728"}

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for comp in COMPARISONS:
        rows = data[comp]
        if not rows:
            continue
        lab = f"{comp[0]} vs {comp[1]}"; col = colors[comp]
        ks = np.arange(1, len(rows) + 1)
        # ΔQALY
        m, lo, hi, _ = running_stats([r["dQALY"] for r in rows])
        axes[0].fill_between(ks, lo, hi, color=col, alpha=0.18)
        axes[0].plot(ks, m, color=col, lw=1.8, label=lab)

        m2, lo2, hi2, _ = running_stats([r["infections_averted"] for r in rows])
        axes[1].fill_between(ks, lo2, hi2, color=col, alpha=0.18)
        axes[1].plot(ks, m2, color=col, lw=1.8, label=lab)

        dom = np.array([r["dominant"] for r in rows], float)
        rate = np.cumsum(dom) / ks
        axes[2].plot(ks, rate, color=col, lw=1.8, marker="o", ms=3, label=lab)
    axes[0].axhline(0, color="k", lw=0.6)
    axes[0].set_title("Running mean ± 95% CI: incremental QALY"); axes[0].set_xlabel("Number of seeds"); axes[0].set_ylabel("ΔQALY")
    axes[1].set_title("Running mean ± 95% CI: infections averted"); axes[1].set_xlabel("Number of seeds"); axes[1].set_ylabel("Infections averted")
    axes[2].set_ylim(-0.05, 1.05); axes[2].set_title("Running dominance rate"); axes[2].set_xlabel("Number of seeds"); axes[2].set_ylabel("Fraction dominant")
    for ax in axes:
        ax.grid(True, alpha=0.3); ax.legend(fontsize=8, loc="lower right")
    fig.suptitle("Monte Carlo convergence diagnostics (running estimates vs number of seeds)")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(OUTPUT_DIR, "convergence.jpg"), dpi=DPI)
    plt.close(fig)

    table = []
    for comp in COMPARISONS:
        rows = data[comp]
        if not rows:
            continue
        N = len(rows)
        for key, lab in [("dQALY", "Incremental QALY"), ("infections_averted", "Infections averted"),
                         ("dCost", "Incremental cost (JPY)"), ("NMB", "NMB (JPY)")]:
            x = np.array([r[key] for r in rows], float)
            mean = x.mean(); sd = x.std(ddof=1) if N > 1 else 0.0
            mcse = sd / math.sqrt(N) if N > 1 else 0.0
            rec = {"comparison": f"{comp[0]}_vs_{comp[1]}", "metric": lab, "n_seeds": N,
                   "mean": round(mean, 1), "sd": round(sd, 1), "MCSE": round(mcse, 1),
                   "rel_MCSE_pct": round(100 * mcse / abs(mean), 2) if mean else float("nan")}
            for pct in TARGET_HALFWIDTH_PCT:
                rec[f"n_needed_{int(pct)}pct"] = round(n_needed(sd, mean, pct), 0)
            table.append(rec)
    pd.DataFrame(table).to_csv(os.path.join(OUTPUT_DIR, "mcse_table.csv"), index=False)

    print("=" * 70)
    print(f"Convergence diagnostics  ({sum(len(v) for v in data.values())//len(COMPARISONS)} seeds)")
    print(f"Saved: {os.path.join(OUTPUT_DIR, 'convergence.jpg')}")
    print(f"Saved: {os.path.join(OUTPUT_DIR, 'mcse_table.csv')}")
    print("=" * 70)
    print(pd.DataFrame(table).to_string(index=False))

if __name__ == "__main__":
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Aggregate per-seed outputs into the main projection and outcome figures."""

import os
import glob
import math
from io import StringIO
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
from matplotlib.ticker import MaxNLocator

def _int_year_axis(ax):
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

# =============================================================================

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
OUTPUT_DIR = os.path.join(BASE_169, "144_scenario", "multiseed_figures")

BAND = "ci95"            # "ci95" | "iqr" | "pi95" | "minmax"
MAKE_A = True
MAKE_B = True
MAKE_C = True

BURNIN_END_CAL = 2023
SCEN_COLS = ["0 month", "1 month", "2 months"]

SCEN_COLOR = {"0 month": "#4C72B0", "1 month": "#55A868", "2 months": "#C44E52"}

DPI = 200

OBS_HIV = {2004:780,2005:832,2006:952,2007:1082,2008:1126,2009:1021,2010:1075,2011:1056,
           2012:1002,2013:1106,2014:1091,2015:1006,2016:1011,2017:976,2018:940,2019:903,
           2020:750,2021:742,2022:632,2023:669,2024:662}
OBS_AIDS = {2004:385,2005:367,2006:406,2007:418,2008:431,2009:431,2010:469,2011:473,
            2012:447,2013:484,2014:455,2015:428,2016:437,2017:413,2018:377,2019:333,
            2020:345,2021:315,2022:252,2023:291,2024:332}

OBS_ONART_POINTS = {2020: 26452}
UNDIAG_REF_BAND = (10.0, 20.0)
# =============================================================================

def _read_csv(path):
    try:
        return pd.read_csv(path)
    except Exception:
        with open(path, "rb") as f:
            data = f.read().replace(b"\x00", b"")
        return pd.read_csv(StringIO(data.decode("utf-8", "ignore")))

def _tval(n):
    tt = {2:12.71,3:4.303,4:3.182,5:2.776,6:2.571,7:2.447,8:2.365,9:2.306,10:2.262,
          11:2.228,12:2.201,13:2.179,14:2.160,15:2.145,16:2.131,18:2.110,20:2.093,
          25:2.064,30:2.045,40:2.021,60:2.000}
    if n in tt:
        return tt[n]
    if n > 60:
        return 1.96
    keys = sorted(tt)
    return tt[min(keys, key=lambda k: abs(k - n))]

def band_stats(mat):
    mat = np.asarray(mat, float)
    center = np.full(mat.shape[1], np.nan)
    lo = np.full(mat.shape[1], np.nan)
    hi = np.full(mat.shape[1], np.nan)
    for j in range(mat.shape[1]):
        col = mat[:, j]
        col = col[~np.isnan(col)]
        n = len(col)
        if n == 0:
            continue
        if BAND == "ci95":
            m = col.mean(); sd = col.std(ddof=1) if n > 1 else 0.0; se = sd / math.sqrt(n) if n > 1 else 0.0
            t = _tval(n); center[j] = m; lo[j] = m - t * se; hi[j] = m + t * se
        elif BAND == "iqr":
            center[j] = np.median(col); lo[j] = np.percentile(col, 25); hi[j] = np.percentile(col, 75)
        elif BAND == "pi95":
            center[j] = np.median(col); lo[j] = np.percentile(col, 2.5); hi[j] = np.percentile(col, 97.5)
        else:  # minmax
            center[j] = np.median(col); lo[j] = col.min(); hi[j] = col.max()
    return center, lo, hi

def band_label(n):
    return {"ci95": f"mean ± 95% CI (n={n})", "iqr": f"median + IQR (n={n})",
            "pi95": f"median + 95% PI (n={n})", "minmax": f"median + min–max (n={n})"}[BAND]

def seed_dirs():
    out = []
    for base in SCAN_DIRS:
        out += sorted(glob.glob(os.path.join(base, "seed_*")))
    return out

def gather_burnin(dirs):
    metrics = ["new_diagnoses_year", "on_art", "new_infections_year", "hiv_positive",
               "diagnosed", "undiagnosed_hiv", "aids"]
    stacks = {m: [] for m in metrics}
    undiag_pct = []
    years = None
    for d in dirs:
        p = os.path.join(d, "burnin_summary.csv")
        if not os.path.exists(p):
            continue
        df = _read_csv(p).sort_values("burn_in_year")
        if years is None:
            years = df["burn_in_year"].to_numpy()
        for m in metrics:
            stacks[m].append(df[m].to_numpy(dtype=float))
        undiag_pct.append((df["undiagnosed_hiv"] / df["hiv_positive"] * 100.0).to_numpy())
    mats = {m: np.vstack(stacks[m]) for m in metrics if stacks[m]}
    mats["undiag_pct"] = np.vstack(undiag_pct) if undiag_pct else None
    cal = np.array([BURNIN_END_CAL - 20 + int(y) for y in years]) if years is not None else None
    return mats, cal

def fig_calibration(dirs, n):
    mats, cal = gather_burnin(dirs)
    if cal is None:
        print("  [A] no burnin data"); return
    obs_total = {y: OBS_HIV.get(y, np.nan) + OBS_AIDS.get(y, np.nan) for y in OBS_HIV}
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    rows = []

    def panel(ax, metric, title, ylab, obs_points=None, ref_band=None):
        c, lo, hi = band_stats(mats[metric])
        ax.fill_between(cal, lo, hi, color="#1f77b4", alpha=0.25, label=band_label(n))
        ax.plot(cal, c, color="#1f77b4", lw=2, label="Model")
        if obs_points:
            xs = sorted(obs_points); ys = [obs_points[x] for x in xs]
            ax.plot(xs, ys, "o", color="#d62728", ms=5, label="Japan observed")
        if ref_band:
            ax.axhspan(ref_band[0], ref_band[1], color="green", alpha=0.10, label="Japan estimate range")
        ax.set_title(title); ax.set_xlabel("Calendar year"); ax.set_ylabel(ylab)
        ax.grid(True, alpha=0.3); _int_year_axis(ax)
        for i, yy in enumerate(cal):
            rows.append({"calendar_year": int(yy), "metric": metric,
                         "center": c[i], "lo": lo[i], "hi": hi[i]})

    panel(axes[0, 0], "new_diagnoses_year", "New diagnoses per year",
          "New HIV+AIDS reports / year", obs_points=obs_total)
    panel(axes[0, 1], "on_art", "Currently on ART", "Persons on ART", obs_points=OBS_ONART_POINTS)
    panel(axes[1, 0], "new_infections_year", "New infections per year (modeled)", "New infections / year")

    c, lo, hi = band_stats(mats["undiag_pct"])
    ax = axes[1, 1]
    ax.fill_between(cal, lo, hi, color="#1f77b4", alpha=0.25, label=band_label(n))
    ax.plot(cal, c, color="#1f77b4", lw=2, label="Model")
    ax.axhspan(UNDIAG_REF_BAND[0], UNDIAG_REF_BAND[1], color="green", alpha=0.10, label="Japan estimate 10–20%")
    ax.set_title("Undiagnosed fraction"); ax.set_xlabel("Calendar year"); ax.set_ylabel("Undiagnosed (%)")
    ax.grid(True, alpha=0.3); _int_year_axis(ax)
    for i, yy in enumerate(cal):
        rows.append({"calendar_year": int(yy), "metric": "undiag_pct", "center": c[i], "lo": lo[i], "hi": hi[i]})

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="#1f77b4", alpha=0.25, label=band_label(n)),
        Line2D([0], [0], color="#1f77b4", lw=2, label="Model"),
        Line2D([0], [0], marker="o", color="#d62728", lw=0, ms=5, label="Japan observed"),
        Patch(facecolor="green", alpha=0.10, label="Japan estimate (10–20%)"),
    ]
    fig.suptitle(f"Calibration: model vs Japan surveillance, 2004–2023  [{band_label(n)}]")
    fig.legend(handles=legend_handles, loc="lower center", ncol=4, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=[0, 0.05, 1, 0.97])
    fig.savefig(os.path.join(OUTPUT_DIR, "A_calibration.jpg"), dpi=DPI)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(os.path.join(OUTPUT_DIR, "compiled_calibration.csv"), index=False)
    print("  [A] saved A_calibration.jpg / compiled_calibration.csv")

B_METRICS = [
    ("yearly_new_infections.csv",            "New infections / year"),
    ("yearly_new_diagnoses.csv",             "New diagnoses / year"),
    ("yearly_on_art_current.csv",            "Currently on ART"),
    ("yearly_cumulative_infected.csv",       "Cumulative infections"),
    ("yearly_deaths_among_hiv_positives.csv", "Deaths among HIV+ / year"),
    ("yearly_new_art_starts.csv",            "New ART starts / year"),
]

def gather_main(dirs, fname):
    per = {s: [] for s in SCEN_COLS}
    cal = None
    for d in dirs:
        p = os.path.join(d, fname)
        if not os.path.exists(p):
            continue
        df = _read_csv(p)

        cols_map = {}
        ok = True
        for s in SCEN_COLS:
            if s in df.columns:
                cols_map[s] = [s]
            elif (s + "_over200") in df.columns and (s + "_under200") in df.columns:
                cols_map[s] = [s + "_over200", s + "_under200"]
            else:
                ok = False
                break
        if not ok:
            continue
        if cal is None:
            cal = np.array([BURNIN_END_CAL + int(y) for y in df["Year"].to_numpy()])
        for s in SCEN_COLS:
            per[s].append(df[cols_map[s]].sum(axis=1).to_numpy(dtype=float))
    mats = {s: np.vstack(per[s]) for s in SCEN_COLS if per[s]}
    return mats, cal

def fig_projection(dirs, n):
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.ravel()
    rows = []
    for k, (fname, ylab) in enumerate(B_METRICS):
        ax = axes[k]
        mats, cal = gather_main(dirs, fname)
        if cal is None:
            ax.set_visible(False); continue
        for s in SCEN_COLS:
            if s not in mats:
                continue
            c, lo, hi = band_stats(mats[s])
            col = SCEN_COLOR[s]
            ax.fill_between(cal, lo, hi, color=col, alpha=0.18)
            ax.plot(cal, c, color=col, lw=1.8, label=s)
            for i, yy in enumerate(cal):
                rows.append({"calendar_year": int(yy), "scenario": s, "metric": fname.replace("yearly_", "").replace(".csv", ""),
                             "center": c[i], "lo": lo[i], "hi": hi[i]})
        ax.set_title(ylab.split(" (")[0]); ax.set_xlabel("Calendar year"); ax.set_ylabel(ylab)
        ax.grid(True, alpha=0.3); ax.legend(fontsize=8, loc="upper left"); _int_year_axis(ax)
    if len(B_METRICS) < len(axes):
        for j in range(len(B_METRICS), len(axes)):
            axes[j].set_visible(False)
    fig.suptitle(f"Forward projection by ART-timing scenario, 2024–2063  [{band_label(n)}]")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(os.path.join(OUTPUT_DIR, "B_projection.jpg"), dpi=DPI)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(os.path.join(OUTPUT_DIR, "compiled_projection.csv"), index=False)
    print("  [B] saved B_projection.jpg / compiled_projection.csv")

def pooled_delay(dirs, kind):
    suffix = "diag_to_art_delay_individual.csv" if kind == "diag" else "eligible_to_art_delay_individual.csv"
    pooled = {s: [] for s in ["scenario0", "scenario1", "scenario2"]}
    for d in dirs:
        for sc in pooled:
            p = os.path.join(d, f"{sc}_{suffix}")
            if os.path.exists(p):
                try:
                    df = _read_csv(p)
                    pooled[sc].append(df["delay_months"].to_numpy(dtype=float))
                except Exception:
                    pass
    return {sc: (np.concatenate(v) if v else np.array([])) for sc, v in pooled.items()}

def fig_delay_and_subgroup(dirs, n):
    label = {"scenario0": "0 month", "scenario1": "1 month", "scenario2": "2 months"}
    rows = []
    for kind, title in [("diag", "Diagnosis to ART"), ("eligible", "Eligibility to ART")]:
        pooled = pooled_delay(dirs, kind)
        fig, ax = plt.subplots(figsize=(7, 5))
        for sc, arr in pooled.items():
            if arr.size == 0:
                continue
            xs = np.sort(arr); ys = np.arange(1, len(xs) + 1) / len(xs)
            ax.step(xs, ys, where="post", color=SCEN_COLOR[label[sc]], lw=1.8, label=label[sc])
            rows.append({"source": kind, "scenario": label[sc], "n_pooled": int(arr.size),
                         "mean": float(arr.mean()), "median": float(np.median(arr)),
                         "p25": float(np.percentile(arr, 25)), "p75": float(np.percentile(arr, 75)),
                         "p90": float(np.percentile(arr, 90))})
        ax.set_xlim(0, 12); ax.set_xlabel("Delay (months)"); ax.set_ylabel("Cumulative probability")
        ax.set_title(f"{title} delay (pooled individuals across {n} seeds)")
        ax.grid(True, alpha=0.3); ax.legend(loc="lower right")
        fig.tight_layout()
        fig.savefig(os.path.join(OUTPUT_DIR, f"C_delay_cdf_{kind}.jpg"), dpi=DPI)
        plt.close(fig)
    pd.DataFrame(rows).to_csv(os.path.join(OUTPUT_DIR, "compiled_delay.csv"), index=False)

    sizes = {s: [] for s in SCEN_COLS}
    for d in dirs:
        p = os.path.join(d, "high_cd4_low_vl_longterm_noART_risk_summary.csv")
        if not os.path.exists(p):
            continue
        df = _read_csv(p)
        for s in SCEN_COLS:
            r = df[df["scenario"] == s]
            if len(r):
                sizes[s].append(float(r["group_size"].iloc[0]))
    fig, ax = plt.subplots(figsize=(7.5, 5)); srows = []
    ax.set_axisbelow(True)
    for i, s in enumerate(SCEN_COLS):
        if not sizes[s]:
            continue
        c, lo, hi = band_stats(np.array(sizes[s]).reshape(-1, 1))
        ax.bar(i, c[0], color=SCEN_COLOR[s], alpha=0.9, zorder=3,
               yerr=[[c[0] - lo[0]], [hi[0] - c[0]]], capsize=5,
               ecolor="black", error_kw={"zorder": 5})
        srows.append({"scenario": s, "metric": "high_cd4_lowvl_group_size",
                      "center": c[0], "lo": lo[0], "hi": hi[0]})
    ax.set_xticks(range(len(SCEN_COLS))); ax.set_xticklabels(SCEN_COLS)
    ax.set_ylabel("Group size (persons)")
    ax.set_title(f"High-CD4 / low-VL, long-term untreated group size\n({band_label(n)})", fontsize=11)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "C_high_cd4_group.jpg"), dpi=DPI)
    plt.close(fig)
    pd.DataFrame(srows).to_csv(os.path.join(OUTPUT_DIR, "compiled_subgroup.csv"), index=False)
    print("  [C] saved C_delay_cdf_*.jpg / C_high_cd4_group.jpg / compiled_delay.csv / compiled_subgroup.csv")

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    dirs = seed_dirs()
    n = len(dirs)
    print("=" * 70)
    print(f"multiseed_figures: {n} seeds, BAND={BAND}")
    print(f"output: {OUTPUT_DIR}")
    print("=" * 70)
    if n == 0:
        print("No seed_* folders found."); return
    if MAKE_A:
        fig_calibration(dirs, n)
    if MAKE_B:
        fig_projection(dirs, n)
    if MAKE_C:
        fig_delay_and_subgroup(dirs, n)
    print("Done.")

if __name__ == "__main__":
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Prevention-decomposition QALY estimator and dominance check."""

import os
import glob
import json
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

# =============================================================================

# =============================================================================
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# --- Input/output locations are repository-relative and match run_multiseed.py ---
HERE = os.path.dirname(os.path.abspath(__file__))
RUN_MODE = "production"   # published 200-seed analysis; set to "validation" for the 30-seed check
_MS_DIR = "core_multiseed_validation" if RUN_MODE == "validation" else "core_multiseed"

SCAN_DIRS = [
    os.path.join(HERE, "..", "outputs", _MS_DIR),
]

COMPARISONS = [("0m", "1m"), ("0m", "2m")]
SCEN_FILE = {"0m": "scenario0", "1m": "scenario1", "2m": "scenario2"}
BE_KEY    = {("0m", "1m"): "0m vs 1m", ("0m", "2m"): "0m vs 2m"}
PRETTY_COMP = {"0m_vs_1m": "0 vs 1 month", "0m_vs_2m": "0 vs 2 months", "1m_vs_2m": "1 vs 2 months"}
def pretty(c): return PRETTY_COMP.get(c, c)

N_YEARS      = 40
DISCOUNT_RATE = 0.02
U_NEG        = 0.94
SURVIVAL_ADJ = 0.648

WTP = 5_000_000

OUTPUT_DIR = os.path.join(HERE, "..", "outputs", "prevention_dominance_figures")
OUTPUT_CSV_PER_SEED = os.path.join(OUTPUT_DIR, "prevention_dominance_per_seed.csv")
OUTPUT_CSV_SUMMARY  = os.path.join(OUTPUT_DIR, "prevention_dominance_summary_stats.csv")
OUTPUT_MD  = os.path.join(OUTPUT_DIR, "prevention_dominance_summary.md")

MAKE_FIGURES   = True
OUTPUT_FIG_DIR = OUTPUT_DIR
FIG_DPI        = 200

# =============================================================================

def residual_discounted_years(y, n_years=N_YEARS, r=DISCOUNT_RATE):
    return float(sum((1.0 + r) ** (-(t - 1)) for t in range(int(y), n_years + 1)))

def _read_incident(seed_dir, scen_prefix):
    path = os.path.join(seed_dir, f"{scen_prefix}_incident_infection_cohort_summary.csv")
    if not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path)
    except Exception:
        with open(path, "rb") as f:
            data = f.read().replace(b"\x00", b"")
        from io import StringIO
        df = pd.read_csv(StringIO(data.decode("utf-8", "ignore")))
    return df[["cohort_year", "n_incident", "mean_discounted_qalys"]].copy()

def _read_breakeven_dcost(seed_dir, be_key):
    path = os.path.join(seed_dir, "breakeven", "breakeven_report_ready.json")
    if not os.path.exists(path):
        return None
    rep = json.load(open(path, encoding="utf-8"))
    comp = rep.get("breakeven_comparisons", {}).get(be_key)
    if comp is None:
        return None

    return float(comp["final_cumulative_billion_jpy"]) * 1e9

def compute_seed_comparison(seed_dir, intervention, comparator):
    inf_i = _read_incident(seed_dir, SCEN_FILE[intervention])     # 0m
    inf_c = _read_incident(seed_dir, SCEN_FILE[comparator])       # 1m or 2m
    if inf_i is None or inf_c is None:
        return None

    m = pd.merge(inf_i, inf_c, on="cohort_year", how="outer",
                 suffixes=("_i", "_c")).fillna(0.0).sort_values("cohort_year")

    prevention = 0.0
    treatment = 0.0
    averted_total = 0.0
    qaly_loss_weight_num = 0.0
    for _, row in m.iterrows():
        y = int(row["cohort_year"])
        n_i = float(row["n_incident_i"]); n_c = float(row["n_incident_c"])
        mq_i = float(row["mean_discounted_qalys_i"]); mq_c = float(row["mean_discounted_qalys_c"])
        averted = n_c - n_i
        common = min(n_i, n_c)
        resid = residual_discounted_years(y)
        qaly_uninf = U_NEG * resid * SURVIVAL_ADJ

        per_inf_loss = qaly_uninf - mq_c
        prevention += averted * per_inf_loss

        treatment += common * (mq_i - mq_c)
        averted_total += averted
        qaly_loss_weight_num += max(averted, 0.0) * per_inf_loss

    dqaly = prevention + treatment
    dcost = _read_breakeven_dcost(seed_dir, BE_KEY[(intervention, comparator)])
    avg_per_inf_loss = (qaly_loss_weight_num / averted_total) if averted_total > 0 else float("nan")
    return {
        "dQALY": dqaly,
        "dQALY_prevention": prevention,
        "dQALY_treatment": treatment,
        "infections_averted": averted_total,
        "per_infection_QALY_loss": avg_per_inf_loss,
        "dCost": dcost,
    }

def _ci(vals):
    x = np.array([v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))], float)
    n = len(x)
    if n == 0:
        return dict(n=0, mean=float("nan"), sd=float("nan"), lo=float("nan"), hi=float("nan"))
    m = x.mean(); sd = x.std(ddof=1) if n > 1 else 0.0; se = sd / math.sqrt(n) if n > 1 else 0.0
    tt = {2:12.71,3:4.303,4:3.182,5:2.776,6:2.571,7:2.447,8:2.365,9:2.306,10:2.262,
          11:2.228,12:2.201,13:2.179,14:2.160,15:2.145,20:2.093,30:2.045}
    t = tt.get(n, 1.96 if n > 30 else 2.262)
    return dict(n=n, mean=m, sd=sd, lo=m - t * se, hi=m + t * se)

def _fmt(metric, x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    if metric in ("dCost_JPY", "NMB_JPY"):
        return f"{x/1e9:,.1f}"
    if metric == "per_infection_QALY_loss":
        return f"{x:,.1f}"
    return f"{x:,.0f}"

MD_METRICS = [
    ("infections_averted",       "Infections averted",                    "count"),
    ("per_infection_QALY_loss",  "QALY loss per infection",               "QALY/infection"),
    ("dQALY_prevention",         "Incremental QALY: prevention",          "QALY"),
    ("dQALY_treatment",          "Incremental QALY: treatment timing",    "QALY"),
    ("dQALY",                    "Incremental QALY: total",               "QALY"),
    ("dCost_JPY",                "Incremental cost",                      "billion JPY"),
    ("NMB_JPY",                  "Net monetary benefit (WTP 5M/QALY)",    "billion JPY"),
]

COMP_LABEL = {
    "0m_vs_1m": "0 month vs 1 month (Rapid ART vs 1-month delay)",
    "0m_vs_2m": "0 month vs 2 months (Rapid ART vs 2-month delay)",
}

def _ci_desc(n):
    tt = {2:12.71,3:4.303,4:3.182,5:2.776,6:2.571,7:2.447,8:2.365,9:2.306,10:2.262,
          11:2.228,12:2.201,13:2.179,14:2.160,15:2.145,20:2.093,30:2.045}
    if n > 30:
        return "normal approximation, 1.96 x SEM"
    return f"t distribution, {tt.get(n, 2.045)} x SEM"

def write_markdown(df, path):
    import datetime
    L = []
    L.append("# Cost-effectiveness of Rapid ART: multi-seed prevention-decomposition summary")
    L.append("")
    L.append(f"Generated: {datetime.datetime.now():%Y-%m-%d %H:%M}")
    L.append("")
    n_seeds = df["seed"].nunique()
    L.append("## Methods summary")
    L.append("")
    L.append("Incremental outcomes were estimated from a stochastic individual-based "
             "microsimulation (hiv_rapid_art_model). To avoid the high variance of whole-population "
             "total-difference QALYs, incremental QALYs were obtained by a prevention "
             "decomposition: averted infections multiplied by the difference between the "
             "counterfactual HIV-negative QALYs and the observed post-infection QALYs of the "
             "incident-infection cohort, plus a treatment-timing component among persons infected "
             "in both arms. Incremental costs were taken from the break-even analysis (cumulative "
             "discounted cost difference over the 40-year horizon). Results were pooled across "
             f"{n_seeds} stochastic replicates (seeds) and are reported as mean with 95% "
             f"confidence intervals ({_ci_desc(n_seeds)}).")
    L.append("")
    L.append(f"Key parameters: HIV-negative utility {U_NEG}, discount rate {DISCOUNT_RATE:.0%}, "
             f"horizon {N_YEARS} years, survival adjustment {SURVIVAL_ADJ}, "
             f"willingness-to-pay JPY {WTP:,} per QALY. "
             "Utility source: Shiroiwa et al. 2021, Value in Health, Japanese EQ-5D-5L norms.")
    L.append("")
    L.append(f"Replicates analyzed: N = {n_seeds} seeds.")
    L.append("")

    for comp, g in df.groupby("comparison"):
        L.append(f"## Results: {COMP_LABEL.get(comp, comp)}")
        L.append("")
        L.append("| Outcome | Unit | Mean | 95% CI | CI excludes 0 |")
        L.append("|---|---|---|---|---|")
        for key, name, unit in MD_METRICS:
            c = _ci(g[key].tolist())
            mean = _fmt(key, c["mean"]); lo = _fmt(key, c["lo"]); hi = _fmt(key, c["hi"])
            sig = "Yes" if (c["n"] > 1 and (c["lo"] > 0 or c["hi"] < 0)) else "No"
            L.append(f"| {name} | {unit} | {mean} | {lo} to {hi} | {sig} |")
        L.append("")
        nd = int(g["dominant_SE"].sum()); nn = int(g["NMB_positive"].sum()); ns = len(g)
        L.append(f"Dominance: Rapid ART was dominant (lower cost and higher QALY) in "
                 f"**{nd}/{ns}** seeds; net monetary benefit positive in **{nn}/{ns}** seeds.")
        L.append("")

    L.append("## Robustness note")
    L.append("")
    _parts = []
    for _comp, _g in df.groupby("comparison"):
        _n = len(_g)
        _parts.append(f"{_comp}: QALY gain positive in {int((_g['dQALY'] > 0).sum())}/{_n}, "
                      f"dominant in {int(_g['dominant_SE'].sum())}/{_n}, "
                      f"NMB positive in {int(_g['NMB_positive'].sum())}/{_n}")
    L.append("Seed-level results: " + "; ".join(_parts) + ". "
             "The magnitude of the incremental QALYs depends on the assumed HIV-negative "
             "counterfactual (U_NEG and the survival adjustment); a sensitivity analysis varying "
             "these is recommended. Infections averted derive from paired simulation outputs, "
             "whereas cost results also reflect the specified treatment and costing assumptions.")
    L.append("")

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L))

def _read_breakeven_yearly(seed_dir, be_key):
    path = os.path.join(seed_dir, "breakeven", "breakeven_report_ready.json")
    if not os.path.exists(path):
        return None
    try:
        rep = json.load(open(path, encoding="utf-8"))
    except Exception:
        return None
    comp = rep.get("breakeven_comparisons", {}).get(be_key)
    if comp is None:
        return None
    yd = comp.get("yearly_data", [])
    if not yd:
        return None
    return np.array([d["cumulative_diff_billion_jpy"] for d in yd], float)

def _yearly_band(mat):
    mat = np.asarray(mat, float)
    n = mat.shape[0]
    m = np.nanmean(mat, axis=0)
    sd = np.nanstd(mat, axis=0, ddof=1) if n > 1 else np.zeros(mat.shape[1])
    se = sd / math.sqrt(n) if n > 1 else np.zeros(mat.shape[1])
    tt = {2:12.71,3:4.303,4:3.182,5:2.776,6:2.571,8:2.365,10:2.262,12:2.201,14:2.16,20:2.093}
    t = 1.96 if n > 30 else tt.get(n, 2.16)
    return m, m - t * se, m + t * se

def make_figures(df):
    os.makedirs(OUTPUT_FIG_DIR, exist_ok=True)
    comps = list(df["comparison"].unique())
    colors = {"0m_vs_1m": "#1f77b4", "0m_vs_2m": "#d62728"}

    fig, ax = plt.subplots(figsize=(7, 6))
    for comp in comps:
        g = df[df["comparison"] == comp]
        x = g["dQALY"].to_numpy(float); y = g["dCost_JPY"].to_numpy(float) / 1e9
        col = colors.get(comp)
        ax.scatter(x, y, s=25, alpha=0.55, color=col, label=f"{pretty(comp)} (per seed)")
        ax.scatter([x.mean()], [y.mean()], s=200, marker="*", color=col,
                   edgecolor="k", zorder=5, label=f"{pretty(comp)} mean")
    ax.axhline(0, color="k", lw=0.8); ax.axvline(0, color="k", lw=0.8)
    xl = ax.get_xlim()
    xs = np.linspace(min(0, xl[0]), max(0, xl[1]), 10)
    ax.plot(xs, (WTP / 1e9) * xs, "g--", lw=1, label=f"WTP {WTP/1e6:.0f}M/QALY")
    ax.set_xlabel("Incremental QALY (vs comparator)")
    ax.set_ylabel("Incremental cost (billion JPY)")
    ax.set_title("Cost-effectiveness plane (per-seed; lower-right = dominant)")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout(); fig.savefig(os.path.join(OUTPUT_FIG_DIR, "D_ce_plane.jpg"), dpi=FIG_DPI); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.set_axisbelow(True)
    for i, comp in enumerate(comps):
        c = _ci((df[df["comparison"] == comp]["NMB_JPY"] / 1e9).tolist())
        ax.bar(i, c["mean"], color=colors.get(comp), alpha=0.9, zorder=3,
               yerr=[[c["mean"] - c["lo"]], [c["hi"] - c["mean"]]], capsize=5)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(comps))); ax.set_xticklabels([pretty(c) for c in comps])
    ax.set_ylabel("Net monetary benefit (billion JPY)")
    ax.set_title(f"NMB (WTP {WTP/1e6:.0f}M/QALY), mean ± 95% CI")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUTPUT_FIG_DIR, "D_nmb_ci.jpg"), dpi=FIG_DPI); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, key, lab in [(axes[0], "dQALY", "Incremental QALY"),
                         (axes[1], "infections_averted", "Infections averted")]:
        data = [df[df["comparison"] == comp][key].to_numpy(float) for comp in comps]
        ax.boxplot(data, showmeans=True)
        ax.set_xticks(range(1, len(comps) + 1)); ax.set_xticklabels([pretty(c) for c in comps])
        ax.set_ylabel(lab); ax.set_title(lab + " across seeds"); ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUTPUT_FIG_DIR, "D_distributions.jpg"), dpi=FIG_DPI); plt.close(fig)

    seed_ds = []
    for base in SCAN_DIRS:
        seed_ds += sorted(glob.glob(os.path.join(base, "seed_*")))
    fig, ax = plt.subplots(figsize=(8, 5))
    for comp in comps:
        be_key = BE_KEY.get(tuple(comp.split("_vs_")))
        if be_key is None:
            continue
        series = [a[:40] for a in (_read_breakeven_yearly(d, be_key) for d in seed_ds)
                  if a is not None and len(a) >= 40]
        if not series:
            continue
        mat = np.vstack(series); m, lo, hi = _yearly_band(mat)
        cal = np.arange(2024, 2024 + mat.shape[1]); col = colors.get(comp)
        ax.fill_between(cal, lo, hi, color=col, alpha=0.2)
        ax.plot(cal, m, color=col, lw=1.8, label=pretty(comp))
    ax.axhline(0, color="k", lw=0.9)
    ax.set_xlabel("Calendar year")
    ax.set_ylabel("Cumulative cost difference, billion JPY\n(0 month − comparator)")
    ax.set_title("Break-even curve (mean ± 95% CI); crossing 0 = break-even")
    ax.grid(True, alpha=0.3); ax.legend(loc="upper right")
    fig.tight_layout(); fig.savefig(os.path.join(OUTPUT_FIG_DIR, "D_breakeven_curve.jpg"), dpi=FIG_DPI); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax in axes:
        ax.set_axisbelow(True)
    for i, comp in enumerate(comps):
        gq = _ci(df[df["comparison"] == comp]["dQALY"].tolist())
        cs = _ci((-df[df["comparison"] == comp]["dCost_JPY"] / 1e9).tolist())
        col = colors.get(comp)
        axes[0].bar(i, gq["mean"], color=col, alpha=0.9, zorder=3,
                    yerr=[[gq["mean"] - gq["lo"]], [gq["hi"] - gq["mean"]]], capsize=5,
                    ecolor="black", error_kw={"zorder": 5})
        axes[1].bar(i, cs["mean"], color=col, alpha=0.9, zorder=3,
                    yerr=[[cs["mean"] - cs["lo"]], [cs["hi"] - cs["mean"]]], capsize=5,
                    ecolor="black", error_kw={"zorder": 5})
    axes[0].set_title("Health benefit: QALY gains"); axes[0].set_ylabel("Incremental QALY (mean ± 95% CI)")
    axes[1].set_title("Budget impact: cost savings"); axes[1].set_ylabel("Cumulative cost saved, billion JPY (mean ± 95% CI)")
    for ax in axes:
        ax.set_xticks(range(len(comps))); ax.set_xticklabels([pretty(c) for c in comps])
        ax.axhline(0, color="k", lw=0.8); ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUTPUT_FIG_DIR, "D_qaly_cost_bars.jpg"), dpi=FIG_DPI); plt.close(fig)

    print(f"Saved figures to: {OUTPUT_FIG_DIR}")

def main():

    seed_dirs = []
    for base in SCAN_DIRS:
        seed_dirs += sorted(glob.glob(os.path.join(base, "seed_*")))

    rows = []
    for sd in seed_dirs:
        seed = os.path.basename(sd).replace("seed_", "")
        for intervention, comparator in COMPARISONS:
            res = compute_seed_comparison(sd, intervention, comparator)
            if res is None or res["dCost"] is None:
                continue
            nmb = res["dQALY"] * WTP - res["dCost"]
            dominant = (res["dQALY"] > 0) and (res["dCost"] < 0)
            rows.append({
                "seed": seed,
                "comparison": f"{intervention}_vs_{comparator}",
                "infections_averted": round(res["infections_averted"], 0),
                "per_infection_QALY_loss": round(res["per_infection_QALY_loss"], 3),
                "dQALY": round(res["dQALY"], 1),
                "dQALY_prevention": round(res["dQALY_prevention"], 1),
                "dQALY_treatment": round(res["dQALY_treatment"], 1),
                "dCost_JPY": round(res["dCost"], 0),
                "NMB_JPY": round(nmb, 0),
                "dominant_SE": dominant,
                "NMB_positive": nmb > 0,
            })

    if not rows:
        print("No completed seeds found (need incident_infection_cohort_summary + breakeven json).")
        return

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = pd.DataFrame(rows)

    summary_rows = []
    for comp, g in df.groupby("comparison"):
        for metric in ["infections_averted", "dQALY", "dQALY_prevention", "dQALY_treatment",
                       "dCost_JPY", "NMB_JPY", "per_infection_QALY_loss"]:
            c = _ci(g[metric].tolist())
            summary_rows.append({
                "comparison": comp, "metric": metric, "n": c["n"],
                "mean": round(c["mean"], 1), "ci_low": round(c["lo"], 1), "ci_high": round(c["hi"], 1),
                "sig_excl_0": (c["lo"] > 0 or c["hi"] < 0) if c["n"] > 1 else False,
            })
        n = len(g)
        summary_rows.append({"comparison": comp, "metric": "DOMINANT_SE_count",
                             "n": n, "mean": int(g["dominant_SE"].sum()),
                             "ci_low": "", "ci_high": "", "sig_excl_0": f"{int(g['dominant_SE'].sum())}/{n}"})
        summary_rows.append({"comparison": comp, "metric": "NMB_positive_count",
                             "n": n, "mean": int(g["NMB_positive"].sum()),
                             "ci_low": "", "ci_high": "", "sig_excl_0": f"{int(g['NMB_positive'].sum())}/{n}"})
    summary = pd.DataFrame(summary_rows)

    df.to_csv(OUTPUT_CSV_PER_SEED, index=False, encoding="utf-8")
    summary.to_csv(OUTPUT_CSV_SUMMARY, index=False, encoding="utf-8")

    write_markdown(df, OUTPUT_MD)

    if MAKE_FIGURES:
        make_figures(df)

    print("=" * 78)
    print(f"Prevention-decomposition ΔQALY + breakeven ΔCost → dominance")
    print(f"scanned seeds dirs: {SCAN_DIRS}")
    print(f"U_NEG={U_NEG}, discount={DISCOUNT_RATE}, N_YEARS={N_YEARS}, SURVIVAL_ADJ={SURVIVAL_ADJ}, WTP={WTP:,}")
    print("=" * 78)
    for comp, g in df.groupby("comparison"):
        n = len(g)
        print(f"\n[{comp}]  n={n} seeds")
        for metric, unit in [("infections_averted", "cases"), ("per_infection_QALY_loss", "QALY/infection"),
                             ("dQALY_prevention", "QALY"), ("dQALY_treatment", "QALY"),
                             ("dQALY", "QALY"), ("dCost_JPY", "JPY"), ("NMB_JPY", "JPY")]:
            c = _ci(g[metric].tolist())
            sig = "*excludes 0" if (c["n"] > 1 and (c["lo"] > 0 or c["hi"] < 0)) else "includes 0"
            print(f"  {metric:24s}: mean {c['mean']:,.1f} {unit}  95%CI[{c['lo']:,.1f}, {c['hi']:,.1f}]  {sig}")
        print(f"  -> dominant(cost saving & QALY gain): {int(g['dominant_SE'].sum())}/{n} seeds, "
              f"NMB>0: {int(g['NMB_positive'].sum())}/{n} seeds")
    print(f"\nSaved: {OUTPUT_CSV_PER_SEED}")
    print(f"Saved: {OUTPUT_CSV_SUMMARY}")
    print(f"Saved: {OUTPUT_MD}")
    print("=" * 78)

if __name__ == "__main__":
    main()

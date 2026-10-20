#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Summary tables for the revision (R1) analyses.

Part 1 (new simulation runs, from run_revision_scenarios.py):
  For every <analysis>/<setting> folder under revision_R1_runs/, computes the
  per-seed incremental outcomes of the "0 month" slot versus the "1 month" slot
  (infections and deaths averted, incremental QALYs by the prespecified
  prevention decomposition and by the whole-population direct estimate,
  incremental discounted cost under the CEA cost convention, net monetary
  benefit, dominance) and aggregates across seeds (mean, 95% CI, dominant
  seeds). The published base case (outputs/core_multiseed) restricted to the same
  seeds is reported alongside for comparison.
  A1: adds the factorial contrasts (waiting-interval effect, eligibility-
      criteria effect, interaction).
  A7: adds the "cost of inaction" relative to implementation in 2024.

Part 2 (post-processing of the published 200-seed run, no new simulation):
  Discount-rate sensitivity (0, 2, 3, 5 per cent per year for costs and QALYs)
  and a 20-year horizon, both using the whole-population direct estimator
  (per-seed yearly outputs in cea_report_ready.json), with 95% CIs across seeds.

Run as-is in VS Code; settings are in the CONFIG block. Outputs go to
outputs/revision_R1_runs/summary/ (and are also copied to PACKAGE_COPY_DIR if set).
"""

import os
import sys
import glob
import json
import math
import shutil
import numpy as np
import pandas as pd

# =============================================================================
# CONFIG
# =============================================================================
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import BASE_MULTISEED, REVISION_RUNS

BASE_MULTISEED_DIR = BASE_MULTISEED     # published 200-seed run (outputs/core_multiseed)
RUNS_ROOT = REVISION_RUNS
OUTPUT_DIR = os.path.join(REVISION_RUNS, "summary")
PACKAGE_COPY_DIR = None   # optional: another folder to which the summary files are copied

WTP = 5_000_000
DISCOUNT_RATES = [0.00, 0.02, 0.03, 0.05]     # Part 2, applied to both costs and QALYs
SHORT_HORIZON_YEARS = 20                      # Part 2
COMPARISONS_PART2 = [("0m", "1m"), ("0m", "2m")]

RUN_PART1 = True
RUN_PART2 = True

# Per-seed metrics are cached in summary/cache/ so that re-runs only read seeds
# that are new since the last run (reading from a synchronised or network folder is slow). Set
# MAX_NEW_SEEDS_PER_RUN to limit the work of one run (the tables are then
# written from the seeds available so far; run again until "cache complete").
CACHE_DIR = os.path.join(OUTPUT_DIR, "cache")
MAX_NEW_SEEDS_PER_RUN = None
# =============================================================================
_new_seed_budget = [None]

import prevention_qaly_dominance as P   # compute_seed_comparison (prevention decomposition), SURVIVAL_ADJ


def _ci(vals):
    x = np.array([v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))], float)
    n = len(x)
    if n == 0:
        return float("nan"), float("nan"), float("nan"), 0
    m = x.mean()
    se = x.std(ddof=1) / math.sqrt(n) if n > 1 else 0.0
    t = 1.96 if n > 30 else 2.045
    return m, m - t * se, m + t * se, n


def _fmt_ci(m, lo, hi, scale=1.0, nd=0):
    if isinstance(m, float) and math.isnan(m):
        return "NA"
    f = f"{{:.{nd}f}}"
    return f"{f.format(m/scale)} ({f.format(lo/scale)} to {f.format(hi/scale)})"


# -----------------------------------------------------------------------------
# Per-seed readers
# -----------------------------------------------------------------------------
def _read_cea(seed_dir):
    p = os.path.join(seed_dir, "cea", "cea_report_ready.json")
    if not os.path.exists(p):
        return None
    return json.load(open(p, encoding="utf-8"))


def _deaths(seed_dir, label):
    p = os.path.join(seed_dir, "heterogeneity_summary.csv")
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p)
    return float(df[df["Scenario"] == label]["Cum_Deaths_HIV_Pos"].sum())


def _breakeven_year(seed_dir, key):
    p = os.path.join(seed_dir, "breakeven", "breakeven_report_ready.json")
    if not os.path.exists(p):
        return None
    rep = json.load(open(p, encoding="utf-8"))
    comp = rep.get("breakeven_comparisons", {}).get(key)
    if comp is None:
        return None
    y = comp.get("break_even_year")
    return float(y) if y is not None else float("nan")


def seed_metrics(seed_dir, comp="1m"):
    """Incremental outcomes of the 0 month slot versus the comparator slot for one seed."""
    r = P.compute_seed_comparison(seed_dir, "0m", comp)
    cea = _read_cea(seed_dir)
    if r is None or cea is None:
        return None
    bc = cea["base_case_results"]
    if "0m" not in bc or comp not in bc:
        return None
    comp_label = {"1m": "1 month", "2m": "2 months"}[comp]
    dcost_cea = float(bc["0m"]["total_cost_discounted"]) - float(bc[comp]["total_cost_discounted"])
    dqaly_direct = float(bc["0m"]["total_qaly_discounted"]) - float(bc[comp]["total_qaly_discounted"])
    d0 = _deaths(seed_dir, "0 month")
    dc = _deaths(seed_dir, comp_label)
    deaths_averted = (dc - d0) if (d0 is not None and dc is not None) else float("nan")
    nmb = r["dQALY"] * WTP - dcost_cea
    return {
        "infections_averted": r["infections_averted"],
        "deaths_averted": deaths_averted,
        "dQALY_decomposition": r["dQALY"],
        "dQALY_direct": dqaly_direct,
        "dCost_CEA": dcost_cea,
        "dCost_BIA": r["dCost"] if r["dCost"] is not None else float("nan"),
        "NMB": nmb,
        "NMB_direct": dqaly_direct * WTP - dcost_cea,
        "dominant": bool(dcost_cea < 0 and r["dQALY"] > 0),
        "NMB_positive": bool(nmb > 0),
        "breakeven_year": _breakeven_year(seed_dir, P.BE_KEY[("0m", comp)]),
    }


METRICS = ["infections_averted", "deaths_averted", "dQALY_decomposition", "dQALY_direct",
           "dCost_CEA", "dCost_BIA", "NMB", "NMB_direct"]

DAYS_PER_MONTH = 365.25 / 12.0


def _delay_stats(seed_dir):
    """Diagnosis-to-ART interval of the two scenario slots (months) from the model's summary."""
    p = os.path.join(seed_dir, "diagnosis_to_ART_delay_summary.csv")
    out = {}
    if not os.path.exists(p):
        return out
    df = pd.read_csv(p)
    for slot, tag in (("0 month", "int"), ("1 month", "comp")):
        r = df[df["scenario"] == slot]
        if len(r) == 0:
            continue
        r = r.iloc[0]
        out[f"delay_{tag}_median_mo"] = float(r["median_delay_months"])
        out[f"delay_{tag}_p25_mo"] = float(r["p25_delay_months"])
        out[f"delay_{tag}_p75_mo"] = float(r["p75_delay_months"])
        out[f"delay_{tag}_mean_mo"] = float(r["mean_delay_months"])
    return out


def collect(setting_dir, comp="1m", seeds=None, cache_name=None):
    """Per-seed metrics for a setting folder, using (and extending) the cache."""
    cache_path = os.path.join(CACHE_DIR, f"{cache_name}.csv") if cache_name else None
    cached = pd.DataFrame()
    if cache_path and os.path.exists(cache_path):
        cached = pd.read_csv(cache_path)
    have = set(cached["seed"].astype(int).tolist()) if len(cached) else set()
    rows = []
    for d in sorted(glob.glob(os.path.join(setting_dir, "seed_*"))):
        seed = int(os.path.basename(d).replace("seed_", ""))
        if seeds is not None and seed not in seeds:
            continue
        if seed in have:
            continue
        if _new_seed_budget[0] is not None and _new_seed_budget[0] <= 0:
            continue
        m = seed_metrics(d, comp)
        if m is None:
            continue
        m["seed"] = seed
        m.update(_delay_stats(d))
        rows.append(m)
        if _new_seed_budget[0] is not None:
            _new_seed_budget[0] -= 1
    new = pd.DataFrame(rows)
    df = pd.concat([cached, new], ignore_index=True) if len(cached) else new
    if cache_path and len(new):
        os.makedirs(CACHE_DIR, exist_ok=True)
        df.to_csv(cache_path, index=False)
    if len(df) and seeds is not None:
        df = df[df["seed"].isin(list(seeds))].reset_index(drop=True)
    return df


def aggregate(df, label_cols):
    if df is None or len(df) == 0:
        return None
    out = dict(label_cols)
    out["n_seeds"] = int(len(df))
    for k in METRICS:
        m, lo, hi, n = _ci(df[k].tolist())
        out[k + "_mean"] = m
        out[k + "_lo"] = lo
        out[k + "_hi"] = hi
    out["dominant_seeds"] = f"{int(df['dominant'].sum())}/{len(df)}"
    out["NMB_positive_seeds"] = f"{int(df['NMB_positive'].sum())}/{len(df)}"
    be = df["breakeven_year"].dropna()
    out["breakeven_year_median"] = float(be.median()) if len(be) else float("nan")
    for tag in ("comp", "int"):
        if f"delay_{tag}_median_mo" in df.columns:
            med = df[f"delay_{tag}_median_mo"].mean() * DAYS_PER_MONTH
            p25 = df[f"delay_{tag}_p25_mo"].mean() * DAYS_PER_MONTH
            p75 = df[f"delay_{tag}_p75_mo"].mean() * DAYS_PER_MONTH
            mean = df[f"delay_{tag}_mean_mo"].mean() * DAYS_PER_MONTH
            out[f"diag_to_ART_{tag}_days"] = f"median {med:.0f} (IQR {p25:.0f} to {p75:.0f}); mean {mean:.0f}"
    return out


def pretty_table(rows):
    """Human-readable table (billion yen, rounded) for Markdown."""
    recs = []
    for r in rows:
        recs.append({
            **{k: v for k, v in r.items() if not any(k.endswith(s) for s in ("_mean", "_lo", "_hi"))},
            "Infections averted (95% CI)": _fmt_ci(r["infections_averted_mean"], r["infections_averted_lo"], r["infections_averted_hi"]),
            "Deaths averted (95% CI)": _fmt_ci(r["deaths_averted_mean"], r["deaths_averted_lo"], r["deaths_averted_hi"]),
            "Incremental QALYs, decomposition (95% CI)": _fmt_ci(r["dQALY_decomposition_mean"], r["dQALY_decomposition_lo"], r["dQALY_decomposition_hi"]),
            "Incremental QALYs, direct (95% CI)": _fmt_ci(r["dQALY_direct_mean"], r["dQALY_direct_lo"], r["dQALY_direct_hi"]),
            "Cost saving, billion yen (95% CI)": _fmt_ci(-r["dCost_CEA_mean"], -r["dCost_CEA_hi"], -r["dCost_CEA_lo"], 1e9, 1),
            "NMB, billion yen (95% CI)": _fmt_ci(r["NMB_mean"], r["NMB_lo"], r["NMB_hi"], 1e9, 1),
        })
    return pd.DataFrame(recs)


def _md_table(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join("" if (isinstance(v, float) and math.isnan(v)) else str(v) for v in row.tolist()) + " |")
    return "\n".join(lines)


# -----------------------------------------------------------------------------
# Part 1: new scenario runs
# -----------------------------------------------------------------------------
def part1():
    all_rows, md = [], []
    per_seed_frames = {}
    analyses = sorted(d for d in os.listdir(RUNS_ROOT)
                      if os.path.isdir(os.path.join(RUNS_ROOT, d)) and d.split("_")[0] in ("A1", "A2", "A4", "A5", "A7", "A8", "A9", "A10", "A11", "A12"))
    for an in analyses:
        an_dir = os.path.join(RUNS_ROOT, an)
        settings = sorted(d for d in os.listdir(an_dir) if os.path.isdir(os.path.join(an_dir, d)))
        seeds_union = set()
        setting_frames = {}
        for st in settings:
            df = collect(os.path.join(an_dir, st), "1m", cache_name=f"{an}__{st}")
            n_dirs = len(glob.glob(os.path.join(an_dir, st, "seed_*")))
            if len(df) == 0:
                print(f"  [{an}/{st}] no completed seeds yet")
                continue
            if len(df) < n_dirs:
                print(f"  [{an}/{st}] {len(df)}/{n_dirs} seed folders summarised so far (cache incomplete)")
            setting_frames[st] = df
            seeds_union |= set(df["seed"].tolist())
        if not setting_frames:
            continue
        # Published base case restricted to the same seeds
        base = collect(BASE_MULTISEED_DIR, "1m", seeds=seeds_union, cache_name="base_core_multiseed_144")
        per_seed_frames[(an, "base_same_seeds")] = base
        rows = []
        b = aggregate(base, {"analysis": an, "setting": "base (published model, same seeds)"})
        if b:
            rows.append(b)
        for st, df in setting_frames.items():
            per_seed_frames[(an, st)] = df
            note = ""
            cfgp = os.path.join(an_dir, st, "setting_config.json")
            if os.path.exists(cfgp):
                note = json.load(open(cfgp, encoding="utf-8")).get("note", "")
            a = aggregate(df, {"analysis": an, "setting": st})
            a["note"] = note
            rows.append(a)
        all_rows.extend(rows)
        md.append(f"## {an}\n\n" + _md_table(pretty_table(rows)) + "\n")

        if an == "A1_factorial" and {"criteria_retained_no_wait", "criteria_removed_one_month"} <= set(setting_frames) \
                and base is not None and len(base):
            md.append(factorial_contrasts(setting_frames["criteria_retained_no_wait"],
                                          setting_frames["criteria_removed_one_month"], base, all_rows))
        if an.startswith("A7") and base is not None and len(base):
            md.append(cost_of_inaction(setting_frames, base, all_rows))
        if an == "A9_observed_ineligible_wait" and {"one_month_plus_observed_tail",
                                                    "criteria_retained_no_wait_plus_observed_tail"} <= set(setting_frames):
            md.append(observed_tail_decomposition(setting_frames["one_month_plus_observed_tail"],
                                                  setting_frames["criteria_retained_no_wait_plus_observed_tail"], all_rows))

    if all_rows:
        pd.DataFrame(all_rows).to_csv(os.path.join(OUTPUT_DIR, "revision_scenarios_summary.csv"), index=False)
        for (an, st), df in per_seed_frames.items():
            df.to_csv(os.path.join(OUTPUT_DIR, f"per_seed_{an}_{st}.csv"), index=False)
        with open(os.path.join(OUTPUT_DIR, "revision_scenarios_summary.md"), "w", encoding="utf-8") as f:
            f.write("# Revision R1: additional scenario analyses (0 month slot versus 1 month slot)\n\n"
                    f"SURVIVAL_ADJ = {P.SURVIVAL_ADJ}; WTP = {WTP:,} yen per QALY; costs under the CEA convention "
                    "(off-ART cost ratio 0.45); 95% CI across seeds (normal approximation).\n\n" + "\n".join(md))
    return all_rows


def factorial_contrasts(df_AB, df_CD, df_AD, all_rows):
    """Factorial decomposition. Cells: A = current pathway (criteria, one-month wait);
    B = criteria retained, no wait; C = criteria removed, one-month wait; D = full reform.
    Setting 1 gives A - B, setting 2 gives C - D, the base run gives A - D (same seeds)."""
    m = df_AB.merge(df_CD, on="seed", suffixes=("_AB", "_CD")).merge(df_AD, on="seed")
    contrasts = {
        "Waiting-interval removal, criteria retained (A - B)": {k: m[k + "_AB"] for k in METRICS},
        "Waiting-interval removal, criteria removed (C - D)": {k: m[k + "_CD"] for k in METRICS},
        "Eligibility-criteria removal, one-month wait retained (A - C) = (A - D) - (C - D)": {k: m[k] - m[k + "_CD"] for k in METRICS},
        "Eligibility-criteria removal, no wait (B - D) = (A - D) - (A - B)": {k: m[k] - m[k + "_AB"] for k in METRICS},
        "Full reform (A - D), same seeds": {k: m[k] for k in METRICS},
        "Interaction (A - B) - (C - D)": {k: m[k + "_AB"] - m[k + "_CD"] for k in METRICS},
    }
    rows = []
    for name, series in contrasts.items():
        out = {"analysis": "A1_factorial", "setting": name, "n_seeds": int(len(m))}
        for k, s in series.items():
            mu, lo, hi, n = _ci(s.tolist())
            out[k + "_mean"], out[k + "_lo"], out[k + "_hi"] = mu, lo, hi
        dcost = series["dCost_CEA"]; dq = series["dQALY_decomposition"]
        out["dominant_seeds"] = f"{int(((dcost < 0) & (dq > 0)).sum())}/{len(m)}"
        out["NMB_positive_seeds"] = f"{int((series['NMB'] > 0).sum())}/{len(m)}"
        out["breakeven_year_median"] = float("nan")
        rows.append(out)
    all_rows.extend(rows)
    share = []
    full = rows[4]
    for r in rows[:4]:
        share.append({"contrast": r["setting"],
                      "share of full-reform infections averted (%)": round(100 * r["infections_averted_mean"] / full["infections_averted_mean"], 1),
                      "share of full-reform NMB (%)": round(100 * r["NMB_mean"] / full["NMB_mean"], 1)})
    return ("### A1 factorial contrasts (per-seed paired differences)\n\n" + _md_table(pretty_table(rows)) +
            "\n\n" + _md_table(pd.DataFrame(share)) + "\n")


def observed_tail_decomposition(df_AD, df_AB, all_rows):
    """Decomposition under the observed ineligible wait: A' - D (total, setting i),
    A' - B' (waiting-interval effect, setting iii), eligibility effect = difference."""
    m = df_AD.merge(df_AB, on="seed", suffixes=("_AD", "_AB"))
    contrasts = {
        "Observed tail: full reform (A' - D)": {k: m[k + "_AD"] for k in METRICS},
        "Observed tail: waiting-interval removal, criteria retained (A' - B')": {k: m[k + "_AB"] for k in METRICS},
        "Observed tail: eligibility-criteria removal (B' - D) = (A' - D) - (A' - B')": {k: m[k + "_AD"] - m[k + "_AB"] for k in METRICS},
    }
    rows = []
    for name, series in contrasts.items():
        out = {"analysis": "A9_observed_ineligible_wait", "setting": name, "n_seeds": int(len(m))}
        for k, sr in series.items():
            mu, lo, hi, n = _ci(sr.tolist())
            out[k + "_mean"], out[k + "_lo"], out[k + "_hi"] = mu, lo, hi
        dcost = series["dCost_CEA"]; dq = series["dQALY_decomposition"]
        out["dominant_seeds"] = f"{int(((dcost < 0) & (dq > 0)).sum())}/{len(m)}"
        out["NMB_positive_seeds"] = f"{int((series['NMB'] > 0).sum())}/{len(m)}"
        out["breakeven_year_median"] = float("nan")
        rows.append(out)
    all_rows.extend(rows)
    full = rows[0]
    share = [{"contrast": r["setting"],
              "share of full-reform infections averted (%)": round(100 * r["infections_averted_mean"] / full["infections_averted_mean"], 1),
              "share of full-reform NMB (%)": round(100 * r["NMB_mean"] / full["NMB_mean"], 1)} for r in rows[1:]]
    return ("### A9 decomposition under the observed ineligible wait (per-seed paired differences)\n\n"
            + _md_table(pretty_table(rows)) + "\n\n" + _md_table(pd.DataFrame(share)) + "\n")


def cost_of_inaction(setting_frames, base, all_rows):
    rows = []
    for st, df in setting_frames.items():
        m = df.merge(base, on="seed", suffixes=("", "_base"))
        out = {"analysis": "A7_delayed_start", "setting": f"cost of inaction: {st} versus 2024", "n_seeds": int(len(m))}
        for k in METRICS:
            mu, lo, hi, n = _ci((m[k + "_base"] - m[k]).tolist())
            out[k + "_mean"], out[k + "_lo"], out[k + "_hi"] = mu, lo, hi
        out["dominant_seeds"] = ""; out["NMB_positive_seeds"] = ""; out["breakeven_year_median"] = float("nan")
        rows.append(out)
    all_rows.extend(rows)
    return ("### A7 cost of inaction (2024 implementation minus delayed implementation, per seed)\n\n"
            + _md_table(pretty_table(rows)) + "\n")


# -----------------------------------------------------------------------------
# Part 2: discount-rate and 20-year horizon sensitivity from the published run
# -----------------------------------------------------------------------------
def _yearly_from_cea(cea, key):
    bc = cea["base_case_results"][key]
    return (np.array(bc["yearly_qaly"], float), np.array(bc["yearly_cost"], float),
            np.array(bc["yearly_qaly_discounted"], float), np.array(bc["yearly_cost_discounted"], float))


def _cum_infections_deaths(seed_dir, prefix, years):
    p = os.path.join(seed_dir, f"{prefix}_yearly_summary.csv")
    if not os.path.exists(p):
        return float("nan"), float("nan")
    df = pd.read_csv(p)
    df = df[df["Year"] <= years]
    return float(df["NewlyInfected"].sum()), float(df["Deaths"].sum())


def part2():
    rows_dr, rows_h = [], []
    seed_dirs = sorted(glob.glob(os.path.join(BASE_MULTISEED_DIR, "seed_*")))
    for inter, comp in COMPARISONS_PART2:
        per = {r: {"dQALY": [], "dCost": [], "NMB": []} for r in DISCOUNT_RATES}
        h20 = {"dQALY": [], "dCost": [], "NMB": [], "inf": [], "deaths": []}
        for d in seed_dirs:
            cea = _read_cea(d)
            if cea is None or inter not in cea["base_case_results"] or comp not in cea["base_case_results"]:
                continue
            qi, ci, qdi, cdi = _yearly_from_cea(cea, inter)
            qc, cc, qdc, cdc = _yearly_from_cea(cea, comp)
            t = np.arange(1, len(qi) + 1)
            for r in DISCOUNT_RATES:
                f = (1.0 + r) ** (-(t - 1))
                dq = float(np.sum((qi - qc) * f)); dc = float(np.sum((ci - cc) * f))
                per[r]["dQALY"].append(dq); per[r]["dCost"].append(dc); per[r]["NMB"].append(dq * WTP - dc)
            n = SHORT_HORIZON_YEARS
            dq20 = float(np.sum(qdi[:n] - qdc[:n])); dc20 = float(np.sum(cdi[:n] - cdc[:n]))
            h20["dQALY"].append(dq20); h20["dCost"].append(dc20); h20["NMB"].append(dq20 * WTP - dc20)
            pi = P.SCEN_FILE[inter]; pc = P.SCEN_FILE[comp]
            ii, di = _cum_infections_deaths(d, pi, n); ic, dcc = _cum_infections_deaths(d, pc, n)
            h20["inf"].append(ic - ii); h20["deaths"].append(dcc - di)
        for r in DISCOUNT_RATES:
            mq, lq, hq, n = _ci(per[r]["dQALY"]); mc, lc, hc, _ = _ci(per[r]["dCost"]); mn, ln, hn, _ = _ci(per[r]["NMB"])
            dom = int(np.sum((np.array(per[r]["dCost"]) < 0) & (np.array(per[r]["dQALY"]) > 0)))
            rows_dr.append({"comparison": f"{inter} vs {comp}", "discount_rate": r, "n_seeds": n,
                            "Incremental QALYs, direct (95% CI)": _fmt_ci(mq, lq, hq),
                            "Cost saving, billion yen (95% CI)": _fmt_ci(-mc, -hc, -lc, 1e9, 1),
                            "NMB, billion yen (95% CI)": _fmt_ci(mn, ln, hn, 1e9, 1),
                            "dominant_seeds": f"{dom}/{n}",
                            "NMB_positive_seeds": f"{int(np.sum(np.array(per[r]['NMB']) > 0))}/{n}"})
        mq, lq, hq, n = _ci(h20["dQALY"]); mc, lc, hc, _ = _ci(h20["dCost"]); mn, ln, hn, _ = _ci(h20["NMB"])
        mi, li, hi, _ = _ci(h20["inf"]); md_, ld, hd, _ = _ci(h20["deaths"])
        dom = int(np.sum((np.array(h20["dCost"]) < 0) & (np.array(h20["dQALY"]) > 0)))
        rows_h.append({"comparison": f"{inter} vs {comp}", "horizon_years": SHORT_HORIZON_YEARS, "n_seeds": n,
                       "Infections averted (95% CI)": _fmt_ci(mi, li, hi),
                       "Deaths averted (95% CI)": _fmt_ci(md_, ld, hd),
                       "Incremental QALYs, direct (95% CI)": _fmt_ci(mq, lq, hq),
                       "Cost saving, billion yen (95% CI)": _fmt_ci(-mc, -hc, -lc, 1e9, 1),
                       "NMB, billion yen (95% CI)": _fmt_ci(mn, ln, hn, 1e9, 1),
                       "dominant_seeds": f"{dom}/{n}"})
    df_dr = pd.DataFrame(rows_dr); df_h = pd.DataFrame(rows_h)
    df_dr.to_csv(os.path.join(OUTPUT_DIR, "discount_rate_sensitivity_direct.csv"), index=False)
    df_h.to_csv(os.path.join(OUTPUT_DIR, "horizon_20y_direct.csv"), index=False)
    with open(os.path.join(OUTPUT_DIR, "base_postprocessing_summary.md"), "w", encoding="utf-8") as f:
        f.write("# Post-processing of the published 200-seed run (whole-population direct estimator)\n\n"
                "## Discount-rate sensitivity (same rate for costs and QALYs; year 1 undiscounted)\n\n"
                + _md_table(df_dr) + "\n\n## 20-year horizon (discounted at 2% per year)\n\n" + _md_table(df_h) + "\n")
    print(df_dr.to_string(index=False)); print(df_h.to_string(index=False))


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    _new_seed_budget[0] = MAX_NEW_SEEDS_PER_RUN
    if RUN_PART1:
        rows = part1()
        if rows:
            print(pd.DataFrame(rows)[["analysis", "setting", "n_seeds", "dominant_seeds"]].to_string(index=False))
    if RUN_PART2:
        part2()
    if PACKAGE_COPY_DIR:
        os.makedirs(PACKAGE_COPY_DIR, exist_ok=True)
        for fn in os.listdir(OUTPUT_DIR):
            src = os.path.join(OUTPUT_DIR, fn)
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join(PACKAGE_COPY_DIR, fn))
        print(f"\nSaved to {OUTPUT_DIR} and copied to {PACKAGE_COPY_DIR}")
    else:
        print(f"\nSaved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()

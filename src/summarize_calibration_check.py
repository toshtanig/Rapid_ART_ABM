#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Summary of the calibration check for the timing convention of ART starts.

Reads the burn-in summaries written by run_calibration_check.py

    outputs/revision_R2_runs/C0_calibration_check/<setting>/seed_<n>/burnin_summary.csv

and the burn-in summaries of versions 1.0.0 and 1.1.0 for the same seeds

    outputs/core_multiseed/seed_<n>/burnin_summary.csv

and reports, for each setting,
  1. the simulated state in calendar year 2023 (the calibration indicators of Supplementary Table S2), with the
     paired difference from the earlier versions on the same seeds and the documented acceptance ranges;
  2. the goodness of fit of the seed-median HIV plus AIDS notification series for 2004 to 2023 against the
     observed national series (the statistics of Supplementary Table S16);
  3. the year-by-year seed-median notification series.

Run as-is in VS Code (no arguments); all choices are in the CONFIG block. It can be run while seeds are still
being completed: only completed seeds are read. Output: C0_calibration_check/summary/ (Markdown and CSV).
"""

import os
import sys
import glob
import datetime
import numpy as np
import pandas as pd

# =============================================================================
# CONFIG
# =============================================================================
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import BASE_MULTISEED, REVISION_RUNS

CHECK_ROOT = os.path.join(REVISION_RUNS, "C0_calibration_check")
PUBLISHED_ROOT = BASE_MULTISEED      # versions 1.0.0 / 1.1.0, 200 seeds (read only)
OUT_DIR = os.path.join(CHECK_ROOT, "summary")
SETTINGS_TO_SUMMARISE = None     # None = every setting folder found under CHECK_ROOT, or a list of names

BURNIN_START_YEAR = 2004         # burn_in_year 1 = calendar year 2004, burn_in_year 20 = 2023
# Observed national notifications (AIDS Surveillance Committee annual reports), as in multiseed_figures.py
OBS_HIV = {2004: 780, 2005: 832, 2006: 952, 2007: 1082, 2008: 1126, 2009: 1021, 2010: 1075, 2011: 1056,
           2012: 1002, 2013: 1106, 2014: 1091, 2015: 1006, 2016: 1011, 2017: 976, 2018: 940, 2019: 903,
           2020: 750, 2021: 742, 2022: 632, 2023: 669}
OBS_AIDS = {2004: 385, 2005: 367, 2006: 406, 2007: 418, 2008: 431, 2009: 431, 2010: 469, 2011: 473,
            2012: 447, 2013: 484, 2014: 455, 2015: 428, 2016: 437, 2017: 413, 2018: 377, 2019: 333,
            2020: 345, 2021: 315, 2022: 252, 2023: 291}
FIT_PERIODS = [("2004 to 2023 (all)", 2004, 2023),
               ("2004 to 2013 (rise to the peak)", 2004, 2013),
               ("2014 to 2023 (decline)", 2014, 2023)]
# Calibration targets and acceptance ranges at 2023 (Supplementary Material section 2 and Table S2).
# (indicator, target text, lower bound or None, upper bound or None)
TARGETS = [
    ("Annual HIV and AIDS notifications", "about 946 (approximately 950 in 2022 to 2024)", None, None),
    ("People currently on ART", "about 29,680 (acceptance range 28,500 to 30,500)", 28500, 30500),
    ("Undiagnosed people with HIV", "about 3,500 (range 3,200 to 3,800)", 3200, 3800),
    ("Undiagnosed fraction (%)", "10 to 20 (national estimate)", 10.0, 20.0),
    ("People with HIV", "33,000 to 36,000 (soft target)", 33000, 36000),
    ("Viral suppression on ART (%)", "above 99", 99.0, None),
    ("New infections per year", "no target", None, None),
]
# =============================================================================


def read_burnin(path):
    df = pd.read_csv(path)
    df["calendar_year"] = BURNIN_START_YEAR - 1 + df["burn_in_year"].astype(int)
    return df.set_index("calendar_year")


def indicators_2023(df):
    r = df.loc[2023]
    hiv = float(r["hiv_positive"])
    return {
        "Annual HIV and AIDS notifications": float(r["new_diagnoses_year"]),
        "People currently on ART": float(r["on_art"]),
        "Undiagnosed people with HIV": float(r["undiagnosed_hiv"]),
        "Undiagnosed fraction (%)": 100.0 * float(r["undiagnosed_hiv"]) / hiv if hiv > 0 else np.nan,
        "People with HIV": hiv,
        "Viral suppression on ART (%)": 100.0 * float(r["vl_suppressed"]) / float(r["on_art"]) if r["on_art"] > 0 else np.nan,
        "New infections per year": float(r["new_infections_year"]),
    }


def fit_statistics(median_series):
    """RMSE, MAPE, and mean bias of a {year: modelled} series against the observed HIV plus AIDS notifications."""
    rows = []
    for name, y0, y1 in FIT_PERIODS:
        years = [y for y in range(y0, y1 + 1) if y in median_series.index]
        obs = np.array([OBS_HIV[y] + OBS_AIDS[y] for y in years], float)
        mod = np.array([median_series.loc[y] for y in years], float)
        err = mod - obs
        rows.append({"period": name,
                     "RMSE (per year)": float(np.sqrt(np.mean(err ** 2))),
                     "MAPE (%)": float(100.0 * np.mean(np.abs(err) / obs)),
                     "Mean bias, modelled minus observed (per year)": float(np.mean(err)),
                     "Observed mean (per year)": float(np.mean(obs))})
    return pd.DataFrame(rows)


def mean_sd(x):
    x = np.asarray(x, float)
    return x.mean(), (x.std(ddof=1) if len(x) > 1 else np.nan)


def fmt(v, dec=0):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ""
    return f"{v:,.{dec}f}"


def summarise_setting(name, lines):
    setting_dir = os.path.join(CHECK_ROOT, name)
    new, old = {}, {}
    for path in sorted(glob.glob(os.path.join(setting_dir, "seed_*", "burnin_summary.csv"))):
        seed = int(os.path.basename(os.path.dirname(path)).split("_")[1])
        if not os.path.exists(os.path.join(os.path.dirname(path), "calibration_check_config.json")):
            continue                                   # seed still running
        new[seed] = read_burnin(path)
        pub = os.path.join(PUBLISHED_ROOT, f"seed_{seed}", "burnin_summary.csv")
        if os.path.exists(pub):
            old[seed] = read_burnin(pub)
    seeds = sorted(new)
    paired = [s for s in seeds if s in old]
    lines.append(f"## Setting: {name}")
    lines.append("")
    lines.append(f"Seeds completed: {len(seeds)} ({', '.join(str(s) for s in seeds)}); "
                 f"paired with the earlier versions: {len(paired)}.")
    lines.append("")
    if not seeds:
        return None

    # ---- 1. state in 2023 -------------------------------------------------------------------------------------
    ind_new = pd.DataFrame({s: indicators_2023(new[s]) for s in seeds}).T
    ind_old = pd.DataFrame({s: indicators_2023(old[s]) for s in paired}).T if paired else None
    lines.append("### 1. Simulated state in 2023 (mean across seeds; SD of seed values)")
    lines.append("")
    lines.append("| Indicator | Versions 1.0.0 / 1.1.0, same seeds | This setting | Paired difference (%) | Target | Within range |")
    lines.append("|---|---|---|---|---|---|")
    table_rows = []
    for ind, target, lo, hi in TARGETS:
        dec = 1 if "%" in ind else 0
        m_new, sd_new = mean_sd(ind_new[ind])
        if ind_old is not None:
            m_old, sd_old = mean_sd(ind_old[ind])
            d = (ind_new.loc[paired, ind] - ind_old.loc[paired, ind]).mean()
            rel = 100.0 * d / m_old if m_old else np.nan
            old_txt = f"{fmt(m_old, dec)} ({fmt(sd_old, dec)})"
            diff_txt = f"{d:+,.{dec}f} ({rel:+.1f}%)"
        else:
            m_old = np.nan; old_txt = ""; diff_txt = ""
        within = ""
        if lo is not None or hi is not None:
            ok = (lo is None or m_new >= lo) and (hi is None or m_new <= hi)
            within = "yes" if ok else "NO"
            if ind_old is not None:
                ok_old = (lo is None or m_old >= lo) and (hi is None or m_old <= hi)
                within += " (earlier versions: " + ("yes" if ok_old else "NO") + ")"
        lines.append(f"| {ind} | {old_txt} | {fmt(m_new, dec)} ({fmt(sd_new, dec)}) | {diff_txt} | {target} | {within} |")
        table_rows.append({"setting": name, "indicator": ind, "published_same_seeds_mean": m_old,
                           "this_setting_mean": m_new, "this_setting_sd": sd_new, "target": target,
                           "within_range": within, "n_seeds": len(seeds), "n_paired": len(paired)})
    lines.append("")

    # ---- 2. fit of the notification series ------------------------------------------------------------------
    med_new = pd.DataFrame({s: new[s]["new_diagnoses_year"] for s in seeds}).median(axis=1)
    fits = [("This setting", fit_statistics(med_new))]
    if paired:
        med_old = pd.DataFrame({s: old[s]["new_diagnoses_year"] for s in paired}).median(axis=1)
        fits.insert(0, ("Versions 1.0.0 / 1.1.0, same seeds", fit_statistics(med_old)))
    all_pub = sorted(glob.glob(os.path.join(PUBLISHED_ROOT, "seed_*", "burnin_summary.csv")))
    if all_pub:
        med_all = pd.DataFrame({i: read_burnin(p)["new_diagnoses_year"] for i, p in enumerate(all_pub)}).median(axis=1)
        fits.insert(0, (f"Versions 1.0.0 / 1.1.0, all {len(all_pub)} seeds (Table S16)", fit_statistics(med_all)))
    lines.append("### 2. Fit of the seed-median HIV plus AIDS notifications to the observed series")
    lines.append("")
    lines.append("| Series | Period | RMSE (per year) | MAPE (%) | Mean bias (per year) | Observed mean (per year) |")
    lines.append("|---|---|---|---|---|---|")
    fit_rows = []
    for label, df in fits:
        for _, r in df.iterrows():
            lines.append(f"| {label} | {r['period']} | {r['RMSE (per year)']:.0f} | {r['MAPE (%)']:.1f} | "
                         f"{r['Mean bias, modelled minus observed (per year)']:+.0f} | {r['Observed mean (per year)']:.0f} |")
            fit_rows.append({"setting": name, "series": label, **r.to_dict()})
    lines.append("")

    # ---- 3. year-by-year series -----------------------------------------------------------------------------
    lines.append("### 3. Seed-median series by calendar year")
    lines.append("")
    lines.append("| Year | Observed notifications | Notifications, earlier versions | Notifications, this setting | Ratio | "
                 "On ART, earlier versions | On ART, this setting | New infections, earlier versions | New infections, this setting |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    year_rows = []
    cols = ["new_diagnoses_year", "on_art", "new_infections_year"]
    mn = {c: pd.DataFrame({s: new[s][c] for s in seeds}).median(axis=1) for c in cols}
    mo = {c: (pd.DataFrame({s: old[s][c] for s in paired}).median(axis=1) if paired else None) for c in cols}
    for y in range(BURNIN_START_YEAR, 2024):
        if y not in mn["new_diagnoses_year"].index:
            continue
        obs = OBS_HIV[y] + OBS_AIDS[y]
        o = {c: (float(mo[c].loc[y]) if mo[c] is not None else np.nan) for c in cols}
        n = {c: float(mn[c].loc[y]) for c in cols}
        ratio = n["new_diagnoses_year"] / o["new_diagnoses_year"] if o["new_diagnoses_year"] else np.nan
        lines.append(f"| {y} | {obs:,} | {fmt(o['new_diagnoses_year'])} | {fmt(n['new_diagnoses_year'])} | {ratio:.3f} | "
                     f"{fmt(o['on_art'])} | {fmt(n['on_art'])} | {fmt(o['new_infections_year'])} | {fmt(n['new_infections_year'])} |")
        year_rows.append({"setting": name, "calendar_year": y, "observed_notifications": obs,
                          **{f"published_{c}": o[c] for c in cols}, **{f"this_setting_{c}": n[c] for c in cols}})
    lines.append("")
    return pd.DataFrame(table_rows), pd.DataFrame(fit_rows), pd.DataFrame(year_rows)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    if SETTINGS_TO_SUMMARISE is None:
        names = sorted(d for d in os.listdir(CHECK_ROOT)
                       if os.path.isdir(os.path.join(CHECK_ROOT, d)) and not d.startswith("_") and d != "summary")
    else:
        names = list(SETTINGS_TO_SUMMARISE)
    lines = ["# Calibration check: timing convention of ART starts (burn-in, 2004 to 2023)", "",
             f"Generated {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}. Earlier versions = outputs/core_multiseed "
             "(first viral-load update of a scheduled start in the start cycle). Paired difference = this setting minus "
             "the earlier versions on the same seeds.", ""]
    parts = [[], [], []]
    for name in names:
        out = summarise_setting(name, lines)
        if out is None:
            continue
        for k in range(3):
            parts[k].append(out[k])
    md = os.path.join(OUT_DIR, "calibration_check_summary.md")
    with open(md, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    for k, fn in enumerate(["calibration_check_state_2023.csv", "calibration_check_fit_statistics.csv",
                            "calibration_check_series_by_year.csv"]):
        if parts[k]:
            pd.concat(parts[k], ignore_index=True).to_csv(os.path.join(OUT_DIR, fn), index=False)
    print("\n".join(lines))
    print("saved", md)


if __name__ == "__main__":
    main()

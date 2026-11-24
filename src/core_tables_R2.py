#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tables of the Supplementary Material that are computed from the per-seed cache of summarize_R2.py.

Nothing is read from the seed folders here: run summarize_R2.py (step "cache") first. Outputs (CSV and one
Markdown overview) go to post/tables/ of the target:

  calibration            modelled state in 2023, notifications by year, fit to the observed notifications,
                         first analytic year (2024), for the base analysis and for every setting whose burn-in
                         differs from the base
  risk_group             infections and deaths among people with HIV averted, by risk group
  transmission_source    domestic infections by the state of the transmitting person, by strategy, and the
                         decomposition of the additional infections under delayed ART
  monte_carlo            standard deviation across seeds, relative Monte Carlo standard error and the number of
                         seeds for a 95% CI half-width of 10% and of 5% of the mean
  absolute_by_arm        40-year totals of each strategy (infections, deaths, ART initiations, diagnoses)
  delays                 diagnosis-to-ART and eligibility-to-ART intervals (pooled over seeds)
  ineligible_at_diagnosis  share of ART starters who did not meet the eligibility criteria at diagnosis, and
                         their time from diagnosis to ART
  breakeven_curve        cumulative discounted cost difference by year and the first five years undiscounted
  mortality_validation   person-years, mortality rate and observed-to-expected ratio by care state and CD4
  survival_weight        one-way sensitivity to the counterfactual survival weight (per-seed costs)
  subgroup               size of the untreated subgroup with high CD4 count and low viral load
  figure_data            yearly series (median and 95% prediction interval across seeds)
  convention             (TARGET = "R2" only) the base analysis against the published outputs on the same
                         seeds: the results under the earlier timing convention for ART starts

Run as-is in VS Code; all settings are in the CONFIG block, no command-line arguments.
"""

import os
import sys
import json
import math

import numpy as np
import pandas as pd

# =============================================================================
# CONFIG
# =============================================================================
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import summarize_R2 as SR

TARGET = "R2"                    # "R2" or "published"
SURVIVAL_WEIGHT_GRID = [0.45, 0.55, 0.60, 0.648, 0.70, 0.75, 0.80, 0.85]
# Observed HIV plus AIDS notifications (AIDS Surveillance Committee annual report 2024)
OBSERVED = {2004: 1165, 2005: 1199, 2006: 1358, 2007: 1500, 2008: 1557, 2009: 1452, 2010: 1544, 2011: 1529,
            2012: 1449, 2013: 1590, 2014: 1546, 2015: 1434, 2016: 1448, 2017: 1389, 2018: 1317, 2019: 1236,
            2020: 1095, 2021: 1057, 2022: 884, 2023: 960}
OBSERVED_2024 = 994
FIRST_BURNIN_YEAR = 2004
PERIODS = [("2004 to 2023", 2004, 2023), ("2004 to 2013", 2004, 2013), ("2014 to 2023", 2014, 2023)]
# =============================================================================

DAYS = 365.25 / 12.0
SOURCES = ["undiagnosed", "diagnosed_pre_art", "on_art_unsuppressed", "disengaged"]


def out_dir():
    return os.path.join(os.path.dirname(SR.OUT_ROOT[TARGET]), "tables")


def ci(vals):
    m, lo, hi, n = SR._ci(list(vals))
    return m, lo, hi


def fmt(m, lo, hi, nd=0, scale=1.0):
    f = f"{{:,.{nd}f}}"
    return f"{f.format(m / scale)} ({f.format(lo / scale)} to {f.format(hi / scale)})"


def md(df):
    return SR.md_table(df)


# -----------------------------------------------------------------------------
def calibration(cache, label):
    """State in 2023, notifications by year and fit, from the burn-in summaries of the seeds."""
    recs = [r for r in cache.values() if "burnin" in r]
    if not recs:
        return None
    b = [r["burnin"] for r in recs]
    n = len(b)
    def col(k):
        return np.array([x[k] for x in b], float)
    hiv, art, und, sup, notif, inf = col("hiv"), col("on_art"), col("undiag"), col("supp"), col("notif"), col("inf")
    state = []
    for name, v in [("Annual HIV and AIDS notifications", notif), ("People currently on ART", art),
                    ("People with HIV", hiv), ("Undiagnosed people with HIV", und),
                    ("Undiagnosed fraction (%)", 100 * und / hiv), ("Viral suppression on ART (%)", 100 * sup / art),
                    ("New infections per year", inf)]:
        sd = v.std(ddof=1) if n > 1 else 0.0
        state.append({"setting": label, "indicator": name, "n_seeds": n, "mean": v.mean(), "sd": sd,
                      "range95_lo": v.mean() - 1.96 * sd, "range95_hi": v.mean() + 1.96 * sd,
                      "ci95_halfwidth": 1.96 * sd / math.sqrt(n)})
    notif_y = np.array([x["notif_y"] for x in b], float)            # seeds x years
    years = [FIRST_BURNIN_YEAR + i for i in range(notif_y.shape[1])]
    med, mean = np.median(notif_y, axis=0), notif_y.mean(axis=0)
    lo, hi = np.percentile(notif_y, 2.5, axis=0), np.percentile(notif_y, 97.5, axis=0)
    series = pd.DataFrame({"setting": label, "year": years, "observed": [OBSERVED.get(y) for y in years],
                           "modelled_mean": mean, "modelled_median": med, "pi95_lo": lo, "pi95_hi": hi})
    fit = []
    for name, y0, y1 in PERIODS:
        idx = [i for i, y in enumerate(years) if y0 <= y <= y1 and y in OBSERVED]
        obs = np.array([OBSERVED[years[i]] for i in idx], float)
        mod = med[idx]
        fit.append({"setting": label, "period": name, "n_seeds": n,
                    "RMSE": float(np.sqrt(np.mean((mod - obs) ** 2))),
                    "MAPE_percent": float(100 * np.mean(np.abs(mod - obs) / obs)),
                    "mean_bias": float(np.mean(mod - obs)), "observed_mean": float(obs.mean()),
                    "modelled_mean_of_medians": float(mod.mean())})
    return pd.DataFrame(state), series, pd.DataFrame(fit)


def first_year(cache, label):
    rows = []
    for arm in SR.ARMS:
        v = np.array([r["arms"][arm]["dx_y"][0] for r in cache.values() if arm in r["arms"] and r["arms"][arm].get("dx_y")], float)
        if len(v):
            rows.append({"setting": label, "arm": arm, "n_seeds": len(v), "notifications_2024_mean": v.mean(),
                         "pi95_lo": np.percentile(v, 2.5), "pi95_hi": np.percentile(v, 97.5), "observed_2024": OBSERVED_2024})
    return pd.DataFrame(rows)


def het_share(cache, label):
    rows = []
    for arm in SR.ARMS:
        tot = het = 0.0
        k = 0
        for r in cache.values():
            g = r["arms"].get(arm, {}).get("by_group")
            if g:
                tot += sum(v[0] for v in g.values())
                het += g.get("Heterosexual", [0.0])[0]
                k += 1
        if k:
            rows.append({"setting": label, "arm": arm, "n_seeds": k, "heterosexual_share_of_infections_percent": 100 * het / tot})
    return pd.DataFrame(rows)


def risk_group(cache):
    rows = []
    for comp in ("1m", "2m"):
        recs = [r for r in cache.values() if comp in r["arms"]]
        if not recs:
            continue
        groups = sorted(recs[0]["arms"]["0m"]["by_group"])
        for g in groups:
            di = np.array([r["arms"][comp]["by_group"][g][0] - r["arms"]["0m"]["by_group"][g][0] for r in recs])
            dd = np.array([r["arms"][comp]["by_group"][g][1] - r["arms"]["0m"]["by_group"][g][1] for r in recs])
            ci_i, ci_d = ci(di), ci(dd)
            base_i = np.mean([r["arms"][comp]["by_group"][g][0] for r in recs])
            rows.append({"comparison": f"0m vs {comp}", "risk_group": g, "n_seeds": len(recs),
                         "infections_averted": fmt(*ci_i), "deaths_hiv_averted": fmt(*ci_d),
                         "infections_averted_mean": ci_i[0], "deaths_hiv_averted_mean": ci_d[0],
                         "comparator_infections_mean": base_i,
                         "relative_reduction_percent": 100 * ci_i[0] / base_i if base_i else float("nan")})
    return pd.DataFrame(rows)


def transmission_source(cache):
    recs = [r for r in cache.values() if all("source" in r["arms"].get(a, {}) for a in ("0m", "1m"))]
    if not recs:
        return None, None
    arms = [a for a in SR.ARMS if all(a in r["arms"] and "source" in r["arms"][a] for r in recs)]
    allc = set().union(*[set(r["arms"]["0m"]["source"]) for r in recs])
    external = sorted(c for c in allc if c.startswith("exogenous"))       # inflow of people with HIV: not domestic transmission
    cats = sorted((c for c in allc if c not in external), key=lambda c: (SOURCES.index(c) if c in SOURCES else 99, c))
    t7, t8 = [], []
    mean = {a: {c: float(np.mean([r["arms"][a]["source"].get(c, 0.0) for r in recs])) for c in cats + external} for a in arms}
    for a in arms:
        tot = sum(mean[a][c] for c in cats)
        row = {"strategy": SR.LABEL[a], "n_seeds": len(recs), "domestic_total": tot}
        for c in cats:
            row[c] = mean[a][c]
            row[c + "_percent"] = 100 * mean[a][c] / tot
        for c in external:
            row[c + " (not in the domestic total)"] = mean[a][c]
        t7.append(row)
    for a in arms:
        if a == "0m":
            continue
        add = {c: mean[a][c] - mean["0m"][c] for c in cats}
        tot = sum(add.values())
        row = {"comparison": f"{SR.LABEL[a]} delay minus immediate ART", "n_seeds": len(recs), "additional_total": tot}
        for c in cats:
            row[c] = add[c]
            row[c + "_percent"] = 100 * add[c] / tot
        t8.append(row)
    return pd.DataFrame(t7), pd.DataFrame(t8)


def _seeds_needed(sd, mean, share):
    """Number of seeds for which the half-width of the 95% CI is the given share of the mean (rounded up)."""
    if not (math.isfinite(sd) and math.isfinite(mean)) or mean == 0:
        return float("nan")
    return math.ceil((1.96 * sd / (share * abs(mean))) ** 2)


def monte_carlo(cache):
    rows = []
    for comp in ("1m", "2m"):
        df = SR.frame(cache, "0m", cache, comp, same_folder_pair=f"0m vs {comp}")
        if not len(df):
            continue
        n = len(df)
        for k in ["infections_averted", "deaths_averted", "deaths_averted_all", "dQALY_decomposition", "dQALY_direct",
                  "dCost_CEA", "NMB", "NMB_direct"]:
            v = df[k].to_numpy(float)
            sd = v.std(ddof=1)
            rows.append({"comparison": f"0m vs {comp}", "outcome": k, "n_seeds": n, "mean": v.mean(), "sd": sd,
                         "se": sd / math.sqrt(n), "relative_mcse_percent": 100 * (sd / math.sqrt(n)) / abs(v.mean()),
                         "seeds_for_ci95_halfwidth_10pct": _seeds_needed(sd, v.mean(), 0.10),
                         "seeds_for_ci95_halfwidth_5pct": _seeds_needed(sd, v.mean(), 0.05),
                         "share_positive_percent": 100 * float(np.mean(v > 0))})
        q = np.array([r["arms"]["0m"]["qaly_d"] for r in cache.values() if "0m" in r["arms"]], float)
        rows.append({"comparison": f"0m vs {comp}", "outcome": "total discounted QALYs of the immediate-ART arm", "n_seeds": len(q),
                     "mean": q.mean(), "sd": q.std(ddof=1), "se": q.std(ddof=1) / math.sqrt(len(q)),
                     "relative_mcse_percent": float("nan"), "seeds_for_ci95_halfwidth_10pct": float("nan"),
                     "seeds_for_ci95_halfwidth_5pct": float("nan"), "share_positive_percent": float("nan")})
    return pd.DataFrame(rows)


def absolute_by_arm(cache):
    """Means per seed of the 40-year totals of each strategy."""
    rows = []
    for arm in SR.ARMS:
        recs = [r["arms"][arm] for r in cache.values() if arm in r["arms"]]
        if not recs:
            continue

        def mean(f):
            v = [f(a) for a in recs]
            v = [x for x in v if x is not None]
            return float(np.mean(v)) if v else float("nan")

        rows.append({"slot": SR.LABEL[arm], "n_seeds": len(recs),
                     "cumulative_new_infections": mean(lambda a: sum(v[0] for v in a["by_group"].values())),
                     "deaths_among_people_with_hiv": mean(lambda a: a["deaths_hiv"]),
                     "art_initiations": mean(lambda a: sum(v[2] for v in a["by_group"].values())),
                     "new_diagnoses": mean(lambda a: sum(a["dx_y"]) if a.get("dx_y") else None),
                     "deaths_whole_population": mean(lambda a: sum(a["deaths_all_y"]) if a.get("deaths_all_y") else None),
                     "total_discounted_qalys": mean(lambda a: a["qaly_d"]),
                     "total_discounted_cost_bn": mean(lambda a: a["cost_d"]) / 1e9})
    return pd.DataFrame(rows)


def _pooled_quantiles(hists):
    """Median, quartiles and mean (months) of the pooled distribution from per-seed histograms."""
    tot = {}
    for h in hists:
        for m, c in h:
            tot[m] = tot.get(m, 0.0) + c
    if not tot:
        return None
    months = np.array(sorted(tot))
    counts = np.array([tot[m] for m in months])
    cum = np.cumsum(counts) / counts.sum()
    q = lambda p: float(months[np.searchsorted(cum, p, side="left")])
    return {"n": float(counts.sum()), "median": q(0.5), "p25": q(0.25), "p75": q(0.75),
            "mean": float(np.sum(months * counts) / counts.sum())}


def delays(caches):
    rows = []
    for (an, st), cache in caches.items():
        for arm in SR.ARMS:
            for kind, key in (("diagnosis to ART", "hist_diag"), ("eligibility to ART", "hist_elig")):
                hs = [r["arms"][arm][key] for r in cache.values() if arm in r["arms"] and key in r["arms"][arm]]
                if not hs:
                    continue
                s = _pooled_quantiles(hs)
                rows.append({"analysis": an, "setting": st, "slot": SR.LABEL[arm], "interval": kind, "n_seeds": len(hs),
                             "people_pooled": s["n"],
                             "median_days": s["median"] * DAYS, "p25_days": s["p25"] * DAYS, "p75_days": s["p75"] * DAYS,
                             "mean_days": s["mean"] * DAYS, "mean_months": s["mean"],
                             "text": f"median {s['median'] * DAYS:.0f} (IQR {s['p25'] * DAYS:.0f} to {s['p75'] * DAYS:.0f}) days; mean {s['mean'] * DAYS:.0f} days"})
    return pd.DataFrame(rows)


def ineligible_wait(caches):
    """People who did not meet the eligibility criteria at diagnosis, among those who started ART: share and
    time from diagnosis to ART (pooled over seeds)."""
    rows = []
    for (an, st), cache in caches.items():
        for arm in SR.ARMS:
            recs = [r["arms"][arm]["inelig"] for r in cache.values() if "inelig" in r["arms"].get(arm, {})]
            if not recs:
                continue
            n, k = sum(x["n"] for x in recs), sum(x["n_inelig"] for x in recs)
            s = _pooled_quantiles([x["hist"] for x in recs])
            row = {"analysis": an, "setting": st, "slot": SR.LABEL[arm], "n_seeds": len(recs),
                   "art_starters_per_seed": n / len(recs), "ineligible_at_diagnosis_per_seed": k / len(recs),
                   "ineligible_at_diagnosis_percent": 100 * k / n if n else float("nan")}
            if s is not None:
                row.update({"diagnosis_to_art_mean_months": s["mean"], "diagnosis_to_art_median_months": s["median"],
                            "diagnosis_to_art_p25_months": s["p25"], "diagnosis_to_art_p75_months": s["p75"]})
            e = _pooled_quantiles([x["hist_until_eligible"] for x in recs if "hist_until_eligible" in x])
            if e is not None:
                row.update({"diagnosis_to_eligibility_mean_months": e["mean"], "diagnosis_to_eligibility_median_months": e["median"],
                            "diagnosis_to_eligibility_p25_months": e["p25"], "diagnosis_to_eligibility_p75_months": e["p75"]})
            rows.append(row)
    return pd.DataFrame(rows)


def breakeven_curves(caches):
    """For every setting: cumulative discounted cost difference by year between the two slots of the setting
    (budget-impact cost convention; mean and 95% CI across seeds) and the undiscounted annual difference in the
    first five years."""
    curve, short = [], []
    for (an, st), cache in caches.items():
        c, s = _breakeven_of_setting(cache)
        for rows, out in ((c, curve), (s, short)):
            for r in rows:
                out.append({"analysis": an, "setting": st, **r})
    return pd.DataFrame(curve), pd.DataFrame(short)


def _breakeven_of_setting(cache):
    curve, short = [], []
    for comp in ("1m", "2m"):
        key = f"0m vs {comp}"
        cum = [r["pairs"][key]["be_cum"] for r in cache.values() if r["pairs"].get(key, {}).get("be_cum")]
        ann = [r["pairs"][key]["be_ann"] for r in cache.values() if r["pairs"].get(key, {}).get("be_ann")]
        if not cum:
            continue
        m = np.array(cum, float)
        se = m.std(axis=0, ddof=1) / math.sqrt(m.shape[0]) if m.shape[0] > 1 else np.zeros(m.shape[1])
        mean = m.mean(axis=0)
        first = next((i + 1 for i, v in enumerate(mean) if v <= 0), None)
        for i in range(m.shape[1]):
            curve.append({"comparison": key, "year": i + 1, "calendar_year": 2024 + i, "n_seeds": m.shape[0],
                          "cumulative_bn_mean": mean[i], "ci95_lo": mean[i] - 1.96 * se[i], "ci95_hi": mean[i] + 1.96 * se[i],
                          "break_even_year_of_the_mean_curve": first})
        if ann:
            a = np.array(ann, float)
            und = a * (1.0 + SR.DISCOUNT_RATE) ** np.arange(a.shape[1])
            for label, v in [(str(i + 1), und[:, i]) for i in range(5)] + [("1 to 5", und[:, :5].sum(axis=1))]:
                mm, lo, hi = ci(v)
                short.append({"comparison": key, "year": label, "n_seeds": a.shape[0],
                              "undiscounted_difference_bn": fmt(mm, lo, hi, nd=2), "mean": mm, "ci95_lo": lo, "ci95_hi": hi})
    return curve, short


def delay_histograms(caches):
    """Pooled counts of the diagnosis-to-ART and eligibility-to-ART intervals (whole months) by setting and slot."""
    rows = []
    for (an, st), cache in caches.items():
        for arm in SR.ARMS:
            for kind, key in (("diagnosis to ART", "hist_diag"), ("eligibility to ART", "hist_elig")):
                tot, n = {}, 0
                for r in cache.values():
                    h = r["arms"].get(arm, {}).get(key)
                    if not h:
                        continue
                    n += 1
                    for m, c in h:
                        tot[m] = tot.get(m, 0.0) + c
                total = sum(tot.values())
                for m in sorted(tot):
                    rows.append({"analysis": an, "setting": st, "slot": SR.LABEL[arm], "interval": kind, "n_seeds": n,
                                 "delay_months": m, "count": tot[m], "proportion": tot[m] / total})
    return pd.DataFrame(rows)


def mortality_validation(cache, arm="0m"):
    recs = [r for r in cache.values() if "mort" in r["arms"].get(arm, {})]
    if not recs:
        return None
    keys = sorted(set().union(*[set(r["arms"][arm]["mort"]) for r in recs]))
    rows = []
    for k in keys:
        d = sum(r["arms"][arm]["mort"].get(k, [0, 0, 0])[0] for r in recs)
        py = sum(r["arms"][arm]["mort"].get(k, [0, 0, 0])[1] for r in recs)
        ex = sum(r["arms"][arm]["mort"].get(k, [0, 0, 0])[2] for r in recs)
        cd4, care = k.split("|")
        rows.append({"cd4_stratum": cd4, "care_state": care, "n_seeds": len(recs), "person_years_per_seed": py / len(recs),
                     "deaths_per_seed": d / len(recs), "mortality_rate_per_1000py": 1000 * d / py if py else float("nan"),
                     "expected_deaths_per_seed": ex / len(recs), "OE_ratio": d / ex if ex else float("nan")})
    return pd.DataFrame(rows)


def survival_weight(cache):
    rows = []
    base_sa = SR.SURVIVAL_ADJ
    try:
        for comp in ("1m", "2m"):
            recs = [r for r in cache.values() if comp in r["arms"]]
            if not recs:
                continue
            dcost = np.array([r["arms"]["0m"]["cost_d"] - r["arms"][comp]["cost_d"] for r in recs])
            for sa in SURVIVAL_WEIGHT_GRID:
                SR.SURVIVAL_ADJ = sa
                dq, loss = [], []
                for r in recs:
                    prev, treat, averted = SR.decomposition(r["arms"]["0m"], r["arms"][comp])
                    dq.append(prev + treat)
                    loss.append(prev / averted if averted else float("nan"))
                dq = np.array(dq)
                nmb = dq * SR.WTP - dcost
                rows.append({"comparison": f"0m vs {comp}", "survival_weight": sa, "n_seeds": len(recs),
                             "qaly_loss_per_averted_infection": float(np.nanmean(loss)),
                             "dQALY": fmt(*ci(dq)), "dQALY_mean": dq.mean(), "dQALY_lo": ci(dq)[1], "dQALY_hi": ci(dq)[2],
                             "NMB_bn": fmt(*ci(nmb), nd=1, scale=1e9), "NMB_mean_bn": nmb.mean() / 1e9,
                             "NMB_lo_bn": ci(nmb)[1] / 1e9, "NMB_hi_bn": ci(nmb)[2] / 1e9,
                             "dominant_seeds": f"{int(np.sum((dcost < 0) & (dq > 0)))}/{len(recs)}",
                             "NMB_positive_seeds": f"{int(np.sum(nmb > 0))}/{len(recs)}"})
    finally:
        SR.SURVIVAL_ADJ = base_sa
    return pd.DataFrame(rows)


def subgroup(cache):
    rows = []
    for arm in SR.ARMS:
        v = np.array([r["arms"][arm]["highcd4"]["group_size"] for r in cache.values() if "highcd4" in r["arms"].get(arm, {})], float)
        f = np.array([r["arms"][arm]["highcd4"].get("fraction_of_hiv_positive", np.nan) for r in cache.values() if "highcd4" in r["arms"].get(arm, {})], float)
        if len(v):
            m, lo, hi = ci(v)
            rows.append({"slot": SR.LABEL[arm], "n_seeds": len(v), "group_size": fmt(m, lo, hi), "group_size_mean": m,
                         "group_size_lo": lo, "group_size_hi": hi, "fraction_of_people_with_hiv_percent": 100 * float(np.nanmean(f))})
    return pd.DataFrame(rows)


def inflow(cache):
    rows = []
    for arm in SR.ARMS:
        v = [r["arms"][arm]["inflow"] for r in cache.values() if "inflow" in r["arms"].get(arm, {})]
        if v:
            rows.append({"slot": SR.LABEL[arm], "n_seeds": len(v), "cumulative_foreign_born_entries_mean": float(np.mean(v))})
    return pd.DataFrame(rows)


def figure_data(caches):
    """Yearly series of every setting (analytic period) and of the burn-in of the base analysis."""
    rows = []
    for (an, st), cache in caches.items():
        for r in _figure_series(cache):
            rows.append({"analysis": an, "setting": st, **r})
    cache = caches[SR.BASE]
    burn = []
    recs = [r["burnin"] for r in cache.values() if "burnin" in r and "hiv_y" in r["burnin"]]
    for key, name in (("notif_y", "notifications"), ("inf_y", "new infections"), ("on_art_y", "on ART"),
                      ("hiv_y", "people with HIV"), ("undiag_y", "undiagnosed")):
        if not recs:
            break
        a = np.array([x[key] for x in recs], float)
        for i in range(a.shape[1]):
            burn.append({"series": name, "year": FIRST_BURNIN_YEAR + i, "n_seeds": a.shape[0], "mean": a[:, i].mean(),
                         "median": np.median(a[:, i]), "pi95_lo": np.percentile(a[:, i], 2.5), "pi95_hi": np.percentile(a[:, i], 97.5),
                         "observed": OBSERVED.get(FIRST_BURNIN_YEAR + i) if key == "notif_y" else None})
    if recs:
        # undiagnosed people as a percentage of people with HIV, per seed (panel d of Supplementary Figure S2)
        a = 100.0 * np.array([x["undiag_y"] for x in recs], float) / np.array([x["hiv_y"] for x in recs], float)
        for i in range(a.shape[1]):
            burn.append({"series": "undiagnosed fraction (%)", "year": FIRST_BURNIN_YEAR + i, "n_seeds": a.shape[0], "mean": a[:, i].mean(),
                         "median": np.median(a[:, i]), "pi95_lo": np.percentile(a[:, i], 2.5), "pi95_hi": np.percentile(a[:, i], 97.5),
                         "observed": None})
    return pd.DataFrame(rows), pd.DataFrame(burn)


def _figure_series(cache):
    rows = []
    for arm in SR.ARMS:
        for key, name, cumulative in (("inf_y", "new infections", False), ("dx_y", "new diagnoses", False),
                                      ("on_art_y", "on ART", False), ("inf_y", "cumulative new infections since 2024", True),
                                      ("deaths_hiv_y", "deaths among people with HIV", False),
                                      ("art_starts_y", "ART initiations", False),
                                      ("deaths_all_y", "deaths, whole population", False)):
            m = [r["arms"][arm][key] for r in cache.values() if arm in r["arms"] and r["arms"][arm].get(key)]
            if not m:
                continue
            a = np.array(m, float)
            if cumulative:
                a = np.cumsum(a, axis=1)
            for i in range(a.shape[1]):
                rows.append({"slot": SR.LABEL[arm], "series": name, "year": 2024 + i, "n_seeds": a.shape[0],
                             "mean": a[:, i].mean(), "median": np.median(a[:, i]),
                             "pi95_lo": np.percentile(a[:, i], 2.5), "pi95_hi": np.percentile(a[:, i], 97.5)})
    return rows


def convention_comparison(new_caches):
    """R2 base analysis and observed practice against the published outputs on the same seeds."""
    saved = SR.TARGET
    rows = []
    try:
        SR.TARGET = "published"
        pub = {k: SR.load_cache(*k) for k in (SR.BASE, SR.OBSERVED_PRACTICE)}
    finally:
        SR.TARGET = saved
    for key, comps, name in ((SR.BASE, ("1m", "2m"), "statutory minimum"), (SR.OBSERVED_PRACTICE, ("1m",), "observed practice")):
        if key not in new_caches or not pub.get(key):
            continue
        for comp in comps:
            fn = SR.frame(new_caches[key], "0m", new_caches[key], comp)
            fp = SR.frame(pub[key], "0m", pub[key], comp)
            m = fn.merge(fp, on="seed", suffixes=("_new", "_old"))
            if not len(m):
                continue
            row = {"comparator": f"{name}, {SR.LABEL[comp]} slot", "n_seeds": len(m),
                   "dominant_uniform": f"{int(m['dominant_new'].sum())}/{len(m)}", "dominant_earlier": f"{int(m['dominant_old'].sum())}/{len(m)}"}
            for k, nd, sc in (("infections_averted", 0, 1.0), ("deaths_averted", 0, 1.0), ("dQALY_decomposition", 0, 1.0),
                              ("dCost_CEA", 1, -1e9), ("NMB", 1, 1e9)):
                a, b = m[k + "_new"].to_numpy(float) / sc, m[k + "_old"].to_numpy(float) / sc
                row[k + " uniform convention"] = fmt(*ci(a), nd=nd)
                row[k + " earlier convention"] = fmt(*ci(b), nd=nd)
                row[k + " paired difference"] = fmt(*ci(a - b), nd=nd)
            rows.append(row)
    return pd.DataFrame(rows)


def main():
    SR.TARGET = TARGET
    od = out_dir()
    os.makedirs(od, exist_ok=True)
    settings = SR.list_settings()
    caches = {(an, st): SR.load_cache(an, st) for an, st, _ in settings}
    caches = {k: v for k, v in caches.items() if v}
    if SR.BASE not in caches:
        raise SystemExit("No cached seeds of the base analysis; run summarize_R2.py first.")
    base = caches[SR.BASE]
    doc = [f"# Supplementary tables from the per-seed cache ({'runs of model version 1.2.0' if TARGET == 'R2' else 'published outputs'})\n",
           f"Base analysis: {len(base)} seeds.\n"]

    # calibration: base and every setting whose burn-in differs from the base
    states, series, fits, firsts, shares = [], [], [], [], []
    for (an, st), cache in caches.items():
        label = f"{an} / {st}"
        if (an, st) != SR.BASE:
            common = sorted(set(cache) & set(base))
            if common and all(cache[s].get("burnin", {}).get("md5") == base[s].get("burnin", {}).get("md5") for s in common):
                continue                                   # same burn-in as the base
        c = calibration(cache, label)
        if c is None:
            continue
        states.append(c[0]); series.append(c[1]); fits.append(c[2])
        firsts.append(first_year(cache, label)); shares.append(het_share(cache, label))
    pd.concat(states).to_csv(os.path.join(od, "calibration_state_2023.csv"), index=False)
    pd.concat(series).to_csv(os.path.join(od, "calibration_notifications_by_year.csv"), index=False)
    pd.concat(fits).to_csv(os.path.join(od, "calibration_fit.csv"), index=False)
    pd.concat(firsts).to_csv(os.path.join(od, "first_analytic_year_notifications.csv"), index=False)
    pd.concat(shares).to_csv(os.path.join(od, "heterosexual_share.csv"), index=False)
    s0 = states[0].copy()
    s0["value"] = [f"{m:,.1f} ({a:,.1f} to {b:,.1f})" if "%" in i else f"{m:,.0f} ({a:,.0f} to {b:,.0f})"
                   for i, m, a, b in zip(s0["indicator"], s0["mean"], s0["range95_lo"], s0["range95_hi"])]
    doc.append("## State in 2023, base analysis (mean and 95% range across seeds, mean ± 1.96 SD)\n\n" + md(s0[["indicator", "n_seeds", "value"]]) + "\n")
    doc.append("## Fit to the observed notifications (seed median)\n\n" + md(pd.concat(fits).round(1)) + "\n")
    doc.append("## First analytic year (2024)\n\n" + md(pd.concat(firsts).round(0)) + "\n")
    doc.append("## Heterosexual share of incident infections over the analytic horizon\n\n" + md(pd.concat(shares).round(1)) + "\n")
    wide = pd.concat(states).pivot_table(index="setting", columns="indicator", values="mean", sort=False).round(1).reset_index()
    doc.append("## State in 2023 by setting (settings whose burn-in differs from the base)\n\n" + md(wide) + "\n")

    rg = risk_group(base)
    rg.to_csv(os.path.join(od, "risk_group.csv"), index=False)
    doc.append("## By risk group\n\n" + md(rg[["comparison", "risk_group", "n_seeds", "infections_averted", "deaths_hiv_averted"]]) + "\n")

    t7, t8 = transmission_source(base)
    if t7 is not None:
        t7.to_csv(os.path.join(od, "transmission_source_by_strategy.csv"), index=False)
        t8.to_csv(os.path.join(od, "transmission_source_additional.csv"), index=False)
        doc.append("## Domestic infections by source state\n\n" + md(t7.round(1)) + "\n\n" + md(t8.round(1)) + "\n")

    mc = monte_carlo(base)
    mc.to_csv(os.path.join(od, "monte_carlo.csv"), index=False)
    doc.append("## Monte Carlo properties\n\n" + md(mc.round(2)) + "\n")

    ab = absolute_by_arm(base)
    ab.to_csv(os.path.join(od, "absolute_by_arm.csv"), index=False)
    doc.append("## 40-year totals by strategy (means per seed)\n\n" + md(ab.round(1)) + "\n")

    dl = delays(caches)
    dl.to_csv(os.path.join(od, "delays.csv"), index=False)
    doc.append("## Waiting times (pooled over seeds)\n\n" + md(dl[["analysis", "setting", "slot", "interval", "n_seeds", "text"]]) + "\n")

    iw = ineligible_wait(caches)
    if len(iw):
        iw.to_csv(os.path.join(od, "ineligible_at_diagnosis.csv"), index=False)
        doc.append("## People not meeting the eligibility criteria at diagnosis, among ART starters (pooled over seeds)\n\n" + md(iw.round(1)) + "\n")

    bc, bs = breakeven_curves(caches)
    if len(bc):
        bc.to_csv(os.path.join(od, "breakeven_curve.csv"), index=False)
        bs.to_csv(os.path.join(od, "short_term_budget_impact.csv"), index=False)
        doc.append("## Undiscounted cost difference in the first five years, budget-impact convention (billion yen; mean and 95% CI across seeds)\n\n"
                   + md(bs[["analysis", "setting", "comparison", "year", "n_seeds", "undiscounted_difference_bn"]]) + "\n")
        be = bc[bc["year"].isin([20, 40])].pivot_table(index=["analysis", "setting", "comparison", "n_seeds", "break_even_year_of_the_mean_curve"],
                                                       columns="year", values="cumulative_bn_mean").round(1).reset_index()
        be.columns = [str(c) if not isinstance(c, str) else c for c in be.columns]
        doc.append("## Break-even year of the mean curve and cumulative discounted cost difference at years 20 and 40 (billion yen)\n\n" + md(be) + "\n")
    dh = delay_histograms(caches)
    if len(dh):
        dh.to_csv(os.path.join(od, "delay_histograms.csv"), index=False)

    mv = mortality_validation(base)
    if mv is not None:
        mv.to_csv(os.path.join(od, "mortality_validation.csv"), index=False)
        doc.append("## Mortality by care state and CD4 stratum, immediate-ART arm\n\n" + md(mv.round(2)) + "\n")

    sw = survival_weight(base)
    sw.to_csv(os.path.join(od, "survival_weight.csv"), index=False)
    doc.append("## Counterfactual survival weight\n\n" + md(sw[["comparison", "survival_weight", "n_seeds", "qaly_loss_per_averted_infection",
                                                              "dQALY", "NMB_bn", "dominant_seeds", "NMB_positive_seeds"]].round(2)) + "\n")

    sg = subgroup(base)
    if len(sg):
        sg.to_csv(os.path.join(od, "subgroup_high_cd4.csv"), index=False)
        doc.append("## Untreated subgroup with high CD4 count and low viral load\n\n" + md(sg[["slot", "n_seeds", "group_size", "fraction_of_people_with_hiv_percent"]].round(2)) + "\n")
    fl = inflow(base)
    if len(fl):
        fl.to_csv(os.path.join(od, "foreign_born_inflow.csv"), index=False)
        doc.append("## Foreign-born entries with HIV over the analytic horizon\n\n" + md(fl.round(0)) + "\n")

    fa, fb = figure_data(caches)
    fa.to_csv(os.path.join(od, "figure_data_analytic_period.csv"), index=False)
    fb.to_csv(os.path.join(od, "figure_data_burnin.csv"), index=False)

    if TARGET == "R2":
        cv = convention_comparison(caches)
        if len(cv):
            cv.to_csv(os.path.join(od, "convention_comparison.csv"), index=False)
            doc.append("## Uniform convention (this analysis) and earlier convention (published outputs), same seeds\n\n" + md(cv) + "\n")

    with open(os.path.join(od, "tables.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(doc))
    print(f"tables written to {od}")


if __name__ == "__main__":
    main()

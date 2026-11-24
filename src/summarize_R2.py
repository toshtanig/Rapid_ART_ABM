#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Summary of the simulation outputs: per-seed cache and contrast tables.

Step "cache": every completed seed of every setting is read once and stored as one line of
    post/<target>/cache/<analysis>__<setting>.jsonl
with, for each scenario slot (arm) of the seed, the absolute outcomes over the 40-year horizon (infections,
deaths among people with HIV, deaths in the whole simulated population, yearly QALYs and costs of the
cost-effectiveness analysis, incident-infection cohort sizes and mean QALYs, waiting times, transmission by
source state), the pairwise break-even results, and the state at the end of the burn-in (2023). Reading
through the synchronised folder is slow, so a run can be limited in time (MAX_SECONDS_READING) and repeated;
seeds already in the cache are not read again.

Step "tables": from the cache, the per-seed incremental outcomes of immediate ART (0 month slot) versus each
comparator slot are computed (infections and deaths averted, incremental QALYs by the prespecified prevention
decomposition and by the whole-population direct estimate, incremental discounted cost under the CEA cost
convention, net monetary benefit, dominance) and aggregated across seeds (mean, 95% CI across seeds). Because
the cache holds absolute outcomes per arm, contrasts between arms of different settings on the same seeds
(factorial decompositions) use the same functions. Also: discount-rate and time-horizon sensitivity with the
direct estimator, and consistency checks of the runs.

TARGET = "R2" summarises the runs of this folder (revision_R2_runs). TARGET = "published" applies the same
code to the published outputs (outputs/core_multiseed and outputs/revision_R1_runs, read only), to verify that it
reproduces the published numbers; its output goes to post/validation_published/ in this folder.

Run as-is in VS Code; all settings are in the CONFIG block, no command-line arguments.
"""

import os
import sys
import glob
import json
import math
import time
import hashlib
from io import StringIO

import numpy as np
import pandas as pd

# =============================================================================
# CONFIG
# =============================================================================
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _paths import BASE_144, BASE_MULTISEED, REVISION_R1_RUNS, REVISION_RUNS

TARGET = "R2"                    # "R2" or "published"
STEPS = ["cache", "tables"]      # which steps to run
SETTINGS_TO_READ = None          # None = every setting found; or a list of (analysis, setting) for the cache step
MAX_SECONDS_READING = None       # stop reading new seeds after this many seconds (None = no limit)

WTP = 5_000_000                  # yen per QALY
N_YEARS = 40
DISCOUNT_RATE = 0.02
U_NEG = 0.94                     # utility without HIV (prevention decomposition)
SURVIVAL_ADJ = 0.648             # counterfactual survival weight (prevention decomposition)
DISCOUNT_RATES = [0.00, 0.02, 0.03, 0.05]     # sensitivity, applied to costs and QALYs
SHORT_HORIZON_YEARS = 20

BASE = ("B0_core", "default")
OBSERVED_PRACTICE = ("A9_observed_ineligible_wait", "empirical_pathways_plus_observed_tail")
OBSERVED_CRITERIA_NO_WAIT = ("A9_observed_ineligible_wait", "criteria_retained_no_wait_plus_observed_tail")
STATUTORY_OBSERVED_TAIL = ("A9_observed_ineligible_wait", "one_month_plus_observed_tail")

R1_RUNS = REVISION_R1_RUNS
OUT_ROOT = {"R2": os.path.join(REVISION_RUNS, "post", "summary"),
            "published": os.path.join(REVISION_RUNS, "post", "validation_published", "summary")}
# =============================================================================

CACHE_VERSION = 5                # raised when read_seed stores more; older records are read again
ARMS = ["0m", "1m", "2m"]
LABEL = {"0m": "0 month", "1m": "1 month", "2m": "2 months"}
PREFIX = {"0m": "scenario0", "1m": "scenario1", "2m": "scenario2"}
DAYS_PER_MONTH = 365.25 / 12.0
MORTALITY_PERIOD = f"Y1-{N_YEARS}"   # row label of the whole horizon in *_mortality_by_cd4_care_state.csv
METRICS = ["infections_averted", "deaths_averted", "deaths_averted_all", "dQALY_decomposition", "dQALY_prevention",
           "dQALY_treatment", "dQALY_direct", "dCost_CEA", "dCost_BIA", "NMB", "NMB_direct"]


# -----------------------------------------------------------------------------
# Settings of the target
# -----------------------------------------------------------------------------
def list_settings():
    """[(analysis, setting, folder)] of the target, base first."""
    out = []
    if TARGET == "R2":
        for an in sorted(os.listdir(REVISION_RUNS)):
            d = os.path.join(REVISION_RUNS, an)
            if not os.path.isdir(d) or not (an[:1] in ("A", "B")) or an.startswith("A0"):
                continue
            for st in sorted(os.listdir(d)):
                if os.path.isdir(os.path.join(d, st)):
                    out.append((an, st, os.path.join(d, st)))
    elif TARGET == "published":
        out.append((BASE[0], BASE[1], BASE_MULTISEED))
        if os.path.isdir(R1_RUNS):
            for an in sorted(os.listdir(R1_RUNS)):
                d = os.path.join(R1_RUNS, an)
                if not os.path.isdir(d) or an.split("_")[0] not in ("A1", "A2", "A4", "A5", "A7", "A8", "A9", "A10", "A11", "A12"):
                    continue
                for st in sorted(os.listdir(d)):
                    if os.path.isdir(os.path.join(d, st)):
                        out.append((an, st, os.path.join(d, st)))
        for mult in (100, 110, 118, 125, 130):
            d = os.path.join(BASE_144, f"core_144_foi{mult}")
            if os.path.isdir(d):
                out.append(("A13_foi_multiplier", f"multiplier_{mult / 100:.2f}", d))
    else:
        raise SystemExit('TARGET must be "R2" or "published"')
    out.sort(key=lambda x: (x[:2] != BASE, x[0], x[1]))
    return out


def cache_path(an, st):
    return os.path.join(OUT_ROOT[TARGET], "cache", f"{an}__{st}.jsonl")


# -----------------------------------------------------------------------------
# Step "cache": read one seed
# -----------------------------------------------------------------------------
def _read_csv(path):
    try:
        return pd.read_csv(path)
    except Exception:
        with open(path, "rb") as f:
            data = f.read().replace(b"\x00", b"")
        return pd.read_csv(StringIO(data.decode("utf-8", "ignore")))


def _col(df, label, n=N_YEARS):
    """Yearly column of a Year-by-scenario table as a list of length n (None if the slot is absent)."""
    if df is None or label not in df.columns:
        return None
    s = df.set_index("Year")[label].reindex(range(1, n + 1)).fillna(0.0)
    return [float(v) for v in s.tolist()]


def _col_total(df, label, n=N_YEARS):
    """As _col; a slot that is split into columns "<label>_<stratum>" is summed over the strata."""
    if df is None:
        return None
    if label in df.columns:
        return _col(df, label, n)
    cols = [c for c in df.columns if c.startswith(label + "_")]
    if not cols:
        return None
    s = df.set_index("Year")[cols].sum(axis=1).reindex(range(1, n + 1)).fillna(0.0)
    return [float(v) for v in s.tolist()]


def _maybe(seed_dir, name):
    p = os.path.join(seed_dir, name)
    return _read_csv(p) if os.path.exists(p) else None


def seed_is_complete(seed_dir):
    return (os.path.exists(os.path.join(seed_dir, "heterogeneity_summary.csv"))
            and os.path.exists(os.path.join(seed_dir, "cea", "cea_report_ready.json"))
            and os.path.exists(os.path.join(seed_dir, "breakeven", "breakeven_report_ready.json")))


def read_seed(seed_dir):
    """All quantities of one seed that the tables need (see the module docstring)."""
    seed = int(os.path.basename(seed_dir).replace("seed_", ""))
    rec = {"seed": seed, "v": CACHE_VERSION, "arms": {}, "pairs": {}}
    with open(os.path.join(seed_dir, "cea", "cea_report_ready.json"), encoding="utf-8") as f:
        cea = json.load(f)["base_case_results"]
    het = _read_csv(os.path.join(seed_dir, "heterogeneity_summary.csv"))
    yearly = {k: _maybe(seed_dir, f"yearly_{k}.csv") for k in
              ("new_deaths", "deaths_among_hiv_positives", "new_infections", "new_diagnoses", "on_art_current",
               "undiagnosed_population", "new_art_starts", "cumulative_infected")}
    delay = {"diag": _maybe(seed_dir, "diagnosis_to_ART_delay_summary.csv"),
             "elig": _maybe(seed_dir, "eligible_to_ART_delay_summary.csv")}
    pop = _maybe(seed_dir, "yearly_population.csv")
    highcd4 = _maybe(seed_dir, "high_cd4_low_vl_longterm_noART_risk_summary.csv")
    inflow = _maybe(seed_dir, "foreign_inflow_summary.csv")
    mort = None
    mort_files = glob.glob(os.path.join(seed_dir, "*_mortality_by_cd4_care_state.csv"))
    if mort_files:
        mort = _read_csv(sorted(mort_files)[-1])
    for arm in ARMS:
        if arm not in cea:
            continue
        lab, pre = LABEL[arm], PREFIX[arm]
        a = {}
        coh = _read_csv(os.path.join(seed_dir, f"{pre}_incident_infection_cohort_summary.csv"))
        coh = coh.set_index("cohort_year").reindex(range(1, N_YEARS + 1)).fillna(0.0)
        a["n_inc"] = [float(v) for v in coh["n_incident"].tolist()]
        a["mq"] = [float(v) for v in coh["mean_discounted_qalys"].tolist()]
        h = het[het["Scenario"] == lab]
        a["deaths_hiv"] = float(h["Cum_Deaths_HIV_Pos"].sum())
        a["by_group"] = {str(r["Risk_Group"]): [float(r["Cum_New_Infections"]), float(r["Cum_Deaths_HIV_Pos"]),
                                                float(r["Cum_ART_Starts"])] for _, r in h.iterrows()}
        a["deaths_all_y"] = _col(yearly["new_deaths"], lab)
        a["deaths_hiv_y"] = _col(yearly["deaths_among_hiv_positives"], lab)
        a["inf_y"] = _col(yearly["new_infections"], lab)
        a["dx_y"] = _col(yearly["new_diagnoses"], lab)
        a["on_art_y"] = _col(yearly["on_art_current"], lab)
        # yearly_undiagnosed_population.csv counts people never tested in the whole population, not undiagnosed HIV
        a["untested_y"] = _col(yearly["undiagnosed_population"], lab)
        a["art_starts_y"] = _col_total(yearly["new_art_starts"], lab)
        a["cum_inf_y"] = _col(yearly["cumulative_infected"], lab)
        c = cea[arm]
        a["qaly_y"] = [float(v) for v in c["yearly_qaly"]]
        a["cost_y"] = [float(v) for v in c["yearly_cost"]]
        a["qaly_dy"] = [float(v) for v in c["yearly_qaly_discounted"]]
        a["cost_dy"] = [float(v) for v in c["yearly_cost_discounted"]]
        a["qaly_d"] = float(c["total_qaly_discounted"])
        a["cost_d"] = float(c["total_cost_discounted"])
        ys = _maybe(seed_dir, f"{pre}_yearly_summary.csv")       # people with HIV only
        if ys is not None:
            g = ys.groupby("Year")[["NewlyInfected", "NewlyDiagnosed", "Deaths"]].sum().reindex(range(1, N_YEARS + 1)).fillna(0.0)
            a["ys_inf_y"] = [float(v) for v in g["NewlyInfected"].tolist()]
            a["ys_dx_y"] = [float(v) for v in g["NewlyDiagnosed"].tolist()]
            a["ys_deaths_y"] = [float(v) for v in g["Deaths"].tolist()]
        for key, df in delay.items():
            if df is not None:
                r = df[df["scenario"] == lab]
                if len(r):
                    r = r.iloc[0]
                    a[f"delay_{key}"] = [float(r["mean_delay_months"]), float(r["median_delay_months"]),
                                         float(r["p25_delay_months"]), float(r["p75_delay_months"]),
                                         float(r["n_art_starters"])]
        ts = _maybe(seed_dir, f"{pre}_transmission_source.csv")
        if ts is not None and "source_category" in ts.columns:
            if "risk_group" in ts.columns and (ts["risk_group"] == "Total").any():
                ts = ts[ts["risk_group"] == "Total"]            # the file also has one block of rows per risk group
            a["source"] = {str(k): float(v) for k, v in
                           ts.groupby("source_category")["allocated_domestic_infections"].sum().items()}
        for key, name in (("hist_diag", f"{pre}_diag_to_art_delay_hist.csv"), ("hist_elig", f"{pre}_eligible_to_art_delay_hist.csv")):
            hd = _maybe(seed_dir, name)
            if hd is not None and "delay_months" in hd.columns:
                a[key] = [[float(m), float(c)] for m, c in zip(hd["delay_months"].tolist(), hd["count"].tolist())]
        # people who did not meet the eligibility criteria at diagnosis: number and time from diagnosis to ART
        di = _maybe(seed_dir, f"{pre}_diag_to_art_delay_individual.csv")
        ei = _maybe(seed_dir, f"{pre}_eligible_to_art_delay_individual.csv")
        if (di is not None and ei is not None and {"person_id", "diagnosis_month", "art_start_month"} <= set(di.columns)
                and {"person_id", "eligible_month"} <= set(ei.columns)):
            j = di.drop_duplicates("person_id").merge(ei.drop_duplicates("person_id")[["person_id", "eligible_month"]],
                                                      on="person_id", how="inner")
            late = j[j["eligible_month"] > j["diagnosis_month"]]
            wait = (late["art_start_month"] - late["diagnosis_month"]).round().astype(int)        # diagnosis to ART
            until = (late["eligible_month"] - late["diagnosis_month"]).round().astype(int)        # diagnosis to eligibility
            a["inelig"] = {"n": float(len(j)), "n_inelig": float(len(late)),
                           "hist": [[float(m), float(c)] for m, c in sorted(wait.value_counts().items())],
                           "hist_until_eligible": [[float(m), float(c)] for m, c in sorted(until.value_counts().items())]}
        a["pop_y"] = _col(pop, f"{lab}_total")
        if highcd4 is not None:
            r = highcd4[highcd4["scenario"] == lab]
            if len(r):
                r = r.iloc[0]
                a["highcd4"] = {k: float(r[k]) for k in ("group_size", "fraction_of_hiv_positive",
                                                         "expected_transmissions_per_year_from_group",
                                                         "approx_yearly_mortality_risk") if k in r.index}
        if inflow is not None:
            r = inflow[inflow["Scenario"] == lab]
            if len(r):
                a["inflow"] = float(r.iloc[0]["Cum_Inflow_Total"])
        if mort is not None:
            # the file has single years, five-year blocks and the whole horizon; only the whole horizon is used
            mh = mort[(mort["scenario"] == lab) & (mort["period"].astype(str) == MORTALITY_PERIOD)]
            g = mh.groupby(["cd4_stratum", "care_state"])[["deaths", "person_years", "expected_deaths"]].sum()
            a["mort"] = {f"{c}|{k}": [float(v["deaths"]), float(v["person_years"]), float(v["expected_deaths"])]
                         for (c, k), v in g.iterrows()}
        rec["arms"][arm] = a
    with open(os.path.join(seed_dir, "breakeven", "breakeven_report_ready.json"), encoding="utf-8") as f:
        be = json.load(f).get("breakeven_comparisons", {})
    for key, comp in be.items():
        y = comp.get("break_even_year")
        cum = ann = None
        yd = comp.get("yearly_data")
        if isinstance(yd, list) and yd and isinstance(yd[0], dict):
            if "cumulative_diff_billion_jpy" in yd[0]:
                cum = [float(d["cumulative_diff_billion_jpy"]) for d in yd]
            if "annual_diff_billion_jpy" in yd[0]:
                ann = [float(d["annual_diff_billion_jpy"]) for d in yd]
        # be_cum, be_ann: discounted cost difference by year (cumulative, annual), budget-impact cost convention
        rec["pairs"][key] = {"be_year": (float(y) if y is not None else None),
                             "be_final": float(comp["final_cumulative_billion_jpy"]), "be_cum": cum, "be_ann": ann}
    bp = os.path.join(seed_dir, "burnin_summary.csv")
    if os.path.exists(bp):
        b = _read_csv(bp)
        r = b.iloc[-1]
        with open(bp, "rb") as f:
            md5 = hashlib.md5(f.read()).hexdigest()
        rec["burnin"] = {"md5": md5,
                         "notif_y": [float(v) for v in b["new_diagnoses_year"].tolist()],
                         "inf_y": [float(v) for v in b["new_infections_year"].tolist()],
                         "on_art_y": [float(v) for v in b["on_art"].tolist()],
                         "hiv_y": [float(v) for v in b["hiv_positive"].tolist()],
                         "undiag_y": [float(v) for v in b["undiagnosed_hiv"].tolist()],
                         "supp_y": [float(v) for v in b["vl_suppressed"].tolist()],
                         "deaths_y": [float(v) for v in b["deaths_year"].tolist()],
                         "alive_y": [float(v) for v in b["total_alive"].tolist()],
                         "hiv": float(r["hiv_positive"]), "on_art": float(r["on_art"]),
                         "undiag": float(r["undiagnosed_hiv"]), "supp": float(r["vl_suppressed"]),
                         "notif": float(r["new_diagnoses_year"]), "inf": float(r["new_infections_year"])}
    ip = os.path.join(seed_dir, "run_info.json")
    if os.path.exists(ip):
        with open(ip, encoding="utf-8") as f:
            info = json.load(f)
        rec["run"] = {"terminal": info.get("terminal"), "minutes": info.get("minutes"),
                      "model_md5": info.get("model_file_md5"), "model_values": info.get("model_values")}
    return rec


def load_cache(an, st):
    """{seed: record} of one setting."""
    out = {}
    p = cache_path(an, st)
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except ValueError:
                    continue                      # line cut off by an interrupted run; the seed is read again
                if rec.get("v") == CACHE_VERSION:
                    out[int(rec["seed"])] = rec
    return out


def _drop_old_records(an, st):
    """Rewrite the cache file without records of an older CACHE_VERSION or cut-off lines."""
    p = cache_path(an, st)
    if not os.path.exists(p):
        return
    good = load_cache(an, st)
    with open(p, encoding="utf-8") as f:
        n_lines = sum(1 for line in f if line.strip())
    if n_lines != len(good):
        with open(p, "w", encoding="utf-8") as f:
            for seed in sorted(good):
                f.write(json.dumps(good[seed], separators=(",", ":")) + "\n")


def step_cache():
    t0 = time.time()
    os.makedirs(os.path.join(OUT_ROOT[TARGET], "cache"), exist_ok=True)
    n_new = n_left = 0
    out_of_time = False
    for an, st, d in list_settings():
        if SETTINGS_TO_READ is not None and (an, st) not in [tuple(x) for x in SETTINGS_TO_READ]:
            continue
        if MAX_SECONDS_READING is not None and time.time() - t0 > MAX_SECONDS_READING:
            out_of_time = True               # the remaining settings are not even listed (listing is slow too)
            break
        _drop_old_records(an, st)
        have = set(load_cache(an, st))
        todo = []
        for sd in glob.glob(os.path.join(d, "seed_*")):
            try:
                seed = int(os.path.basename(sd).replace("seed_", ""))
            except ValueError:
                continue
            if seed not in have and seed_is_complete(sd):
                todo.append((seed, sd))
        todo.sort()
        done_here = 0
        with open(cache_path(an, st), "a", encoding="utf-8") as f:
            for seed, sd in todo:
                if MAX_SECONDS_READING is not None and time.time() - t0 > MAX_SECONDS_READING:
                    break
                try:
                    rec = read_seed(sd)
                except Exception as e:
                    print(f"  [{an}/{st}] seed {seed}: not read ({type(e).__name__}: {e})")
                    continue
                f.write(json.dumps(rec, separators=(",", ":")) + "\n")
                f.flush()
                done_here += 1
        n_new += done_here
        n_left += len(todo) - done_here
        if todo:
            print(f"  [{an}/{st}] cached {len(have) + done_here} seeds ({done_here} new, {len(todo) - done_here} still to read)",
                  flush=True)
    if out_of_time:
        print(f"cache step: {n_new} seeds read in {time.time() - t0:.0f} s; time limit reached before all settings "
              "were looked at (run again)", flush=True)
        return max(n_left, 1)
    print(f"cache step: {n_new} seeds read in {time.time() - t0:.0f} s; {n_left} complete seeds still to read"
          + ("" if n_left == 0 else " (run again)"), flush=True)
    return n_left


# -----------------------------------------------------------------------------
# Step "tables": contrasts between arms
# -----------------------------------------------------------------------------
def _ci(vals):
    x = np.array([v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))], float)
    n = len(x)
    if n == 0:
        return float("nan"), float("nan"), float("nan"), 0
    m = x.mean()
    se = x.std(ddof=1) / math.sqrt(n) if n > 1 else 0.0
    t = 1.96 if n > 30 else 2.045
    return m, m - t * se, m + t * se, n


_RESID = np.array([sum((1.0 + DISCOUNT_RATE) ** (-(t - 1)) for t in range(y, N_YEARS + 1)) for y in range(1, N_YEARS + 1)])


def decomposition(arm_i, arm_c):
    """Prevention decomposition of the incremental QALYs (as prevention_qaly_dominance.compute_seed_comparison):
    averted infections valued at the QALY loss per infection of the comparator cohort of the same year, plus
    the change in mean QALYs of the infections common to both arms."""
    n_i, n_c = np.array(arm_i["n_inc"]), np.array(arm_c["n_inc"])
    mq_i, mq_c = np.array(arm_i["mq"]), np.array(arm_c["mq"])
    averted = n_c - n_i
    per_inf_loss = U_NEG * _RESID * SURVIVAL_ADJ - mq_c
    prevention = float(np.sum(averted * per_inf_loss))
    treatment = float(np.sum(np.minimum(n_i, n_c) * (mq_i - mq_c)))
    return prevention, treatment, float(averted.sum())


def contrast(arm_i, arm_c, pair=None):
    """Incremental outcomes of arm_i (intervention) versus arm_c (comparator) for one seed."""
    prev, treat, averted = decomposition(arm_i, arm_c)
    dq = prev + treat
    dcost = arm_i["cost_d"] - arm_c["cost_d"]
    dq_direct = arm_i["qaly_d"] - arm_c["qaly_d"]
    all_i, all_c = arm_i.get("deaths_all_y"), arm_c.get("deaths_all_y")
    out = {
        "infections_averted": averted,
        "deaths_averted": arm_c["deaths_hiv"] - arm_i["deaths_hiv"],
        "deaths_averted_all": (float(np.sum(all_c)) - float(np.sum(all_i))) if (all_i and all_c) else float("nan"),
        "dQALY_decomposition": dq, "dQALY_prevention": prev, "dQALY_treatment": treat,
        "dQALY_direct": dq_direct,
        "dCost_CEA": dcost,
        "dCost_BIA": (pair["be_final"] * 1e9) if pair else float("nan"),
        "NMB": dq * WTP - dcost,
        "NMB_direct": dq_direct * WTP - dcost,
        "dominant": bool(dcost < 0 and dq > 0),
        "NMB_positive": bool(dq * WTP - dcost > 0),
        "breakeven_year": (pair["be_year"] if pair and pair["be_year"] is not None else float("nan")),
        "comparator_infections": float(np.sum(arm_c["n_inc"])),
    }
    return out


def frame(cache_i, arm_i, cache_c, arm_c, seeds=None, same_folder_pair=None):
    """Per-seed contrast table of arm_i in cache_i versus arm_c in cache_c on the seeds common to both."""
    rows = []
    common = sorted(set(cache_i) & set(cache_c))
    if seeds is not None:
        common = [s for s in common if s in seeds]
    for s in common:
        ri, rc = cache_i[s], cache_c[s]
        if arm_i not in ri["arms"] or arm_c not in rc["arms"]:
            continue
        pair = ri["pairs"].get(same_folder_pair) if same_folder_pair else None
        row = contrast(ri["arms"][arm_i], rc["arms"][arm_c], pair)
        row["seed"] = s
        rows.append(row)
    return pd.DataFrame(rows)


def combine(frames, signs):
    """Per-seed linear combination of contrast tables (for differences of contrasts)."""
    base = frames[0][["seed"]].copy()
    for f in frames[1:]:
        base = base.merge(f[["seed"]], on="seed")
    out = base.copy()
    for k in METRICS:
        tot = np.zeros(len(base))
        for f, sg in zip(frames, signs):
            tot = tot + sg * base.merge(f[["seed", k]], on="seed")[k].to_numpy()
        out[k] = tot
    out["dominant"] = (out["dCost_CEA"] < 0) & (out["dQALY_decomposition"] > 0)
    out["NMB_positive"] = out["NMB"] > 0
    out["breakeven_year"] = float("nan")
    return out


def aggregate(df, labels):
    if df is None or len(df) == 0:
        return None
    out = dict(labels)
    out["n_seeds"] = int(len(df))
    for k in METRICS:
        m, lo, hi, n = _ci(df[k].tolist())
        out[k + "_mean"], out[k + "_lo"], out[k + "_hi"] = m, lo, hi
    out["dominant_seeds"] = f"{int(df['dominant'].sum())}/{len(df)}"
    out["NMB_positive_seeds"] = f"{int(df['NMB_positive'].sum())}/{len(df)}"
    be = df["breakeven_year"].dropna() if "breakeven_year" in df.columns else pd.Series(dtype=float)
    out["breakeven_year_median"] = float(be.median()) if len(be) else float("nan")
    if "comparator_infections" in df.columns:
        out["comparator_infections_mean"] = float(df["comparator_infections"].mean())
    return out


def _fmt(m, lo, hi, scale=1.0, nd=0):
    if isinstance(m, float) and math.isnan(m):
        return "NA"
    f = f"{{:,.{nd}f}}"
    return f"{f.format(m / scale)} ({f.format(lo / scale)} to {f.format(hi / scale)})"


def pretty(rows):
    recs = []
    for r in rows:
        recs.append({
            "analysis": r.get("analysis", ""), "contrast": r.get("contrast", ""), "n": r["n_seeds"],
            "Infections averted": _fmt(r["infections_averted_mean"], r["infections_averted_lo"], r["infections_averted_hi"]),
            "Deaths averted, people with HIV": _fmt(r["deaths_averted_mean"], r["deaths_averted_lo"], r["deaths_averted_hi"]),
            "Deaths averted, whole population": _fmt(r["deaths_averted_all_mean"], r["deaths_averted_all_lo"], r["deaths_averted_all_hi"]),
            "QALYs, decomposition": _fmt(r["dQALY_decomposition_mean"], r["dQALY_decomposition_lo"], r["dQALY_decomposition_hi"]),
            "QALYs, direct": _fmt(r["dQALY_direct_mean"], r["dQALY_direct_lo"], r["dQALY_direct_hi"]),
            "Cost saving, bn yen": _fmt(-r["dCost_CEA_mean"], -r["dCost_CEA_hi"], -r["dCost_CEA_lo"], 1e9, 1),
            "NMB, bn yen": _fmt(r["NMB_mean"], r["NMB_lo"], r["NMB_hi"], 1e9, 1),
            "Dominant": r.get("dominant_seeds", ""), "NMB > 0": r.get("NMB_positive_seeds", ""),
            "Break-even year (median)": ("" if math.isnan(r.get("breakeven_year_median", float("nan"))) else f"{r['breakeven_year_median']:.0f}"),
        })
    return pd.DataFrame(recs)


def md_table(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join("" if (isinstance(v, float) and math.isnan(v)) else str(v) for v in row.tolist()) + " |")
    return "\n".join(lines)


def sensitivity_direct(cache, arm_i, arm_c, label):
    """Discount-rate sensitivity and shorter horizon with the whole-population direct estimator."""
    rows_dr, rows_h = [], []
    recs = [r for r in cache.values() if arm_i in r["arms"] and arm_c in r["arms"]]
    if not recs:
        return rows_dr, rows_h
    t = np.arange(1, N_YEARS + 1)
    for rate in DISCOUNT_RATES:
        f = (1.0 + rate) ** (-(t - 1))
        dq = [float(np.sum((np.array(r["arms"][arm_i]["qaly_y"]) - np.array(r["arms"][arm_c]["qaly_y"])) * f)) for r in recs]
        dc = [float(np.sum((np.array(r["arms"][arm_i]["cost_y"]) - np.array(r["arms"][arm_c]["cost_y"])) * f)) for r in recs]
        nmb = [q * WTP - c for q, c in zip(dq, dc)]
        mq, lq, hq, n = _ci(dq); mc, lc, hc, _ = _ci(dc); mn, ln, hn, _ = _ci(nmb)
        rows_dr.append({"comparison": label, "discount_rate": rate, "n_seeds": n,
                        "Incremental QALYs, direct (95% CI)": _fmt(mq, lq, hq),
                        "Cost saving, billion yen (95% CI)": _fmt(-mc, -hc, -lc, 1e9, 1),
                        "NMB, billion yen (95% CI)": _fmt(mn, ln, hn, 1e9, 1),
                        "dominant_seeds": f"{int(np.sum((np.array(dc) < 0) & (np.array(dq) > 0)))}/{n}",
                        "NMB_positive_seeds": f"{int(np.sum(np.array(nmb) > 0))}/{n}"})
    h = SHORT_HORIZON_YEARS
    dq = [float(np.sum(np.array(r["arms"][arm_i]["qaly_dy"][:h]) - np.array(r["arms"][arm_c]["qaly_dy"][:h]))) for r in recs]
    dc = [float(np.sum(np.array(r["arms"][arm_i]["cost_dy"][:h]) - np.array(r["arms"][arm_c]["cost_dy"][:h]))) for r in recs]
    nmb = [q * WTP - c for q, c in zip(dq, dc)]
    have_ys = all("ys_inf_y" in r["arms"][arm_i] and "ys_inf_y" in r["arms"][arm_c] for r in recs)
    inf = [float(np.sum(r["arms"][arm_c]["ys_inf_y"][:h]) - np.sum(r["arms"][arm_i]["ys_inf_y"][:h])) for r in recs] if have_ys else []
    dth = [float(np.sum(r["arms"][arm_c]["ys_deaths_y"][:h]) - np.sum(r["arms"][arm_i]["ys_deaths_y"][:h])) for r in recs] if have_ys else []
    mq, lq, hq, n = _ci(dq); mc, lc, hc, _ = _ci(dc); mn, ln, hn, _ = _ci(nmb)
    mi, li, hi, _ = _ci(inf); md_, ld, hd, _ = _ci(dth)
    rows_h.append({"comparison": label, "horizon_years": h, "n_seeds": n,
                   "Infections averted (95% CI)": _fmt(mi, li, hi),
                   "Deaths averted, people with HIV (95% CI)": _fmt(md_, ld, hd),
                   "Incremental QALYs, direct (95% CI)": _fmt(mq, lq, hq),
                   "Cost saving, billion yen (95% CI)": _fmt(-mc, -hc, -lc, 1e9, 1),
                   "NMB, billion yen (95% CI)": _fmt(mn, ln, hn, 1e9, 1),
                   "dominant_seeds": f"{int(np.sum((np.array(dc) < 0) & (np.array(dq) > 0)))}/{n}"})
    return rows_dr, rows_h


def step_tables():
    out_dir = OUT_ROOT[TARGET]
    os.makedirs(out_dir, exist_ok=True)
    settings = list_settings()
    caches = {(an, st): load_cache(an, st) for an, st, _ in settings}
    caches = {k: v for k, v in caches.items() if v}
    if BASE not in caches:
        print("No cached seeds of the base analysis yet; nothing to tabulate.")
        return
    base = caches[BASE]
    rows, md, per_seed = [], [], {}

    def add(df, analysis, name, key):
        a = aggregate(df, {"analysis": analysis, "contrast": name})
        if a:
            rows.append(a)
            per_seed[key] = df
        return a

    # ---- base analysis ------------------------------------------------------------------
    block = []
    for comp in ("1m", "2m"):
        df = frame(base, "0m", base, comp, same_folder_pair=f"0m vs {comp}")
        a = add(df, BASE[0], f"immediate ART vs {LABEL[comp]} delay", f"{BASE[0]}__{BASE[1]}__0m_vs_{comp}")
        if a:
            block.append(a)
    md.append("## Base analysis (statutory-minimum comparators)\n\n" + md_table(pretty(block)) + "\n")

    # ---- every other setting: 0 month slot versus its comparator slot(s), with the base on the same seeds
    for (an, st), cache in caches.items():
        if (an, st) == BASE:
            continue
        block = []
        seeds = set(cache)
        arms_here = set().union(*[set(r["arms"]) for r in cache.values()])
        for comp in ("1m", "2m"):
            if comp not in arms_here:
                continue
            df = frame(cache, "0m", cache, comp, same_folder_pair=f"0m vs {comp}")
            a = add(df, an, f"{st}: 0 month slot vs {LABEL[comp]} slot", f"{an}__{st}__0m_vs_{comp}")
            if a:
                block.append(a)
            dfb = frame(base, "0m", base, comp, seeds=seeds, same_folder_pair=f"0m vs {comp}")
            b = add(dfb, an, f"{st}: base analysis on the same seeds, vs {LABEL[comp]} delay", f"{an}__{st}__base_same_seeds_{comp}")
            if b:
                block.append(b)
        if block:
            md.append(f"## {an} / {st}\n\n" + md_table(pretty(block)) + "\n")

    # ---- factorial decomposition, statutory-minimum assumptions (A1) -------------------------------
    kAB, kCD = ("A1_factorial", "criteria_retained_no_wait"), ("A1_factorial", "criteria_removed_one_month")
    if kAB in caches and kCD in caches:
        fAB = frame(caches[kAB], "0m", caches[kAB], "1m")
        fCD = frame(caches[kCD], "0m", caches[kCD], "1m")
        seeds = set(fAB["seed"]) & set(fCD["seed"])
        fAD = frame(base, "0m", base, "1m", seeds=seeds)
        fAB, fCD = fAB[fAB["seed"].isin(fAD["seed"])], fCD[fCD["seed"].isin(fAD["seed"])]
        block = [
            add(fAD, "A1_factorial", "Full reform (A - D), same seeds", "A1__full"),
            add(fAB, "A1_factorial", "Waiting interval removed, criteria retained (A - B)", "A1__wait_given_criteria"),
            add(fCD, "A1_factorial", "Waiting interval removed, criteria removed (C - D)", "A1__wait_given_no_criteria"),
            add(combine([fAD, fCD], [1, -1]), "A1_factorial", "Criteria removed, one-month wait retained (A - C)", "A1__criteria_given_wait"),
            add(combine([fAD, fAB], [1, -1]), "A1_factorial", "Criteria removed, no wait (B - D)", "A1__criteria_given_no_wait"),
            add(combine([fAB, fCD], [1, -1]), "A1_factorial", "Interaction (A - B) - (C - D)", "A1__interaction"),
        ]
        md.append("## Factorial decomposition, statutory-minimum assumptions (A1)\n\n" + md_table(pretty([b for b in block if b])) + "\n")

    # ---- decomposition with the observed wait of people not meeting the criteria (statutory minimum + tail)
    # Cell A' (statutory minimum with the observed ineligible wait) is the 1 month slot of the setting kB. The
    # dedicated setting STATUTORY_OBSERVED_TAIL holds the same arm on fewer seeds and is used for the economics.
    # Cell C (criteria removed, one-month wait retained) does not involve the ineligible wait and is taken from
    # the A1 factorial on the same seeds.
    kB = OBSERVED_CRITERIA_NO_WAIT
    if kB in caches:
        fAB = frame(caches[kB], "0m", caches[kB], "1m")                           # B' vs A'
        fAD = frame(base, "0m", caches[kB], "1m", seeds=set(fAB["seed"]))         # D vs A'
        fAB = fAB[fAB["seed"].isin(fAD["seed"])]
        if len(fAD):
            block = [
                add(fAD, "A9_statutory_tail", "Full reform vs statutory minimum with observed ineligible wait (A' - D)", "A9tail__full"),
                add(fAB, "A9_statutory_tail", "Waiting interval removed, criteria retained (A' - B')", "A9tail__wait"),
                add(combine([fAD, fAB], [1, -1]), "A9_statutory_tail", "Criteria removed, no wait (B' - D)", "A9tail__criteria"),
            ]
            if kCD in caches:
                fCD = frame(caches[kCD], "0m", caches[kCD], "1m", seeds=set(fAD["seed"]))
                if len(fCD) == len(fAD):
                    block += [
                        add(fCD, "A9_statutory_tail", "Waiting interval removed, criteria removed (C - D)", "A9tail__wait_given_no_criteria"),
                        add(combine([fAD, fCD], [1, -1]), "A9_statutory_tail", "Criteria removed, one-month wait retained (A' - C)", "A9tail__criteria_given_wait"),
                        add(combine([fAB, fCD], [1, -1]), "A9_statutory_tail", "Interaction (A' - B') - (C - D)", "A9tail__interaction"),
                    ]
            md.append("## Decomposition against the statutory minimum with the observed ineligible wait\n\n"
                      + md_table(pretty([b for b in block if b])) + "\n")

    # ---- decomposition against the observed-practice comparator ------------------------------------
    kO = OBSERVED_PRACTICE
    if kO in caches and kB in caches:
        seeds = set(caches[kO]) & set(caches[kB])
        f_total = frame(caches[kO], "0m", caches[kO], "1m", seeds=seeds)                      # D vs A''
        f_wait = frame(caches[kB], "0m", caches[kO], "1m", seeds=seeds)                       # B' vs A''
        f_crit = frame(caches[kO], "0m", caches[kB], "0m", seeds=seeds)                       # D vs B' (direct)
        if len(f_total):
            block = [
                add(f_total, "A9_observed_practice", "Full reform vs observed practice (A'' - D)", "A9obs__full"),
                add(f_wait, "A9_observed_practice", "Waiting interval removed, criteria retained (A'' - B')", "A9obs__wait"),
                add(combine([f_total, f_wait], [1, -1]), "A9_observed_practice",
                    "Criteria removed (B' - D) = (A'' - D) - (A'' - B')", "A9obs__criteria"),
                add(f_crit, "A9_observed_practice", "Criteria removed, direct contrast of the two arms (check)", "A9obs__criteria_direct"),
            ]
            md.append("## Decomposition against the observed-practice comparator\n\n" + md_table(pretty([b for b in block if b])) + "\n")

    # ---- cost of inaction (delayed implementation) -----------------------------------------------------
    block = []
    for (an, st), cache in caches.items():
        if not an.startswith("A7"):
            continue
        fd = frame(cache, "0m", cache, "1m")
        fb = frame(base, "0m", base, "1m", seeds=set(fd["seed"]))
        fd = fd[fd["seed"].isin(fb["seed"])]
        if len(fd):
            block.append(add(combine([fb, fd], [1, -1]), an, f"cost of inaction: {st} versus 2024", f"{an}__{st}__cost_of_inaction"))
    if block:
        md.append("## Cost of inaction (implementation in 2024 minus delayed implementation)\n\n"
                  + md_table(pretty([b for b in block if b])) + "\n")

    # ---- discount rate and horizon ---------------------------------------------------------------------
    rows_dr, rows_h = [], []
    for label, cache, comp in [("statutory minimum: 0m vs 1m", base, "1m"), ("statutory minimum: 0m vs 2m", base, "2m"),
                               ("observed practice: 0m vs comparator", caches.get(kO, {}), "1m")]:
        a, b = sensitivity_direct(cache, "0m", comp, label)
        rows_dr += a
        rows_h += b
    if rows_dr:
        pd.DataFrame(rows_dr).to_csv(os.path.join(out_dir, "discount_rate_sensitivity_direct.csv"), index=False)
        pd.DataFrame(rows_h).to_csv(os.path.join(out_dir, "horizon_sensitivity_direct.csv"), index=False)
        md.append("## Discount-rate sensitivity (direct estimator; same rate for costs and QALYs)\n\n" + md_table(pd.DataFrame(rows_dr))
                  + f"\n\n## {SHORT_HORIZON_YEARS}-year horizon (direct estimator; discounted at 2% per year)\n\n" + md_table(pd.DataFrame(rows_h)) + "\n")

    # ---- consistency checks ----------------------------------------------------------------------------------
    checks = []
    md5s = {}
    for (an, st), cache in caches.items():
        for r in cache.values():
            m = (r.get("run") or {}).get("model_md5")
            if m:
                md5s[m] = md5s.get(m, 0) + 1
    if md5s:
        checks.append(f"- model file versions among cached seeds: {md5s}")
    for (an, st), cache in caches.items():
        if (an, st) == BASE:
            continue
        common = sorted(set(cache) & set(base))
        if not common:
            continue
        same_burn = sum(1 for s in common if cache[s].get("burnin", {}).get("md5") == base[s].get("burnin", {}).get("md5"))
        same_d = sum(1 for s in common if "0m" in cache[s]["arms"] and cache[s]["arms"]["0m"]["n_inc"] == base[s]["arms"]["0m"]["n_inc"]
                     and abs(cache[s]["arms"]["0m"]["qaly_d"] - base[s]["arms"]["0m"]["qaly_d"]) < 1e-6)
        checks.append(f"- {an} / {st}: {len(cache)} seeds; burn-in identical to the base on {same_burn}/{len(common)} common seeds; "
                      f"0 month slot identical to the base on {same_d}/{len(common)}")
    checks.insert(0, f"- {BASE[0]} / {BASE[1]}: {len(base)} seeds")

    pd.DataFrame(rows).to_csv(os.path.join(out_dir, "summary_contrasts.csv"), index=False)
    os.makedirs(os.path.join(out_dir, "per_seed"), exist_ok=True)
    for key, df in per_seed.items():
        df.to_csv(os.path.join(out_dir, "per_seed", f"{key}.csv"), index=False)
    title = "runs of model version 1.2.0 (revision_R2_runs)" if TARGET == "R2" else "published outputs (validation of this script)"
    with open(os.path.join(out_dir, "summary.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# Summary of the {title}\n\n"
                f"Generated {time.strftime('%Y-%m-%d %H:%M')}. Immediate ART (0 month slot) versus the comparator slot; 40-year horizon; "
                f"discounted at 2% per year; costs under the CEA cost convention; WTP {WTP:,} yen per QALY; SURVIVAL_ADJ = {SURVIVAL_ADJ}; "
                "mean with 95% CI across seeds (normal approximation; stochastic variation between seeds only).\n\n"
                + "\n".join(md) + "\n## Consistency checks\n\n" + "\n".join(checks) + "\n")
    print(f"tables step: {len(rows)} contrasts written to {out_dir}")
    for c in checks[:60]:
        print(c)


def main():
    if "cache" in STEPS:
        step_cache()
    if "tables" in STEPS:
        step_tables()


if __name__ == "__main__":
    main()

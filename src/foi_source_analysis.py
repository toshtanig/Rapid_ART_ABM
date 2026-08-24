#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Transmission-source (care-state) decomposition of averted infections."""

import os
import sys
import glob
import math
import json
import numpy as np
import pandas as pd

# =============================================================================
# CONFIG
# =============================================================================
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# --- Data location: set BASE_144 to your local 144_scenario folder ---
BASE_144 = r"/path/to/20250315_169/144_scenario"
PROJ_169 = os.path.dirname(BASE_144)
BASE_169 = PROJ_169
BASE_104  = os.path.join(BASE_169, "144_scenario")
SCRIPTS   = os.path.join(BASE_104, "src")

SOURCE_SCAN_DIR = os.path.join(BASE_104, "core_foi118")

FOI_DIRS = {
    1.00: os.path.join(BASE_104, "core_foi100"),
    1.10: os.path.join(BASE_104, "core_foi110"),
    1.18: os.path.join(BASE_104, "core_foi118"),
    1.25: os.path.join(BASE_104, "core_foi125"),
    1.30: os.path.join(BASE_104, "core_foi130"),
}

OUTPUT_DIR = os.path.join(BASE_104, "foi_source_analysis")
WTP = 5_000_000
SCEN_MAP = {"scenario0": "0 month", "scenario1": "1 month", "scenario2": "2 months"}
DOMESTIC_CATS = ["undiagnosed", "diagnosed_pre_art", "on_art_unsuppressed", "disengaged"]
# =============================================================================

sys.path.insert(0, SCRIPTS)
import prevention_qaly_dominance as P   # reuse compute_seed_comparison + SURVIVAL_ADJ=0.645

def _ci(x):
    x = np.array([v for v in x if v is not None and not (isinstance(v, float) and math.isnan(v))], float)
    n = len(x)
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    m = x.mean()
    se = x.std(ddof=1) / math.sqrt(n) if n > 1 else 0.0
    t = 1.96 if n > 30 else 2.045
    return m, m - t * se, m + t * se

# -----------------------------------------------------------------------------
# (A) Transmission-source decomposition
# -----------------------------------------------------------------------------
def source_decomposition():
    seed_dirs = sorted(glob.glob(os.path.join(SOURCE_SCAN_DIR, "seed_*")))
    # per-seed cumulative allocated infections by (scenario, category) using Total rows
    per = {sc: {c: [] for c in DOMESTIC_CATS + ["exogenous_foreign_inflow"]} for sc in SCEN_MAP.values()}
    n_used = 0
    for d in seed_dirs:
        ok = True
        for scf, sclabel in SCEN_MAP.items():
            path = os.path.join(d, f"{scf}_transmission_source.csv")
            if not os.path.exists(path):
                ok = False
                break
        if not ok:
            continue
        n_used += 1
        for scf, sclabel in SCEN_MAP.items():
            df = pd.read_csv(os.path.join(d, f"{scf}_transmission_source.csv"))
            tot = df[df["risk_group"] == "Total"]
            for c in DOMESTIC_CATS + ["exogenous_foreign_inflow"]:
                per[sclabel][c].append(
                    float(tot[tot["source_category"] == c]["allocated_domestic_infections"].sum()))

    # Table 1: composition by scenario (mean cumulative + % of domestic total)
    rows1, rows1pct = [], []
    means = {}
    for sclabel in SCEN_MAP.values():
        means[sclabel] = {c: np.mean(per[sclabel][c]) if per[sclabel][c] else float("nan")
                          for c in DOMESTIC_CATS + ["exogenous_foreign_inflow"]}
        dom_tot = sum(means[sclabel][c] for c in DOMESTIC_CATS)
        rows1.append({"scenario": sclabel, **{c: round(means[sclabel][c]) for c in DOMESTIC_CATS},
                      "domestic_total": round(dom_tot),
                      "exogenous_foreign_inflow": round(means[sclabel]["exogenous_foreign_inflow"])})
        rows1pct.append({"scenario": sclabel,
                         **{c: (round(100 * means[sclabel][c] / dom_tot, 1) if dom_tot > 0 else float("nan"))
                            for c in DOMESTIC_CATS}})
    df1 = pd.DataFrame(rows1); df1pct = pd.DataFrame(rows1pct)

    # Table 2: additional domestic infections under delay, by source state
    rows2 = []
    for comp, sclabel in [("1m minus 0m", "1 month"), ("2m minus 0m", "2 months")]:
        add = {c: means[sclabel][c] - means["0 month"][c] for c in DOMESTIC_CATS}
        add_tot = sum(add.values())
        rows2.append({"comparison": comp, "additional_domestic_infections": round(add_tot),
                      **{c: round(add[c]) for c in DOMESTIC_CATS},
                      **{c + "_pct": (round(100 * add[c] / add_tot, 1) if add_tot != 0 else float("nan"))
                         for c in DOMESTIC_CATS}})
    df2 = pd.DataFrame(rows2)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df1.to_csv(os.path.join(OUTPUT_DIR, "source_composition_by_scenario.csv"), index=False)
    df1pct.to_csv(os.path.join(OUTPUT_DIR, "source_composition_pct.csv"), index=False)
    df2.to_csv(os.path.join(OUTPUT_DIR, "source_additional_under_delay.csv"), index=False)
    print(f"[source] seeds used: {n_used}  (from {SOURCE_SCAN_DIR})")
    print(df1.to_string(index=False)); print(); print(df1pct.to_string(index=False))
    print(); print(df2.to_string(index=False))
    return df1, df1pct, df2

# -----------------------------------------------------------------------------
# (B) FOI multiplier sensitivity
# -----------------------------------------------------------------------------
def _total_infections(seed_dir, scen_file):
    """Cumulative total new infections for a scenario from heterogeneity_summary.csv."""
    p = os.path.join(seed_dir, "heterogeneity_summary.csv")
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p)
    label = {"scenario0": "0 month", "scenario1": "1 month", "scenario2": "2 months"}[scen_file]
    return float(df[df["Scenario"] == label]["Cum_New_Infections"].sum())

def _cea_dcost(seed_dir, inter, comp):
    """Per-seed cost-effectiveness incremental cost (off-ART cost ratio 0.45),
    i.e. the whole-population discounted cost difference inter - comp, read from
    the seed's cea_report_ready.json. This keeps the FOI sensitivity dominance
    assessment on the same cost definition as the main analysis."""
    p = os.path.join(seed_dir, "cea", "cea_report_ready.json")
    if not os.path.exists(p):
        return None
    try:
        bc = json.load(open(p, encoding="utf-8"))["base_case_results"]
        return float(bc[inter]["total_cost_discounted"]) - float(bc[comp]["total_cost_discounted"])
    except Exception:
        return None

def _missing_inputs(d):
    miss = []
    for scf in ("scenario0", "scenario1", "scenario2"):
        if not os.path.exists(os.path.join(d, f"{scf}_incident_infection_cohort_summary.csv")):
            miss.append(scf + "_cohort")
    if not os.path.exists(os.path.join(d, "breakeven", "breakeven_report_ready.json")):
        miss.append("breakeven_json")
    return miss

def foi_sensitivity():
    rows = []
    per_seed_rows = []
    skipped = []
    for mult, mdir in FOI_DIRS.items():
        seed_dirs = sorted(glob.glob(os.path.join(mdir, "seed_*")))
        for comp_lbl, comp in [("0m vs 1m", "1m"), ("0m vs 2m", "2m")]:
            av, dq, dc, nmb, dom = [], [], [], [], []
            inf0, infc = [], []
            for d in seed_dirs:
                r = P.compute_seed_comparison(d, "0m", comp)
                _dcost = _cea_dcost(d, "0m", comp)   # CEA cost definition (off-ART 0.45)
                if r is None or _dcost is None:
                    miss = _missing_inputs(d)
                    if _dcost is None:
                        miss.append("cea_json")
                    skipped.append({"analytic_FOI_multiplier": mult, "comparison": comp_lbl,
                                    "seed": os.path.basename(d).replace("seed_", ""),
                                    "missing": ";".join(miss) if miss else "incomplete/other"})
                    continue
                _nmb = r["dQALY"] * WTP - _dcost
                _dom = _dcost < 0 and r["dQALY"] > 0
                av.append(r["infections_averted"])
                dq.append(r["dQALY"])
                dc.append(_dcost)
                nmb.append(_nmb)
                dom.append(_dom)
                per_seed_rows.append({
                    "analytic_FOI_multiplier": mult, "comparison": comp_lbl,
                    "seed": os.path.basename(d).replace("seed_", ""),
                    "infections_averted": round(r["infections_averted"]),
                    "dCost_JPY": round(_dcost), "dCost_billionJPY": round(_dcost/1e9, 2),
                    "dQALY": round(r["dQALY"], 1),
                    "NMB_JPY": round(_nmb), "NMB_billionJPY": round(_nmb/1e9, 2),
                    "dominant": bool(_dom), "NMB_positive": bool(_nmb > 0),
                })
                i0 = _total_infections(d, "scenario0")
                ic = _total_infections(d, "scenario1" if comp == "1m" else "scenario2")
                if i0 is not None:
                    inf0.append(i0); infc.append(ic)
            m_av, lo_av, hi_av = _ci(av)
            m_dq, _, _ = _ci(dq)
            m_dc, _, _ = _ci(dc)
            m_nmb, lo_nmb, hi_nmb = _ci(nmb)
            rows.append({
                "analytic_FOI_multiplier": mult, "comparison": comp_lbl, "n_seeds": len(av),
                "infections_0m": round(np.mean(inf0)) if inf0 else None,
                "infections_comparator": round(np.mean(infc)) if infc else None,
                "infections_averted": round(m_av), "infections_averted_CI": f"{round(lo_av)}-{round(hi_av)}",
                "incremental_QALY": round(m_dq),
                "cost_saving_billionJPY": round(-m_dc / 1e9, 1),
                "NMB_billionJPY": round(m_nmb / 1e9, 1), "NMB_CI": f"{round(lo_nmb/1e9,1)}-{round(hi_nmb/1e9,1)}",
                "dominant_seeds": f"{int(np.sum(dom))}/{len(dom)}",
            })
    df = pd.DataFrame(rows)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df.to_csv(os.path.join(OUTPUT_DIR, "foi_sensitivity.csv"), index=False)

    df_ps = pd.DataFrame(per_seed_rows).sort_values(
        ["analytic_FOI_multiplier", "comparison", "seed"]).reset_index(drop=True)
    df_ps.to_csv(os.path.join(OUTPUT_DIR, "foi_sensitivity_per_seed.csv"), index=False)

    short = {"0m vs 1m": "0v1", "0m vs 2m": "0v2"}
    dpw = df_ps.copy()
    dpw["c"] = dpw["comparison"].map(short)
    wide = dpw.pivot_table(
        index=["analytic_FOI_multiplier", "seed"], columns="c",
        values=["infections_averted", "dCost_billionJPY", "dQALY",
                "NMB_billionJPY", "dominant", "NMB_positive"], aggfunc="first")
    wide.columns = [f"{v}_{c}" for v, c in wide.columns]
    wide = wide.reset_index().sort_values(["analytic_FOI_multiplier", "seed"])
    wide.to_csv(os.path.join(OUTPUT_DIR, "foi_sensitivity_per_seed_wide.csv"), index=False)

    if skipped:
        df_sk = pd.DataFrame(skipped).sort_values(
            ["analytic_FOI_multiplier", "comparison", "seed"])
        df_sk.to_csv(os.path.join(OUTPUT_DIR, "foi_skipped_seeds.csv"), index=False)
        print("\n[skipped seeds: folder present but inputs incomplete -> excluded]")
        print(df_sk.to_string(index=False))
    print("\n[FOI sensitivity]")
    print(df.to_string(index=False))

    bad = df_ps[(~df_ps["dominant"]) | (~df_ps["NMB_positive"])]
    if len(bad):
        print("\n[non-dominant or NMB<=0 seeds]")
        print(bad[["analytic_FOI_multiplier", "comparison", "seed",
                   "dCost_billionJPY", "dQALY", "NMB_billionJPY",
                   "dominant", "NMB_positive"]].to_string(index=False))
    return df

def main():
    print("SURVIVAL_ADJ (from prevention module):", P.SURVIVAL_ADJ)
    source_decomposition()
    foi_sensitivity()
    print(f"\nSaved to: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()

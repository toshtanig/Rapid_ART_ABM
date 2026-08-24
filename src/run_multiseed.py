#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Single-machine driver: run the microsimulation for the fixed publication seeds.

Runs hiv_rapid_art_model.py once per seed (SEEDS in CONFIG; 200 for the published
analysis), and for each seed runs the cost-effectiveness and break-even analyses,
appending one row per seed to multiseed_comparison.csv (interruptible). A common
random-number design re-seeds the three ART-timing scenarios identically within
each replicate so incremental outcomes are paired per-seed differences.

All run conditions are set in the CONFIG section below; no command-line arguments.
"""

import os
import sys
import gc
import json
import time
import datetime
from io import StringIO

# =============================================================================

# =============================================================================

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = HERE

MULTISEED_DIR = os.path.join(HERE, "..", "outputs", "core_multiseed")

SEEDS = [
    14, 18, 24, 42, 55, 58, 61, 66, 69, 71,
    75, 77, 87, 97, 101, 102, 104, 107, 110, 119,
    123, 128, 145, 150, 157, 161, 176, 178, 188, 194,
    202, 203, 206, 207, 209, 216, 219, 225, 228, 234,
    244, 247, 249, 276, 286, 303, 306, 308, 309, 313,
    317, 327, 330, 340, 341, 345, 347, 355, 356, 363,
    364, 380, 389, 404, 417, 424, 427, 441, 448, 456,
    465, 466, 469, 472, 477, 484, 486, 488, 492, 495,
    504, 505, 508, 515, 516, 527, 536, 543, 547, 567,
    581, 582, 585, 586, 590, 592, 596, 606, 608, 613,
    615, 619, 624, 637, 641, 642, 650, 652, 663, 678,
    684, 698, 700, 707, 713, 720, 721, 731, 747, 754,
    762, 769, 771, 777, 780, 789, 790, 798, 803, 804,
    808, 811, 814, 820, 824, 826, 842, 844, 845, 857,
    865, 871, 880, 887, 890, 892, 895, 901, 909, 911,
    912, 925, 936, 937, 940, 943, 944, 956, 968, 973,
    989, 1001, 1002, 1009, 1010, 1013, 1014, 1025, 1026, 1032,
    1044, 1088, 1098, 1106, 1111, 1123, 1124, 1133, 1159, 1167,
    1185, 1190, 1191, 1196, 1198, 1210, 1211, 1212, 1226, 1234,
    1237, 1245, 1247, 1259, 1267, 1271, 1282, 1293, 1296, 1298,
]   # 200 publication seeds

RUN_CEA       = True
RUN_BREAKEVEN = True

PARALLEL_SCENARIOS = True

SKIP_PLOTS     = True
KEEP_PLOT_FILES = ["summary_plots.jpg"]

CEA_PSA_ITERATIONS = 200

BREAKEVEN_OK_RANGE = (5, 15)

OUTPUT_CSV = os.path.join(MULTISEED_DIR, "multiseed_comparison.csv")
# =============================================================================

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

def _keep(fname):
    if not SKIP_PLOTS:
        return True
    base = os.path.basename(fname) if isinstance(fname, str) else ""
    return base in KEEP_PLOT_FILES

_ORIG_FIG_SAVEFIG = Figure.savefig
def _selective_fig_savefig(self, fname=None, *a, **k):
    if fname is not None and _keep(fname):
        return _ORIG_FIG_SAVEFIG(self, fname, *a, **k)
    return None
Figure.savefig = _selective_fig_savefig

_ORIG_PLT_SAVEFIG = plt.savefig
def _selective_plt_savefig(fname=None, *a, **k):
    if fname is not None and _keep(fname):
        return _ORIG_PLT_SAVEFIG(fname, *a, **k)
    return None
plt.savefig = _selective_plt_savefig

import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, SCRIPTS_DIR)
import hiv_rapid_art_model as M
if RUN_CEA:
    import cea_analysis as CEA
if RUN_BREAKEVEN:
    import breakeven_analysis_detailed_json as BE

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
def _read_csv_robust(path):
    try:
        return pd.read_csv(path)
    except Exception:
        with open(path, "rb") as f:
            data = f.read().replace(b"\x00", b"")
        return pd.read_csv(StringIO(data.decode("utf-8", "ignore")))

def _read_burnin_end(seed_dir):
    df = _read_csv_robust(os.path.join(seed_dir, "burnin_summary.csv"))
    r = df.iloc[-1]
    hiv = float(r["hiv_positive"]); onart = float(r["on_art"])
    return {
        "HIVpos_end":  int(hiv),
        "onART_end":   int(onart),
        "undiag_pct":  round(100.0 * float(r["undiagnosed_hiv"]) / hiv, 1) if hiv else float("nan"),
        "VLsupp_pct":  round(100.0 * float(r["vl_suppressed"]) / onart, 1) if onart else float("nan"),
        "newInf_yr":   int(r["new_infections_year"]),
        "newDiag_yr":  int(r["new_diagnoses_year"]),
    }

def _read_infections_and_het(seed_dir):
    df = _read_csv_robust(os.path.join(seed_dir, "heterogeneity_summary.csv"))
    out = {}
    label_map = {"0 month": "0m", "1 month": "1m", "2 months": "2m"}
    for sc_label, short in label_map.items():
        sub = df[df["Scenario"] == sc_label]
        total = float(sub["Cum_New_Infections"].sum())
        het = float(sub[sub["Risk_Group"] == "Heterosexual"]["Cum_New_Infections"].sum())
        out[f"cumInf_{short}"] = int(total)
        out[f"hetPct_{short}"] = round(100.0 * het / total, 1) if total else float("nan")
    return out

def _run_cea_and_read_qaly(seed_dir):
    CEA.Config.INPUT_DIR  = seed_dir
    CEA.Config.OUTPUT_DIR = os.path.join(seed_dir, "cea")
    CEA.Config.MAKE_FIGURES = False
    if CEA_PSA_ITERATIONS is not None:
        CEA.Config.PSA_N_ITERATIONS = CEA_PSA_ITERATIONS
    os.makedirs(CEA.Config.OUTPUT_DIR, exist_ok=True)
    CEA.main()
    with open(os.path.join(CEA.Config.OUTPUT_DIR, "cea_report_ready.json"), encoding="utf-8") as f:
        rep = json.load(f)
    bc = rep["base_case_results"]
    return {f"QALY_{s}": round(float(bc[s]["total_qaly_discounted"]), 1) for s in ("0m", "1m", "2m")}

def _run_breakeven_and_read(seed_dir):
    BE.INPUT_DIR  = seed_dir
    BE.OUTPUT_DIR = os.path.join(seed_dir, "breakeven")
    os.makedirs(BE.OUTPUT_DIR, exist_ok=True)
    BE.main()
    with open(os.path.join(BE.OUTPUT_DIR, "breakeven_report_ready.json"), encoding="utf-8") as f:
        rep = json.load(f)
    comps = rep["breakeven_comparisons"]
    return {
        "BE_0m_vs_1m": comps.get("0m vs 1m", {}).get("break_even_year", -1),
        "BE_0m_vs_2m": comps.get("0m vs 2m", {}).get("break_even_year", -1),
    }

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
def run_one_seed(seed):
    seed_dir = os.path.join(MULTISEED_DIR, f"seed_{seed}")
    os.makedirs(seed_dir, exist_ok=True)

    M.SEED = seed
    M.OUTPUT_DIR = seed_dir
    burnin_pop, _ = M.run_burnin(seed=seed)
    scenarios = M.get_enabled_scenarios()                # ["0 month","1 month","2 months"]
    all_results, all_pops = {}, {}
    if PARALLEL_SCENARIOS and len(scenarios) >= 2:

        worker_args = [(sc, burnin_pop, seed) for sc in scenarios]
        with ProcessPoolExecutor(max_workers=len(scenarios)) as ex:
            for sc, res, pop in ex.map(M._run_scenario_worker, worker_args):
                all_results[sc] = res
                all_pops[sc] = pop
    else:
        for sc in scenarios:
            res, pop = M.run_simulation(label=sc, initial_pop=burnin_pop, seed=seed)
            all_results[sc] = res
            all_pops[sc] = pop

    M.save_all_outputs(all_results, all_pops, scenarios)
    plt.close("all")

    row = {"seed": seed}
    row.update(_read_burnin_end(seed_dir))
    row.update(_read_infections_and_het(seed_dir))

    del burnin_pop, all_results, all_pops
    gc.collect()

    if RUN_CEA:
        row.update(_run_cea_and_read_qaly(seed_dir))
    if RUN_BREAKEVEN:
        row.update(_run_breakeven_and_read(seed_dir))

    row["INF_MONO_PASS"] = bool(row["cumInf_0m"] < row["cumInf_1m"] < row["cumInf_2m"])
    if RUN_CEA:
        row["QALY_ORDER_PASS"] = bool(row["QALY_0m"] > row["QALY_1m"] > row["QALY_2m"])
    if RUN_BREAKEVEN:
        lo, hi = BREAKEVEN_OK_RANGE
        row["BE_IN_RANGE_PASS"] = bool(all(
            (row.get(k, -1) is not None and row.get(k, -1) > 0 and lo <= row.get(k, -1) <= hi)
            for k in ("BE_0m_vs_1m", "BE_0m_vs_2m")
        ))
    flags = [row.get("INF_MONO_PASS", False)]
    if RUN_CEA:       flags.append(row.get("QALY_ORDER_PASS", False))
    if RUN_BREAKEVEN: flags.append(row.get("BE_IN_RANGE_PASS", False))
    row["ALL_PASS"] = bool(all(flags))
    return row

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
def main():
    os.makedirs(MULTISEED_DIR, exist_ok=True)
    print("=" * 70)
    print(f"Multi-seed runner: {len(SEEDS)} seeds")
    print(f"Output dir : {MULTISEED_DIR}")
    print(f"Seeds      : {SEEDS}")
    print(f"RUN_CEA={RUN_CEA}  RUN_BREAKEVEN={RUN_BREAKEVEN}  SKIP_PLOTS={SKIP_PLOTS}")
    print("=" * 70, flush=True)

    rows = []
    t0 = time.time()
    for i, seed in enumerate(SEEDS, 1):
        ts = time.time()
        print(f"\n[{i}/{len(SEEDS)}] SEED={seed}  start {datetime.datetime.now():%H:%M:%S}", flush=True)
        try:
            row = run_one_seed(seed)
            rows.append(row)
            print(f"   done in {(time.time()-ts)/60:.1f} min | "
                  f"HIV+={row.get('HIVpos_end')}, onART={row.get('onART_end')}, "
                  f"INF_MONO={row.get('INF_MONO_PASS')}, QALY_ORDER={row.get('QALY_ORDER_PASS')}, "
                  f"BE={row.get('BE_0m_vs_1m')}/{row.get('BE_0m_vs_2m')}, ALL={row.get('ALL_PASS')}",
                  flush=True)
        except Exception as e:
            import traceback; traceback.print_exc()
            print(f"   [ERROR] seed {seed} failed: {e}", flush=True)
            rows.append({"seed": seed, "error": str(e)})

        pd.DataFrame(rows).to_csv(OUTPUT_CSV, index=False)

    df = pd.DataFrame(rows)
    ok = df[df["error"].isna()] if "error" in df.columns else df
    flag_cols = [c for c in ("INF_MONO_PASS", "QALY_ORDER_PASS", "BE_IN_RANGE_PASS", "ALL_PASS") if c in df.columns]
    num_cols = [c for c in df.columns
                if c not in ("seed", "error") and c not in flag_cols
                and pd.api.types.is_numeric_dtype(df[c])]

    mean_row = {"seed": "MEAN"}
    range_row = {"seed": "MIN-MAX"}
    for c in num_cols:
        try:
            mean_row[c] = round(float(ok[c].mean()), 1)
            range_row[c] = f"{ok[c].min()}-{ok[c].max()}"
        except Exception:
            pass
    pass_row = {"seed": "PASS(n/total)"}
    for c in flag_cols:
        try:
            pass_row[c] = f"{int(ok[c].sum())}/{len(ok)}"
        except Exception:
            pass

    df_out = pd.concat([df, pd.DataFrame([mean_row, range_row, pass_row])], ignore_index=True)
    df_out.to_csv(OUTPUT_CSV, index=False)

    print("\n" + "=" * 70)
    print(f"TOTAL: {(time.time()-t0)/60:.1f} min")
    print(f"Saved: {OUTPUT_CSV}")
    for c in flag_cols:
        try:
            print(f"  {c}: {int(ok[c].sum())}/{len(ok)} seeds passed")
        except Exception:
            pass
    print("=" * 70, flush=True)

if __name__ == "__main__":
    main()

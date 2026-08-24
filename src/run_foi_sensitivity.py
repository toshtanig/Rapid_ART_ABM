#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Driver for the analytic-period force-of-infection multiplier sensitivity."""

import os
import sys
import gc
import json
import time
import datetime

# =============================================================================

# =============================================================================
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# --- Data location: set BASE_144 to your local 144_scenario folder ---
BASE_144 = r"/path/to/20250315_169/144_scenario"
PROJ_169 = os.path.dirname(BASE_144)
BASE_169 = PROJ_169
SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_BASE    = os.path.join(BASE_169, "144_scenario")

ANALYTIC_MULTIPLIERS = [1.00, 1.10, 1.18, 1.25, 1.30]

FOI_TERMINAL = "Raptorlake"   # "Raptorlake" / "Blackwell" / "MSI" / "Mouse" / "Dell"
FOI_SEED_SETS = {
    "Raptorlake": [   # 27 seeds (22min) ~49.5h
        14, 18, 24, 42, 55, 58, 61, 66, 69, 71,
        75, 77, 87, 97, 101, 102, 104, 107, 110, 119,
        123, 128, 145, 150, 157, 161, 176,
    ],
    "Blackwell": [    # 23 seeds (25min) ~47.9h
        178, 188, 194, 202, 203, 206, 207, 209, 216, 219,
        225, 228, 234, 244, 247, 249, 276, 286, 303, 306,
        308, 309, 313,
    ],
    "MSI": [          # 21 seeds (28min) ~49.0h
        317, 327, 330, 340, 341, 345, 347, 355, 356, 363,
        364, 380, 389, 404, 417, 424, 427, 441, 448, 456,
        465,
    ],
    "Mouse": [        # 15 seeds (38min) ~47.5h
        466, 469, 472, 477, 484, 486, 488, 492, 495, 504,
        505, 508, 515, 516, 527,
    ],
    "Dell": [         # 14 seeds (42min) ~49.0h
        536, 543, 547, 567, 581, 582, 585, 586, 590, 592,
        596, 606, 608, 613,
    ],
}
if FOI_TERMINAL not in FOI_SEED_SETS:
    raise ValueError("FOI_TERMINAL must be one of: " + ", ".join(FOI_SEED_SETS))
SEEDS = FOI_SEED_SETS[FOI_TERMINAL]

SEED_ORDER = "forward"
if SEED_ORDER == "reverse":
    SEEDS = list(reversed(SEEDS))
elif SEED_ORDER == "shuffle":
    import random as _random
    SEEDS = list(SEEDS)
    _random.Random(12345).shuffle(SEEDS)
elif SEED_ORDER != "forward":
    raise ValueError('SEED_ORDER must be "forward", "reverse", or "shuffle"')

MULT_ORDER = "forward"
if MULT_ORDER == "reverse":
    ANALYTIC_MULTIPLIERS = list(reversed(ANALYTIC_MULTIPLIERS))
elif MULT_ORDER != "forward":
    raise ValueError('MULT_ORDER must be "forward" or "reverse"')
# ==========================================================================

RUN_CEA       = True
RUN_BREAKEVEN = True
PARALLEL_SCENARIOS = True
SKIP_PLOTS    = True
KEEP_PLOT_FILES = ["summary_plots.jpg"]
CEA_PSA_ITERATIONS = 200

def out_dir_for(mult):
    return os.path.join(OUT_BASE, f"core_foi{int(round(mult*100))}")
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

def _run_cea_and_read_qaly(seed_dir):
    CEA.Config.INPUT_DIR  = seed_dir
    CEA.Config.OUTPUT_DIR = os.path.join(seed_dir, "cea")
    CEA.Config.MAKE_FIGURES = False
    if CEA_PSA_ITERATIONS is not None:
        CEA.Config.PSA_N_ITERATIONS = CEA_PSA_ITERATIONS
    os.makedirs(CEA.Config.OUTPUT_DIR, exist_ok=True)
    CEA.main()

def _run_breakeven(seed_dir):
    BE.INPUT_DIR  = seed_dir
    BE.OUTPUT_DIR = os.path.join(seed_dir, "breakeven")
    os.makedirs(BE.OUTPUT_DIR, exist_ok=True)
    BE.main()

def run_one_seed(seed, multiseed_dir):
    seed_dir = os.path.join(multiseed_dir, f"seed_{seed}")
    os.makedirs(seed_dir, exist_ok=True)

    M.SEED = seed
    M.OUTPUT_DIR = seed_dir
    burnin_pop, _ = M.run_burnin(seed=seed)
    scenarios = M.get_enabled_scenarios()
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

    del burnin_pop, all_results, all_pops
    gc.collect()

    if RUN_CEA:
        _run_cea_and_read_qaly(seed_dir)
    if RUN_BREAKEVEN:
        _run_breakeven(seed_dir)

def _seed_done(seed_dir):
    if not os.path.exists(os.path.join(seed_dir, "heterogeneity_summary.csv")):
        return False
    if RUN_CEA and not os.path.exists(os.path.join(seed_dir, "cea", "cea_report_ready.json")):
        return False
    if RUN_BREAKEVEN and not os.path.exists(os.path.join(seed_dir, "breakeven", "breakeven_report_ready.json")):
        return False
    return True

def main():
    print("=" * 70)
    print(f"FOI sensitivity [{FOI_TERMINAL}]: multipliers={ANALYTIC_MULTIPLIERS}  seeds={len(SEEDS)}")
    print("=" * 70, flush=True)
    t0 = time.time()
    for mult in ANALYTIC_MULTIPLIERS:

        M.MAIN_TRANSMISSION_MULTIPLIER = float(mult)
        mdir = out_dir_for(mult)
        os.makedirs(mdir, exist_ok=True)
        print(f"\n### ANALYTIC FOI MULTIPLIER = {mult}  ->  {mdir}", flush=True)
        for i, seed in enumerate(SEEDS, 1):
            seed_dir = os.path.join(mdir, f"seed_{seed}")
            if _seed_done(seed_dir):
                print(f"  [{i}/{len(SEEDS)}] mult={mult} SEED={seed} skip (already done)", flush=True)
                continue
            ts = time.time()
            print(f"  [{i}/{len(SEEDS)}] mult={mult} SEED={seed} start {datetime.datetime.now():%H:%M:%S}", flush=True)
            try:
                run_one_seed(seed, mdir)
                print(f"     done in {(time.time()-ts)/60:.1f} min", flush=True)
            except Exception as e:
                import traceback; traceback.print_exc()
                print(f"     [ERROR] mult={mult} seed={seed}: {e}", flush=True)
    print("\n" + "=" * 70)
    print(f"TOTAL: {(time.time()-t0)/60:.1f} min")
    print("Next: run foi_source_analysis.py to build the sensitivity and source tables.")
    print("=" * 70, flush=True)

if __name__ == "__main__":
    main()

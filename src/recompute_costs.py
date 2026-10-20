#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Recompute cost summaries from per-seed outputs."""

import os
import sys

# ============================ CONFIG ============================

import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# --- Data location: set BASE_144 to your local 144_scenario folder ---
BASE_144 = r"/path/to/project/144_scenario"
PROJ_169 = os.path.dirname(BASE_144)
_PROJ_144 = BASE_144

# Parent folder(s) that contain the per-seed output folders (seed_1, seed_2, ...).
# List every multi-seed folder you used (e.g. one per machine).
SCAN_DIRS = [
    os.path.join(_PROJ_144, "core_multiseed"),
]

RUN_CEA            = True    # recompute per-seed cost-effectiveness outputs
RUN_BREAKEVEN      = True    # recompute per-seed break-even outputs (needed by prevention_qaly_dominance.py)
CEA_PSA_ITERATIONS = 200     # per-seed PSA iterations; None keeps the engine default
# ===============================================================

# Disable figure output during the sweep (faster, less disk); CSV/JSON are unaffected.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
Figure.savefig = lambda self, *a, **k: None
plt.savefig = lambda *a, **k: None

# Import the economic engines from this script's directory.
SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)
if RUN_CEA:
    import cea_analysis as CEA
if RUN_BREAKEVEN:
    import breakeven_analysis_detailed_json as BE

def iter_seed_dirs():
    import glob
    for base in SCAN_DIRS:
        for d in sorted(glob.glob(os.path.join(base, "seed_*"))):
            if os.path.isdir(d):
                yield d

def recompute_one_seed(seed_dir):
    """Recompute the cost-effectiveness and break-even outputs for one existing seed folder."""
    if RUN_CEA:
        CEA.Config.INPUT_DIR = seed_dir
        CEA.Config.OUTPUT_DIR = os.path.join(seed_dir, "cea")
        CEA.Config.MAKE_FIGURES = False
        if CEA_PSA_ITERATIONS is not None:
            CEA.Config.PSA_N_ITERATIONS = CEA_PSA_ITERATIONS
        os.makedirs(CEA.Config.OUTPUT_DIR, exist_ok=True)
        CEA.main()

    if RUN_BREAKEVEN:
        BE.INPUT_DIR = seed_dir
        BE.OUTPUT_DIR = os.path.join(seed_dir, "breakeven")
        os.makedirs(BE.OUTPUT_DIR, exist_ok=True)
        BE.main()

def main():
    seed_dirs = list(iter_seed_dirs())
    if not seed_dirs:
        print("No seed_* folders found. Check SCAN_DIRS in the CONFIG block.")
        return
    print(f"Recomputing economic layer for {len(seed_dirs)} seed folders "
          f"(microsimulation is NOT re-run).")
    for i, seed_dir in enumerate(seed_dirs, 1):
        print(f"[{i}/{len(seed_dirs)}] {seed_dir}", flush=True)
        try:
            recompute_one_seed(seed_dir)
            plt.close("all")
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"  [ERROR] {seed_dir}: {e}", flush=True)
    print("Done. Next, run: prevention_qaly_dominance.py, cea_param_uncertainty.py, "
          "survival_adj_sensitivity.py")

if __name__ == "__main__":
    main()

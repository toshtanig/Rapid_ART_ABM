#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Probabilistic and one-way parameter uncertainty (PSA, DSA, CEAC)."""

import os
import glob
import matplotlib
matplotlib.use("Agg")
import pandas as pd

# --- Lancet TIFF export hook (300 dpi, LZW; multipanel 180 mm / single 107 mm) ---
EXPORT_TIFF = True
if EXPORT_TIFF:
    try:
        import sys as _sys
        _sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import lancet_export
        lancet_export.enable_tiff_export()
    except Exception as _e:
        print("[lancet_export] disabled:", _e)

# =============================================================================
# CONFIG
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
OUTPUT_DIR     = os.path.join(BASE_169, "144_scenario", "cea_param_uncertainty")
MEAN_INPUT_DIR = os.path.join(OUTPUT_DIR, "mean_input")
RUN_BREAKEVEN  = True
# =============================================================================

SCEN_FILES = ["scenario0_yearly_summary.csv", "scenario1_yearly_summary.csv", "scenario2_yearly_summary.csv"]
GROUP_KEYS = ["Year", "AgeGroup", "CD4Category", "OnART", "OnPrEP"]
NUM_COLS   = ["PersonMonthsAlive", "NewlyInfected", "NewlyDiagnosed", "Deaths"]

def _read_csv(path):
    try:
        return pd.read_csv(path)
    except Exception:
        with open(path, "rb") as f:
            data = f.read().replace(b"\x00", b"")
        from io import StringIO
        return pd.read_csv(StringIO(data.decode("utf-8", "ignore")))

def build_mean_input():
    os.makedirs(MEAN_INPUT_DIR, exist_ok=True)
    seed_dirs = []
    for base in SCAN_DIRS:
        seed_dirs += sorted(glob.glob(os.path.join(base, "seed_*")))

    n_used = 0
    for scen_file in SCEN_FILES:
        frames = []
        for d in seed_dirs:
            p = os.path.join(d, scen_file)
            if os.path.exists(p):
                frames.append(_read_csv(p))
        if not frames:
            print(f"  [skip] no seeds have {scen_file}")
            continue
        n = len(frames)
        n_used = max(n_used, n)
        allc = pd.concat(frames, ignore_index=True)

        summed = allc.groupby(GROUP_KEYS, as_index=False)[NUM_COLS].sum()
        summed[NUM_COLS] = summed[NUM_COLS] / n
        summed.to_csv(os.path.join(MEAN_INPUT_DIR, scen_file), index=False)
        print(f"  averaged {scen_file} over {n} seeds → {len(summed)} rows")
    return n_used

def run_cea_on_mean():
    import cea_analysis as CEA
    CEA.Config.INPUT_DIR = MEAN_INPUT_DIR
    CEA.Config.OUTPUT_DIR = OUTPUT_DIR
    CEA.Config.MAKE_FIGURES = True
    CEA.main()

def run_breakeven_on_mean():
    import breakeven_analysis_detailed_json as BE
    BE.INPUT_DIR = MEAN_INPUT_DIR
    BE.OUTPUT_DIR = OUTPUT_DIR
    os.makedirs(BE.OUTPUT_DIR, exist_ok=True)
    BE.main()

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print("=" * 70)
    print("CEA parameter-uncertainty (DSA Tornado / PSA CEAC) on seed-averaged input")
    print(f"output: {OUTPUT_DIR}")
    print("=" * 70)
    print("[1] building mean input ...")
    n = build_mean_input()
    if n == 0:
        print("No seed scenario summaries found."); return
    print(f"[2] running cea_analysis_rapid_v3 on mean input (n={n} seeds averaged) ...")
    run_cea_on_mean()
    if RUN_BREAKEVEN:
        print("[3] running breakeven on mean input ...")
        run_breakeven_on_mean()
    print("=" * 70)
    print(f"Done. Tornado/CEAC/budget figures & CSVs in: {OUTPUT_DIR}")
    print("Note: this is the 'parameter uncertainty' layer; for 'stochastic (seed) uncertainty' see prevention_qaly_dominance.py.")
    print("=" * 70)

if __name__ == "__main__":
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Calibration check for the timing convention of ART starts (burn-in only).

Runs the 20-year burn-in (calendar years 2004 to 2023) for a few seeds on this terminal under each setting in
SETTINGS and writes burnin_summary.csv per seed to

    outputs/revision_R2_runs/C0_calibration_check/<setting>/seed_<n>/

No analytic-period simulation and no cost-effectiveness post-processing are run, so one seed takes a few minutes.
The seeds are the first seeds of each terminal's list within the 200 publication seeds, so the burn-in of
versions 1.0.0 and 1.1.0 for the same seeds is available in outputs/core_multiseed/seed_<n>/burnin_summary.csv;
summarize_calibration_check.py compares the two and reports the calibration indicators and fit statistics.

Run as-is in VS Code: all choices are in the CONFIG block below, no command-line arguments. Set TERMINAL, then run.
Seeds already completed are skipped, so the script can be restarted safely.
"""

import os
import sys
import gc
import json
import time
import datetime

# =============================================================================
# CONFIG
# =============================================================================
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import REVISION_RUNS

TERMINAL = "node1"          # "node1" / "node2" / "node3" / "node4" / "node5"
SEEDS_PER_TERMINAL = 4           # seeds taken from the head of this terminal's list (20 seeds over five terminals)
SEED_SETS = {                    # head of each terminal's list in run_revision_scenarios.py (publication seeds)
    "node1": [14, 18, 24, 42, 55, 58, 61, 66],
    "node2":  [178, 188, 194, 202, 203, 206, 207, 209],
    "node3":        [317, 327, 330, 340, 341, 345, 347, 355],
    "node4":      [466, 469, 472, 477, 484, 486, 488, 492],
    "node5":       [536, 543, 547, 567, 581, 582, 585, 586],
}

# Each setting: model-global overrides applied on top of the defaults of hiv_rapid_art_model.py in this folder.
# Add a setting here to test a candidate recalibration (for example a different MORTALITY_RECALIB_FOI_SCALE).
SETTINGS = {
    "uniform_convention": {
        "model_overrides": {},
        "note": "Defaults of this folder: for every ART start the first viral-load update is in the cycle "
                "after the start cycle, in the burn-in as well as in the analytic period.",
    },
    # Candidate recalibrations (2026-11-17). With the uniform convention and the calibration scalar of the
    # earlier versions (1.05), the 2023 stock indicators exceed their acceptance ranges (people on ART 31,311
    # against 28,500 to 30,500), so the scalar is re-derived. It multiplies the force of infection in the
    # burn-in and in the analytic period alike. The settings run in the order listed.
    "uniform_foi_scale_1p02": {
        "model_overrides": {"MORTALITY_RECALIB_FOI_SCALE": 1.02},
        "note": "Uniform convention with the calibration scalar set to 1.02.",
    },
    "uniform_foi_scale_1p03": {
        "model_overrides": {"MORTALITY_RECALIB_FOI_SCALE": 1.03},
        "note": "Uniform convention with the calibration scalar set to 1.03.",
    },
    "uniform_foi_scale_1p01": {
        "model_overrides": {"MORTALITY_RECALIB_FOI_SCALE": 1.01},
        "note": "Uniform convention with the calibration scalar set to 1.01.",
    },
    # Recalibration for the sensitivity analysis of reduced transmission after diagnosis (2026-11-17). The
    # relative infectivity of diagnosed, unsuppressed people applies in the burn-in as well as in the
    # analytic period, so for each value the calibration scalar is re-derived to reproduce the 2023 targets.
    # Two scalars per value bracket the expected solution (interpolation on the logarithmic scale).
    "factor_0p5_foi_scale_1p12": {
        "model_overrides": {"DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR": 0.5, "MORTALITY_RECALIB_FOI_SCALE": 1.12},
        "note": "Relative infectivity of diagnosed, unsuppressed people 0.5; calibration scalar 1.12.",
    },
    "factor_0p5_foi_scale_1p16": {
        "model_overrides": {"DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR": 0.5, "MORTALITY_RECALIB_FOI_SCALE": 1.16},
        "note": "Relative infectivity of diagnosed, unsuppressed people 0.5; calibration scalar 1.16.",
    },
    "factor_0p3_foi_scale_1p16": {
        "model_overrides": {"DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR": 0.3, "MORTALITY_RECALIB_FOI_SCALE": 1.16},
        "note": "Relative infectivity of diagnosed, unsuppressed people 0.3; calibration scalar 1.16.",
    },
    "factor_0p3_foi_scale_1p21": {
        "model_overrides": {"DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR": 0.3, "MORTALITY_RECALIB_FOI_SCALE": 1.21},
        "note": "Relative infectivity of diagnosed, unsuppressed people 0.3; calibration scalar 1.21.",
    },
}
SETTINGS_TO_RUN = None           # None = all settings, or a list of setting names

OUT_ROOT = os.path.join(REVISION_RUNS, "C0_calibration_check")
# =============================================================================

os.environ["OUTPUT_BASE_DIR"] = os.path.join(OUT_ROOT, "_model_default_output")
os.makedirs(os.environ["OUTPUT_BASE_DIR"], exist_ok=True)

import matplotlib
matplotlib.use("Agg")

import hiv_rapid_art_model as M

_REPORTED_KEYS = ["ART_START_VL_UPDATE_NEXT_MONTH", "ART_START_VL_UPDATE_NEXT_MONTH_IN_BURNIN",
                  "MORTALITY_RECALIB_FOI_SCALE", "RECENT_FOI_BOOST", "RECENT_FOI_BOOST_FROM_YEAR",
                  "MAIN_TRANSMISSION_MULTIPLIER", "AGE_MORTALITY_OFFSET",
                  "DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR"]
_BASE_MODEL_VALUES = {}


def apply_setting(setting):
    """Reset the overridden model globals to their defaults, then apply this setting's overrides."""
    for k, v in _BASE_MODEL_VALUES.items():
        setattr(M, k, v)
    for k, v in setting.get("model_overrides", {}).items():
        if k not in _BASE_MODEL_VALUES:
            _BASE_MODEL_VALUES[k] = getattr(M, k)
        setattr(M, k, v)


def run_one_seed(seed, setting_name, setting, setting_dir):
    seed_dir = os.path.join(setting_dir, f"seed_{seed}")
    os.makedirs(seed_dir, exist_ok=True)
    apply_setting(setting)
    M.SEED = seed
    M.OUTPUT_DIR = seed_dir
    t0 = time.time()
    pop, _ = M.run_burnin(seed=seed)
    elapsed = time.time() - t0
    del pop
    gc.collect()
    info = {
        "setting": setting_name,
        "seed": int(seed),
        "terminal": TERMINAL,
        "finished": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_seconds": round(elapsed, 1),
        "model_values": {k: getattr(M, k) for k in _REPORTED_KEYS if hasattr(M, k)},
        "burnin_transmission_multiplier_by_year": {str(k): v for k, v in M.BURNIN_TRANSMISSION_MULTIPLIER_BY_YEAR.items()},
        "model_overrides": setting.get("model_overrides", {}),
        "note": setting.get("note", ""),
    }
    with open(os.path.join(seed_dir, "calibration_check_config.json"), "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2, ensure_ascii=False)
    return elapsed


def main():
    if TERMINAL not in SEED_SETS:
        raise SystemExit(f"TERMINAL must be one of {list(SEED_SETS)}; got {TERMINAL!r}")
    seeds = SEED_SETS[TERMINAL][:SEEDS_PER_TERMINAL]
    names = list(SETTINGS) if SETTINGS_TO_RUN is None else list(SETTINGS_TO_RUN)
    print(f"Calibration check (burn-in only) on {TERMINAL}: settings {names}, seeds {seeds}")
    print(f"Output: {OUT_ROOT}")
    for name in names:
        setting = SETTINGS[name]
        setting_dir = os.path.join(OUT_ROOT, name)
        os.makedirs(setting_dir, exist_ok=True)
        for i, seed in enumerate(seeds, 1):
            done = os.path.join(setting_dir, f"seed_{seed}", "burnin_summary.csv")
            if os.path.exists(done) and os.path.exists(os.path.join(setting_dir, f"seed_{seed}", "calibration_check_config.json")):
                print(f"[{name}] seed {seed} ({i}/{len(seeds)}): already done, skipped")
                continue
            print(f"[{name}] seed {seed} ({i}/{len(seeds)}): running burn-in ...", flush=True)
            el = run_one_seed(seed, name, setting, setting_dir)
            print(f"[{name}] seed {seed}: finished in {el / 60:.1f} min", flush=True)
    apply_setting({})
    print("All requested burn-in runs are complete. Next: run summarize_calibration_check.py on any one terminal.")


if __name__ == "__main__":
    main()

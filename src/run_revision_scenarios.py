#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Driver for the main and supplementary simulations (model version 1.2.0).

Runs the model (hiv_rapid_art_model.py in this folder) for a queue of
(analysis, setting, number of seeds) items, from the top of the queue to the
bottom. For every seed the burn-in is run once and the scenario slots of the
setting are run from that burn-in (common random numbers), followed by the
cost-effectiveness (cea_analysis) and break-even
(breakeven_analysis_detailed_json) post-processing. Output layout:

    outputs/revision_R2_runs/<analysis>/<setting>/seed_<n>/...

Several terminals share one queue through the synchronised output folder:

  * every seed has an owner terminal (shares proportional to the speed of the
    terminals), and a terminal first runs its own seeds of a queue item;
  * a terminal that has finished its own seeds of an item takes seeds of that
    item that no other terminal has started, from the end of the list of the
    terminal with the most seeds left (HELP_OTHER_TERMINALS);
  * a terminal claims a seed by writing a small file into its folder before
    it starts, so that other terminals leave the seed alone; a claim older
    than CLAIM_TIMEOUT_MINUTES on a seed that is not done is treated as
    abandoned;
  * completed seeds are skipped, so the script can be stopped and restarted
    at any time, on any terminal.

Run as-is in VS Code: all choices are in the CONFIG block below, no
command-line arguments. Set TERMINAL, then run.

Analyses:
  B0_core                     main analysis: immediate ART, one-month and two-month delay
  A1_factorial                eligibility criteria versus waiting interval
  A2_foi_trajectory           time-varying analytic-period force of infection
  A4_time_to_suppression      faster / slower time to viral suppression
  A5_post_diagnosis_behaviour reduced infectivity of diagnosed, unsuppressed people
  A7_delayed_start            immediate ART implemented from 2027 / 2030
  A8_status_quo_empirical     comparator with the observed waiting-time distribution
  A9_observed_ineligible_wait comparators with the observed wait of people not meeting the criteria
  A10_conservative_combined   observed-practice comparator combined with reduced infectivity after diagnosis
  A12_age_offset              background mortality read at the actual age
  A13_foi_multiplier          analytic-period force-of-infection multiplier
  A0_regression_check         earlier timing convention and calibration scalar (must reproduce version 1.1.0)
"""

import os
import sys
import gc
import glob
import json
import time
import atexit
import socket
import signal
import hashlib
import datetime
import traceback
import faulthandler
import concurrent.futures

# =============================================================================
# CONFIG
# =============================================================================
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import REVISION_RUNS

RUN = "main"                     # "main" = MAIN_QUEUE; "regression" = REGRESSION_QUEUE (one seed, on one terminal,
                                 # before the main queue; then run compare_regression.py)
TERMINAL = ""                    # set to the label of this machine: "node1" / "node2" / "node3" / "node4" / "node5"

# Queue of (analysis, setting, number of seeds), run from top to bottom. The seeds of an item are the first
# n of PUBLICATION_SEEDS. An item may be listed again further down with a larger number of seeds.
MAIN_QUEUE = [
    ("B0_core", "default", 200),
    ("A9_observed_ineligible_wait", "empirical_pathways_plus_observed_tail", 200),
    ("A9_observed_ineligible_wait", "criteria_retained_no_wait_plus_observed_tail", 100),
    ("A5_post_diagnosis_behaviour", "factor_0.5", 100),
    ("A5_post_diagnosis_behaviour", "factor_0.3", 100),
    ("A10_conservative_combined", "empirical_tail_behaviour_0.3", 100),
    ("A10_conservative_combined", "empirical_tail_behaviour_0.5", 50),
    ("A1_factorial", "criteria_retained_no_wait", 100),
    ("A1_factorial", "criteria_removed_one_month", 100),
    ("A9_observed_ineligible_wait", "one_month_plus_observed_tail", 50),
    ("A8_status_quo_empirical", "empirical_pathways", 50),
    ("A7_delayed_start", "start_2027", 50),
    ("A7_delayed_start", "start_2030", 50),
    ("A2_foi_trajectory", "decline", 50),
    ("A2_foi_trajectory", "resurgence", 50),
    ("A4_time_to_suppression", "fast", 30),
    ("A4_time_to_suppression", "slow", 30),
    ("A12_age_offset", "age_offset_0", 30),
    ("A13_foi_multiplier", "multiplier_1.00", 50),
    ("A13_foi_multiplier", "multiplier_1.10", 50),
    ("A13_foi_multiplier", "multiplier_1.25", 50),
    ("A13_foi_multiplier", "multiplier_1.30", 50),
    ("A13_foi_multiplier", "multiplier_1.00", 100),
    ("A13_foi_multiplier", "multiplier_1.10", 100),
    ("A13_foi_multiplier", "multiplier_1.25", 100),
    ("A13_foi_multiplier", "multiplier_1.30", 100),
]
REGRESSION_QUEUE = [
    ("A0_regression_check", "earlier_convention", 1),
]
REGRESSION_SEED = 14             # compared with outputs/core_multiseed/seed_14 by compare_regression.py

PUBLICATION_SEEDS = [
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
]   # the 200 publication seeds, in their fixed order

# Minutes per seed-run of each terminal (measured in the earlier runs). They set the share of the seeds
# that each terminal owns; the terminals also balance the work among themselves while running (see below),
# so the values need not be exact. Do not change them while a run is in progress.
MINUTES_PER_SEED_RUN = {"node1": 17, "node2": 20, "node3": 10, "node4": 29, "node5": 23}

HELP_OTHER_TERMINALS = True      # after its own seeds of a queue item, take unstarted seeds of other terminals
CLAIM_TIMEOUT_MINUTES = 180      # a seed started by another terminal and not finished within this time is abandoned
WAIT_FOR_OTHERS = True           # at the end of the queue, wait for seeds still running elsewhere (and take
                                 # them over if they are abandoned) instead of stopping
WAIT_POLL_MINUTES = 10
MAX_SEED_RUNS = None             # stop after this many seed-runs on this terminal (for example 1 for a pilot)

# Largest number of scenario slots run at the same time on each terminal. Three slots at once need about
# one and a half times the memory of two. A terminal on which that exceeds the available memory runs the
# slots of a three-scenario setting two at a time and then the third; the results are identical.
MAX_PARALLEL_SCENARIOS = {"node1": 3, "node2": 3, "node3": 2, "node4": 3, "node5": 2}
SCENARIO_TIMEOUT_MINUTES = 75    # a group of scenario slots still running after this time is stopped (error)
PROGRESS_NOTE_MINUTES = 5        # while scenario slots run, print a line this often
MAX_CONSECUTIVE_ERRORS = 2       # stop this terminal when this many seeds in a row end with an error

# Records for the diagnosis of stops (they do not affect the simulation). A crash of the interpreter
# (segmentation fault, abort) writes the Python stack to _logs/fault_<TERMINAL>.log. The free memory of this
# machine is read while a seed runs: _logs/memory_<TERMINAL>.log gets a line every MEMORY_LOG_MINUTES, and
# run_info.json of each seed records the lowest free memory and the largest memory use of the processes.
MEMORY_SAMPLE_SECONDS = 15       # while scenario slots run, the free memory is read this often
MEMORY_LOG_MINUTES = 1           # a line in the memory log this often (and at each stage of a seed)

# A terminal name can run only once at a time. The start is refused if the name is already running on this
# machine (a second start), or was active on another machine within TERMINAL_IN_USE_MINUTES (TERMINAL not
# set to the name of this machine). Set ALLOW_TERMINAL_IN_USE = True only to take over the name of a
# machine that has been stopped.
TERMINAL_IN_USE_MINUTES = 30
ALLOW_TERMINAL_IN_USE = False

# The main queue starts only if the defaults of hiv_rapid_art_model.py in this folder equal these values
# (protection against a stale copy of the model file on one terminal). None = not yet fixed: the main
# queue does not start.
EXPECTED_MODEL_DEFAULTS = {
    "ART_START_VL_UPDATE_NEXT_MONTH": True,
    "ART_START_VL_UPDATE_NEXT_MONTH_IN_BURNIN": True,
    "MORTALITY_RECALIB_FOI_SCALE": 1.02,
    "RECENT_FOI_BOOST": 1.04,
    "MAIN_TRANSMISSION_MULTIPLIER": 1.18,
    "AGE_MORTALITY_OFFSET": 5.0,
    "DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR": 1.0,
    "TTS_MIN_MONTHS": 1.0,
    "TTS_RANGE_MONTHS": 5.0,
}

# Calibration scalar for each relative infectivity of diagnosed, unsuppressed people (analyses A5 and A10).
# The relative infectivity applies in the burn-in as well as in the analytic period, so the scalar is
# re-derived for each value to reproduce the 2023 calibration targets. None = not yet fixed: the main
# queue does not start.
FOI_SCALE_BY_INFECTIVITY_FACTOR = {0.5: 1.12, 0.3: 1.16}

RUN_CEA = True
RUN_BREAKEVEN = True
PARALLEL_SCENARIOS = True        # run the scenario slots of one seed in parallel processes
SKIP_PLOTS = True                # suppress the model's per-scenario figure files
KEEP_PLOT_FILES = ["summary_plots.jpg"]
CEA_PSA_ITERATIONS = 200

OUT_ROOT = REVISION_RUNS
LOG_DIR = os.path.join(OUT_ROOT, "_logs")
# =============================================================================


# -----------------------------------------------------------------------------
# Crash log. It is switched on here, before matplotlib and the model are imported, so that a crash of the
# interpreter while the modules load is recorded as well.
# -----------------------------------------------------------------------------
_FAULT_LOG = [None]


def _fault_note(msg):
    """A dated line in the crash log, so that a crash record (which carries no time) can be attributed."""
    f = _FAULT_LOG[0]
    if f is None:
        return
    try:
        f.write(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} [{TERMINAL}] {msg}\n")
        f.flush()
    except Exception:
        pass


def _enable_fault_log():
    """If the interpreter crashes (segmentation fault, abort, illegal instruction), the Python stack of every
    thread is written to _logs/fault_<TERMINAL>.log (to standard error while TERMINAL is not set). Scenario
    processes inherit the setting."""
    if _FAULT_LOG[0] is not None:
        return
    try:
        if not TERMINAL:
            raise ValueError("TERMINAL is not set")
        os.makedirs(LOG_DIR, exist_ok=True)
        f = open(os.path.join(LOG_DIR, f"fault_{TERMINAL}.log"), "a", encoding="utf-8")
        faulthandler.enable(file=f, all_threads=True)
        _FAULT_LOG[0] = f                       # the file must stay open for the lifetime of the process
        _fault_note(f"process {os.getpid()} on {socket.gethostname()} started and loads the modules; "
                    "a crash of the interpreter is recorded below the last dated line")
    except Exception:
        try:
            faulthandler.enable()               # standard error
        except Exception:
            pass


_enable_fault_log()


def _linear_schedule(y0, v0, y1, v1, y_end):
    """{year: value} rising/falling linearly from (y0, v0) to (y1, v1), then constant to y_end."""
    sched = {}
    for y in range(y0, y_end + 1):
        if y >= y1:
            sched[y] = float(v1)
        else:
            sched[y] = float(v0 + (v1 - v0) * (y - y0) / (y1 - y0))
    return sched


_THREE_ARMS = ["0 month", "1 month", "2 months"]
_TWO_ARMS = ["0 month", "1 month"]

# Observed-practice comparator (1 month slot): eligibility-to-ART delay matched to the multicentre cohort,
# AIDS exemption pathway, and the observed wait of people who do not meet the certification criteria.
_OBSERVED_PRACTICE = {
    "immediate": False, "apply_criteria": True,
    "delay_mean": 2.3, "delay_sd": 1.2,
    "fast_track": {"delay_mean": 0.5, "delay_sd": 0.25, "share": 0.65, "cd4_below": 200.0},
    "ineligible_fallback": {"delay_mean": 14.0, "delay_sd": 11.0},
}

# Each setting: scenario slots to run, per-slot scenario-config overrides, and model-global overrides.
# The "0 month" slot is always compared against the "1 month" slot (and the "2 months" slot, when
# present) by the downstream CEA / break-even scripts.
ANALYSES = {
    "B0_core": {
        "description": "Main analysis: immediate ART (0 month) versus the one-month and two-month delay "
                       "comparators, defaults of the model file in this folder.",
        "settings": {
            "default": {
                "scenarios": _THREE_ARMS,
                "scenario_overrides": {},
                "model_overrides": {},
                "note": "Base case (no overrides).",
            },
        },
    },
    "A0_regression_check": {
        "description": "Timing convention and calibration scalar of versions 1.0.0 and 1.1.0. The outputs "
                       "must match outputs/core_multiseed/seed_<n>/ exactly (compare_regression.py) before "
                       "the main queue is launched.",
        "settings": {
            "earlier_convention": {
                "scenarios": _THREE_ARMS,
                "scenario_overrides": {},
                "model_overrides": {"ART_START_VL_UPDATE_NEXT_MONTH": False,
                                    "MORTALITY_RECALIB_FOI_SCALE": 1.05},
                "note": "First viral-load update of a scheduled start in the start cycle "
                        "(ART_START_VL_UPDATE_NEXT_MONTH = False) and calibration scalar 1.05: the "
                        "behaviour of versions 1.0.0 and 1.1.0.",
            },
        },
    },
    "A1_factorial": {
        "description": "Factorial decomposition of the certification reform into removal of the "
                       "eligibility criteria and removal of the waiting interval.",
        "settings": {
            "criteria_retained_no_wait": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {"0 month": {"immediate": True, "apply_criteria": True}},
                "model_overrides": {},
                "note": "0 month slot = cell B (criteria retained, waiting interval removed); "
                        "1 month slot = cell A (current pathway). Contrast A - B = waiting removal given criteria.",
            },
            "criteria_removed_one_month": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {"1 month": {"immediate": False, "apply_criteria": False}},
                "model_overrides": {},
                "note": "0 month slot = cell D (full reform); 1 month slot = cell C (criteria removed, "
                        "one-month delay retained). Contrast C - D = waiting removal given no criteria.",
            },
        },
    },
    "A2_foi_trajectory": {
        "description": "Time-varying analytic-period force-of-infection multiplier "
                       "(base case: constant 1.18).",
        "settings": {
            "decline": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {},
                "model_overrides": {"ANALYTIC_FOI_MULTIPLIER_BY_YEAR": _linear_schedule(2024, 1.18, 2043, 1.00, 2063)},
                "note": "Multiplier falls linearly from 1.18 (2024) to 1.00 (2043), then constant.",
            },
            "resurgence": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {},
                "model_overrides": {"ANALYTIC_FOI_MULTIPLIER_BY_YEAR": _linear_schedule(2024, 1.18, 2033, 1.30, 2063)},
                "note": "Multiplier rises linearly from 1.18 (2024) to 1.30 (2033), then constant.",
            },
        },
    },
    "A4_time_to_suppression": {
        "description": "Time from ART initiation to viral suppression "
                       "(base case: 1 + 5 x Beta(2,5) months, median about 2.3 months).",
        "settings": {
            "fast": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {},
                "model_overrides": {"TTS_RANGE_MONTHS": 2.0},
                "note": "1 + 2 x Beta(2,5) months (median about 1.5, maximum 3).",
            },
            "slow": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {},
                "model_overrides": {"TTS_RANGE_MONTHS": 8.0},
                "note": "1 + 8 x Beta(2,5) months (median about 3.1, maximum 9).",
            },
        },
    },
    "A5_post_diagnosis_behaviour": {
        "description": "Reduced transmission after diagnosis: relative infectivity of diagnosed, "
                       "virally unsuppressed sources (base case 1.0), in the burn-in and in the analytic "
                       "period, with the calibration scalar re-derived for each value so that the 2023 "
                       "calibration targets are reproduced.",
        "settings": {
            "factor_0.5": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {},
                "model_overrides": {"DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR": 0.5,
                                    "MORTALITY_RECALIB_FOI_SCALE": FOI_SCALE_BY_INFECTIVITY_FACTOR[0.5]},
                "note": "50% reduction in transmission from diagnosed unsuppressed sources; recalibrated.",
            },
            "factor_0.3": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {},
                "model_overrides": {"DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR": 0.3,
                                    "MORTALITY_RECALIB_FOI_SCALE": FOI_SCALE_BY_INFECTIVITY_FACTOR[0.3]},
                "note": "70% reduction in transmission from diagnosed unsuppressed sources; recalibrated.",
            },
        },
    },
    "A7_delayed_start": {
        "description": "Immediate ART implemented later than 2024; before the start year the "
                       "0 month slot follows the current pathway (cost of inaction).",
        "settings": {
            "start_2027": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {"0 month": {"policy_start_year": 2027}},
                "model_overrides": {},
                "note": "Immediate ART from calendar year 2027 (three-year implementation delay).",
            },
            "start_2030": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {"0 month": {"policy_start_year": 2030}},
                "model_overrides": {},
                "note": "Immediate ART from calendar year 2030 (six-year implementation delay).",
            },
        },
    },
    "A8_status_quo_empirical": {
        "description": "Comparator with a delay distribution matched to the multicentre J-HIV RWD cohort "
                       "(11 hospitals, 2015-2024): two-test pathway median 70 days (IQR 42 to 106) after "
                       "the first visit, and an AIDS exemption pathway (about 20% of people, median 11 days).",
        "settings": {
            "empirical_pathways": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {"1 month": {
                    "immediate": False, "apply_criteria": True,
                    "delay_mean": 2.3, "delay_sd": 1.2,
                    "fast_track": {"delay_mean": 0.5, "delay_sd": 0.25, "share": 0.65, "cd4_below": 200.0},
                }},
                "model_overrides": {},
                "note": "1 month slot = empirical status quo: eligibility-to-ART delay gamma(mean 2.3, SD 1.2 "
                        "months) for the two-test pathway; 65% of individuals with CD4 below 200 (AIDS "
                        "exemption) start after gamma(0.5, 0.25) months. Monthly steps round delays up to "
                        "the next cycle.",
            },
        },
    },
    "A9_observed_ineligible_wait": {
        "description": "Comparators in which people who do not meet the certification criteria start ART "
                       "after an observed-like wait instead of waiting for laboratory progression: fallback "
                       "gamma(mean 14, SD 11 months), median about 11 months, matching the multicentre cohort "
                       "(not certifiable: median 334 days, IQR 115 to 708; half started within one year).",
        "settings": {
            "one_month_plus_observed_tail": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {"1 month": {
                    "immediate": False, "apply_criteria": True, "delay_mean": 1.0, "delay_sd": 0.5,
                    "ineligible_fallback": {"delay_mean": 14.0, "delay_sd": 11.0},
                }},
                "model_overrides": {},
                "note": "1 month slot = statutory-minimum comparator (criteria + one-month delay) with the "
                        "ineligible wait capped by the observed fallback distribution.",
            },
            "empirical_pathways_plus_observed_tail": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {"1 month": dict(_OBSERVED_PRACTICE)},
                "model_overrides": {},
                "note": "1 month slot = observed-practice comparator: empirical pathways (A8) with the "
                        "ineligible wait capped by the observed fallback distribution.",
            },
            "criteria_retained_no_wait_plus_observed_tail": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {
                    "0 month": {"immediate": True, "apply_criteria": True,
                                "ineligible_fallback": {"delay_mean": 14.0, "delay_sd": 11.0}},
                    "1 month": {"immediate": False, "apply_criteria": True, "delay_mean": 1.0, "delay_sd": 0.5,
                                "ineligible_fallback": {"delay_mean": 14.0, "delay_sd": 11.0}},
                },
                "model_overrides": {},
                "note": "Factorial cell B under the observed ineligible wait: 0 month slot = criteria retained, "
                        "waiting interval removed; 1 month slot = statutory-minimum comparator; both with the "
                        "fallback. With empirical_pathways_plus_observed_tail and B0_core on the same seeds this "
                        "gives the decomposition of the reform against the observed-practice comparator.",
            },
        },
    },
    "A10_conservative_combined": {
        "description": "Conservative estimates: the observed-practice comparator (empirical pathways plus "
                       "observed ineligible wait) combined with reduced transmission from diagnosed, "
                       "unsuppressed people (recalibrated as in A5).",
        "settings": {
            "empirical_tail_behaviour_0.3": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {"1 month": dict(_OBSERVED_PRACTICE)},
                "model_overrides": {"DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR": 0.3,
                                    "MORTALITY_RECALIB_FOI_SCALE": FOI_SCALE_BY_INFECTIVITY_FACTOR[0.3]},
                "note": "Observed-practice comparator with DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR = 0.3 "
                        "(70% reduction); recalibrated.",
            },
            "empirical_tail_behaviour_0.5": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {"1 month": dict(_OBSERVED_PRACTICE)},
                "model_overrides": {"DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR": 0.5,
                                    "MORTALITY_RECALIB_FOI_SCALE": FOI_SCALE_BY_INFECTIVITY_FACTOR[0.5]},
                "note": "Observed-practice comparator with DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR = 0.5 "
                        "(50% reduction); recalibrated.",
            },
        },
    },
    "A12_age_offset": {
        "description": "Background mortality age offset. The base case reads the sex-specific life table at "
                       "the person's age plus five years (AGE_MORTALITY_OFFSET = 5.0). This setting uses the "
                       "actual age (offset 0.0) for HIV-negative and HIV-positive people alike, without "
                       "recalibration.",
        "settings": {
            "age_offset_0": {
                "scenarios": _TWO_ARMS,
                "scenario_overrides": {},
                "model_overrides": {"AGE_MORTALITY_OFFSET": 0.0},
                "note": "AGE_MORTALITY_OFFSET = 0.0 (base case 5.0); all other inputs unchanged.",
            },
        },
    },
    "A13_foi_multiplier": {
        "description": "Analytic-period force-of-infection multiplier (base case 1.18; B0_core on the same "
                       "seeds is the 1.18 row).",
        "settings": {
            "multiplier_1.00": {
                "scenarios": _THREE_ARMS,
                "scenario_overrides": {},
                "model_overrides": {"MAIN_TRANSMISSION_MULTIPLIER": 1.00},
                "note": "MAIN_TRANSMISSION_MULTIPLIER = 1.00.",
            },
            "multiplier_1.10": {
                "scenarios": _THREE_ARMS,
                "scenario_overrides": {},
                "model_overrides": {"MAIN_TRANSMISSION_MULTIPLIER": 1.10},
                "note": "MAIN_TRANSMISSION_MULTIPLIER = 1.10.",
            },
            "multiplier_1.25": {
                "scenarios": _THREE_ARMS,
                "scenario_overrides": {},
                "model_overrides": {"MAIN_TRANSMISSION_MULTIPLIER": 1.25},
                "note": "MAIN_TRANSMISSION_MULTIPLIER = 1.25.",
            },
            "multiplier_1.30": {
                "scenarios": _THREE_ARMS,
                "scenario_overrides": {},
                "model_overrides": {"MAIN_TRANSMISSION_MULTIPLIER": 1.30},
                "note": "MAIN_TRANSMISSION_MULTIPLIER = 1.30.",
            },
        },
    },
}

# =============================================================================
# Model import (output base must be set before import) and plotting control
# =============================================================================
os.environ["OUTPUT_BASE_DIR"] = os.path.join(OUT_ROOT, "_model_default_output")
os.makedirs(os.environ["OUTPUT_BASE_DIR"], exist_ok=True)

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

from concurrent.futures import ProcessPoolExecutor

import hiv_rapid_art_model as M
if RUN_CEA:
    import cea_analysis as CEA
if RUN_BREAKEVEN:
    import breakeven_analysis_detailed_json as BE

# Pristine copies of the scenario configuration, restored before every setting.
_BASE_SCENARIO_CONFIGS = json.loads(json.dumps(M.SCENARIO_CONFIGS))
_MODEL_OVERRIDE_KEYS = ["ANALYTIC_FOI_MULTIPLIER_BY_YEAR", "TTS_RANGE_MONTHS", "TTS_MIN_MONTHS",
                        "DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR", "MAIN_TRANSMISSION_MULTIPLIER",
                        "ART_START_VL_UPDATE_NEXT_MONTH", "ART_START_VL_UPDATE_NEXT_MONTH_IN_BURNIN",
                        "MORTALITY_RECALIB_FOI_SCALE", "RECENT_FOI_BOOST", "AGE_MORTALITY_OFFSET"]
_BASE_MODEL_VALUES = {k: getattr(M, k) for k in _MODEL_OVERRIDE_KEYS}


def apply_setting(setting):
    """Reset the model to its defaults, then apply this setting's overrides."""
    M.SCENARIO_CONFIGS = json.loads(json.dumps(_BASE_SCENARIO_CONFIGS))
    for k, v in _BASE_MODEL_VALUES.items():
        setattr(M, k, v)
    for label, over in setting.get("scenario_overrides", {}).items():
        M.SCENARIO_CONFIGS[label].update(json.loads(json.dumps(over)))
    for k, v in setting.get("model_overrides", {}).items():
        if k not in _BASE_MODEL_VALUES:
            raise KeyError(f"model override '{k}' is not listed in _MODEL_OVERRIDE_KEYS")
        if k == "ANALYTIC_FOI_MULTIPLIER_BY_YEAR" and v is not None:
            v = {int(y): float(m) for y, m in v.items()}
        setattr(M, k, v)
    M.ENABLED_SCENARIO_LABELS = list(setting["scenarios"])


def _scenario_worker(args):
    """Worker for parallel scenario runs; re-applies the setting so that the
    overrides are in force even if the process pool uses spawn/forkserver."""
    scenario, burnin_pop, seed, setting, output_dir = args
    apply_setting(setting)
    M.SEED = seed
    M.OUTPUT_DIR = output_dir
    results, pop = M.run_simulation(label=scenario, initial_pop=burnin_pop, seed=seed)
    return scenario, results, pop


def _run_cea(seed_dir):
    CEA.Config.INPUT_DIR = seed_dir
    CEA.Config.OUTPUT_DIR = os.path.join(seed_dir, "cea")
    CEA.Config.MAKE_FIGURES = False
    if CEA_PSA_ITERATIONS is not None:
        CEA.Config.PSA_N_ITERATIONS = CEA_PSA_ITERATIONS
    os.makedirs(CEA.Config.OUTPUT_DIR, exist_ok=True)
    CEA.main()


def _run_breakeven(seed_dir):
    BE.INPUT_DIR = seed_dir
    BE.OUTPUT_DIR = os.path.join(seed_dir, "breakeven")
    os.makedirs(BE.OUTPUT_DIR, exist_ok=True)
    BE.main()


def _worker_init():
    """Linux: the kernel stops a scenario process when the process that started it dies, so that no
    scenario process is left behind (holding memory) if the driver is killed."""
    try:
        import ctypes
        ctypes.CDLL("libc.so.6", use_errno=True).prctl(1, signal.SIGKILL)      # PR_SET_PDEATHSIG
    except Exception:
        pass


def _kill_pool(ex):
    """Stop the processes of a pool at once (after an error, a timeout or an interruption)."""
    try:
        procs = list((getattr(ex, "_processes", None) or {}).values())
    except Exception:
        procs = []
    for p in procs:
        try:
            p.kill()
        except Exception:
            pass
    try:
        ex.shutdown(wait=False, cancel_futures=True)
    except Exception:
        pass


def _run_scenarios_in_processes(scenarios, burnin_pop, seed, setting, seed_dir):
    """Run the scenario slots in separate processes, at most MAX_PARALLEL_SCENARIOS[TERMINAL] at a time
    (group after group, each group in fresh processes). Every slot is an independent call of
    run_simulation from the same burn-in with the same seed, so the grouping does not affect the results.
    Returns (results, populations) keyed in the order of `scenarios`. Raises if a process fails or a
    group runs longer than SCENARIO_TIMEOUT_MINUTES."""
    limit = max(1, int(MAX_PARALLEL_SCENARIOS.get(TERMINAL, len(scenarios))))
    results, pops = {}, {}
    for i in range(0, len(scenarios), limit):
        group = scenarios[i:i + limit]
        ex = ProcessPoolExecutor(max_workers=len(group), initializer=_worker_init)
        finished = False
        try:
            futures = {ex.submit(_scenario_worker, (sc, burnin_pop, seed, setting, seed_dir)): sc for sc in group}
            pending = set(futures)
            t0 = time.time()
            next_note = PROGRESS_NOTE_MINUTES * 60.0
            step = max(1.0, min(float(MEMORY_SAMPLE_SECONDS), PROGRESS_NOTE_MINUTES * 60.0))
            while pending:
                done, pending = concurrent.futures.wait(pending, timeout=step,
                                                        return_when=concurrent.futures.FIRST_EXCEPTION)
                _note_memory("scenario slots [" + ", ".join(group) + "]")
                for fut in done:
                    sc, res, pop = fut.result()          # raises here if the process failed
                    results[sc] = res
                    pops[sc] = pop
                if pending:
                    elapsed = (time.time() - t0) / 60.0
                    names = ", ".join(futures[f] for f in futures if f in pending)
                    if elapsed > SCENARIO_TIMEOUT_MINUTES:
                        raise TimeoutError(f"scenario slots [{names}] still running after {elapsed:.0f} min "
                                           f"(SCENARIO_TIMEOUT_MINUTES = {SCENARIO_TIMEOUT_MINUTES})")
                    if elapsed * 60.0 >= next_note:
                        print(f"      ... [{names}] running ({elapsed:.0f} min)", flush=True)
                        _touch_instance()
                        next_note += PROGRESS_NOTE_MINUTES * 60.0
            _note_memory("scenario slots [" + ", ".join(group) + "] finished", force=True)
            finished = True
        finally:
            if finished:
                ex.shutdown(wait=True)
            else:
                _kill_pool(ex)
    return {sc: results[sc] for sc in scenarios}, {sc: pops[sc] for sc in scenarios}


def run_one_seed(seed, setting, seed_dir, should_yield=None):
    """Burn-in, scenario slots, outputs and post-processing of one seed. Returns False (nothing completed)
    if should_yield() reports after the burn-in that another terminal started this seed first."""
    os.makedirs(seed_dir, exist_ok=True)

    apply_setting(setting)
    M.SEED = seed
    M.OUTPUT_DIR = seed_dir
    _note_memory("burn-in", force=True)
    burnin_pop, _ = M.run_burnin(seed=seed)
    _note_memory("burn-in finished", force=True)
    if should_yield is not None and should_yield():
        del burnin_pop
        gc.collect()
        return False
    scenarios = M.get_enabled_scenarios()
    all_results, all_pops = {}, {}
    if PARALLEL_SCENARIOS and len(scenarios) >= 2:
        all_results, all_pops = _run_scenarios_in_processes(scenarios, burnin_pop, seed, setting, seed_dir)
    else:
        for sc in scenarios:
            res, pop = M.run_simulation(label=sc, initial_pop=burnin_pop, seed=seed)
            all_results[sc] = res
            all_pops[sc] = pop
    apply_setting(setting)          # the parent process state is used by save_all_outputs
    M.SEED = seed
    M.OUTPUT_DIR = seed_dir
    _note_memory("saving the outputs", force=True)
    M.save_all_outputs(all_results, all_pops, scenarios)
    plt.close("all")
    _note_memory("outputs saved", force=True)

    del burnin_pop, all_results, all_pops
    gc.collect()

    if RUN_CEA:
        _note_memory("cost-effectiveness post-processing", force=True)
        _run_cea(seed_dir)
    if RUN_BREAKEVEN:
        _note_memory("break-even post-processing", force=True)
        _run_breakeven(seed_dir)
    _note_memory("post-processing finished", force=True)
    return True


# =============================================================================
# Seed ownership, claims, and logging
# =============================================================================
def build_owner_map():
    """{seed: terminal}. The seeds are dealt out in their fixed order, each to the terminal whose queue
    would finish earliest, so that for any leading subset of the seeds the number owned by a terminal
    is proportional to its speed (1 / MINUTES_PER_SEED_RUN)."""
    minutes = {t: float(m) for t, m in MINUTES_PER_SEED_RUN.items()}
    count = {t: 0 for t in minutes}
    owner = {}
    for s in PUBLICATION_SEEDS:
        t = min(minutes, key=lambda k: ((count[k] + 1) * minutes[k], k))
        count[t] += 1
        owner[s] = t
    return owner


def _log(msg, also_print=True):
    line = f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} [{TERMINAL}] {msg}"
    if also_print:
        print(line, flush=True)
    try:
        os.makedirs(LOG_DIR, exist_ok=True)
        with open(os.path.join(LOG_DIR, f"run_{TERMINAL}.log"), "a", encoding="utf-8") as f:
            f.write(line + "\n")
    except OSError:
        pass


# -----------------------------------------------------------------------------
# Records for the diagnosis of stops: memory (the crash log is switched on after the CONFIG block)
# -----------------------------------------------------------------------------
_MEMORY = {"what": "", "lowest_available": None, "lowest_swap_free": None, "low_since_line": None, "last_line": 0.0}


def _meminfo_gb():
    """MemTotal, MemAvailable, SwapTotal, SwapFree of this machine in GB (Linux); None elsewhere."""
    try:
        v = {}
        with open("/proc/meminfo", encoding="ascii") as f:
            for line in f:
                k, _, rest = line.partition(":")
                if k in ("MemTotal", "MemAvailable", "SwapTotal", "SwapFree"):
                    v[k] = float(rest.split()[0]) / 1048576.0
        return v if len(v) == 4 else None
    except Exception:
        return None


def _peak_rss_gb():
    """(largest memory use so far of this process, of the largest finished scenario process) in GB; Linux."""
    try:
        import resource
        return (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576.0,
                resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1048576.0)
    except Exception:
        return None, None


def _reset_memory(what):
    _MEMORY.update({"what": what, "lowest_available": None, "lowest_swap_free": None})


def _note_memory(stage, force=False):
    """Read the free memory; keep the lowest value of the seed; write a line to _logs/memory_<TERMINAL>.log
    at most every MEMORY_LOG_MINUTES, or at once with force=True (stages of a seed)."""
    m = _meminfo_gb()
    if m is None:
        return
    for key, val in (("lowest_available", m["MemAvailable"]), ("lowest_swap_free", m["SwapFree"]),
                     ("low_since_line", m["MemAvailable"])):
        if _MEMORY[key] is None or val < _MEMORY[key]:
            _MEMORY[key] = val
    now = time.time()
    if not force and now - _MEMORY["last_line"] < MEMORY_LOG_MINUTES * 60.0:
        return
    line = (f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} [{TERMINAL}] {_MEMORY['what']} {stage}: available "
            f"{m['MemAvailable']:.1f} of {m['MemTotal']:.1f} GB (lowest since the last line "
            f"{_MEMORY['low_since_line']:.1f}); swap free {m['SwapFree']:.1f} of {m['SwapTotal']:.1f} GB")
    _MEMORY["last_line"] = now
    _MEMORY["low_since_line"] = None
    try:
        os.makedirs(LOG_DIR, exist_ok=True)
        with open(os.path.join(LOG_DIR, f"memory_{TERMINAL}.log"), "a", encoding="utf-8") as f:
            f.write(line + "\n")
    except OSError:
        pass


def _memory_summary():
    """Memory record of the seed just finished, for run_info.json."""
    m = _meminfo_gb()
    own, child = _peak_rss_gb()
    r = lambda x: None if x is None else round(x, 2)
    return {"total_gb": r(m["MemTotal"]) if m else None, "swap_total_gb": r(m["SwapTotal"]) if m else None,
            "lowest_available_gb": r(_MEMORY["lowest_available"]), "lowest_swap_free_gb": r(_MEMORY["lowest_swap_free"]),
            "driver_peak_rss_gb": r(own), "largest_scenario_process_peak_rss_gb": r(child),
            "note": "lowest values during this seed; peak values since the start of the driver process"}


def _seed_done(seed_dir):
    if not os.path.exists(os.path.join(seed_dir, "heterogeneity_summary.csv")):
        return False
    if RUN_CEA and not os.path.exists(os.path.join(seed_dir, "cea", "cea_report_ready.json")):
        return False
    if RUN_BREAKEVEN and not os.path.exists(os.path.join(seed_dir, "breakeven", "breakeven_report_ready.json")):
        return False
    return True


def _claim_path(seed_dir, terminal=None):
    return os.path.join(seed_dir, f"_claim_{terminal or TERMINAL}.json")


def _read_claims(seed_dir):
    """[(terminal, epoch)] of the claims of other terminals that are not yet timed out."""
    out = []
    now = time.time()
    for path in glob.glob(os.path.join(seed_dir, "_claim_*.json")):
        term = os.path.basename(path)[len("_claim_"):-len(".json")]
        if term == TERMINAL:
            continue
        try:
            with open(path, encoding="utf-8") as f:
                epoch = float(json.load(f)["epoch"])
        except Exception:
            try:
                epoch = os.path.getmtime(path)      # file not yet fully synchronised
            except OSError:
                continue
        if (now - epoch) / 60.0 < CLAIM_TIMEOUT_MINUTES:
            out.append((term, epoch))
    return out


def _write_claim(seed_dir):
    os.makedirs(seed_dir, exist_ok=True)
    epoch = time.time()
    with open(_claim_path(seed_dir), "w", encoding="utf-8") as f:
        json.dump({"terminal": TERMINAL, "epoch": epoch,
                   "started": datetime.datetime.fromtimestamp(epoch).isoformat(timespec="seconds")}, f)
    return epoch


def _remove_claim(seed_dir):
    """Withdraw the claim of this terminal (after an error, an interruption, or when the seed is left to
    another terminal). The claim of a completed seed is kept: its outputs may reach the other terminals
    later than the removal of the claim would, and the seed would then look unstarted to them."""
    try:
        os.remove(_claim_path(seed_dir))
    except OSError:
        pass


def seed_state(seed_dir):
    """('done', None) / ('elsewhere', terminal) / ('free', None)."""
    if _seed_done(seed_dir):
        return "done", None
    claims = _read_claims(seed_dir)
    if claims:
        return "elsewhere", min(claims, key=lambda c: c[1])[0]
    return "free", None


def _model_file_md5():
    with open(M.__file__, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def _jsonable(v):
    return {str(a): b for a, b in v.items()} if isinstance(v, dict) else v


def _write_setting_manifest(analysis, setting_name, setting, setting_dir):
    manifest = {
        "analysis": analysis,
        "setting": setting_name,
        "description": ANALYSES[analysis]["description"],
        "note": setting.get("note", ""),
        "scenarios": setting["scenarios"],
        "scenario_overrides": setting.get("scenario_overrides", {}),
        "model_overrides": {k: _jsonable(v) for k, v in setting.get("model_overrides", {}).items()},
        "model_defaults": {k: _jsonable(v) for k, v in _BASE_MODEL_VALUES.items()},
        "model_file_md5": _MODEL_MD5,
    }
    path = os.path.join(setting_dir, "setting_config.json")
    try:
        with open(path, encoding="utf-8") as f:
            if json.load(f) == manifest:
                return                       # already written (by this or another terminal)
    except Exception:
        pass
    with open(path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)


def _check_model_defaults():
    problems = []
    for k, expected in EXPECTED_MODEL_DEFAULTS.items():
        actual = getattr(M, k, "<missing>")
        if expected is None:
            problems.append(f"{k}: the expected value is not yet fixed in EXPECTED_MODEL_DEFAULTS (model file: {actual})")
        elif actual != expected:
            problems.append(f"{k}: model file has {actual}, expected {expected}")
    if problems:
        raise SystemExit("The main queue was not started because the model file in this folder does not have "
                         "the expected defaults:\n  - " + "\n  - ".join(problems) +
                         f"\nModel file: {M.__file__}")


_MODEL_MD5 = _model_file_md5()
_FAILED = set()          # (analysis, setting, seed) that raised an error in this process
_NOTED = set()           # log lines already written once by this process (seeds running elsewhere, item headers)
_N_RUNS = [0]            # seed-runs completed by this process
_ERRORS_IN_A_ROW = [0]   # seeds that ended with an error since the last completed seed
_STARTED = time.time()


def _limit_reached():
    return MAX_SEED_RUNS is not None and _N_RUNS[0] >= int(MAX_SEED_RUNS)


def _too_many_errors():
    return MAX_CONSECUTIVE_ERRORS is not None and _ERRORS_IN_A_ROW[0] >= int(MAX_CONSECUTIVE_ERRORS)


def _stop_now():
    return _limit_reached() or _too_many_errors()


def _instance_path():
    return os.path.join(LOG_DIR, f"instance_{TERMINAL}.json")


def _touch_instance():
    """Record that this terminal name is running here (host, process, time of the last sign of life)."""
    try:
        os.makedirs(LOG_DIR, exist_ok=True)
        with open(_instance_path(), "w", encoding="utf-8") as f:
            json.dump({"terminal": TERMINAL, "host": socket.gethostname(), "pid": os.getpid(),
                       "started": datetime.datetime.fromtimestamp(_STARTED).isoformat(timespec="seconds"),
                       "heartbeat": time.time()}, f)
    except OSError:
        pass


def _release_terminal_name():
    try:
        with open(_instance_path(), encoding="utf-8") as f:
            mine = json.load(f)
        if mine.get("host") == socket.gethostname() and int(mine.get("pid", -1)) == os.getpid():
            os.remove(_instance_path())
    except Exception:
        pass


def _pid_is_this_driver(pid):
    """Linux: True if process `pid` exists and is a run of this script."""
    try:
        with open(f"/proc/{int(pid)}/cmdline", "rb") as f:
            return b"run_revision_scenarios" in f.read()
    except Exception:
        return False


def _claim_terminal_name():
    """Refuse to start if this terminal name is already running (see TERMINAL_IN_USE_MINUTES in CONFIG)."""
    host = socket.gethostname()
    try:
        with open(_instance_path(), encoding="utf-8") as f:
            other = json.load(f)
    except Exception:
        other = None
    if other and not ALLOW_TERMINAL_IN_USE:
        age = (time.time() - float(other.get("heartbeat", 0.0))) / 60.0
        if other.get("host") == host:
            if int(other.get("pid", -1)) != os.getpid() and _pid_is_this_driver(other.get("pid")):
                raise SystemExit(f'TERMINAL = "{TERMINAL}" is already running on this machine (process '
                                 f'{other.get("pid")}, started {other.get("started")}). Stop that run first; '
                                 "a second run would repeat the same seeds.")
        elif age < TERMINAL_IN_USE_MINUTES:
            raise SystemExit(f'TERMINAL = "{TERMINAL}" is in use on another machine ({other.get("host")}, active '
                             f"{age:.0f} min ago). Set TERMINAL to the name of this machine. If that machine has "
                             f"been stopped, wait {TERMINAL_IN_USE_MINUTES} min or set ALLOW_TERMINAL_IN_USE = True.")
    _touch_instance()
    atexit.register(_release_terminal_name)


def _log_once(msg):
    if msg not in _NOTED:
        _NOTED.add(msg)
        _log(msg)


def try_seed(analysis, setting_name, setting, setting_dir, seed, owner_terminal, tag, header=None):
    """Run one seed unless it is done, running elsewhere, or failed earlier in this process.
    Returns True if this terminal completed it."""
    key = (analysis, setting_name, seed)
    seed_dir = os.path.join(setting_dir, f"seed_{seed}")
    state, who = seed_state(seed_dir)
    if state == "done" or key in _FAILED:
        return False
    if state == "elsewhere":
        _log_once(f"  {tag} seed {seed}: running on {who}, left alone")
        return False
    if header:
        _log_once(header)
    my_epoch = _write_claim(seed_dir)

    def should_yield():
        earlier = [c for c in _read_claims(seed_dir) if (c[1], c[0]) < (my_epoch, TERMINAL)]
        return bool(earlier) or _seed_done(seed_dir)

    taken = "" if owner_terminal in (None, TERMINAL) else f" (taken over from {owner_terminal})"
    _log(f"  {tag} seed {seed}: start{taken}")
    _touch_instance()
    _reset_memory(f"{tag} seed {seed}")
    _fault_note(f"{tag} seed {seed}: start")
    ts = time.time()
    try:
        completed = run_one_seed(seed, setting, seed_dir, should_yield)
    except KeyboardInterrupt:
        _remove_claim(seed_dir)
        _log(f"  {tag} seed {seed}: interrupted by the user")
        raise
    except Exception as e:
        _FAILED.add(key)
        _ERRORS_IN_A_ROW[0] += 1
        _remove_claim(seed_dir)
        _log(f"  {tag} seed {seed}: ERROR {type(e).__name__}: {e}\n{traceback.format_exc()}")
        return False
    minutes = (time.time() - ts) / 60.0
    if not completed:
        _remove_claim(seed_dir)
        _log(f"  {tag} seed {seed}: another terminal started it first, left to that terminal "
             f"({minutes:.1f} min spent on the burn-in)")
        return False
    info = {
        "analysis": analysis, "setting": setting_name, "seed": int(seed),
        "terminal": TERMINAL, "owner_terminal": owner_terminal,
        "started": datetime.datetime.fromtimestamp(ts).isoformat(timespec="seconds"),
        "finished": datetime.datetime.now().isoformat(timespec="seconds"),
        "finished_epoch": time.time(),
        "minutes": round(minutes, 2),
        "model_file_md5": _MODEL_MD5,
        "model_values": {k: _jsonable(getattr(M, k)) for k in _MODEL_OVERRIDE_KEYS},
        "scenarios": list(setting["scenarios"]),
        "max_parallel_scenarios": int(MAX_PARALLEL_SCENARIOS.get(TERMINAL, len(setting["scenarios"]))),
        "memory": _memory_summary(),
    }
    with open(os.path.join(seed_dir, "run_info.json"), "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2, ensure_ascii=False)
    _N_RUNS[0] += 1
    _ERRORS_IN_A_ROW[0] = 0
    _touch_instance()
    mem = info["memory"]
    low = ("" if mem["lowest_available_gb"] is None else
           f"; lowest available memory {mem['lowest_available_gb']:.1f} of {mem['total_gb']:.1f} GB")
    _log(f"  {tag} seed {seed}: done in {minutes:.1f} min{low}")
    _fault_note(f"{tag} seed {seed}: done")
    return True


def _seeds_of_item(n_seeds):
    if RUN == "regression":
        return [REGRESSION_SEED]
    return list(PUBLICATION_SEEDS[: int(n_seeds)])


def _pick_seed_to_take(analysis, setting_name, setting_dir, seeds, owner):
    """A free seed near the end of the list of the other terminal that has the most free seeds left in
    this item, or None. Each terminal takes a different position from the end (its rank among the
    terminal names), so that two terminals that become free at the same moment do not take the same seed."""
    free = {}
    for s in seeds:
        o = owner[s]
        if o == TERMINAL or (analysis, setting_name, s) in _FAILED:
            continue
        if seed_state(os.path.join(setting_dir, f"seed_{s}"))[0] == "free":
            free.setdefault(o, []).append(s)
    if not free:
        return None
    o = max(free, key=lambda t: (len(free[t]), t))
    rank = sorted(MINUTES_PER_SEED_RUN).index(TERMINAL)
    return free[o][-1 - (rank % len(free[o]))]


def run_item(analysis, setting_name, n_seeds, owner):
    """Run this terminal's part of one queue item. Returns the number of seeds it completed."""
    setting = ANALYSES[analysis]["settings"][setting_name]
    setting_dir = os.path.join(OUT_ROOT, analysis, setting_name)
    os.makedirs(setting_dir, exist_ok=True)
    _write_setting_manifest(analysis, setting_name, setting, setting_dir)
    seeds = _seeds_of_item(n_seeds)
    tag = f"[{analysis}/{setting_name}]"
    if RUN == "regression":
        own = list(seeds)
    else:
        own = [s for s in seeds if owner[s] == TERMINAL]
    n_done_before = sum(1 for s in seeds if _seed_done(os.path.join(setting_dir, f"seed_{s}")))
    if n_done_before == len(seeds):
        return 0
    header = f"{tag} {len(seeds)} seeds, {n_done_before} done; own seeds on this terminal: {len(own)}"
    n = 0
    for s in own:
        if _stop_now():
            return n
        n += bool(try_seed(analysis, setting_name, setting, setting_dir, s, TERMINAL, tag, header))
    if HELP_OTHER_TERMINALS and RUN == "main":
        while not _stop_now():
            s = _pick_seed_to_take(analysis, setting_name, setting_dir, seeds, owner)
            if s is None:
                break
            before = len(_FAILED)
            ok = try_seed(analysis, setting_name, setting, setting_dir, s, owner[s], tag, header)
            n += bool(ok)
            if not ok and len(_FAILED) == before and seed_state(os.path.join(setting_dir, f"seed_{s}"))[0] == "free":
                # neither completed, failed, nor claimed by another terminal: avoid an endless loop
                _FAILED.add((analysis, setting_name, s))
    return n


def _queue_status(queue):
    """(number of seeds not done, number of those marked as running on another terminal)."""
    not_done = elsewhere = 0
    seen = set()
    for analysis, setting_name, n_seeds in queue:
        setting_dir = os.path.join(OUT_ROOT, analysis, setting_name)
        for s in _seeds_of_item(n_seeds):
            if (analysis, setting_name, s) in seen:
                continue
            seen.add((analysis, setting_name, s))
            state, _ = seed_state(os.path.join(setting_dir, f"seed_{s}"))
            if state != "done":
                not_done += 1
                elsewhere += (state == "elsewhere")
    return not_done, elsewhere


def main():
    if RUN not in ("main", "regression"):
        raise SystemExit('RUN must be "main" or "regression"')
    if TERMINAL not in MINUTES_PER_SEED_RUN:
        raise SystemExit('Set TERMINAL in the CONFIG block to the name of this machine, one of: '
                         + ", ".join(MINUTES_PER_SEED_RUN) + f' (it is "{TERMINAL}" now).')
    if len(set(PUBLICATION_SEEDS)) != len(PUBLICATION_SEEDS):
        raise SystemExit("PUBLICATION_SEEDS contains duplicates")
    queue = MAIN_QUEUE if RUN == "main" else REGRESSION_QUEUE
    for analysis, setting_name, n_seeds in queue:
        if analysis not in ANALYSES or setting_name not in ANALYSES[analysis]["settings"]:
            raise SystemExit(f"Unknown queue item: {analysis} / {setting_name}")
        if not 1 <= int(n_seeds) <= len(PUBLICATION_SEEDS):
            raise SystemExit(f"Number of seeds out of range for {analysis} / {setting_name}: {n_seeds}")
    if RUN == "main":
        _check_model_defaults()
        unset = [f"{a} / {st}: {k}" for a, st, _ in queue
                 for k, v in ANALYSES[a]["settings"][st].get("model_overrides", {}).items() if v is None]
        if unset:
            raise SystemExit("The main queue was not started because these model overrides are not yet fixed "
                             "(see FOI_SCALE_BY_INFECTIVITY_FACTOR in CONFIG):\n  - " + "\n  - ".join(sorted(set(unset))))
    _claim_terminal_name()
    _enable_fault_log()                         # already on from the import, unless TERMINAL was set afterwards
    _fault_note(f"driver started (process {os.getpid()} on {socket.gethostname()})")
    owner = build_owner_map()
    share = {t: sum(1 for s in PUBLICATION_SEEDS if owner[s] == t) for t in MINUTES_PER_SEED_RUN}

    _log("=" * 70)
    _log(f"RUN = {RUN}; {len(queue)} queue items; model file md5 {_MODEL_MD5}")
    _log(f"output root: {OUT_ROOT}")
    _log(f"seeds owned per terminal (of {len(PUBLICATION_SEEDS)}): {share}")
    _log(f"scenario slots at a time on this terminal: {MAX_PARALLEL_SCENARIOS.get(TERMINAL, 'all')}; "
         f"host {socket.gethostname()}, process {os.getpid()}")
    _log("=" * 70)
    t0 = time.time()
    pass_no = 0
    while True:
        pass_no += 1
        n_pass = 0
        for analysis, setting_name, n_seeds in queue:
            if _stop_now():
                break
            n_pass += run_item(analysis, setting_name, n_seeds, owner)
        if _too_many_errors():
            _log(f"STOPPED: {_ERRORS_IN_A_ROW[0]} seeds in a row ended with an error on this terminal "
                 f"(MAX_CONSECUTIVE_ERRORS = {MAX_CONSECUTIVE_ERRORS}). The seeds are left to the other terminals. "
                 f"See {os.path.join(LOG_DIR, 'run_' + TERMINAL + '.log')} before restarting.")
            break
        if _limit_reached():
            _log(f"MAX_SEED_RUNS = {MAX_SEED_RUNS} reached; stopped")
            break
        not_done, elsewhere = _queue_status(queue)
        if not_done == 0:
            _log("queue complete: every seed of every item is done")
            break
        if n_pass > 0:
            continue                     # look again from the top: seeds may have been freed meanwhile
        if elsewhere > 0 and WAIT_FOR_OTHERS and RUN == "main":
            _log(f"{not_done} seeds not done, {elsewhere} of them running on other terminals; "
                 f"checking again in {WAIT_POLL_MINUTES} min")
            _touch_instance()
            time.sleep(WAIT_POLL_MINUTES * 60.0)
            continue
        _log(f"stopped with {not_done} seeds not done ({elsewhere} running on other terminals, "
             f"{len(_FAILED)} failed on this terminal; see {os.path.join(LOG_DIR, 'run_' + TERMINAL + '.log')})")
        break
    _log(f"TOTAL: {(time.time() - t0) / 60:.1f} min, {_N_RUNS[0]} seed-runs completed on this terminal")
    if RUN == "regression":
        _log("Next: run compare_regression.py")
    else:
        _log("Next: run check_progress.py for the overall state")


if __name__ == "__main__":
    main()

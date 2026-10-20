#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Driver for the revision (R1) scenario analyses.

Runs the additional analyses requested at peer review with the modified model
(hiv_rapid_art_model.py in this folder), one "setting" at a time, for a set of
random-number seeds. Every setting runs the "0 month" and "1 month" scenario
slots from the same burn-in (common random numbers), then the cost-
effectiveness (cea_analysis) and break-even (breakeven_analysis_detailed_json)
post-processing, exactly as in the published pipeline. Output layout:

    outputs/revision_R1_runs/<ANALYSIS>/<setting>/seed_<n>/...

Run as-is in VS Code: all choices are in the CONFIG block below, no
command-line arguments. Set ANALYSIS and TERMINAL, then run. Seeds already
completed are skipped, so the script can be restarted safely.

Analyses (see README.md, section 'Version 1.1.0'):
  A1_factorial               eligibility criteria vs waiting interval (2 settings)
  A2_foi_trajectory          time-varying analytic-period force of infection (2)
  A4_time_to_suppression     faster / slower time to viral suppression (2)
  A5_post_diagnosis_behaviour reduced infectivity of diagnosed unsuppressed people (2)
  A7_delayed_start           immediate ART implemented from 2027 / 2030 (2)
  A11_timing_convention      scheduled ART starts: first viral-load update deferred one cycle (1)
  A12_age_offset             background mortality read at the actual age, offset 0 (1)
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
from _paths import BASE_144, REVISION_RUNS

ANALYSIS = "A1_factorial"        # one of the keys of ANALYSES below
SETTINGS_TO_RUN = None           # None = all settings of the analysis, or a list of setting names

TERMINAL = "node1"          # "node1" / "node2" / "node3" / "node4" / "node5"
# node1 to node5 label the five-way split of the seed list used to share the runs across
# workstations; on a single machine, run each label in turn (or all seeds of one label).
SEED_SETS = {                    # same 100-seed split as the FOI sensitivity run
    "node1": [14, 18, 24, 42, 55, 58, 61, 66, 69, 71,
                   75, 77, 87, 97, 101, 102, 104, 107, 110, 119,
                   123, 128, 145, 150, 157, 161, 176],
    "node2":  [178, 188, 194, 202, 203, 206, 207, 209, 216, 219,
                   225, 228, 234, 244, 247, 249, 276, 286, 303, 306,
                   308, 309, 313],
    "node3":        [317, 327, 330, 340, 341, 345, 347, 355, 356, 363,
                   364, 380, 389, 404, 417, 424, 427, 441, 448, 456,
                   465],
    "node4":      [466, 469, 472, 477, 484, 486, 488, 492, 495, 504,
                   505, 508, 515, 516, 527],
    "node5":       [536, 543, 547, 567, 581, 582, 585, 586, 590, 592,
                   596, 606, 608, 613],
}
# Number of seeds per analysis (union over the five terminals). Each terminal runs
# the first k seeds of its list, with k proportional to its share of the 100-seed
# split, so the machines coordinate without further configuration. Rationale
# (convergence_diagnostics/mcse_table.csv, published run): relative Monte Carlo SE
# at 200 seeds is 2.1% (infections averted) and 2.6% (NMB); 10% relative MCSE needs
# about 32 and 51 seeds respectively, so 50 seeds give about 4 to 5% and 30 seeds
# about 5 to 7%. The factorial decomposition keeps 100 because its contrasts are
# differences of differences.
SEED_TARGET_BY_ANALYSIS = {
    "A0_regression_check": 1,
    "A1_factorial": 100,
    "A2_foi_trajectory": 50,
    "A4_time_to_suppression": 30,
    "A5_post_diagnosis_behaviour": 100,   # key sensitivity analysis
    "A7_delayed_start": 50,
    "A8_status_quo_empirical": 50,
    "A9_observed_ineligible_wait": None,   # per-setting targets below (SEED_TARGET_BY_SETTING)
    "A10_conservative_combined": 50,
    "A11_timing_convention": 50,     # first viral-load update of scheduled starters deferred one cycle
    "A12_age_offset": 30,            # background mortality read at the actual age (offset 0)
}
# Per-setting override of the seed target (analysis, setting) -> seeds.
SEED_TARGET_BY_SETTING = {
    ("A9_observed_ineligible_wait", "one_month_plus_observed_tail"): 50,
    ("A9_observed_ineligible_wait", "empirical_pathways_plus_observed_tail"): 100,
    ("A9_observed_ineligible_wait", "criteria_retained_no_wait_plus_observed_tail"): 50,
    ("A11_timing_convention", "scheduled_start_vl_update_next_month"): 50,
    ("A11_timing_convention", "scheduled_start_vl_update_next_month_observed_practice"): 50,
}
SEED_LIMIT = None                # overrides the allocation above for THIS terminal (e.g. 1 for a pilot)
SEED_ORDER = "forward"           # "forward" / "reverse" / "shuffle"

RUN_CEA = True
RUN_BREAKEVEN = True
PARALLEL_SCENARIOS = True        # run the scenario slots of one seed in parallel processes
SKIP_PLOTS = True                # suppress the model's per-scenario figure files
KEEP_PLOT_FILES = ["summary_plots.jpg"]
CEA_PSA_ITERATIONS = 200

OUT_ROOT = REVISION_RUNS
# =============================================================================


def _linear_schedule(y0, v0, y1, v1, y_end):
    """{year: value} rising/falling linearly from (y0, v0) to (y1, v1), then constant to y_end."""
    sched = {}
    for y in range(y0, y_end + 1):
        if y >= y1:
            sched[y] = float(v1)
        else:
            sched[y] = float(v0 + (v1 - v0) * (y - y0) / (y1 - y0))
    return sched


# Each setting: scenario slots to run, per-slot scenario-config overrides, and
# model-global overrides. The "0 month" slot is always compared against the
# "1 month" slot by the downstream CEA / break-even scripts.
ANALYSES = {
    "A0_regression_check": {
        "description": "Default configuration (no overrides): the three published scenarios. Outputs must "
                       "match outputs/core_multiseed/seed_<n>/ exactly (run with SEED_LIMIT = 1, then compare "
                       "with compare_regression.py) before any other analysis is launched.",
        "settings": {
            "default": {
                "scenarios": ["0 month", "1 month", "2 months"],
                "scenario_overrides": {},
                "model_overrides": {},
                "note": "Published model behaviour; used only to verify that the modified code reproduces it.",
            },
        },
    },
    "A1_factorial": {
        "description": "Factorial decomposition of the certification reform into removal of the "
                       "eligibility criteria and removal of the waiting interval.",
        "settings": {
            "criteria_retained_no_wait": {
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {"0 month": {"immediate": True, "apply_criteria": True}},
                "model_overrides": {},
                "note": "0 month slot = cell B (criteria retained, waiting interval removed); "
                        "1 month slot = cell A (current pathway). Contrast A - B = waiting removal given criteria.",
            },
            "criteria_removed_one_month": {
                "scenarios": ["0 month", "1 month"],
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
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {},
                "model_overrides": {"ANALYTIC_FOI_MULTIPLIER_BY_YEAR": _linear_schedule(2024, 1.18, 2043, 1.00, 2063)},
                "note": "Multiplier falls linearly from 1.18 (2024) to 1.00 (2043), then constant.",
            },
            "resurgence": {
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {},
                "model_overrides": {"ANALYTIC_FOI_MULTIPLIER_BY_YEAR": _linear_schedule(2024, 1.18, 2033, 1.30, 2063)},
                "note": "Multiplier rises linearly from 1.18 (2024) to 1.30 (2033), then constant.",
            },
        },
    },
    "A11_timing_convention": {
        "description": "Timing convention for scheduled ART starts (sensitivity analysis A11). In the "
                       "published model a scheduled starter receives the first viral-load update in the start "
                       "month, whereas an immediate starter (who starts after the viral-load update of the "
                       "diagnosis month) receives it one month later, so the realised difference in the onset of "
                       "viral-load decline between immediate ART and a wait of k months is k - 1 months. This "
                       "setting defers the first update of scheduled starters by one cycle so that both pathways "
                       "share the same convention (difference k months). Compare with the published base case on "
                       "the same seeds (summarize_revision_scenarios.py adds it automatically).",
        "settings": {
            "scheduled_start_vl_update_next_month": {
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {},
                "model_overrides": {"ART_START_VL_UPDATE_NEXT_MONTH": True},
                "note": "Statutory-minimum comparison. ART_START_VL_UPDATE_NEXT_MONTH = True: the first "
                        "viral-load update of every start processed at a scheduled start month (eligible "
                        "people after their wait, restarts after re-engagement) is deferred to the next cycle, "
                        "the convention already applied to starts at diagnosis; analytic period only.",
            },
            "scheduled_start_vl_update_next_month_observed_practice": {
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {"1 month": {
                    "immediate": False, "apply_criteria": True,
                    "delay_mean": 2.3, "delay_sd": 1.2,
                    "fast_track": {"delay_mean": 0.5, "delay_sd": 0.25, "share": 0.65, "cd4_below": 200.0},
                    "ineligible_fallback": {"delay_mean": 14.0, "delay_sd": 11.0},
                }},
                "model_overrides": {"ART_START_VL_UPDATE_NEXT_MONTH": True},
                "note": "Observed-practice comparison (A9 empirical_pathways_plus_observed_tail) under the "
                        "same timing convention; compare with A9 empirical_pathways_plus_observed_tail on the "
                        "same seeds (the published convention). Optional second setting.",
            },
        },
    },
    "A12_age_offset": {
        "description": "Background mortality age offset (sensitivity analysis A12). The published model reads "
                       "the sex-specific life table at the person's age plus five years (AGE_MORTALITY_OFFSET = 5.0, "
                       "inherited from the calibration of earlier versions). This setting uses the actual age "
                       "(offset 0.0) for HIV-negative and HIV-positive people alike, without recalibration.",
        "settings": {
            "age_offset_0": {
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {},
                "model_overrides": {"AGE_MORTALITY_OFFSET": 0.0},
                "note": "AGE_MORTALITY_OFFSET = 0.0 (published model 5.0); all other inputs unchanged.",
            },
        },
    },
    "A4_time_to_suppression": {
        "description": "Time from ART initiation to viral suppression "
                       "(base case: 1 + 5 x Beta(2,5) months, median about 2.3 months).",
        "settings": {
            "fast": {
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {},
                "model_overrides": {"TTS_RANGE_MONTHS": 2.0},
                "note": "1 + 2 x Beta(2,5) months (median about 1.5, maximum 3).",
            },
            "slow": {
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {},
                "model_overrides": {"TTS_RANGE_MONTHS": 8.0},
                "note": "1 + 8 x Beta(2,5) months (median about 3.1, maximum 9).",
            },
        },
    },
    "A5_post_diagnosis_behaviour": {
        "description": "Reduced transmission behaviour after diagnosis: relative infectivity of "
                       "diagnosed, virally unsuppressed sources (base case 1.0).",
        "settings": {
            "factor_0.7": {
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {},
                "model_overrides": {"DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR": 0.7},
                "note": "30% reduction in transmission from diagnosed unsuppressed sources.",
            },
            "factor_0.5": {
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {},
                "model_overrides": {"DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR": 0.5},
                "note": "50% reduction in transmission from diagnosed unsuppressed sources.",
            },
        },
    },
    "A8_status_quo_empirical": {
        "description": "Comparator with a delay distribution matched to the multicentre J-HIV RWD cohort "
                       "(11 hospitals, 2015-2024): two-test pathway median 70 days (IQR 42 to 106) after "
                       "the first visit, and an AIDS exemption pathway (about 20% of people, median 11 days).",
        "settings": {
            "empirical_pathways": {
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {"1 month": {
                    "immediate": False, "apply_criteria": True,
                    "delay_mean": 2.3, "delay_sd": 1.2,
                    "fast_track": {"delay_mean": 0.5, "delay_sd": 0.25, "share": 0.65, "cd4_below": 200.0},
                }},
                "model_overrides": {},
                "note": "1 month slot = empirical status quo: eligibility-to-ART delay gamma(mean 2.3, SD 1.2 "
                        "months) for the two-test pathway; 65% of individuals with CD4 below 200 (AIDS "
                        "exemption) start after gamma(0.5, 0.25) months. Monthly steps round delays up to "
                        "the next cycle, so the modelled medians are 61 to 91 days versus 70 observed.",
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
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {"1 month": {
                    "immediate": False, "apply_criteria": True, "delay_mean": 1.0, "delay_sd": 0.5,
                    "ineligible_fallback": {"delay_mean": 14.0, "delay_sd": 11.0},
                }},
                "model_overrides": {},
                "note": "1 month slot = statutory-minimum comparator (criteria + one-month delay) with the "
                        "ineligible wait capped by the observed fallback distribution.",
            },
            "empirical_pathways_plus_observed_tail": {
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {"1 month": {
                    "immediate": False, "apply_criteria": True,
                    "delay_mean": 2.3, "delay_sd": 1.2,
                    "fast_track": {"delay_mean": 0.5, "delay_sd": 0.25, "share": 0.65, "cd4_below": 200.0},
                    "ineligible_fallback": {"delay_mean": 14.0, "delay_sd": 11.0},
                }},
                "model_overrides": {},
                "note": "1 month slot = empirical status quo (A8 pathways) with the ineligible wait capped by "
                        "the observed fallback distribution: the most realistic current-practice comparator.",
            },
            "criteria_retained_no_wait_plus_observed_tail": {
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {
                    "0 month": {"immediate": True, "apply_criteria": True,
                                "ineligible_fallback": {"delay_mean": 14.0, "delay_sd": 11.0}},
                    "1 month": {"immediate": False, "apply_criteria": True, "delay_mean": 1.0, "delay_sd": 0.5,
                                "ineligible_fallback": {"delay_mean": 14.0, "delay_sd": 11.0}},
                },
                "model_overrides": {},
                "note": "Factorial cell B under the observed ineligible wait: 0 month slot = criteria retained, "
                        "waiting interval removed; 1 month slot = statutory-minimum comparator; both with the "
                        "fallback. Contrast = waiting-interval effect under observed assumptions; the "
                        "eligibility effect is one_month_plus_observed_tail minus this contrast.",
            },
        },
    },
    "A10_conservative_combined": {
        "description": "Conservative estimate: the realistic current-practice comparator (empirical pathways "
                       "plus observed ineligible wait) combined with a 30% reduction in transmission from "
                       "diagnosed, unsuppressed people (plausible lower bound of the benefit).",
        "settings": {
            "empirical_tail_behaviour_0.7": {
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {"1 month": {
                    "immediate": False, "apply_criteria": True,
                    "delay_mean": 2.3, "delay_sd": 1.2,
                    "fast_track": {"delay_mean": 0.5, "delay_sd": 0.25, "share": 0.65, "cd4_below": 200.0},
                    "ineligible_fallback": {"delay_mean": 14.0, "delay_sd": 11.0},
                }},
                "model_overrides": {"DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR": 0.7},
                "note": "A9 empirical_pathways_plus_observed_tail with DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR = 0.7.",
            },
        },
    },
    "A7_delayed_start": {
        "description": "Immediate ART implemented later than 2024; before the start year the "
                       "0 month slot follows the current pathway (cost of inaction).",
        "settings": {
            "start_2027": {
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {"0 month": {"policy_start_year": 2027}},
                "model_overrides": {},
                "note": "Immediate ART from calendar year 2027 (three-year implementation delay).",
            },
            "start_2030": {
                "scenarios": ["0 month", "1 month"],
                "scenario_overrides": {"0 month": {"policy_start_year": 2030}},
                "model_overrides": {},
                "note": "Immediate ART from calendar year 2030 (six-year implementation delay).",
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
                        "ART_START_VL_UPDATE_NEXT_MONTH", "AGE_MORTALITY_OFFSET"]
_BASE_MODEL_VALUES = {k: getattr(M, k) for k in _MODEL_OVERRIDE_KEYS}


def apply_setting(setting):
    """Reset the model to its defaults, then apply this setting's overrides."""
    M.SCENARIO_CONFIGS = json.loads(json.dumps(_BASE_SCENARIO_CONFIGS))
    for k, v in _BASE_MODEL_VALUES.items():
        setattr(M, k, v)
    for label, over in setting.get("scenario_overrides", {}).items():
        M.SCENARIO_CONFIGS[label].update(over)
    for k, v in setting.get("model_overrides", {}).items():
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


def run_one_seed(seed, setting, setting_dir):
    seed_dir = os.path.join(setting_dir, f"seed_{seed}")
    os.makedirs(seed_dir, exist_ok=True)

    apply_setting(setting)
    M.SEED = seed
    M.OUTPUT_DIR = seed_dir
    burnin_pop, _ = M.run_burnin(seed=seed)
    scenarios = M.get_enabled_scenarios()
    all_results, all_pops = {}, {}
    if PARALLEL_SCENARIOS and len(scenarios) >= 2:
        worker_args = [(sc, burnin_pop, seed, setting, seed_dir) for sc in scenarios]
        with ProcessPoolExecutor(max_workers=len(scenarios)) as ex:
            for sc, res, pop in ex.map(_scenario_worker, worker_args):
                all_results[sc] = res
                all_pops[sc] = pop
    else:
        for sc in scenarios:
            res, pop = M.run_simulation(label=sc, initial_pop=burnin_pop, seed=seed)
            all_results[sc] = res
            all_pops[sc] = pop
    apply_setting(setting)          # the parent process state is used by save_all_outputs
    M.SEED = seed
    M.OUTPUT_DIR = seed_dir
    M.save_all_outputs(all_results, all_pops, scenarios)
    plt.close("all")

    del burnin_pop, all_results, all_pops
    gc.collect()

    if RUN_CEA:
        _run_cea(seed_dir)
    if RUN_BREAKEVEN:
        _run_breakeven(seed_dir)


def allocate_seeds(analysis, terminal, setting=None):
    """First k seeds of this terminal's list, k proportional to the terminal's share of
    the 100-seed split (largest-remainder rounding so that the union has the target size)."""
    target = SEED_TARGET_BY_SETTING.get((analysis, setting)) if setting is not None else None
    if target is None:
        target = SEED_TARGET_BY_ANALYSIS.get(analysis)
    if target is None:
        target = max([v for (a, _), v in SEED_TARGET_BY_SETTING.items() if a == analysis] or [100])
    target = int(target)
    names = list(SEED_SETS)
    sizes = {n: len(SEED_SETS[n]) for n in names}
    total = sum(sizes.values())
    raw = {n: target * sizes[n] / total for n in names}
    k = {n: int(raw[n]) for n in names}
    for n in sorted(names, key=lambda n: raw[n] - k[n], reverse=True)[: target - sum(k.values())]:
        k[n] += 1
    k = {n: min(k[n], sizes[n]) for n in names}
    return list(SEED_SETS[terminal][: k[terminal]]), k


def _seed_done(seed_dir):
    if not os.path.exists(os.path.join(seed_dir, "heterogeneity_summary.csv")):
        return False
    if RUN_CEA and not os.path.exists(os.path.join(seed_dir, "cea", "cea_report_ready.json")):
        return False
    if RUN_BREAKEVEN and not os.path.exists(os.path.join(seed_dir, "breakeven", "breakeven_report_ready.json")):
        return False
    return True


def _write_setting_manifest(setting_name, setting, setting_dir):
    manifest = {
        "analysis": ANALYSIS,
        "setting": setting_name,
        "description": ANALYSES[ANALYSIS]["description"],
        "note": setting.get("note", ""),
        "scenarios": setting["scenarios"],
        "scenario_overrides": setting.get("scenario_overrides", {}),
        "model_overrides": {k: (v if not isinstance(v, dict) else {str(a): b for a, b in v.items()})
                            for k, v in setting.get("model_overrides", {}).items()},
        "model_defaults": {k: (v if not isinstance(v, dict) else None) for k, v in _BASE_MODEL_VALUES.items()},
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
    }
    with open(os.path.join(setting_dir, "setting_config.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)


def main():
    if ANALYSIS not in ANALYSES:
        raise ValueError(f"ANALYSIS must be one of: {list(ANALYSES)}")
    if TERMINAL not in SEED_SETS:
        raise ValueError("TERMINAL must be one of: " + ", ".join(SEED_SETS))
    def _seeds_for(setting_name):
        seeds, alloc = allocate_seeds(ANALYSIS, TERMINAL, setting_name)
        if SEED_ORDER == "reverse":
            seeds = list(reversed(seeds))
        elif SEED_ORDER == "shuffle":
            import random as _random
            _random.Random(12345).shuffle(seeds)
        elif SEED_ORDER != "forward":
            raise ValueError('SEED_ORDER must be "forward", "reverse", or "shuffle"')
        if SEED_LIMIT is not None:
            seeds = seeds[:int(SEED_LIMIT)]
        return seeds, alloc

    settings = ANALYSES[ANALYSIS]["settings"]
    names = list(settings) if SETTINGS_TO_RUN is None else list(SETTINGS_TO_RUN)
    for nm in names:
        if nm not in settings:
            raise ValueError(f"Unknown setting '{nm}' for {ANALYSIS}; choose from {list(settings)}")

    print("=" * 70)
    print(f"{ANALYSIS} [{TERMINAL}]: settings={names}")
    print(ANALYSES[ANALYSIS]["description"])
    print("=" * 70, flush=True)
    t0 = time.time()
    for nm in names:
        setting = settings[nm]
        seeds, alloc = _seeds_for(nm)
        setting_dir = os.path.join(OUT_ROOT, ANALYSIS, nm)
        os.makedirs(setting_dir, exist_ok=True)
        _write_setting_manifest(nm, setting, setting_dir)
        print(f"\n### SETTING {nm}: {setting.get('note','')}\n    -> {setting_dir}", flush=True)
        print(f"    seeds on this terminal: {len(seeds)}; allocation across terminals: {alloc}", flush=True)
        for i, seed in enumerate(seeds, 1):
            seed_dir = os.path.join(setting_dir, f"seed_{seed}")
            if _seed_done(seed_dir):
                print(f"  [{i}/{len(seeds)}] {nm} SEED={seed} skip (already done)", flush=True)
                continue
            ts = time.time()
            print(f"  [{i}/{len(seeds)}] {nm} SEED={seed} start {datetime.datetime.now():%H:%M:%S}", flush=True)
            try:
                run_one_seed(seed, setting, setting_dir)
                print(f"     done in {(time.time()-ts)/60:.1f} min", flush=True)
            except Exception as e:
                import traceback; traceback.print_exc()
                print(f"     [ERROR] {nm} seed={seed}: {e}", flush=True)
    print("\n" + "=" * 70)
    print(f"TOTAL: {(time.time()-t0)/60:.1f} min")
    print("Next: run summarize_revision_scenarios.py to build the summary tables.")
    print("=" * 70, flush=True)


if __name__ == "__main__":
    main()

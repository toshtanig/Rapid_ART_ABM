#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Parameter uncertainty and budget impact on seed-averaged inputs (PSA, DSA, CEAC, break-even).

For each comparison set (a setting folder of the runs):

Step "mean_input": the yearly summaries of the scenario slots (scenario0/1/2_yearly_summary.csv) are averaged
    over the completed seeds of the setting, exactly as in cea_param_uncertainty.py, and written to
    post/economics/<set>/mean_input/ together with the list of seeds used.

Step "economics": cea_analysis.py (probabilistic sensitivity analysis with the fixed random seed of that
    script, one-way sensitivity analysis, acceptability curve) and breakeven_analysis_detailed_json.py
    (cumulative discounted cost difference under the budget-impact cost convention, break-even year) are run
    on the averaged input; a short table of the undiscounted cost difference in the first five years is added.
    This step needs scipy.

The two steps can be run on different machines: the averaged input is a few hundred kilobytes.

TARGET = "R2" uses the runs of this folder. TARGET = "published" applies the same steps to the published
outputs (read only) to verify that the published numbers are reproduced; output in post/validation_published/.

Run as-is in VS Code; all settings are in the CONFIG block, no command-line arguments.
"""

import os
import sys
import glob
import json
import time
from io import StringIO

import pandas as pd

# =============================================================================
# CONFIG
# =============================================================================
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _paths import BASE_144, BASE_MULTISEED, REVISION_R1_RUNS, REVISION_RUNS

TARGET = "R2"                         # "R2" or "published"
STEPS = ["mean_input", "economics"]   # which steps to run
SETS_TO_RUN = None                    # None = every set whose folder exists; or a list of set names
MIN_SEEDS = 30                        # a set with fewer completed seeds is skipped
PSA_N_ITERATIONS = 1000
MAKE_FIGURES = False                  # figures of cea_analysis / break-even (the manuscript figures are drawn separately)

R1_RUNS = REVISION_R1_RUNS
# name of the set -> folder of the setting (0 month slot = immediate ART; 1 month slot = the comparator)
SETS = {
    "R2": {
        "statutory_minimum": os.path.join(REVISION_RUNS, "B0_core", "default"),
        "observed_practice": os.path.join(REVISION_RUNS, "A9_observed_ineligible_wait", "empirical_pathways_plus_observed_tail"),
        "conservative_0p3": os.path.join(REVISION_RUNS, "A10_conservative_combined", "empirical_tail_behaviour_0.3"),
        "observed_practice_0p5": os.path.join(REVISION_RUNS, "A10_conservative_combined", "empirical_tail_behaviour_0.5"),
        "statutory_observed_tail": os.path.join(REVISION_RUNS, "A9_observed_ineligible_wait", "one_month_plus_observed_tail"),
        "statutory_0p5": os.path.join(REVISION_RUNS, "A5_post_diagnosis_behaviour", "factor_0.5"),
        "statutory_0p3": os.path.join(REVISION_RUNS, "A5_post_diagnosis_behaviour", "factor_0.3"),
    },
    "published": {
        "statutory_minimum": BASE_MULTISEED,
        "observed_practice": os.path.join(R1_RUNS, "A9_observed_ineligible_wait", "empirical_pathways_plus_observed_tail"),
        "conservative_0p7": os.path.join(R1_RUNS, "A10_conservative_combined", "empirical_tail_behaviour_0.7"),
        "statutory_observed_tail": os.path.join(R1_RUNS, "A9_observed_ineligible_wait", "one_month_plus_observed_tail"),
    },
}
OUT_ROOT = {"R2": os.path.join(REVISION_RUNS, "post", "economics"),
            "published": os.path.join(REVISION_RUNS, "post", "validation_published", "economics")}
# =============================================================================

SCEN_FILES = ["scenario0_yearly_summary.csv", "scenario1_yearly_summary.csv", "scenario2_yearly_summary.csv"]
GROUP_KEYS = ["Year", "AgeGroup", "CD4Category", "OnART", "OnPrEP"]
NUM_COLS = ["PersonMonthsAlive", "NewlyInfected", "NewlyDiagnosed", "Deaths"]


def _read_csv(path):
    try:
        return pd.read_csv(path)
    except Exception:
        with open(path, "rb") as f:
            data = f.read().replace(b"\x00", b"")
        return pd.read_csv(StringIO(data.decode("utf-8", "ignore")))


def _complete(seed_dir):
    return (os.path.exists(os.path.join(seed_dir, "heterogeneity_summary.csv"))
            and os.path.exists(os.path.join(seed_dir, "cea", "cea_report_ready.json"))
            and os.path.exists(os.path.join(seed_dir, "breakeven", "breakeven_report_ready.json")))


def sets_of_target():
    names = list(SETS[TARGET]) if SETS_TO_RUN is None else list(SETS_TO_RUN)
    return [(n, SETS[TARGET][n]) for n in names]


def build_mean_input(name, setting_dir):
    """Average the yearly summaries over the completed seeds (sum over seeds divided by the number of seeds,
    so that a stratum absent from a seed counts as zero), as cea_param_uncertainty.build_mean_input."""
    out_dir = os.path.join(OUT_ROOT[TARGET], name, "mean_input")
    seed_dirs = sorted((d for d in glob.glob(os.path.join(setting_dir, "seed_*")) if _complete(d)),
                       key=lambda d: int(os.path.basename(d).replace("seed_", "")))
    if len(seed_dirs) < MIN_SEEDS:
        print(f"  [{name}] {len(seed_dirs)} completed seeds (fewer than {MIN_SEEDS}); skipped")
        return 0
    os.makedirs(out_dir, exist_ok=True)
    files = []
    for scen_file in SCEN_FILES:
        frames = [_read_csv(os.path.join(d, scen_file)) for d in seed_dirs if os.path.exists(os.path.join(d, scen_file))]
        if not frames:
            stale = os.path.join(out_dir, scen_file)
            if os.path.exists(stale):
                os.remove(stale)
            continue
        if len(frames) != len(seed_dirs):
            raise RuntimeError(f"{name}: {scen_file} is present in {len(frames)} of {len(seed_dirs)} completed seeds")
        allc = pd.concat(frames, ignore_index=True)
        summed = allc.groupby(GROUP_KEYS, as_index=False)[NUM_COLS].sum()
        summed[NUM_COLS] = summed[NUM_COLS] / len(frames)
        summed.to_csv(os.path.join(out_dir, scen_file), index=False)
        files.append(scen_file)
    seeds = [int(os.path.basename(d).replace("seed_", "")) for d in seed_dirs]
    with open(os.path.join(out_dir, "mean_input_info.json"), "w", encoding="utf-8") as f:
        json.dump({"set": name, "setting_folder": os.path.relpath(setting_dir, BASE_144), "n_seeds": len(seeds),
                   "seeds": seeds, "files": files, "written": time.strftime("%Y-%m-%d %H:%M:%S")}, f, indent=1)
    print(f"  [{name}] averaged {files} over {len(seeds)} seeds", flush=True)
    return len(seeds)


def run_economics(name):
    set_dir = os.path.join(OUT_ROOT[TARGET], name)
    mean_dir = os.path.join(set_dir, "mean_input")
    if not os.path.exists(os.path.join(mean_dir, "scenario0_yearly_summary.csv")):
        print(f"  [{name}] no averaged input; skipped")
        return False
    import matplotlib
    matplotlib.use("Agg")
    import cea_analysis as CEA
    import breakeven_analysis_detailed_json as BE
    CEA.Config.INPUT_DIR = mean_dir
    CEA.Config.OUTPUT_DIR = set_dir
    CEA.Config.MAKE_FIGURES = bool(MAKE_FIGURES)
    CEA.Config.PSA_N_ITERATIONS = int(PSA_N_ITERATIONS)
    CEA.main()
    BE.INPUT_DIR = mean_dir
    BE.OUTPUT_DIR = set_dir
    if hasattr(BE, "MAKE_FIGURES"):
        BE.MAKE_FIGURES = bool(MAKE_FIGURES)
    BE.main()
    short_term_table(set_dir)
    print(f"  [{name}] PSA, DSA, CEAC and break-even written to {set_dir}", flush=True)
    return True


def short_term_table(set_dir, rate=0.02, years=5):
    """Undiscounted annual cost difference in the first years (budget-impact convention), from the discounted
    yearly differences of the break-even analysis (discount factor (1 + rate) ** -(year - 1))."""
    p = os.path.join(set_dir, "yearly_cost_differences.csv")
    if not os.path.exists(p):
        return
    df = pd.read_csv(p)
    out = {"Year": df["Year"].tolist()[:years]}
    for c in df.columns:
        if c.endswith("_annual_diff_B"):
            und = (df[c] * (1.0 + rate) ** (df["Year"] - 1)).tolist()[:years]
            out[c.replace("_annual_diff_B", " undiscounted annual difference (billion yen)")] = und
    t = pd.DataFrame(out)
    tot = {"Year": f"years 1 to {years}"}
    for c in t.columns[1:]:
        tot[c] = float(t[c].sum())
    pd.concat([t, pd.DataFrame([tot])], ignore_index=True).to_csv(os.path.join(set_dir, "short_term_budget_impact_undiscounted.csv"), index=False)


def main():
    if TARGET not in SETS:
        raise SystemExit('TARGET must be "R2" or "published"')
    for name, setting_dir in sets_of_target():
        if "mean_input" in STEPS:
            if os.path.isdir(setting_dir):
                build_mean_input(name, setting_dir)
            else:
                print(f"  [{name}] setting folder not found; averaged input not rebuilt")
        if "economics" in STEPS:
            run_economics(name)


if __name__ == "__main__":
    main()

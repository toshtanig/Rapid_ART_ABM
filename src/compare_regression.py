#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Regression check: compares one seed of A0_regression_check (modified model,
default configuration) with the same seed of the published run
(outputs/core_multiseed). Every shared CSV must be identical; any difference means
the code changes altered the default behaviour or the random-number stream.

Run as-is in VS Code after run_revision_scenarios.py with ANALYSIS =
"A0_regression_check" and SEED_LIMIT = 1.
"""

import os
import sys
import glob
import pandas as pd

# =============================================================================
# CONFIG
# =============================================================================
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import BASE_MULTISEED, REVISION_RUNS

SEED = 14
NEW_DIR = os.path.join(REVISION_RUNS, "A0_regression_check", "default", f"seed_{SEED}")
REF_DIR = os.path.join(BASE_MULTISEED, f"seed_{SEED}")

# Files that must match exactly (model outputs). CEA/break-even JSONs are
# derived from these and are checked through their CSV tables.
PATTERNS = ["scenario*_yearly_summary.csv", "heterogeneity_summary.csv",
            "scenario*_incident_infection_cohort_summary.csv",
            "scenario*_incident_diagnosis_cohort_summary.csv",
            "scenario*_diag_to_art_delay_hist.csv", "scenario*_eligible_to_art_delay_hist.csv",
            "diagnosis_to_ART_delay_summary.csv", "eligible_to_ART_delay_summary.csv",
            "burnin_summary.csv", "monthly_art_starts.csv", "monthly_tests.csv",
            "cea/base_case_results.csv", "cea/incremental_results.csv",
            "breakeven/*.csv"]
# =============================================================================


def _norm_cols(df):
    return [str(c).replace(" discounted", "") for c in df.columns]


def _frames_equal(a, b):
    if a.shape != b.shape or list(a.columns) != list(b.columns):
        return False
    num = a.select_dtypes("number")
    if len(num.columns) and (num - b[num.columns]).abs().max().max() > 1e-9:
        return False
    obj = a.select_dtypes(exclude="number")
    return obj.equals(b[obj.columns]) if len(obj.columns) else True


def main():
    if not os.path.isdir(NEW_DIR):
        raise SystemExit(f"New run not found: {NEW_DIR}")
    if not os.path.isdir(REF_DIR):
        raise SystemExit(f"Reference run not found: {REF_DIR}")
    n_ok = n_diff = n_missing = 0
    for pat in PATTERNS:
        for ref in sorted(glob.glob(os.path.join(REF_DIR, pat))):
            rel = os.path.relpath(ref, REF_DIR)
            new = os.path.join(NEW_DIR, rel)
            if not os.path.exists(new):
                print(f"  MISSING  {rel}")
                n_missing += 1
                continue
            header_only = False
            try:
                a = pd.read_csv(ref)
                b = pd.read_csv(new)
                same = _frames_equal(a, b)
                if not same and a.shape == b.shape and _norm_cols(a) == _norm_cols(b):
                    # Known difference: the published post-processing renamed the break-even columns
                    # ("Year n cumulative (B JPY)" -> "Year n cumulative discounted (B JPY)") after the
                    # reference seeds were generated. Compare the values with the headers normalised.
                    b2 = b.copy(); b2.columns = a.columns
                    header_only = _frames_equal(a, b2)
                    same = header_only
            except Exception as e:
                same = open(ref, "rb").read() == open(new, "rb").read()
            if same and header_only:
                n_ok += 1
                print(f"  ok       {rel}   (values identical; column headers differ only by the known 'discounted' wording)")
            elif same:
                n_ok += 1
                print(f"  ok       {rel}")
            else:
                n_diff += 1
                print(f"  DIFFERS  {rel}")
    print("\n" + "=" * 60)
    print(f"identical: {n_ok}   different: {n_diff}   missing: {n_missing}")
    print("REGRESSION CHECK PASSED" if (n_diff == 0 and n_missing == 0)
          else "REGRESSION CHECK FAILED: do not launch the analyses until this is resolved")
    print("=" * 60)


if __name__ == "__main__":
    main()

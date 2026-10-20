#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Output locations used by the revision (R1) scripts.

By default every output lives under <repository>/outputs:
    outputs/core_multiseed/        published 200-seed run (src/run_multiseed.py)
    outputs/revision_R1_runs/      revision analyses (src/run_revision_scenarios.py)
    outputs/figures_R1/            revision figures (src/make_figure5.py,
                                   src/make_supp_figure_S1_trends.py)
Set OUTPUT_ROOT_OVERRIDE to an absolute path to keep the outputs elsewhere.
"""
import os

OUTPUT_ROOT_OVERRIDE = None

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_ROOT = OUTPUT_ROOT_OVERRIDE or os.path.join(REPO_ROOT, "outputs")
BASE_MULTISEED = os.path.join(OUTPUT_ROOT, "core_multiseed")
REVISION_RUNS = os.path.join(OUTPUT_ROOT, "revision_R1_runs")
FIGURES_R1 = os.path.join(OUTPUT_ROOT, "figures_R1")

# Aliases used by some of the revision scripts.
BASE_144 = OUTPUT_ROOT
PROJ_169 = OUTPUT_ROOT

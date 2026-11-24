#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Output locations used by the scripts of versions 1.1.0 and 1.2.0.

By default every output lives under <repository>/outputs:
    outputs/revision_R2_runs/      analyses of version 1.2.0 (src/run_revision_scenarios.py) and their
                                   post-processing (post/summary, post/tables, post/economics)
    outputs/figures_R2/            figures of version 1.2.0 (src/make_figures_R2.py)
    outputs/core_multiseed/        200-seed run of versions 1.0.0 and 1.1.0 (optional, read only: used by the
                                   regression check and by TARGET = "published" of the version 1.2.0 scripts)
    outputs/revision_R1_runs/      analyses of version 1.1.0 (optional, read only: as above)
    outputs/figures_R1/            figures of version 1.1.0 (src/make_figure5.py)
Set OUTPUT_ROOT_OVERRIDE to an absolute path to keep the outputs elsewhere.
"""
import os

OUTPUT_ROOT_OVERRIDE = None

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_ROOT = OUTPUT_ROOT_OVERRIDE or os.path.join(REPO_ROOT, "outputs")
REVISION_RUNS = os.path.join(OUTPUT_ROOT, "revision_R2_runs")
FIGURES_R2 = os.path.join(OUTPUT_ROOT, "figures_R2")
BASE_MULTISEED = os.path.join(OUTPUT_ROOT, "core_multiseed")
REVISION_R1_RUNS = os.path.join(OUTPUT_ROOT, "revision_R1_runs")
FIGURES_R1 = os.path.join(OUTPUT_ROOT, "figures_R1")

# Aliases used by some of the scripts.
BASE_144 = OUTPUT_ROOT
PROJ_169 = OUTPUT_ROOT

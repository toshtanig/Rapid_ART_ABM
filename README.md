# Rapid ART cost-effectiveness microsimulation (HIV, Japan)

Source code and reproducibility materials for:

> Taniguchi T, Imahashi M, Sato D, Noda T. Cost-effectiveness and budget impact
> of immediate versus delayed antiretroviral therapy for HIV in Japan: an
> agent-based modelling study. *The Lancet Regional Health - Western Pacific*
> (under review).

A stochastic individual-based dynamic transmission microsimulation of HIV in
Japan, linked to a cost-effectiveness and budget-impact analysis of total
HIV-related healthcare costs (healthcare payer perspective). It evaluates
reform of the statutory immune-dysfunction certification pathway to permit
immediate ART at diagnosis, against a comparator calibrated to observed
practice in a multicentre cohort, a prespecified statutory-minimum comparator
(one-month delay for eligible people; people not meeting the CD4/HIV-RNA
criteria remain untreated until laboratory progression), and a two-month
comparator, over a 40-year horizon after a 20-year burn-in. HIV-related
mortality is represented by a CD4-stratum by care-state standardised mortality
ratio applied to the age- and sex-specific general-population background
hazard. The full model specification (ODD protocol) and the principal
parameter tables are in the article's Supplementary Material.

## Versions

- **1.2.0 (November 2026)**: one order of events for ART starts within the
  monthly cycle in all strategies and in all periods, with the calibration
  re-derived; every analysis of the article was run with this version (see
  "Version 1.2.0" below and `CHANGELOG.md`). Results differ from versions 1.0.0
  and 1.1.0, which the model reproduces exactly when the earlier settings are
  restored. Archived on Zenodo
  (https://doi.org/10.5281/zenodo.24494835).
- **1.1.0 (October 2026)**: scenario and sensitivity analyses (comparators
  calibrated to observed practice, factorial decomposition, structural
  sensitivity analyses). Archived on Zenodo
  (https://doi.org/10.5281/zenodo.23437249).
- **1.0.0 (August 2026)**: first release (base case and principal sensitivity
  analyses), archived on Zenodo (https://doi.org/10.5281/zenodo.21931762).

## Layout

```
├── requirements.txt / environment.yml      pinned dependencies
├── LICENSE / LICENSE-DATA.txt              MIT (code) / CC BY 4.0 (data)
├── CITATION.cff / CHANGELOG.md             citation metadata / version history
├── data/                                   aggregated model outputs (CC BY 4.0)
│   ├── revision_R2/                        outputs of the version 1.2.0 analyses (the article)
│   ├── revision_R1/                        outputs of the version 1.1.0 analyses
│   └── (files in data/ itself)             outputs of version 1.0.0
└── src/
    ├── hiv_rapid_art_model.py             the microsimulation model (version 1.2.0)
    ├── _paths.py                          output locations
    ├── run_revision_scenarios.py          driver: queue of the main and supplementary analyses
    ├── run_supervised.py                  starts the driver again after an abnormal end
    ├── check_progress.py                  progress of the queue
    ├── run_calibration_check.py / summarize_calibration_check.py
    │                                      burn-in only: calibration indicators and fit
    ├── compare_regression.py              regression check against version 1.1.0
    ├── summarize_R2.py                    per-seed cache and contrast tables
    ├── core_tables_R2.py                  tables of the Supplementary Material
    ├── economics_R2.py                    PSA, one-way DSA, CEAC and break-even on seed-averaged inputs
    ├── make_figures_R2.py                 Figures 2 to 5 and Supplementary Figures S2 to S8
    ├── cea_analysis.py / breakeven_analysis_detailed_json.py
    │                                      cost-effectiveness and break-even routines used by the driver
    │                                      and by economics_R2.py
    ├── make_figure1.py / make_supp_figure_S1_trends.py / check_timing_realisations.py
    └── scripts of versions 1.0.0 and 1.1.0 (see "Earlier versions")
```

## Dependencies

```
conda env create -f environment.yml        # or: pip install -r requirements.txt
```

Published versions: Python 3.11.15, NumPy 2.4.3, SciPy 1.17.1, pandas 3.0.1,
Numba 0.64.0, llvmlite 0.46.0, Matplotlib 3.10.8, seaborn 0.13.2, Pillow 12.1.1.
(python-docx is only for the optional Word build pipeline, not the model;
openpyxl is needed only to redraw Supplementary Figure S1 from the
annual-report Excel tables.)

## Version 1.2.0: the analyses of the article

Each script sets its parameters in a `CONFIG` block at the top and takes no
command-line arguments. Outputs are written under `outputs/` (see
`src/_paths.py`).

**Model changes** (`src/hiv_rapid_art_model.py`).

- Order of events for ART starts (`ART_START_VL_UPDATE_NEXT_MONTH = True`,
  `ART_START_VL_UPDATE_NEXT_MONTH_IN_BURNIN = True`). For every ART start the
  first on-ART person-month is the cycle in which the start is processed, and
  the first viral-load reduction is applied in the following cycle. This holds
  for starts at diagnosis and for starts at a scheduled start month (eligible
  people after their wait, and restarts after re-engagement), in the burn-in as
  well as in the analytic period. The difference in the onset of viral-load
  decline between immediate ART and a wait of k months is therefore k months,
  the same as the difference in on-ART person-time. In versions 1.0.0 and 1.1.0
  a scheduled starter received the first reduction in the start cycle itself,
  so that the difference was k - 1 months.
- Calibration scalar on the force of infection
  (`MORTALITY_RECALIB_FOI_SCALE = 1.02`, previously 1.05), re-derived under the
  new order so that the state of the epidemic at the end of the burn-in (2023)
  reproduces the calibration targets (`src/run_calibration_check.py`,
  `src/summarize_calibration_check.py`).
- With `ART_START_VL_UPDATE_NEXT_MONTH = False` and
  `MORTALITY_RECALIB_FOI_SCALE = 1.05` every output file is identical to
  versions 1.0.0 and 1.1.0 for the same seed (analysis `A0_regression_check` of
  the driver, then `src/compare_regression.py`; the reference outputs of
  release 1.1.0 are expected in `outputs/core_multiseed/`).

**Running the simulations.** `src/run_revision_scenarios.py` runs a queue of
(analysis, setting, number of seeds) items from top to bottom: `B0_core` (main
analysis: immediate ART, one-month and two-month delay, 200 seeds),
`A9_observed_ineligible_wait` (observed-practice comparator and the
comparators with the observed wait of people not meeting the criteria),
`A1_factorial`, `A5_post_diagnosis_behaviour` and `A10_conservative_combined`
(relative transmission probability of diagnosed, unsuppressed people 0.5 and
0.3, with the calibration scalar re-derived for each value: 1.12 and 1.16),
`A8_status_quo_empirical`, `A7_delayed_start`, `A2_foi_trajectory`,
`A4_time_to_suppression`, `A12_age_offset`, and `A13_foi_multiplier`. For every
seed the burn-in is run once and the scenario slots of the setting are run
from it (common random numbers), followed by the cost-effectiveness and
break-even post-processing. Outputs are written to
`outputs/revision_R2_runs/<analysis>/<setting>/seed_<n>/`. The queue can be
shared by several machines through a synchronised output folder: set
`TERMINAL` to one of the labels `node1` to `node5` on each machine (on a single
machine, any one label runs the whole queue, because a machine that has
finished its own seeds takes the seeds that no other machine has started).
Completed seeds are skipped, so a run can be stopped and started again.
`src/run_supervised.py` starts the driver again after an abnormal end, and
`src/check_progress.py` reports progress and the expected finish time.

**Post-processing.** Run in this order:

1. `src/summarize_R2.py` reads every completed seed once into a per-seed cache
   and writes the contrast tables (`post/summary/`): infections and deaths
   averted, deaths in the whole population, incremental QALYs by the
   prespecified prevention decomposition and by the whole-population direct
   estimate, incremental costs, net monetary benefit, dominance, the factorial
   decompositions, and the discount-rate and time-horizon sensitivity.
2. `src/core_tables_R2.py` writes the tables of the Supplementary Material from
   the cache (`post/tables/`): calibration and fit, risk groups, transmission
   sources, Monte Carlo properties, waiting times, break-even curves, mortality
   validation, survival-weight sensitivity, yearly series for the figures, and
   the comparison of the two orders of events on the same seeds.
3. `src/economics_R2.py` runs the probabilistic and one-way sensitivity
   analyses, the acceptability curve and the break-even analysis on the
   seed-averaged yearly outputs of each comparison (`post/economics/<set>/`).
4. `src/make_figures_R2.py` draws Figures 2 to 5 and Supplementary Figures S2
   to S8 from these outputs (`outputs/figures_R2/`). The figures can also be
   drawn without running the model, from the copies of the post-processing
   outputs in `data/revision_R2/` (set `POST` in its CONFIG block).

Setting `TARGET = "published"` in steps 1 to 4 applies the same code to the
outputs of release 1.1.0 (`outputs/core_multiseed/` and
`outputs/revision_R1_runs/`, if present), as a check of the scripts against the
earlier published numbers.

**Other figures and tables.** `src/make_figure1.py` draws Figure 1.
`src/make_supp_figure_S1_trends.py` draws Supplementary Figure S1 from the
national HIV and AIDS annual report (AIDS Surveillance Committee, Ministry of
Health, Labour and Welfare, 2024;
https://api-net.jfap.or.jp/status/japan/nenpo.html), or from the bundled CSV of
the plotted counts if the Excel tables are not present.
`src/check_timing_realisations.py` reproduces Supplementary Tables S1e and
S10a (the time to viral suppression realised under the implemented viral-load
rule, and the waiting times realised at monthly resolution) by direct
calculation, independent of the population model.

Reproducibility rests on the fixed 200-seed list (`PUBLICATION_SEEDS` in the
driver) and a common-random-numbers design: within each replicate the
ART-timing scenarios share re-seeded random streams, so incremental outcomes
are paired per-seed differences and between-scenario noise is minimised.

## Earlier versions

The scripts of versions 1.0.0 and 1.1.0 remain in `src/`. With the model of
this release they run under the new order of events, so they do not reproduce
the outputs of the earlier releases unless the earlier settings are restored
(see "Model changes"); the archived releases 1.0.0 and 1.1.0 contain the model
as it was. `run_multiseed.py` (single-machine driver of 200 seeds) and
`run_foi_sensitivity.py` are superseded by the queue of
`run_revision_scenarios.py`. `summarize_revision_scenarios.py`,
`make_figure5.py`, `make_figures_3_4_S4_S5.py`, `make_figures_S3_S7_S8.py`,
`survival_weight_sensitivity_per_seed.py`, `multiseed_figures.py`, and the
reporting scripts of version 1.0.0 read the outputs of those versions
(`outputs/core_multiseed/`, `outputs/revision_R1_runs/`) and are superseded by
the four post-processing scripts above.

## Data and licensing

`data/revision_R2/` contains the aggregated outputs of the version 1.2.0
analyses, in the layout of `outputs/revision_R2_runs/post/`:

- `summary/summary.md` and `summary/summary_contrasts.csv`: every setting and
  contrast (mean and 95% confidence interval across seeds; the intervals and
  replicate counts describe simulation variability with the epidemiological
  parameters fixed at their calibrated values);
  `summary/discount_rate_sensitivity_direct.csv` and
  `summary/horizon_sensitivity_direct.csv`;
- `summary/per_seed/`: the per-seed incremental outcomes behind each setting,
  the main analysis restricted to the same seeds, and the per-seed factorial
  contrasts;
- `tables/`: the values behind the supplementary tables and figures (one CSV
  per table and `tables.md` as an overview);
- `economics/<set>/`: probabilistic, one-way, acceptability-curve and
  break-even results for the statutory-minimum, observed-practice,
  conservative, and statutory-minimum-with-observed-wait comparisons, with the
  seed-averaged yearly inputs (`mean_input/`).

`data/revision_R1/` and the files in `data/` itself are the outputs of
versions 1.1.0 and 1.0.0 (earlier order of events and calibration scalar).
They are kept because the article reports the results under the earlier order
as a sensitivity analysis (`data/revision_R2/tables/convention_comparison.csv`
compares the two on the same seeds); they are described in the README of
release 1.1.0.

The scripts produce aggregated CSV/JSON summaries and figures; no individual
participant data are used or generated (all calibration targets are public
aggregate surveillance and care-cascade estimates). Code is released under the
MIT License; aggregated outputs under CC BY 4.0. Please cite the article and the
archived Zenodo release (https://doi.org/10.5281/zenodo.24494835; see `CITATION.cff`).

## Funding and contact

Health and Labour Sciences Research Grant, AIDS Control Policy Research Program,
MHLW Japan (21HB1003, 23HB1001, 26HB1001). Toshibumi Taniguchi, Department of
Infectious Diseases, Chiba University Hospital. E-mail: tosh-tanig@chiba-u.jp

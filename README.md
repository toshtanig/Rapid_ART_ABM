# Rapid ART cost-effectiveness microsimulation (HIV, Japan)

Source code and reproducibility materials for:

> Taniguchi T, Imahashi M, Sato D, Noda T. Cost-effectiveness and budget impact
> of immediate versus delayed antiretroviral therapy for HIV in Japan: an
> agent-based modelling study. *The Lancet Regional Health - Western Pacific*
> (under review).

A stochastic individual-based dynamic transmission microsimulation of HIV in
Japan, linked to a cost-effectiveness and budget-impact analysis from the
public-payer perspective. It evaluates reform of the statutory
immune-dysfunction certification pathway to permit immediate ART at diagnosis,
against a prespecified statutory-minimum comparator (one-month delay for
eligible people; people not meeting the CD4/HIV-RNA criteria remain untreated
until laboratory progression), a two-month comparator, and, from version 1.1.0,
comparators calibrated to observed practice in a multicentre cohort, over a
40-year horizon after a 20-year burn-in. HIV-related mortality is represented
by a CD4-stratum by care-state standardised mortality ratio applied to the age-
and sex-specific general-population background hazard. The full model
specification (ODD protocol) and the principal parameter tables are in the
article's Supplementary Material.

## Versions

- **1.1.0 (October 2026)**: scenario and sensitivity analyses reported in the
  article (see "Version 1.1.0" below and `CHANGELOG.md`). With default settings the model reproduces
  version 1.0.0 exactly. Archived on Zenodo
  (https://doi.org/10.5281/zenodo.23437249).
- **1.0.0 (August 2026)**: first release (base case and principal sensitivity
  analyses), archived on Zenodo (https://doi.org/10.5281/zenodo.21931762).

## Layout

```
├── requirements.txt / environment.yml      pinned dependencies
├── LICENSE / LICENSE-DATA.txt              MIT (code) / CC BY 4.0 (data)
├── CITATION.cff / CHANGELOG.md             citation metadata / version history
├── data/                                   aggregated model outputs (CC BY 4.0)
│   └── revision_R1/                        outputs of the version 1.1.0 analyses
└── src/
    ├── hiv_rapid_art_model.py             the microsimulation model
    ├── run_multiseed.py                   single-machine driver: 200 publication seeds
    ├── cea_analysis.py                    cost-effectiveness analysis
    ├── cea_param_uncertainty.py           PSA, one-way DSA, CEAC
    ├── cea_dominant_reporting_json.py     dominance reporting
    ├── breakeven_analysis_detailed_json.py / bia_rapid_enhanced_json.py
    ├── prevention_qaly_dominance.py / survival_adj_sensitivity.py
    ├── run_foi_sensitivity.py / foi_source_analysis.py
    ├── convergence_diagnostics.py / recompute_costs.py
    ├── make_figure1.py / multiseed_figures.py / compose_figures.py /
    │   visual_yearly_cost_json.py / lancet_export.py
    ├── _paths.py                          output locations for the 1.1.0 scripts
    ├── run_revision_scenarios.py          driver for the 1.1.0 scenario analyses
    ├── summarize_revision_scenarios.py    summary tables for the 1.1.0 analyses
    ├── compare_regression.py / check_progress.py
    └── make_figure5.py / make_supp_figure_S1_trends.py / check_timing_realisations.py
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

## Reproducing the analysis

Each script sets its parameters in a `CONFIG` block at the top and takes no
command-line arguments. Set the base/output paths in those blocks to your local
checkout, then:

1. `python src/run_multiseed.py` runs the 200 fixed seeds (`SEEDS` in CONFIG) on
   a single machine, writing one `seed_<n>` folder per seed under
   `outputs/core_multiseed/` and appending to `multiseed_comparison.csv` after
   each seed (interruptible). Shorten `SEEDS` for a quick check.
2. Run the cost-effectiveness, budget-impact, sensitivity, decomposition, and
   figure scripts on the per-seed outputs (each self-configured via its CONFIG).

Reproducibility rests on the fixed 200-seed list and a common-random-numbers
design: within each replicate the ART-timing scenarios share re-seeded random
streams, so incremental outcomes are paired per-seed differences and
between-scenario noise is minimised.

## Version 1.1.0: scenario and sensitivity analyses

**Model changes** (`src/hiv_rapid_art_model.py`). All new options are switched
off by default, and with default settings every output CSV is identical to
version 1.0.0 for the same seed (checked with `src/compare_regression.py`):
scenario options `immediate`, `apply_criteria`, and `policy_start_year`
(factorial decomposition and delayed implementation); a year-specific
analytic-period force-of-infection multiplier
(`ANALYTIC_FOI_MULTIPLIER_BY_YEAR`); the time from ART initiation to viral
suppression (`TTS_MIN_MONTHS`, `TTS_RANGE_MONTHS`); a relative infectivity of
diagnosed, virally unsuppressed people
(`DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR`); an observed-practice waiting
distribution with a shorter wait for people with CD4 below 200 (`fast_track`);
and an observed waiting time for people who do not meet the certification
criteria (`ineligible_fallback`); a timing-convention switch for scheduled ART starts
(`ART_START_VL_UPDATE_NEXT_MONTH`, default off) and a run-time-adjustable background-mortality
age offset (`AGE_MORTALITY_OFFSET`, default 5.0 as in version 1.0.0), both used in the
sensitivity analyses A11 and A12.

**Running the scenario analyses.** `src/run_revision_scenarios.py` runs one
analysis at a time (`ANALYSIS` in CONFIG): A0 regression check, A1 factorial
decomposition, A2 time-varying force of infection, A4 time to suppression, A5
post-diagnosis transmission, A7 delayed implementation, A8 and A9
observed-practice comparators, A10 the conservative combined scenario, A11 the
timing-convention sensitivity analysis, and A12 the age-offset sensitivity analysis. Each
setting runs the immediate and comparator slots from the same burn-in (common
random numbers), followed by the cost-effectiveness and break-even
post-processing. The seed list is split into five labels (`node1` to `node5`)
that were used to share the runs across workstations; on a single machine,
run each label in turn. Completed seeds are skipped, so runs can be resumed.
Outputs are written to `outputs/revision_R1_runs/<ANALYSIS>/<setting>/seed_<n>/`
(see `src/_paths.py`). `src/check_progress.py` reports progress, and
`src/summarize_revision_scenarios.py` writes the summary tables (Part 1) and the
discount-rate and 20-year-horizon post-processing of the published run
(Part 2) to `outputs/revision_R1_runs/summary/`. The probabilistic and one-way
analyses for the new comparators use `src/cea_param_uncertainty.py` on the
seed-averaged yearly outputs of each setting.

**Figures.** `src/make_figure5.py` draws Figure 5; `src/multiseed_figures.py`
draws Figure 2 and Supplementary Figure S2 (set `BAND = "pi95"` and point
`SCAN_DIRS` to `outputs/core_multiseed`); `src/make_supp_figure_S1_trends.py`
draws Supplementary Figure S1 from the national HIV and AIDS annual report
(AIDS Surveillance Committee, Ministry of Health, Labour and Welfare, 2024;
https://api-net.jfap.or.jp/status/japan/nenpo.html), or from the bundled CSV of
the plotted counts if the Excel tables are not present. `src/check_timing_realisations.py`
reproduces Supplementary Tables S1e and S10a (the time to viral suppression realised
under the implemented viral-load rule, and the waiting times realised at monthly
resolution) by direct calculation, independent of the population model.

## Data and licensing

The data/ folder contains selected aggregated outputs underlying the
principal tables and figures of the article (multiseed comparison across the
200 publication seeds, Monte Carlo standard errors, cost-effectiveness and
break-even summaries, survival-weight and force-of-infection sensitivity
results). Per-seed prevention-decomposition results and their summary
statistics are provided as separate files (prevention_dominance_per_seed.csv
and prevention_dominance_summary_stats.csv). In
survival_adj_sensitivity.csv, the SURVIVAL_ADJ = 0.45 lower-bound stress-test
row reproduces Supplementary Table S4; cells not reported there are left
blank. Cumulative values in breakeven_summary.csv are discounted (2% per
year); undiscounted annual cash flows are in breakeven_report_ready.json.
Two files matter for reproducing the replicate counts of the original
submission: `per_seed_cea_incremental_cost.csv` gives the per-seed incremental
costs under the CEA convention (off-ART ratio 0.45) together with the noisy
whole-population direct QALY difference (not used for the replicate counts),
and `prevention_dominance_per_seed.csv` gives the per-seed
prevention-decomposition QALY gains (`dQALY`) and the budget-impact-convention
costs; `table1_replicate_counts_reproduction.csv` provides their join
precomputed, reproducing dominance in 199 of 200 (primary) and 200 of 200
(secondary) replicates with a positive net monetary benefit in all 200.

`data/revision_R1/` (version 1.1.0) contains `revision_scenarios_summary.csv`
and `.md` (all settings: infections and deaths averted, incremental
QALYs by both estimators, incremental costs, net monetary benefit, dominant
seeds), `per_seed/` (the per-seed contrasts behind each setting and the
published run restricted to the same seeds),
`discount_rate_sensitivity_direct.csv`, `horizon_20y_direct.csv`, and
`base_postprocessing_summary.md` (Part 2), `psa_dsa/` (probabilistic, one-way,
and break-even results for the observed-practice, statutory-minimum with
observed wait, and conservative comparators), and `figures/` (the plotted
values of Figure 2 and Supplementary Figures S1 and S2).

The scripts produce aggregated CSV/JSON summaries and figures; no individual
participant data are used or generated (all calibration targets are public
aggregate surveillance and care-cascade estimates). Code is released under the
MIT License; aggregated outputs under CC BY 4.0. Please cite the article and the
archived Zenodo release (https://doi.org/10.5281/zenodo.23437249; see
`CITATION.cff`).

## Funding and contact

Health and Labour Sciences Research Grant, AIDS Control Policy Research Program,
MHLW Japan (21HB1003, 23HB1001, 26HB1001). Toshibumi Taniguchi, Department of
Infectious Diseases, Chiba University Hospital. E-mail: tosh-tanig@chiba-u.jp

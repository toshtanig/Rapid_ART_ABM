# Changelog

## 1.1.0 (2026-10-20)

Scenario and sensitivity analyses reported in the article. With default settings
the model reproduces version 1.0.0 exactly. Zenodo DOI 10.5281/zenodo.23437249.

- Model: optional scenario settings for the factorial decomposition (eligibility
  criteria and waiting interval), delayed implementation, year-specific
  analytic-period force of infection, time to viral suppression, reduced
  infectivity of diagnosed unsuppressed people, an observed-practice waiting
  distribution, and an observed waiting time for people not meeting the
  certification criteria; a timing-convention switch for scheduled ART starts
  (`ART_START_VL_UPDATE_NEXT_MONTH`) and a run-time-adjustable background-mortality age
  offset (`AGE_MORTALITY_OFFSET`; its default of five years, inherited from version 1.0.0,
  is now documented), for the sensitivity analyses A11 and A12.
- New scripts: `run_revision_scenarios.py`, `summarize_revision_scenarios.py`,
  `compare_regression.py`, `check_progress.py`, `_paths.py`,
  `make_figure5.py`, `make_supp_figure_S1_trends.py`, `check_timing_realisations.py`,
  `survival_weight_sensitivity_per_seed.py` (Supplementary Table S4 and Figure S6 with per-seed costs),
  `make_figures_3_4_S4_S5.py` (Figures 3 and 4 and Supplementary Figures S4 and S5 without in-figure titles),
  `make_figures_S3_S7_S8.py` (Supplementary Figures S3, S7 and S8 from the per-seed outputs).
- `multiseed_figures.py`: Figure 2 and Supplementary Figure S2 drawn at the final
  180 mm width with lowercase panel letters, no title in the artwork, and a 95%
  prediction-interval band; panel d of Figure 2 now shows cumulative new infections
  from 2024 (the per-seed cumulative file counts from the start of the burn-in); the
  other plotted values are unchanged.
- Data: `data/revision_R1/` with the summary tables, per-seed contrasts,
  probabilistic and one-way results for the new comparators, and figure data.
- Comments in `hiv_rapid_art_model.py` revised for readability (section
  headings and descriptions of the calibrated parameters); the code is unchanged
  apart from the options listed above.

## 1.0.0 (2026-08-24)

First release (Zenodo DOI 10.5281/zenodo.21931762).

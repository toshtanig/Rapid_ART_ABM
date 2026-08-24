# Rapid ART cost-effectiveness microsimulation (HIV, Japan)

Source code and reproducibility materials for:

> Taniguchi T, Imahashi M, Sato D, Noda T. Treating earlier, spending less:
> cost-effectiveness and budget impact of immediate versus delayed
> antiretroviral therapy for HIV in Japan. *The Lancet Regional Health - Western Pacific* (under review).

A stochastic individual-based dynamic transmission microsimulation of HIV in
Japan, linked to a cost-effectiveness and budget-impact analysis from the
public-payer perspective. It evaluates permitting immediate ART at diagnosis
(relaxing or reforming the statutory immune-dysfunction certification pathway)
versus one-month and two-month delays, over a 40-year horizon after a 20-year
burn-in. HIV-related mortality is represented by a CD4-stratum by care-state
standardised mortality ratio applied to the age- and sex-specific
general-population background hazard. The full model specification (ODD
protocol) and the principal parameter tables are in the article's Supplementary
Material.

## Layout

```
├── requirements.txt / environment.yml      pinned dependencies
├── LICENSE / LICENSE-DATA.txt              MIT (code) / CC BY 4.0 (data)
├── CITATION.cff                            citation metadata
├── data/                                   aggregated model outputs (CC BY 4.0)
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
    └── make_figure1.py / multiseed_figures.py / compose_figures.py /
        visual_yearly_cost_json.py / lancet_export.py
```

## Dependencies

```
conda env create -f environment.yml        # or: pip install -r requirements.txt
```

Published versions: Python 3.11.15, NumPy 2.4.3, SciPy 1.17.1, pandas 3.0.1,
Numba 0.64.0, llvmlite 0.46.0, Matplotlib 3.10.8, seaborn 0.13.2, Pillow 12.1.1.
(python-docx, listed separately, is only for the optional Word build pipeline,
not the model.)

## Reproducing the analysis

Each script sets its parameters in a `CONFIG` block at the top and takes no
command-line arguments. Set the base/output paths in those blocks to your local
checkout, then:

1. `python src/run_multiseed.py` runs the 200 fixed seeds (`SEEDS` in CONFIG) on
   a single machine, writing one `seed_<n>` folder per seed and appending to
   `multiseed_comparison.csv` after each seed (interruptible). Shorten `SEEDS`
   for a quick check.
2. Run the cost-effectiveness, budget-impact, sensitivity, decomposition, and
   figure scripts on the per-seed outputs (each self-configured via its CONFIG).

Reproducibility rests on the fixed 200-seed list and a common-random-numbers
design: within each replicate the three ART-timing scenarios share re-seeded
random streams, so incremental outcomes are paired per-seed differences and
between-scenario noise is minimised.

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
Two files matter for reproducing the Table 1 replicate counts:
`per_seed_cea_incremental_cost.csv` gives the per-seed incremental costs under
the CEA convention (off-ART ratio 0.45) together with the noisy
whole-population direct QALY difference (not used for the replicate counts),
and `prevention_dominance_per_seed.csv` gives the per-seed
prevention-decomposition QALY gains (`dQALY`) and the budget-impact-convention
costs. The replicate counts in Table 1 pair the 0.45 costs with the
prevention-decomposition `dQALY` (join on seed and comparison);
`table1_replicate_counts_reproduction.csv` provides this join precomputed,
reproducing dominance in 199 of 200 (primary) and 200 of 200 (secondary)
replicates with a positive net monetary benefit in all 200. The
scripts produce aggregated CSV/JSON summaries and figures; no individual
participant data are used or generated (all calibration targets are public
aggregate surveillance and care-cascade estimates). Code is released under the
MIT License; aggregated outputs under CC BY 4.0. Please cite the article and the archived
Zenodo release (https://doi.org/10.5281/zenodo.21931762; see `CITATION.cff`).

## Funding and contact

Health and Labour Sciences Research Grant, AIDS Control Policy Research Program,
MHLW Japan (21HB1003, 23HB1001, 26HB1001). Toshibumi Taniguchi, Department of
Infectious Diseases, Chiba University Hospital. E-mail: tosh-tanig@chiba-u.jp

# Supplementary tables from the per-seed cache (runs of model version 1.2.0)

Base analysis: 200 seeds.

## State in 2023, base analysis (mean and 95% range across seeds, mean ± 1.96 SD)

| indicator | n_seeds | value |
|---|---|---|
| Annual HIV and AIDS notifications | 200 | 891 (800 to 982) |
| People currently on ART | 200 | 29,884 (29,100 to 30,669) |
| People with HIV | 200 | 34,824 (33,755 to 35,893) |
| Undiagnosed people with HIV | 200 | 3,824 (3,505 to 4,142) |
| Undiagnosed fraction (%) | 200 | 11.0 (10.3 to 11.6) |
| Viral suppression on ART (%) | 200 | 99.5 (99.4 to 99.6) |
| New infections per year | 200 | 849 (748 to 950) |

## Fit to the observed notifications (seed median)

| setting | period | n_seeds | RMSE | MAPE_percent | mean_bias | observed_mean | modelled_mean_of_medians |
|---|---|---|---|---|---|---|---|
| B0_core / default | 2004 to 2023 | 200 | 198.0 | 12.6 | -163.5 | 1335.4 | 1171.9 |
| B0_core / default | 2004 to 2013 | 200 | 155.2 | 9.1 | -123.2 | 1434.3 | 1311.2 |
| B0_core / default | 2014 to 2023 | 200 | 233.2 | 16.0 | -203.9 | 1236.6 | 1032.7 |
| A10_conservative_combined / empirical_tail_behaviour_0.3 | 2004 to 2023 | 100 | 276.1 | 18.2 | -250.0 | 1335.4 | 1085.4 |
| A10_conservative_combined / empirical_tail_behaviour_0.3 | 2004 to 2013 | 100 | 258.3 | 16.6 | -239.0 | 1434.3 | 1195.3 |
| A10_conservative_combined / empirical_tail_behaviour_0.3 | 2014 to 2023 | 100 | 292.8 | 19.9 | -261.0 | 1236.6 | 975.6 |
| A10_conservative_combined / empirical_tail_behaviour_0.5 | 2004 to 2023 | 50 | 245.4 | 15.9 | -216.8 | 1335.4 | 1118.6 |
| A10_conservative_combined / empirical_tail_behaviour_0.5 | 2004 to 2013 | 50 | 221.9 | 13.9 | -200.2 | 1434.3 | 1234.2 |
| A10_conservative_combined / empirical_tail_behaviour_0.5 | 2014 to 2023 | 50 | 266.9 | 18.0 | -233.5 | 1236.6 | 1003.1 |
| A12_age_offset / age_offset_0 | 2004 to 2023 | 30 | 142.4 | 9.4 | -78.3 | 1335.4 | 1257.1 |
| A12_age_offset / age_offset_0 | 2004 to 2013 | 30 | 136.0 | 8.3 | -72.8 | 1434.3 | 1361.5 |
| A12_age_offset / age_offset_0 | 2014 to 2023 | 30 | 148.6 | 10.6 | -83.8 | 1236.6 | 1152.8 |
| A4_time_to_suppression / fast | 2004 to 2023 | 30 | 205.1 | 13.0 | -172.0 | 1335.4 | 1163.5 |
| A4_time_to_suppression / fast | 2004 to 2013 | 30 | 161.2 | 9.6 | -132.0 | 1434.3 | 1302.4 |
| A4_time_to_suppression / fast | 2014 to 2023 | 30 | 241.2 | 16.3 | -212.0 | 1236.6 | 1024.6 |
| A4_time_to_suppression / slow | 2004 to 2023 | 30 | 192.6 | 12.2 | -155.2 | 1335.4 | 1180.3 |
| A4_time_to_suppression / slow | 2004 to 2013 | 30 | 151.9 | 8.9 | -117.0 | 1434.3 | 1317.4 |
| A4_time_to_suppression / slow | 2014 to 2023 | 30 | 226.0 | 15.5 | -193.4 | 1236.6 | 1043.2 |
| A5_post_diagnosis_behaviour / factor_0.3 | 2004 to 2023 | 100 | 276.1 | 18.2 | -250.0 | 1335.4 | 1085.4 |
| A5_post_diagnosis_behaviour / factor_0.3 | 2004 to 2013 | 100 | 258.3 | 16.6 | -239.0 | 1434.3 | 1195.3 |
| A5_post_diagnosis_behaviour / factor_0.3 | 2014 to 2023 | 100 | 292.8 | 19.9 | -261.0 | 1236.6 | 975.6 |
| A5_post_diagnosis_behaviour / factor_0.5 | 2004 to 2023 | 100 | 243.7 | 15.8 | -215.3 | 1335.4 | 1120.2 |
| A5_post_diagnosis_behaviour / factor_0.5 | 2004 to 2013 | 100 | 220.6 | 13.8 | -199.6 | 1434.3 | 1234.8 |
| A5_post_diagnosis_behaviour / factor_0.5 | 2014 to 2023 | 100 | 264.8 | 17.8 | -231.0 | 1236.6 | 1005.6 |

## First analytic year (2024)

| setting | arm | n_seeds | notifications_2024_mean | pi95_lo | pi95_hi | observed_2024 |
|---|---|---|---|---|---|---|
| B0_core / default | 0m | 200 | 869.0 | 781.0 | 961.0 | 994 |
| B0_core / default | 1m | 200 | 872.0 | 790.0 | 946.0 | 994 |
| B0_core / default | 2m | 200 | 873.0 | 774.0 | 962.0 | 994 |
| A10_conservative_combined / empirical_tail_behaviour_0.3 | 0m | 100 | 872.0 | 781.0 | 999.0 | 994 |
| A10_conservative_combined / empirical_tail_behaviour_0.3 | 1m | 100 | 874.0 | 776.0 | 1013.0 | 994 |
| A10_conservative_combined / empirical_tail_behaviour_0.5 | 0m | 50 | 873.0 | 793.0 | 973.0 | 994 |
| A10_conservative_combined / empirical_tail_behaviour_0.5 | 1m | 50 | 879.0 | 807.0 | 969.0 | 994 |
| A12_age_offset / age_offset_0 | 0m | 30 | 1018.0 | 906.0 | 1096.0 | 994 |
| A12_age_offset / age_offset_0 | 1m | 30 | 1030.0 | 957.0 | 1120.0 | 994 |
| A4_time_to_suppression / fast | 0m | 30 | 866.0 | 767.0 | 931.0 | 994 |
| A4_time_to_suppression / fast | 1m | 30 | 863.0 | 786.0 | 939.0 | 994 |
| A4_time_to_suppression / slow | 0m | 30 | 884.0 | 800.0 | 984.0 | 994 |
| A4_time_to_suppression / slow | 1m | 30 | 890.0 | 810.0 | 985.0 | 994 |
| A5_post_diagnosis_behaviour / factor_0.3 | 0m | 100 | 872.0 | 781.0 | 999.0 | 994 |
| A5_post_diagnosis_behaviour / factor_0.3 | 1m | 100 | 872.0 | 782.0 | 967.0 | 994 |
| A5_post_diagnosis_behaviour / factor_0.5 | 0m | 100 | 876.0 | 791.0 | 977.0 | 994 |
| A5_post_diagnosis_behaviour / factor_0.5 | 1m | 100 | 880.0 | 791.0 | 994.0 | 994 |

## Heterosexual share of incident infections over the analytic horizon

| setting | arm | n_seeds | heterosexual_share_of_infections_percent |
|---|---|---|---|
| B0_core / default | 0m | 200 | 19.1 |
| B0_core / default | 1m | 200 | 17.9 |
| B0_core / default | 2m | 200 | 17.7 |
| A10_conservative_combined / empirical_tail_behaviour_0.3 | 0m | 100 | 19.4 |
| A10_conservative_combined / empirical_tail_behaviour_0.3 | 1m | 100 | 19.3 |
| A10_conservative_combined / empirical_tail_behaviour_0.5 | 0m | 50 | 19.3 |
| A10_conservative_combined / empirical_tail_behaviour_0.5 | 1m | 50 | 19.2 |
| A12_age_offset / age_offset_0 | 0m | 30 | 17.3 |
| A12_age_offset / age_offset_0 | 1m | 30 | 15.4 |
| A4_time_to_suppression / fast | 0m | 30 | 19.6 |
| A4_time_to_suppression / fast | 1m | 30 | 18.1 |
| A4_time_to_suppression / slow | 0m | 30 | 19.4 |
| A4_time_to_suppression / slow | 1m | 30 | 17.7 |
| A5_post_diagnosis_behaviour / factor_0.3 | 0m | 100 | 19.4 |
| A5_post_diagnosis_behaviour / factor_0.3 | 1m | 100 | 19.0 |
| A5_post_diagnosis_behaviour / factor_0.5 | 0m | 100 | 19.3 |
| A5_post_diagnosis_behaviour / factor_0.5 | 1m | 100 | 18.5 |

## State in 2023 by setting (settings whose burn-in differs from the base)

| setting | Annual HIV and AIDS notifications | People currently on ART | People with HIV | Undiagnosed people with HIV | Undiagnosed fraction (%) | Viral suppression on ART (%) | New infections per year |
|---|---|---|---|---|---|---|---|
| B0_core / default | 890.9 | 29884.5 | 34824.4 | 3823.6 | 11.0 | 99.5 | 848.8 |
| A10_conservative_combined / empirical_tail_behaviour_0.3 | 881.6 | 28516.1 | 33406.2 | 3812.2 | 11.4 | 99.5 | 877.9 |
| A10_conservative_combined / empirical_tail_behaviour_0.5 | 899.2 | 29065.9 | 34005.5 | 3840.6 | 11.3 | 99.5 | 872.5 |
| A12_age_offset / age_offset_0 | 1042.5 | 32596.3 | 38339.1 | 4508.7 | 11.8 | 99.5 | 1013.8 |
| A4_time_to_suppression / fast | 885.9 | 29816.9 | 34751.7 | 3825.1 | 11.0 | 99.6 | 849.3 |
| A4_time_to_suppression / slow | 904.8 | 30038.1 | 35076.5 | 3912.2 | 11.1 | 99.4 | 872.5 |
| A5_post_diagnosis_behaviour / factor_0.3 | 881.6 | 28516.1 | 33406.2 | 3812.2 | 11.4 | 99.5 | 877.9 |
| A5_post_diagnosis_behaviour / factor_0.5 | 901.7 | 29073.9 | 34037.5 | 3865.6 | 11.4 | 99.5 | 881.4 |

## By risk group

| comparison | risk_group | n_seeds | infections_averted | deaths_hiv_averted |
|---|---|---|---|---|
| 0m vs 1m | Heterosexual | 200 | 539 (489 to 588) | 280 (253 to 306) |
| 0m vs 1m | MSM | 200 | 3,207 (3,110 to 3,305) | 1,944 (1,882 to 2,006) |
| 0m vs 2m | Heterosexual | 200 | 744 (695 to 793) | 379 (352 to 406) |
| 0m vs 2m | MSM | 200 | 4,396 (4,305 to 4,487) | 2,661 (2,601 to 2,722) |

## Domestic infections by source state

| strategy | n_seeds | domestic_total | undiagnosed | undiagnosed_percent | diagnosed_pre_art | diagnosed_pre_art_percent | on_art_unsuppressed | on_art_unsuppressed_percent | disengaged | disengaged_percent | exogenous_foreign_inflow (not in the domestic total) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 month | 200 | 11950.5 | 11455.1 | 95.9 | 7.6 | 0.1 | 432.4 | 3.6 | 55.4 | 0.5 | 2677.2 |
| 1 month | 200 | 15696.5 | 13801.7 | 87.9 | 1080.6 | 6.9 | 515.7 | 3.3 | 298.5 | 1.9 | 2683.6 |
| 2 months | 200 | 17090.5 | 14653.3 | 85.7 | 1546.6 | 9.0 | 543.3 | 3.2 | 347.4 | 2.0 | 2679.7 |

| comparison | n_seeds | additional_total | undiagnosed | undiagnosed_percent | diagnosed_pre_art | diagnosed_pre_art_percent | on_art_unsuppressed | on_art_unsuppressed_percent | disengaged | disengaged_percent |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 month delay minus immediate ART | 200 | 3745.9 | 2346.6 | 62.6 | 1072.9 | 28.6 | 83.3 | 2.2 | 243.1 | 6.5 |
| 2 months delay minus immediate ART | 200 | 5140.0 | 3198.2 | 62.2 | 1538.9 | 29.9 | 110.9 | 2.2 | 292.0 | 5.7 |

## Monte Carlo properties

| comparison | outcome | n_seeds | mean | sd | se | relative_mcse_percent | seeds_for_ci95_halfwidth_10pct | seeds_for_ci95_halfwidth_5pct | share_positive_percent |
|---|---|---|---|---|---|---|---|---|---|
| 0m vs 1m | infections_averted | 200 | 3745.94 | 776.08 | 54.88 | 1.46 | 17.0 | 66.0 | 100.0 |
| 0m vs 1m | deaths_averted | 200 | 2223.76 | 492.04 | 34.79 | 1.56 | 19.0 | 76.0 | 100.0 |
| 0m vs 1m | deaths_averted_all | 200 | 520.69 | 843.63 | 59.65 | 11.46 | 1009.0 | 4034.0 | 75.5 |
| 0m vs 1m | dQALY_decomposition | 200 | 11723.78 | 2663.56 | 188.34 | 1.61 | 20.0 | 80.0 | 100.0 |
| 0m vs 1m | dQALY_direct | 200 | 13760.26 | 10765.92 | 761.27 | 5.53 | 236.0 | 941.0 | 91.5 |
| 0m vs 1m | dCost_CEA | 200 | -71167877665.17 | 20500151692.53 | 1449579627.71 | 2.04 | 32.0 | 128.0 | 0.0 |
| 0m vs 1m | NMB | 200 | 129786792341.19 | 32478807192.01 | 2296598481.03 | 1.77 | 25.0 | 97.0 | 100.0 |
| 0m vs 1m | NMB_direct | 200 | 139969156050.73 | 59859985470.53 | 4232740164.79 | 3.02 | 71.0 | 282.0 | 98.0 |
| 0m vs 1m | total discounted QALYs of the immediate-ART arm | 200 | 64202031.17 | 25949.56 | 1834.91 |  |  |  |  |
| 0m vs 2m | infections_averted | 200 | 5139.97 | 759.38 | 53.7 | 1.04 | 9.0 | 34.0 | 100.0 |
| 0m vs 2m | deaths_averted | 200 | 3039.96 | 482.43 | 34.11 | 1.12 | 10.0 | 39.0 | 100.0 |
| 0m vs 2m | deaths_averted_all | 200 | 810.22 | 909.19 | 64.29 | 7.93 | 484.0 | 1935.0 | 83.5 |
| 0m vs 2m | dQALY_decomposition | 200 | 16243.79 | 2673.76 | 189.06 | 1.16 | 11.0 | 42.0 | 100.0 |
| 0m vs 2m | dQALY_direct | 200 | 17803.83 | 10771.6 | 761.67 | 4.28 | 141.0 | 563.0 | 96.0 |
| 0m vs 2m | dCost_CEA | 200 | -95813732433.07 | 20006834396.85 | 1414696827.21 | 1.48 | 17.0 | 67.0 | 0.0 |
| 0m vs 2m | NMB | 200 | 177032703431.39 | 32043940283.63 | 2265848747.05 | 1.28 | 13.0 | 51.0 | 100.0 |
| 0m vs 2m | NMB_direct | 200 | 184832896906.43 | 61491879046.82 | 4348132466.19 | 2.35 | 43.0 | 171.0 | 100.0 |
| 0m vs 2m | total discounted QALYs of the immediate-ART arm | 200 | 64202031.17 | 25949.56 | 1834.91 |  |  |  |  |

## 40-year totals by strategy (means per seed)

| slot | n_seeds | cumulative_new_infections | deaths_among_people_with_hiv | art_initiations | new_diagnoses | deaths_whole_population | total_discounted_qalys | total_discounted_cost_bn |
|---|---|---|---|---|---|---|---|---|
| 0 month | 200 | 11950.5 | 31675.2 | 24441.7 | 14611.1 | 1715004.7 | 64202031.2 | 2176.8 |
| 1 month | 200 | 15696.5 | 33898.9 | 27252.2 | 17365.8 | 1715525.4 | 64188270.9 | 2248.0 |
| 2 months | 200 | 17090.5 | 34715.1 | 28321.3 | 18389.0 | 1715814.9 | 64184227.3 | 2272.6 |

## Waiting times (pooled over seeds)

| analysis | setting | slot | interval | n_seeds | text |
|---|---|---|---|---|---|
| B0_core | default | 0 month | diagnosis to ART | 200 | median 0 (IQR 0 to 0) days; mean 0 days |
| B0_core | default | 0 month | eligibility to ART | 200 | median 0 (IQR 0 to 0) days; mean 6 days |
| B0_core | default | 1 month | diagnosis to ART | 200 | median 30 (IQR 30 to 61) days; mean 126 days |
| B0_core | default | 1 month | eligibility to ART | 200 | median 30 (IQR 30 to 61) days; mean 48 days |
| B0_core | default | 2 months | diagnosis to ART | 200 | median 61 (IQR 61 to 91) days; mean 160 days |
| B0_core | default | 2 months | eligibility to ART | 200 | median 61 (IQR 61 to 91) days; mean 76 days |
| A10_conservative_combined | empirical_tail_behaviour_0.3 | 0 month | diagnosis to ART | 100 | median 0 (IQR 0 to 0) days; mean 0 days |
| A10_conservative_combined | empirical_tail_behaviour_0.3 | 0 month | eligibility to ART | 100 | median 0 (IQR 0 to 0) days; mean 5 days |
| A10_conservative_combined | empirical_tail_behaviour_0.3 | 1 month | diagnosis to ART | 100 | median 61 (IQR 30 to 91) days; mean 92 days |
| A10_conservative_combined | empirical_tail_behaviour_0.3 | 1 month | eligibility to ART | 100 | median 61 (IQR 30 to 91) days; mean 60 days |
| A10_conservative_combined | empirical_tail_behaviour_0.5 | 0 month | diagnosis to ART | 50 | median 0 (IQR 0 to 0) days; mean 0 days |
| A10_conservative_combined | empirical_tail_behaviour_0.5 | 0 month | eligibility to ART | 50 | median 0 (IQR 0 to 0) days; mean 5 days |
| A10_conservative_combined | empirical_tail_behaviour_0.5 | 1 month | diagnosis to ART | 50 | median 61 (IQR 30 to 91) days; mean 92 days |
| A10_conservative_combined | empirical_tail_behaviour_0.5 | 1 month | eligibility to ART | 50 | median 61 (IQR 30 to 91) days; mean 60 days |
| A12_age_offset | age_offset_0 | 0 month | diagnosis to ART | 30 | median 0 (IQR 0 to 0) days; mean 0 days |
| A12_age_offset | age_offset_0 | 0 month | eligibility to ART | 30 | median 0 (IQR 0 to 0) days; mean 5 days |
| A12_age_offset | age_offset_0 | 1 month | diagnosis to ART | 30 | median 30 (IQR 30 to 61) days; mean 127 days |
| A12_age_offset | age_offset_0 | 1 month | eligibility to ART | 30 | median 30 (IQR 30 to 61) days; mean 47 days |
| A13_foi_multiplier | multiplier_1.00 | 0 month | diagnosis to ART | 100 | median 0 (IQR 0 to 0) days; mean 0 days |
| A13_foi_multiplier | multiplier_1.00 | 0 month | eligibility to ART | 100 | median 0 (IQR 0 to 0) days; mean 7 days |
| A13_foi_multiplier | multiplier_1.00 | 1 month | diagnosis to ART | 100 | median 30 (IQR 30 to 61) days; mean 121 days |
| A13_foi_multiplier | multiplier_1.00 | 1 month | eligibility to ART | 100 | median 30 (IQR 30 to 61) days; mean 48 days |
| A13_foi_multiplier | multiplier_1.00 | 2 months | diagnosis to ART | 100 | median 61 (IQR 61 to 91) days; mean 156 days |
| A13_foi_multiplier | multiplier_1.00 | 2 months | eligibility to ART | 100 | median 61 (IQR 61 to 91) days; mean 76 days |
| A13_foi_multiplier | multiplier_1.10 | 0 month | diagnosis to ART | 100 | median 0 (IQR 0 to 0) days; mean 0 days |
| A13_foi_multiplier | multiplier_1.10 | 0 month | eligibility to ART | 100 | median 0 (IQR 0 to 0) days; mean 7 days |
| A13_foi_multiplier | multiplier_1.10 | 1 month | diagnosis to ART | 100 | median 30 (IQR 30 to 61) days; mean 124 days |
| A13_foi_multiplier | multiplier_1.10 | 1 month | eligibility to ART | 100 | median 30 (IQR 30 to 61) days; mean 48 days |
| A13_foi_multiplier | multiplier_1.10 | 2 months | diagnosis to ART | 100 | median 61 (IQR 61 to 91) days; mean 159 days |
| A13_foi_multiplier | multiplier_1.10 | 2 months | eligibility to ART | 100 | median 61 (IQR 61 to 91) days; mean 76 days |
| A13_foi_multiplier | multiplier_1.25 | 0 month | diagnosis to ART | 100 | median 0 (IQR 0 to 0) days; mean 0 days |
| A13_foi_multiplier | multiplier_1.25 | 0 month | eligibility to ART | 100 | median 0 (IQR 0 to 0) days; mean 5 days |
| A13_foi_multiplier | multiplier_1.25 | 1 month | diagnosis to ART | 100 | median 30 (IQR 30 to 61) days; mean 128 days |
| A13_foi_multiplier | multiplier_1.25 | 1 month | eligibility to ART | 100 | median 30 (IQR 30 to 61) days; mean 47 days |
| A13_foi_multiplier | multiplier_1.25 | 2 months | diagnosis to ART | 100 | median 61 (IQR 61 to 91) days; mean 163 days |
| A13_foi_multiplier | multiplier_1.25 | 2 months | eligibility to ART | 100 | median 61 (IQR 61 to 91) days; mean 76 days |
| A13_foi_multiplier | multiplier_1.30 | 0 month | diagnosis to ART | 100 | median 0 (IQR 0 to 0) days; mean 0 days |
| A13_foi_multiplier | multiplier_1.30 | 0 month | eligibility to ART | 100 | median 0 (IQR 0 to 0) days; mean 5 days |
| A13_foi_multiplier | multiplier_1.30 | 1 month | diagnosis to ART | 100 | median 30 (IQR 30 to 61) days; mean 129 days |
| A13_foi_multiplier | multiplier_1.30 | 1 month | eligibility to ART | 100 | median 30 (IQR 30 to 61) days; mean 47 days |
| A13_foi_multiplier | multiplier_1.30 | 2 months | diagnosis to ART | 100 | median 61 (IQR 61 to 91) days; mean 163 days |
| A13_foi_multiplier | multiplier_1.30 | 2 months | eligibility to ART | 100 | median 61 (IQR 61 to 91) days; mean 76 days |
| A1_factorial | criteria_removed_one_month | 0 month | diagnosis to ART | 100 | median 0 (IQR 0 to 0) days; mean 0 days |
| A1_factorial | criteria_removed_one_month | 0 month | eligibility to ART | 100 | median 0 (IQR 0 to 0) days; mean 6 days |
| A1_factorial | criteria_removed_one_month | 1 month | diagnosis to ART | 100 | median 30 (IQR 30 to 61) days; mean 47 days |
| A1_factorial | criteria_removed_one_month | 1 month | eligibility to ART | 100 | median 30 (IQR 30 to 61) days; mean 47 days |
| A1_factorial | criteria_retained_no_wait | 0 month | diagnosis to ART | 100 | median 0 (IQR 0 to 0) days; mean 80 days |
| A1_factorial | criteria_retained_no_wait | 0 month | eligibility to ART | 100 | median 0 (IQR 0 to 0) days; mean 7 days |
| A1_factorial | criteria_retained_no_wait | 1 month | diagnosis to ART | 100 | median 30 (IQR 30 to 61) days; mean 126 days |
| A1_factorial | criteria_retained_no_wait | 1 month | eligibility to ART | 100 | median 30 (IQR 30 to 61) days; mean 48 days |
| A2_foi_trajectory | decline | 0 month | diagnosis to ART | 50 | median 0 (IQR 0 to 0) days; mean 0 days |
| A2_foi_trajectory | decline | 0 month | eligibility to ART | 50 | median 0 (IQR 0 to 0) days; mean 6 days |
| A2_foi_trajectory | decline | 1 month | diagnosis to ART | 50 | median 30 (IQR 30 to 61) days; mean 126 days |
| A2_foi_trajectory | decline | 1 month | eligibility to ART | 50 | median 30 (IQR 30 to 61) days; mean 48 days |
| A2_foi_trajectory | resurgence | 0 month | diagnosis to ART | 50 | median 0 (IQR 0 to 0) days; mean 0 days |
| A2_foi_trajectory | resurgence | 0 month | eligibility to ART | 50 | median 0 (IQR 0 to 0) days; mean 5 days |
| A2_foi_trajectory | resurgence | 1 month | diagnosis to ART | 50 | median 30 (IQR 30 to 61) days; mean 128 days |
| A2_foi_trajectory | resurgence | 1 month | eligibility to ART | 50 | median 30 (IQR 30 to 61) days; mean 47 days |
| A4_time_to_suppression | fast | 0 month | diagnosis to ART | 30 | median 0 (IQR 0 to 0) days; mean 0 days |
| A4_time_to_suppression | fast | 0 month | eligibility to ART | 30 | median 0 (IQR 0 to 0) days; mean 6 days |
| A4_time_to_suppression | fast | 1 month | diagnosis to ART | 30 | median 30 (IQR 30 to 61) days; mean 126 days |
| A4_time_to_suppression | fast | 1 month | eligibility to ART | 30 | median 30 (IQR 30 to 61) days; mean 48 days |
| A4_time_to_suppression | slow | 0 month | diagnosis to ART | 30 | median 0 (IQR 0 to 0) days; mean 0 days |
| A4_time_to_suppression | slow | 0 month | eligibility to ART | 30 | median 0 (IQR 0 to 0) days; mean 6 days |
| A4_time_to_suppression | slow | 1 month | diagnosis to ART | 30 | median 30 (IQR 30 to 61) days; mean 126 days |
| A4_time_to_suppression | slow | 1 month | eligibility to ART | 30 | median 30 (IQR 30 to 61) days; mean 48 days |
| A5_post_diagnosis_behaviour | factor_0.3 | 0 month | diagnosis to ART | 100 | median 0 (IQR 0 to 0) days; mean 0 days |
| A5_post_diagnosis_behaviour | factor_0.3 | 0 month | eligibility to ART | 100 | median 0 (IQR 0 to 0) days; mean 5 days |
| A5_post_diagnosis_behaviour | factor_0.3 | 1 month | diagnosis to ART | 100 | median 30 (IQR 30 to 61) days; mean 128 days |
| A5_post_diagnosis_behaviour | factor_0.3 | 1 month | eligibility to ART | 100 | median 30 (IQR 30 to 61) days; mean 48 days |
| A5_post_diagnosis_behaviour | factor_0.5 | 0 month | diagnosis to ART | 100 | median 0 (IQR 0 to 0) days; mean 0 days |
| A5_post_diagnosis_behaviour | factor_0.5 | 0 month | eligibility to ART | 100 | median 0 (IQR 0 to 0) days; mean 5 days |
| A5_post_diagnosis_behaviour | factor_0.5 | 1 month | diagnosis to ART | 100 | median 30 (IQR 30 to 61) days; mean 128 days |
| A5_post_diagnosis_behaviour | factor_0.5 | 1 month | eligibility to ART | 100 | median 30 (IQR 30 to 61) days; mean 48 days |
| A7_delayed_start | start_2027 | 0 month | diagnosis to ART | 50 | median 0 (IQR 0 to 0) days; mean 14 days |
| A7_delayed_start | start_2027 | 0 month | eligibility to ART | 50 | median 0 (IQR 0 to 0) days; mean 13 days |
| A7_delayed_start | start_2027 | 1 month | diagnosis to ART | 50 | median 30 (IQR 30 to 61) days; mean 126 days |
| A7_delayed_start | start_2027 | 1 month | eligibility to ART | 50 | median 30 (IQR 30 to 61) days; mean 48 days |
| A7_delayed_start | start_2030 | 0 month | diagnosis to ART | 50 | median 0 (IQR 0 to 30) days; mean 33 days |
| A7_delayed_start | start_2030 | 0 month | eligibility to ART | 50 | median 0 (IQR 0 to 30) days; mean 19 days |
| A7_delayed_start | start_2030 | 1 month | diagnosis to ART | 50 | median 30 (IQR 30 to 61) days; mean 126 days |
| A7_delayed_start | start_2030 | 1 month | eligibility to ART | 50 | median 30 (IQR 30 to 61) days; mean 48 days |
| A8_status_quo_empirical | empirical_pathways | 0 month | diagnosis to ART | 50 | median 0 (IQR 0 to 0) days; mean 0 days |
| A8_status_quo_empirical | empirical_pathways | 0 month | eligibility to ART | 50 | median 0 (IQR 0 to 0) days; mean 6 days |
| A8_status_quo_empirical | empirical_pathways | 1 month | diagnosis to ART | 50 | median 61 (IQR 30 to 91) days; mean 146 days |
| A8_status_quo_empirical | empirical_pathways | 1 month | eligibility to ART | 50 | median 61 (IQR 30 to 91) days; mean 64 days |
| A9_observed_ineligible_wait | criteria_retained_no_wait_plus_observed_tail | 0 month | diagnosis to ART | 100 | median 0 (IQR 0 to 0) days; mean 28 days |
| A9_observed_ineligible_wait | criteria_retained_no_wait_plus_observed_tail | 0 month | eligibility to ART | 100 | median 0 (IQR 0 to 0) days; mean 6 days |
| A9_observed_ineligible_wait | criteria_retained_no_wait_plus_observed_tail | 1 month | diagnosis to ART | 100 | median 30 (IQR 30 to 61) days; mean 72 days |
| A9_observed_ineligible_wait | criteria_retained_no_wait_plus_observed_tail | 1 month | eligibility to ART | 100 | median 30 (IQR 30 to 61) days; mean 45 days |
| A9_observed_ineligible_wait | empirical_pathways_plus_observed_tail | 0 month | diagnosis to ART | 200 | median 0 (IQR 0 to 0) days; mean 0 days |
| A9_observed_ineligible_wait | empirical_pathways_plus_observed_tail | 0 month | eligibility to ART | 200 | median 0 (IQR 0 to 0) days; mean 6 days |
| A9_observed_ineligible_wait | empirical_pathways_plus_observed_tail | 1 month | diagnosis to ART | 200 | median 61 (IQR 30 to 91) days; mean 91 days |
| A9_observed_ineligible_wait | empirical_pathways_plus_observed_tail | 1 month | eligibility to ART | 200 | median 61 (IQR 30 to 91) days; mean 61 days |
| A9_observed_ineligible_wait | one_month_plus_observed_tail | 0 month | diagnosis to ART | 50 | median 0 (IQR 0 to 0) days; mean 0 days |
| A9_observed_ineligible_wait | one_month_plus_observed_tail | 0 month | eligibility to ART | 50 | median 0 (IQR 0 to 0) days; mean 6 days |
| A9_observed_ineligible_wait | one_month_plus_observed_tail | 1 month | diagnosis to ART | 50 | median 30 (IQR 30 to 61) days; mean 73 days |
| A9_observed_ineligible_wait | one_month_plus_observed_tail | 1 month | eligibility to ART | 50 | median 30 (IQR 30 to 61) days; mean 45 days |

## People not meeting the eligibility criteria at diagnosis, among ART starters (pooled over seeds)

| analysis | setting | slot | n_seeds | art_starters_per_seed | ineligible_at_diagnosis_per_seed | ineligible_at_diagnosis_percent | diagnosis_to_art_mean_months | diagnosis_to_art_median_months | diagnosis_to_art_p25_months | diagnosis_to_art_p75_months | diagnosis_to_eligibility_mean_months | diagnosis_to_eligibility_median_months | diagnosis_to_eligibility_p25_months | diagnosis_to_eligibility_p75_months |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| B0_core | default | 0 month | 200 | 14684.3 | 5.3 | 0.0 | 3.4 | 1.0 | 1.0 | 1.0 | 2.9 | 1.0 | 1.0 | 1.0 |
| B0_core | default | 1 month | 200 | 17224.7 | 969.0 | 5.6 | 48.9 | 40.0 | 18.0 | 68.0 | 47.0 | 38.0 | 16.0 | 66.0 |
| B0_core | default | 2 months | 200 | 18199.9 | 1192.4 | 6.6 | 44.9 | 34.0 | 14.0 | 64.0 | 42.4 | 32.0 | 12.0 | 61.0 |
| A10_conservative_combined | empirical_tail_behaviour_0.3 | 0 month | 100 | 18064.7 | 5.3 | 0.0 | 3.3 | 1.0 | 1.0 | 1.0 | 2.8 | 1.0 | 1.0 | 1.0 |
| A10_conservative_combined | empirical_tail_behaviour_0.3 | 1 month | 100 | 18740.3 | 1309.5 | 7.0 | 16.1 | 11.0 | 6.0 | 19.0 | 15.2 | 10.0 | 6.0 | 18.0 |
| A10_conservative_combined | empirical_tail_behaviour_0.5 | 0 month | 50 | 16970.0 | 5.3 | 0.0 | 2.7 | 1.0 | 1.0 | 1.0 | 2.4 | 1.0 | 1.0 | 1.0 |
| A10_conservative_combined | empirical_tail_behaviour_0.5 | 1 month | 50 | 18298.3 | 1272.6 | 7.0 | 16.1 | 11.0 | 6.0 | 19.0 | 15.2 | 10.0 | 6.0 | 18.0 |
| A12_age_offset | age_offset_0 | 0 month | 30 | 19197.8 | 6.0 | 0.0 | 3.1 | 1.0 | 1.0 | 1.0 | 2.7 | 1.0 | 1.0 | 1.0 |
| A12_age_offset | age_offset_0 | 1 month | 30 | 23047.1 | 1295.1 | 5.6 | 49.6 | 40.0 | 18.0 | 69.0 | 47.7 | 38.0 | 16.0 | 67.0 |
| A13_foi_multiplier | multiplier_1.00 | 0 month | 100 | 11521.1 | 5.2 | 0.0 | 3.5 | 1.0 | 1.0 | 1.0 | 3.0 | 1.0 | 1.0 | 1.0 |
| A13_foi_multiplier | multiplier_1.00 | 1 month | 100 | 12752.2 | 672.6 | 5.3 | 48.6 | 39.0 | 17.0 | 68.0 | 46.7 | 37.0 | 15.0 | 66.0 |
| A13_foi_multiplier | multiplier_1.00 | 2 months | 100 | 13229.9 | 827.2 | 6.3 | 44.6 | 33.0 | 14.0 | 63.0 | 42.1 | 30.0 | 11.0 | 60.0 |
| A13_foi_multiplier | multiplier_1.10 | 0 month | 100 | 13116.6 | 5.2 | 0.0 | 3.7 | 1.0 | 1.0 | 1.0 | 3.2 | 1.0 | 1.0 | 1.0 |
| A13_foi_multiplier | multiplier_1.10 | 1 month | 100 | 14862.4 | 816.4 | 5.5 | 48.9 | 40.0 | 17.0 | 68.0 | 47.1 | 38.0 | 16.0 | 66.0 |
| A13_foi_multiplier | multiplier_1.10 | 2 months | 100 | 15617.9 | 1002.2 | 6.4 | 45.0 | 34.0 | 14.0 | 64.0 | 42.5 | 31.0 | 12.0 | 61.0 |
| A13_foi_multiplier | multiplier_1.25 | 0 month | 100 | 16437.5 | 5.2 | 0.0 | 3.4 | 1.0 | 1.0 | 1.0 | 2.9 | 1.0 | 1.0 | 1.0 |
| A13_foi_multiplier | multiplier_1.25 | 1 month | 100 | 19812.5 | 1143.5 | 5.8 | 48.6 | 40.0 | 18.0 | 68.0 | 46.7 | 38.0 | 16.0 | 66.0 |
| A13_foi_multiplier | multiplier_1.25 | 2 months | 100 | 21333.5 | 1425.3 | 6.7 | 45.1 | 35.0 | 14.0 | 64.0 | 42.7 | 32.0 | 12.0 | 62.0 |
| A13_foi_multiplier | multiplier_1.30 | 0 month | 100 | 17970.4 | 5.2 | 0.0 | 3.6 | 1.0 | 1.0 | 1.0 | 3.1 | 1.0 | 1.0 | 1.0 |
| A13_foi_multiplier | multiplier_1.30 | 1 month | 100 | 22160.1 | 1299.6 | 5.9 | 48.7 | 40.0 | 18.0 | 68.0 | 46.8 | 38.0 | 16.0 | 66.0 |
| A13_foi_multiplier | multiplier_1.30 | 2 months | 100 | 23982.6 | 1611.0 | 6.7 | 44.9 | 35.0 | 15.0 | 64.0 | 42.5 | 32.0 | 12.0 | 62.0 |
| A1_factorial | criteria_removed_one_month | 0 month | 100 | 14692.8 | 5.2 | 0.0 | 3.6 | 1.0 | 1.0 | 1.0 | 3.2 | 1.0 | 1.0 | 1.0 |
| A1_factorial | criteria_removed_one_month | 1 month | 100 | 15794.0 | 73.9 | 0.5 | 15.1 | 12.0 | 8.0 | 17.0 | 12.6 | 9.0 | 5.0 | 14.0 |
| A1_factorial | criteria_retained_no_wait | 0 month | 100 | 15983.6 | 824.8 | 5.2 | 50.5 | 42.0 | 20.0 | 70.0 | 49.6 | 41.0 | 19.0 | 69.0 |
| A1_factorial | criteria_retained_no_wait | 1 month | 100 | 17190.0 | 966.6 | 5.6 | 48.8 | 40.0 | 18.0 | 68.0 | 46.9 | 38.0 | 16.0 | 66.0 |
| A2_foi_trajectory | decline | 0 month | 50 | 13260.7 | 5.0 | 0.0 | 3.0 | 1.0 | 1.0 | 1.0 | 2.6 | 1.0 | 1.0 | 1.0 |
| A2_foi_trajectory | decline | 1 month | 50 | 15075.0 | 839.4 | 5.6 | 49.2 | 40.0 | 17.0 | 68.0 | 47.3 | 38.0 | 16.0 | 66.0 |
| A2_foi_trajectory | resurgence | 0 month | 50 | 16725.2 | 5.0 | 0.0 | 3.1 | 1.0 | 1.0 | 1.0 | 2.8 | 1.0 | 1.0 | 1.0 |
| A2_foi_trajectory | resurgence | 1 month | 50 | 20526.9 | 1186.7 | 5.8 | 48.9 | 40.0 | 18.0 | 68.0 | 47.0 | 38.0 | 16.0 | 66.0 |
| A4_time_to_suppression | fast | 0 month | 30 | 14580.0 | 5.2 | 0.0 | 3.4 | 1.0 | 1.0 | 1.0 | 2.9 | 1.0 | 1.0 | 1.0 |
| A4_time_to_suppression | fast | 1 month | 30 | 17022.4 | 962.8 | 5.7 | 48.6 | 39.0 | 17.0 | 68.0 | 46.7 | 38.0 | 16.0 | 66.0 |
| A4_time_to_suppression | slow | 0 month | 30 | 15032.8 | 5.1 | 0.0 | 2.8 | 1.0 | 1.0 | 1.0 | 2.4 | 1.0 | 1.0 | 1.0 |
| A4_time_to_suppression | slow | 1 month | 30 | 17486.1 | 990.0 | 5.7 | 48.6 | 39.0 | 17.0 | 68.0 | 46.7 | 38.0 | 16.0 | 66.0 |
| A5_post_diagnosis_behaviour | factor_0.3 | 0 month | 100 | 18064.7 | 5.3 | 0.0 | 3.3 | 1.0 | 1.0 | 1.0 | 2.8 | 1.0 | 1.0 | 1.0 |
| A5_post_diagnosis_behaviour | factor_0.3 | 1 month | 100 | 19039.9 | 1098.6 | 5.8 | 48.9 | 40.0 | 18.0 | 68.0 | 47.0 | 38.0 | 16.0 | 66.0 |
| A5_post_diagnosis_behaviour | factor_0.5 | 0 month | 100 | 17172.4 | 5.2 | 0.0 | 3.2 | 1.0 | 1.0 | 1.0 | 2.8 | 1.0 | 1.0 | 1.0 |
| A5_post_diagnosis_behaviour | factor_0.5 | 1 month | 100 | 18706.4 | 1069.8 | 5.7 | 49.0 | 40.0 | 18.0 | 68.0 | 47.1 | 38.0 | 16.0 | 66.0 |
| A7_delayed_start | start_2027 | 0 month | 50 | 15187.1 | 165.8 | 1.1 | 20.5 | 15.0 | 8.0 | 27.0 | 19.8 | 15.0 | 7.0 | 26.0 |
| A7_delayed_start | start_2027 | 1 month | 50 | 17228.2 | 969.8 | 5.6 | 48.7 | 40.0 | 18.0 | 68.0 | 46.8 | 38.0 | 16.0 | 66.0 |
| A7_delayed_start | start_2030 | 0 month | 50 | 15781.1 | 322.9 | 2.0 | 31.0 | 24.0 | 11.0 | 44.0 | 29.9 | 23.0 | 10.0 | 43.0 |
| A7_delayed_start | start_2030 | 1 month | 50 | 17228.2 | 969.8 | 5.6 | 48.7 | 40.0 | 18.0 | 68.0 | 46.8 | 38.0 | 16.0 | 66.0 |
| A8_status_quo_empirical | empirical_pathways | 0 month | 50 | 14670.5 | 5.0 | 0.0 | 3.1 | 1.0 | 1.0 | 1.0 | 2.7 | 1.0 | 1.0 | 1.0 |
| A8_status_quo_empirical | empirical_pathways | 1 month | 50 | 17678.1 | 1087.3 | 6.2 | 47.1 | 37.0 | 16.0 | 66.0 | 44.4 | 34.0 | 14.0 | 63.0 |
| A9_observed_ineligible_wait | criteria_retained_no_wait_plus_observed_tail | 0 month | 100 | 15066.3 | 871.3 | 5.8 | 15.6 | 10.0 | 5.0 | 18.0 | 15.4 | 10.0 | 5.0 | 18.0 |
| A9_observed_ineligible_wait | criteria_retained_no_wait_plus_observed_tail | 1 month | 100 | 16231.5 | 1025.6 | 6.3 | 15.7 | 11.0 | 6.0 | 19.0 | 15.1 | 10.0 | 5.0 | 18.0 |
| A9_observed_ineligible_wait | empirical_pathways_plus_observed_tail | 0 month | 200 | 14684.3 | 5.3 | 0.0 | 3.4 | 1.0 | 1.0 | 1.0 | 2.9 | 1.0 | 1.0 | 1.0 |
| A9_observed_ineligible_wait | empirical_pathways_plus_observed_tail | 1 month | 200 | 16518.9 | 1131.0 | 6.8 | 16.0 | 11.0 | 6.0 | 19.0 | 15.2 | 10.0 | 6.0 | 18.0 |
| A9_observed_ineligible_wait | one_month_plus_observed_tail | 0 month | 50 | 14670.5 | 5.0 | 0.0 | 3.1 | 1.0 | 1.0 | 1.0 | 2.7 | 1.0 | 1.0 | 1.0 |
| A9_observed_ineligible_wait | one_month_plus_observed_tail | 1 month | 50 | 16261.6 | 1028.2 | 6.3 | 15.8 | 11.0 | 6.0 | 19.0 | 15.2 | 10.0 | 5.0 | 18.0 |

## Undiscounted cost difference in the first five years, budget-impact convention (billion yen; mean and 95% CI across seeds)

| analysis | setting | comparison | year | n_seeds | undiscounted_difference_bn |
|---|---|---|---|---|---|
| B0_core | default | 0m vs 1m | 1 | 200 | 0.67 (0.66 to 0.68) |
| B0_core | default | 0m vs 1m | 2 | 200 | 0.65 (0.63 to 0.67) |
| B0_core | default | 0m vs 1m | 3 | 200 | 0.55 (0.52 to 0.58) |
| B0_core | default | 0m vs 1m | 4 | 200 | 0.41 (0.38 to 0.45) |
| B0_core | default | 0m vs 1m | 5 | 200 | 0.23 (0.19 to 0.27) |
| B0_core | default | 0m vs 1m | 1 to 5 | 200 | 2.51 (2.39 to 2.63) |
| B0_core | default | 0m vs 2m | 1 | 200 | 0.82 (0.81 to 0.83) |
| B0_core | default | 0m vs 2m | 2 | 200 | 0.80 (0.77 to 0.82) |
| B0_core | default | 0m vs 2m | 3 | 200 | 0.68 (0.65 to 0.71) |
| B0_core | default | 0m vs 2m | 4 | 200 | 0.50 (0.46 to 0.54) |
| B0_core | default | 0m vs 2m | 5 | 200 | 0.27 (0.22 to 0.32) |
| B0_core | default | 0m vs 2m | 1 to 5 | 200 | 3.06 (2.94 to 3.19) |
| A10_conservative_combined | empirical_tail_behaviour_0.3 | 0m vs 1m | 1 | 100 | 0.65 (0.64 to 0.67) |
| A10_conservative_combined | empirical_tail_behaviour_0.3 | 0m vs 1m | 2 | 100 | 0.54 (0.51 to 0.56) |
| A10_conservative_combined | empirical_tail_behaviour_0.3 | 0m vs 1m | 3 | 100 | 0.46 (0.42 to 0.49) |
| A10_conservative_combined | empirical_tail_behaviour_0.3 | 0m vs 1m | 4 | 100 | 0.42 (0.37 to 0.47) |
| A10_conservative_combined | empirical_tail_behaviour_0.3 | 0m vs 1m | 5 | 100 | 0.39 (0.32 to 0.46) |
| A10_conservative_combined | empirical_tail_behaviour_0.3 | 0m vs 1m | 1 to 5 | 100 | 2.46 (2.30 to 2.62) |
| A10_conservative_combined | empirical_tail_behaviour_0.5 | 0m vs 1m | 1 | 50 | 0.66 (0.64 to 0.68) |
| A10_conservative_combined | empirical_tail_behaviour_0.5 | 0m vs 1m | 2 | 50 | 0.51 (0.48 to 0.54) |
| A10_conservative_combined | empirical_tail_behaviour_0.5 | 0m vs 1m | 3 | 50 | 0.38 (0.34 to 0.43) |
| A10_conservative_combined | empirical_tail_behaviour_0.5 | 0m vs 1m | 4 | 50 | 0.31 (0.25 to 0.38) |
| A10_conservative_combined | empirical_tail_behaviour_0.5 | 0m vs 1m | 5 | 50 | 0.27 (0.18 to 0.35) |
| A10_conservative_combined | empirical_tail_behaviour_0.5 | 0m vs 1m | 1 to 5 | 50 | 2.13 (1.91 to 2.34) |
| A12_age_offset | age_offset_0 | 0m vs 1m | 1 | 30 | 0.75 (0.73 to 0.78) |
| A12_age_offset | age_offset_0 | 0m vs 1m | 2 | 30 | 0.72 (0.67 to 0.77) |
| A12_age_offset | age_offset_0 | 0m vs 1m | 3 | 30 | 0.60 (0.54 to 0.67) |
| A12_age_offset | age_offset_0 | 0m vs 1m | 4 | 30 | 0.43 (0.35 to 0.51) |
| A12_age_offset | age_offset_0 | 0m vs 1m | 5 | 30 | 0.20 (0.08 to 0.33) |
| A12_age_offset | age_offset_0 | 0m vs 1m | 1 to 5 | 30 | 2.70 (2.41 to 3.00) |
| A13_foi_multiplier | multiplier_1.00 | 0m vs 1m | 1 | 100 | 0.67 (0.65 to 0.68) |
| A13_foi_multiplier | multiplier_1.00 | 0m vs 1m | 2 | 100 | 0.64 (0.61 to 0.67) |
| A13_foi_multiplier | multiplier_1.00 | 0m vs 1m | 3 | 100 | 0.55 (0.51 to 0.59) |
| A13_foi_multiplier | multiplier_1.00 | 0m vs 1m | 4 | 100 | 0.42 (0.37 to 0.47) |
| A13_foi_multiplier | multiplier_1.00 | 0m vs 1m | 5 | 100 | 0.26 (0.21 to 0.32) |
| A13_foi_multiplier | multiplier_1.00 | 0m vs 1m | 1 to 5 | 100 | 2.54 (2.38 to 2.70) |
| A13_foi_multiplier | multiplier_1.00 | 0m vs 2m | 1 | 100 | 0.81 (0.79 to 0.82) |
| A13_foi_multiplier | multiplier_1.00 | 0m vs 2m | 2 | 100 | 0.78 (0.75 to 0.81) |
| A13_foi_multiplier | multiplier_1.00 | 0m vs 2m | 3 | 100 | 0.66 (0.62 to 0.70) |
| A13_foi_multiplier | multiplier_1.00 | 0m vs 2m | 4 | 100 | 0.49 (0.44 to 0.53) |
| A13_foi_multiplier | multiplier_1.00 | 0m vs 2m | 5 | 100 | 0.29 (0.24 to 0.35) |
| A13_foi_multiplier | multiplier_1.00 | 0m vs 2m | 1 to 5 | 100 | 3.02 (2.85 to 3.19) |
| A13_foi_multiplier | multiplier_1.10 | 0m vs 1m | 1 | 100 | 0.67 (0.65 to 0.68) |
| A13_foi_multiplier | multiplier_1.10 | 0m vs 1m | 2 | 100 | 0.62 (0.59 to 0.65) |
| A13_foi_multiplier | multiplier_1.10 | 0m vs 1m | 3 | 100 | 0.54 (0.50 to 0.57) |
| A13_foi_multiplier | multiplier_1.10 | 0m vs 1m | 4 | 100 | 0.41 (0.37 to 0.46) |
| A13_foi_multiplier | multiplier_1.10 | 0m vs 1m | 5 | 100 | 0.26 (0.21 to 0.32) |
| A13_foi_multiplier | multiplier_1.10 | 0m vs 1m | 1 to 5 | 100 | 2.50 (2.35 to 2.66) |
| A13_foi_multiplier | multiplier_1.10 | 0m vs 2m | 1 | 100 | 0.81 (0.79 to 0.82) |
| A13_foi_multiplier | multiplier_1.10 | 0m vs 2m | 2 | 100 | 0.76 (0.73 to 0.80) |
| A13_foi_multiplier | multiplier_1.10 | 0m vs 2m | 3 | 100 | 0.65 (0.61 to 0.69) |
| A13_foi_multiplier | multiplier_1.10 | 0m vs 2m | 4 | 100 | 0.49 (0.43 to 0.54) |
| A13_foi_multiplier | multiplier_1.10 | 0m vs 2m | 5 | 100 | 0.28 (0.22 to 0.35) |
| A13_foi_multiplier | multiplier_1.10 | 0m vs 2m | 1 to 5 | 100 | 2.99 (2.79 to 3.18) |
| A13_foi_multiplier | multiplier_1.25 | 0m vs 1m | 1 | 100 | 0.67 (0.65 to 0.68) |
| A13_foi_multiplier | multiplier_1.25 | 0m vs 1m | 2 | 100 | 0.61 (0.58 to 0.64) |
| A13_foi_multiplier | multiplier_1.25 | 0m vs 1m | 3 | 100 | 0.53 (0.49 to 0.57) |
| A13_foi_multiplier | multiplier_1.25 | 0m vs 1m | 4 | 100 | 0.40 (0.35 to 0.46) |
| A13_foi_multiplier | multiplier_1.25 | 0m vs 1m | 5 | 100 | 0.22 (0.15 to 0.29) |
| A13_foi_multiplier | multiplier_1.25 | 0m vs 1m | 1 to 5 | 100 | 2.43 (2.25 to 2.61) |
| A13_foi_multiplier | multiplier_1.25 | 0m vs 2m | 1 | 100 | 0.81 (0.80 to 0.83) |
| A13_foi_multiplier | multiplier_1.25 | 0m vs 2m | 2 | 100 | 0.77 (0.74 to 0.80) |
| A13_foi_multiplier | multiplier_1.25 | 0m vs 2m | 3 | 100 | 0.64 (0.60 to 0.69) |
| A13_foi_multiplier | multiplier_1.25 | 0m vs 2m | 4 | 100 | 0.43 (0.37 to 0.49) |
| A13_foi_multiplier | multiplier_1.25 | 0m vs 2m | 5 | 100 | 0.16 (0.09 to 0.24) |
| A13_foi_multiplier | multiplier_1.25 | 0m vs 2m | 1 to 5 | 100 | 2.82 (2.63 to 3.01) |
| A13_foi_multiplier | multiplier_1.30 | 0m vs 1m | 1 | 100 | 0.66 (0.65 to 0.67) |
| A13_foi_multiplier | multiplier_1.30 | 0m vs 1m | 2 | 100 | 0.61 (0.59 to 0.64) |
| A13_foi_multiplier | multiplier_1.30 | 0m vs 1m | 3 | 100 | 0.51 (0.47 to 0.55) |
| A13_foi_multiplier | multiplier_1.30 | 0m vs 1m | 4 | 100 | 0.37 (0.32 to 0.42) |
| A13_foi_multiplier | multiplier_1.30 | 0m vs 1m | 5 | 100 | 0.17 (0.10 to 0.23) |
| A13_foi_multiplier | multiplier_1.30 | 0m vs 1m | 1 to 5 | 100 | 2.32 (2.15 to 2.48) |
| A13_foi_multiplier | multiplier_1.30 | 0m vs 2m | 1 | 100 | 0.80 (0.79 to 0.82) |
| A13_foi_multiplier | multiplier_1.30 | 0m vs 2m | 2 | 100 | 0.77 (0.74 to 0.79) |
| A13_foi_multiplier | multiplier_1.30 | 0m vs 2m | 3 | 100 | 0.64 (0.60 to 0.68) |
| A13_foi_multiplier | multiplier_1.30 | 0m vs 2m | 4 | 100 | 0.45 (0.39 to 0.50) |
| A13_foi_multiplier | multiplier_1.30 | 0m vs 2m | 5 | 100 | 0.18 (0.11 to 0.25) |
| A13_foi_multiplier | multiplier_1.30 | 0m vs 2m | 1 to 5 | 100 | 2.84 (2.67 to 3.00) |
| A1_factorial | criteria_removed_one_month | 0m vs 1m | 1 | 100 | 0.28 (0.27 to 0.30) |
| A1_factorial | criteria_removed_one_month | 0m vs 1m | 2 | 100 | 0.23 (0.20 to 0.25) |
| A1_factorial | criteria_removed_one_month | 0m vs 1m | 3 | 100 | 0.19 (0.15 to 0.22) |
| A1_factorial | criteria_removed_one_month | 0m vs 1m | 4 | 100 | 0.10 (0.06 to 0.15) |
| A1_factorial | criteria_removed_one_month | 0m vs 1m | 5 | 100 | 0.04 (-0.01 to 0.10) |
| A1_factorial | criteria_removed_one_month | 0m vs 1m | 1 to 5 | 100 | 0.84 (0.69 to 1.00) |
| A1_factorial | criteria_retained_no_wait | 0m vs 1m | 1 | 100 | 0.22 (0.21 to 0.24) |
| A1_factorial | criteria_retained_no_wait | 0m vs 1m | 2 | 100 | 0.21 (0.18 to 0.24) |
| A1_factorial | criteria_retained_no_wait | 0m vs 1m | 3 | 100 | 0.16 (0.11 to 0.20) |
| A1_factorial | criteria_retained_no_wait | 0m vs 1m | 4 | 100 | 0.08 (0.02 to 0.14) |
| A1_factorial | criteria_retained_no_wait | 0m vs 1m | 5 | 100 | -0.02 (-0.08 to 0.05) |
| A1_factorial | criteria_retained_no_wait | 0m vs 1m | 1 to 5 | 100 | 0.65 (0.47 to 0.84) |
| A2_foi_trajectory | decline | 0m vs 1m | 1 | 50 | 0.68 (0.65 to 0.70) |
| A2_foi_trajectory | decline | 0m vs 1m | 2 | 50 | 0.64 (0.60 to 0.69) |
| A2_foi_trajectory | decline | 0m vs 1m | 3 | 50 | 0.56 (0.50 to 0.63) |
| A2_foi_trajectory | decline | 0m vs 1m | 4 | 50 | 0.41 (0.32 to 0.50) |
| A2_foi_trajectory | decline | 0m vs 1m | 5 | 50 | 0.26 (0.15 to 0.36) |
| A2_foi_trajectory | decline | 0m vs 1m | 1 to 5 | 50 | 2.55 (2.25 to 2.84) |
| A2_foi_trajectory | resurgence | 0m vs 1m | 1 | 50 | 0.68 (0.65 to 0.70) |
| A2_foi_trajectory | resurgence | 0m vs 1m | 2 | 50 | 0.65 (0.60 to 0.69) |
| A2_foi_trajectory | resurgence | 0m vs 1m | 3 | 50 | 0.53 (0.47 to 0.60) |
| A2_foi_trajectory | resurgence | 0m vs 1m | 4 | 50 | 0.39 (0.32 to 0.47) |
| A2_foi_trajectory | resurgence | 0m vs 1m | 5 | 50 | 0.21 (0.12 to 0.30) |
| A2_foi_trajectory | resurgence | 0m vs 1m | 1 to 5 | 50 | 2.46 (2.19 to 2.74) |
| A4_time_to_suppression | fast | 0m vs 1m | 1 | 30 | 0.67 (0.64 to 0.70) |
| A4_time_to_suppression | fast | 0m vs 1m | 2 | 30 | 0.65 (0.59 to 0.71) |
| A4_time_to_suppression | fast | 0m vs 1m | 3 | 30 | 0.58 (0.50 to 0.67) |
| A4_time_to_suppression | fast | 0m vs 1m | 4 | 30 | 0.48 (0.38 to 0.59) |
| A4_time_to_suppression | fast | 0m vs 1m | 5 | 30 | 0.31 (0.18 to 0.45) |
| A4_time_to_suppression | fast | 0m vs 1m | 1 to 5 | 30 | 2.69 (2.33 to 3.06) |
| A4_time_to_suppression | slow | 0m vs 1m | 1 | 30 | 0.68 (0.66 to 0.70) |
| A4_time_to_suppression | slow | 0m vs 1m | 2 | 30 | 0.65 (0.60 to 0.69) |
| A4_time_to_suppression | slow | 0m vs 1m | 3 | 30 | 0.56 (0.49 to 0.64) |
| A4_time_to_suppression | slow | 0m vs 1m | 4 | 30 | 0.41 (0.32 to 0.49) |
| A4_time_to_suppression | slow | 0m vs 1m | 5 | 30 | 0.23 (0.12 to 0.33) |
| A4_time_to_suppression | slow | 0m vs 1m | 1 to 5 | 30 | 2.52 (2.23 to 2.82) |
| A5_post_diagnosis_behaviour | factor_0.3 | 0m vs 1m | 1 | 100 | 0.68 (0.66 to 0.69) |
| A5_post_diagnosis_behaviour | factor_0.3 | 0m vs 1m | 2 | 100 | 0.72 (0.69 to 0.75) |
| A5_post_diagnosis_behaviour | factor_0.3 | 0m vs 1m | 3 | 100 | 0.69 (0.65 to 0.73) |
| A5_post_diagnosis_behaviour | factor_0.3 | 0m vs 1m | 4 | 100 | 0.65 (0.60 to 0.70) |
| A5_post_diagnosis_behaviour | factor_0.3 | 0m vs 1m | 5 | 100 | 0.60 (0.54 to 0.66) |
| A5_post_diagnosis_behaviour | factor_0.3 | 0m vs 1m | 1 to 5 | 100 | 3.34 (3.18 to 3.51) |
| A5_post_diagnosis_behaviour | factor_0.5 | 0m vs 1m | 1 | 100 | 0.68 (0.66 to 0.69) |
| A5_post_diagnosis_behaviour | factor_0.5 | 0m vs 1m | 2 | 100 | 0.67 (0.64 to 0.70) |
| A5_post_diagnosis_behaviour | factor_0.5 | 0m vs 1m | 3 | 100 | 0.63 (0.59 to 0.67) |
| A5_post_diagnosis_behaviour | factor_0.5 | 0m vs 1m | 4 | 100 | 0.59 (0.53 to 0.64) |
| A5_post_diagnosis_behaviour | factor_0.5 | 0m vs 1m | 5 | 100 | 0.52 (0.45 to 0.59) |
| A5_post_diagnosis_behaviour | factor_0.5 | 0m vs 1m | 1 to 5 | 100 | 3.08 (2.90 to 3.27) |
| A7_delayed_start | start_2027 | 0m vs 1m | 1 | 50 | 0.00 (0.00 to 0.00) |
| A7_delayed_start | start_2027 | 0m vs 1m | 2 | 50 | 0.00 (0.00 to 0.00) |
| A7_delayed_start | start_2027 | 0m vs 1m | 3 | 50 | 0.00 (0.00 to 0.00) |
| A7_delayed_start | start_2027 | 0m vs 1m | 4 | 50 | 0.63 (0.60 to 0.65) |
| A7_delayed_start | start_2027 | 0m vs 1m | 5 | 50 | 0.60 (0.57 to 0.64) |
| A7_delayed_start | start_2027 | 0m vs 1m | 1 to 5 | 50 | 1.23 (1.18 to 1.28) |
| A7_delayed_start | start_2030 | 0m vs 1m | 1 | 50 | 0.00 (0.00 to 0.00) |
| A7_delayed_start | start_2030 | 0m vs 1m | 2 | 50 | 0.00 (0.00 to 0.00) |
| A7_delayed_start | start_2030 | 0m vs 1m | 3 | 50 | 0.00 (0.00 to 0.00) |
| A7_delayed_start | start_2030 | 0m vs 1m | 4 | 50 | 0.00 (0.00 to 0.00) |
| A7_delayed_start | start_2030 | 0m vs 1m | 5 | 50 | 0.00 (0.00 to 0.00) |
| A7_delayed_start | start_2030 | 0m vs 1m | 1 to 5 | 50 | 0.00 (0.00 to 0.00) |
| A8_status_quo_empirical | empirical_pathways | 0m vs 1m | 1 | 50 | 0.75 (0.72 to 0.77) |
| A8_status_quo_empirical | empirical_pathways | 0m vs 1m | 2 | 50 | 0.74 (0.70 to 0.78) |
| A8_status_quo_empirical | empirical_pathways | 0m vs 1m | 3 | 50 | 0.66 (0.61 to 0.71) |
| A8_status_quo_empirical | empirical_pathways | 0m vs 1m | 4 | 50 | 0.51 (0.45 to 0.57) |
| A8_status_quo_empirical | empirical_pathways | 0m vs 1m | 5 | 50 | 0.33 (0.26 to 0.40) |
| A8_status_quo_empirical | empirical_pathways | 0m vs 1m | 1 to 5 | 50 | 2.98 (2.78 to 3.19) |
| A9_observed_ineligible_wait | criteria_retained_no_wait_plus_observed_tail | 0m vs 1m | 1 | 100 | 0.23 (0.22 to 0.24) |
| A9_observed_ineligible_wait | criteria_retained_no_wait_plus_observed_tail | 0m vs 1m | 2 | 100 | 0.22 (0.19 to 0.25) |
| A9_observed_ineligible_wait | criteria_retained_no_wait_plus_observed_tail | 0m vs 1m | 3 | 100 | 0.16 (0.12 to 0.20) |
| A9_observed_ineligible_wait | criteria_retained_no_wait_plus_observed_tail | 0m vs 1m | 4 | 100 | 0.08 (0.03 to 0.13) |
| A9_observed_ineligible_wait | criteria_retained_no_wait_plus_observed_tail | 0m vs 1m | 5 | 100 | -0.02 (-0.08 to 0.04) |
| A9_observed_ineligible_wait | criteria_retained_no_wait_plus_observed_tail | 0m vs 1m | 1 to 5 | 100 | 0.67 (0.51 to 0.83) |
| A9_observed_ineligible_wait | empirical_pathways_plus_observed_tail | 0m vs 1m | 1 | 200 | 0.65 (0.64 to 0.66) |
| A9_observed_ineligible_wait | empirical_pathways_plus_observed_tail | 0m vs 1m | 2 | 200 | 0.47 (0.45 to 0.50) |
| A9_observed_ineligible_wait | empirical_pathways_plus_observed_tail | 0m vs 1m | 3 | 200 | 0.32 (0.29 to 0.35) |
| A9_observed_ineligible_wait | empirical_pathways_plus_observed_tail | 0m vs 1m | 4 | 200 | 0.18 (0.14 to 0.21) |
| A9_observed_ineligible_wait | empirical_pathways_plus_observed_tail | 0m vs 1m | 5 | 200 | 0.02 (-0.02 to 0.06) |
| A9_observed_ineligible_wait | empirical_pathways_plus_observed_tail | 0m vs 1m | 1 to 5 | 200 | 1.64 (1.52 to 1.76) |
| A9_observed_ineligible_wait | one_month_plus_observed_tail | 0m vs 1m | 1 | 50 | 0.58 (0.56 to 0.60) |
| A9_observed_ineligible_wait | one_month_plus_observed_tail | 0m vs 1m | 2 | 50 | 0.41 (0.37 to 0.45) |
| A9_observed_ineligible_wait | one_month_plus_observed_tail | 0m vs 1m | 3 | 50 | 0.29 (0.23 to 0.35) |
| A9_observed_ineligible_wait | one_month_plus_observed_tail | 0m vs 1m | 4 | 50 | 0.18 (0.11 to 0.24) |
| A9_observed_ineligible_wait | one_month_plus_observed_tail | 0m vs 1m | 5 | 50 | 0.04 (-0.04 to 0.12) |
| A9_observed_ineligible_wait | one_month_plus_observed_tail | 0m vs 1m | 1 to 5 | 50 | 1.49 (1.27 to 1.72) |

## Break-even year of the mean curve and cumulative discounted cost difference at years 20 and 40 (billion yen)

| analysis | setting | comparison | n_seeds | break_even_year_of_the_mean_curve | 20 | 40 |
|---|---|---|---|---|---|---|
| A10_conservative_combined | empirical_tail_behaviour_0.3 | 0m vs 1m | 100 | 22 | 0.8 | -8.6 |
| A10_conservative_combined | empirical_tail_behaviour_0.5 | 0m vs 1m | 50 | 13 | -7.8 | -29.7 |
| A12_age_offset | age_offset_0 | 0m vs 1m | 30 | 10 | -28.3 | -99.4 |
| A13_foi_multiplier | multiplier_1.00 | 0m vs 1m | 100 | 12 | -10.7 | -32.2 |
| A13_foi_multiplier | multiplier_1.00 | 0m vs 2m | 100 | 12 | -14.7 | -43.8 |
| A13_foi_multiplier | multiplier_1.10 | 0m vs 1m | 100 | 11 | -14.2 | -44.2 |
| A13_foi_multiplier | multiplier_1.10 | 0m vs 2m | 100 | 11 | -20.5 | -62.8 |
| A13_foi_multiplier | multiplier_1.25 | 0m vs 1m | 100 | 10 | -25.1 | -85.0 |
| A13_foi_multiplier | multiplier_1.25 | 0m vs 2m | 100 | 10 | -35.4 | -119.9 |
| A13_foi_multiplier | multiplier_1.30 | 0m vs 1m | 100 | 10 | -30.2 | -103.1 |
| A13_foi_multiplier | multiplier_1.30 | 0m vs 2m | 100 | 9 | -41.6 | -145.0 |
| A1_factorial | criteria_removed_one_month | 0m vs 1m | 100 | 10 | -9.5 | -27.5 |
| A1_factorial | criteria_retained_no_wait | 0m vs 1m | 100 | 9 | -10.9 | -31.5 |
| A2_foi_trajectory | decline | 0m vs 1m | 50 | 11 | -17.3 | -49.3 |
| A2_foi_trajectory | resurgence | 0m vs 1m | 50 | 10 | -25.7 | -92.4 |
| A4_time_to_suppression | fast | 0m vs 1m | 30 | 11 | -17.9 | -60.0 |
| A4_time_to_suppression | slow | 0m vs 1m | 30 | 11 | -19.6 | -63.4 |
| A5_post_diagnosis_behaviour | factor_0.3 | 0m vs 1m | 100 | 21 | 0.2 | -16.3 |
| A5_post_diagnosis_behaviour | factor_0.5 | 0m vs 1m | 100 | 16 | -5.5 | -31.8 |
| A7_delayed_start | start_2027 | 0m vs 1m | 50 | 1 | -10.1 | -46.4 |
| A7_delayed_start | start_2030 | 0m vs 1m | 50 | 1 | -3.3 | -29.7 |
| A8_status_quo_empirical | empirical_pathways | 0m vs 1m | 50 | 11 | -23.3 | -75.4 |
| A9_observed_ineligible_wait | criteria_retained_no_wait_plus_observed_tail | 0m vs 1m | 100 | 8 | -10.8 | -30.1 |
| A9_observed_ineligible_wait | empirical_pathways_plus_observed_tail | 0m vs 1m | 200 | 9 | -16.8 | -47.3 |
| A9_observed_ineligible_wait | one_month_plus_observed_tail | 0m vs 1m | 50 | 9 | -15.3 | -42.5 |
| B0_core | default | 0m vs 1m | 200 | 11 | -20.4 | -65.0 |
| B0_core | default | 0m vs 2m | 200 | 10 | -27.4 | -87.7 |

## Mortality by care state and CD4 stratum, immediate-ART arm

| cd4_stratum | care_state | n_seeds | person_years_per_seed | deaths_per_seed | mortality_rate_per_1000py | expected_deaths_per_seed | OE_ratio |
|---|---|---|---|---|---|---|---|
| 200-349 | disengaged_dropout | 200 | 2585.49 | 333.14 | 128.85 | 332.96 | 1.0 |
| 200-349 | no_art | 200 | 16099.6 | 832.98 | 51.74 | 835.32 | 1.0 |
| 200-349 | on_art_suppressed | 200 | 217359.38 | 6964.7 | 32.04 | 6964.79 | 1.0 |
| 200-349 | on_art_unsuppressed | 200 | 507.42 | 15.84 | 31.23 | 16.11 | 0.98 |
| 350-499 | disengaged_dropout | 200 | 3854.7 | 195.22 | 50.65 | 195.42 | 1.0 |
| 350-499 | no_art | 200 | 12310.1 | 225.72 | 18.34 | 225.59 | 1.0 |
| 350-499 | on_art_suppressed | 200 | 207369.99 | 5674.77 | 27.37 | 5678.98 | 1.0 |
| 350-499 | on_art_unsuppressed | 200 | 221.68 | 5.58 | 25.17 | 5.57 | 1.0 |
| <200 | disengaged_dropout | 200 | 1671.56 | 340.17 | 203.5 | 338.72 | 1.0 |
| <200 | no_art | 200 | 18991.69 | 2091.59 | 110.13 | 2096.84 | 1.0 |
| <200 | on_art_suppressed | 200 | 134371.55 | 5771.44 | 42.95 | 5769.29 | 1.0 |
| <200 | on_art_unsuppressed | 200 | 1580.01 | 90.37 | 57.2 | 90.36 | 1.0 |
| >=500 | disengaged_dropout | 200 | 9070.97 | 379.78 | 41.87 | 378.58 | 1.0 |
| >=500 | no_art | 200 | 10963.17 | 175.01 | 15.96 | 174.58 | 1.0 |
| >=500 | on_art_suppressed | 200 | 299218.13 | 8573.02 | 28.65 | 8576.54 | 1.0 |
| >=500 | on_art_unsuppressed | 200 | 205.95 | 5.8 | 28.19 | 5.98 | 0.97 |

## Counterfactual survival weight

| comparison | survival_weight | n_seeds | qaly_loss_per_averted_infection | dQALY | NMB_bn | dominant_seeds | NMB_positive_seeds |
|---|---|---|---|---|---|---|---|
| 0m vs 1m | 0.45 | 200 | 0.1 | 1,366 (1,241 to 1,491) | 78.0 (75.1 to 80.9) | 187/200 | 200/200 |
| 0m vs 1m | 0.55 | 200 | 1.5 | 6,597 (6,377 to 6,818) | 104.2 (100.5 to 107.8) | 200/200 | 200/200 |
| 0m vs 1m | 0.6 | 200 | 2.19 | 9,213 (8,919 to 9,507) | 117.2 (113.1 to 121.3) | 200/200 | 200/200 |
| 0m vs 1m | 0.65 | 200 | 2.86 | 11,724 (11,355 to 12,093) | 129.8 (125.3 to 134.3) | 200/200 | 200/200 |
| 0m vs 1m | 0.7 | 200 | 3.58 | 14,444 (13,991 to 14,897) | 143.4 (138.5 to 148.3) | 200/200 | 200/200 |
| 0m vs 1m | 0.75 | 200 | 4.28 | 17,060 (16,525 to 17,594) | 156.5 (151.1 to 161.8) | 200/200 | 200/200 |
| 0m vs 1m | 0.8 | 200 | 4.97 | 19,675 (19,059 to 20,292) | 169.5 (163.8 to 175.3) | 200/200 | 200/200 |
| 0m vs 1m | 0.85 | 200 | 5.67 | 22,291 (21,591 to 22,990) | 182.6 (176.4 to 188.8) | 200/200 | 200/200 |
| 0m vs 2m | 0.45 | 200 | 0.14 | 2,211 (2,089 to 2,333) | 106.9 (104.1 to 109.7) | 197/200 | 200/200 |
| 0m vs 2m | 0.55 | 200 | 1.51 | 9,298 (9,076 to 9,520) | 142.3 (138.7 to 145.9) | 200/200 | 200/200 |
| 0m vs 2m | 0.6 | 200 | 2.2 | 12,842 (12,546 to 13,138) | 160.0 (156.0 to 164.1) | 200/200 | 200/200 |
| 0m vs 2m | 0.65 | 200 | 2.86 | 16,244 (15,873 to 16,614) | 177.0 (172.6 to 181.5) | 200/200 | 200/200 |
| 0m vs 2m | 0.7 | 200 | 3.57 | 19,929 (19,475 to 20,383) | 195.5 (190.6 to 200.3) | 200/200 | 200/200 |
| 0m vs 2m | 0.75 | 200 | 4.26 | 23,473 (22,938 to 24,008) | 213.2 (207.9 to 218.5) | 200/200 | 200/200 |
| 0m vs 2m | 0.8 | 200 | 4.95 | 27,017 (26,400 to 27,634) | 230.9 (225.2 to 236.6) | 200/200 | 200/200 |
| 0m vs 2m | 0.85 | 200 | 5.64 | 30,560 (29,861 to 31,260) | 248.6 (242.5 to 254.7) | 200/200 | 200/200 |

## Untreated subgroup with high CD4 count and low viral load

| slot | n_seeds | group_size | fraction_of_people_with_hiv_percent |
|---|---|---|---|
| 0 month | 200 | 58 (57 to 59) | 0.78 |
| 1 month | 200 | 88 (87 to 89) | 1.04 |
| 2 months | 200 | 93 (91 to 94) | 1.04 |

## Foreign-born entries with HIV over the analytic horizon

| slot | n_seeds | cumulative_foreign_born_entries_mean |
|---|---|---|
| 0 month | 200 | 2677.0 |
| 1 month | 200 | 2684.0 |
| 2 months | 200 | 2680.0 |

## Uniform convention (this analysis) and earlier convention (published outputs), same seeds

| comparator | n_seeds | dominant_uniform | dominant_earlier | infections_averted uniform convention | infections_averted earlier convention | infections_averted paired difference | deaths_averted uniform convention | deaths_averted earlier convention | deaths_averted paired difference | dQALY_decomposition uniform convention | dQALY_decomposition earlier convention | dQALY_decomposition paired difference | dCost_CEA uniform convention | dCost_CEA earlier convention | dCost_CEA paired difference | NMB uniform convention | NMB earlier convention | NMB paired difference |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| statutory minimum, 1 month slot | 200 | 200/200 | 199/200 | 3,746 (3,638 to 3,854) | 2,991 (2,871 to 3,112) | 755 (589 to 920) | 2,224 (2,156 to 2,292) | 1,761 (1,686 to 1,837) | 462 (358 to 567) | 11,724 (11,355 to 12,093) | 9,071 (8,671 to 9,470) | 2,653 (2,093 to 3,213) | 71.2 (68.3 to 74.0) | 54.4 (51.3 to 57.4) | 16.8 (12.6 to 21.0) | 129.8 (125.3 to 134.3) | 99.7 (94.8 to 104.6) | 30.1 (23.3 to 36.8) |
| statutory minimum, 2 months slot | 200 | 200/200 | 200/200 | 5,140 (5,035 to 5,245) | 4,535 (4,411 to 4,658) | 605 (437 to 774) | 3,040 (2,973 to 3,107) | 2,688 (2,610 to 2,766) | 352 (244 to 460) | 16,244 (15,873 to 16,614) | 14,095 (13,679 to 14,511) | 2,149 (1,567 to 2,731) | 95.8 (93.0 to 98.6) | 82.5 (79.3 to 85.8) | 13.3 (8.9 to 17.6) | 177.0 (172.6 to 181.5) | 153.0 (147.9 to 158.2) | 24.0 (17.0 to 31.0) |
| observed practice, 1 month slot | 100 | 100/100 | 96/100 | 2,625 (2,480 to 2,770) | 1,839 (1,690 to 1,989) | 786 (567 to 1,005) | 1,610 (1,519 to 1,701) | 1,124 (1,024 to 1,225) | 485 (341 to 629) | 8,728 (8,217 to 9,239) | 6,366 (5,845 to 6,887) | 2,362 (1,604 to 3,120) | 51.9 (48.1 to 55.7) | 35.6 (31.8 to 39.4) | 16.3 (10.7 to 21.8) | 95.5 (89.4 to 101.7) | 67.4 (61.3 to 73.6) | 28.1 (19.1 to 37.1) |

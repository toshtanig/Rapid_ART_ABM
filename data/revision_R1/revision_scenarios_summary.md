# Revision R1: additional scenario analyses (0 month slot versus 1 month slot)

SURVIVAL_ADJ = 0.648; WTP = 5,000,000 yen per QALY; costs under the CEA convention (off-ART cost ratio 0.45); 95% CI across seeds (normal approximation).

## A10_conservative_combined

| analysis | setting | n_seeds | dominant_seeds | NMB_positive_seeds | breakeven_year_median | diag_to_ART_comp_days | diag_to_ART_int_days | Infections averted (95% CI) | Deaths averted (95% CI) | Incremental QALYs, decomposition (95% CI) | Incremental QALYs, direct (95% CI) | Cost saving, billion yen (95% CI) | NMB, billion yen (95% CI) | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A10_conservative_combined | base (published model, same seeds) | 50 | 49/50 | 50/50 | 11.0 | median 30 (IQR 30 to 61); mean 126 | median 0 (IQR 0 to 0); mean 0 | 3129 (2864 to 3393) | 1833 (1662 to 2003) | 9447 (8534 to 10359) | 7364 (4150 to 10577) | 58.1 (51.1 to 65.0) | 105.3 (94.1 to 116.4) |  |
| A10_conservative_combined | empirical_tail_behaviour_0.7 | 50 | 36/50 | 42/50 | 11.0 | median 61 (IQR 30 to 91); mean 91 | median 0 (IQR 0 to 0); mean 0 | 935 (740 to 1131) | 570 (438 to 701) | 3367 (2700 to 4035) | 3553 (744 to 6363) | 15.3 (10.2 to 20.3) | 32.1 (23.9 to 40.3) | A9 empirical_pathways_plus_observed_tail with DIAGNOSED_UNSUPPRESSED_INFECTIVITY_FACTOR = 0.7. |

## A11_timing_convention

| analysis | setting | n_seeds | dominant_seeds | NMB_positive_seeds | breakeven_year_median | diag_to_ART_comp_days | diag_to_ART_int_days | Infections averted (95% CI) | Deaths averted (95% CI) | Incremental QALYs, decomposition (95% CI) | Incremental QALYs, direct (95% CI) | Cost saving, billion yen (95% CI) | NMB, billion yen (95% CI) | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A11_timing_convention | base (published model, same seeds) | 50 | 49/50 | 50/50 | 11.0 | median 30 (IQR 30 to 61); mean 126 | median 0 (IQR 0 to 0); mean 0 | 3129 (2864 to 3393) | 1833 (1662 to 2003) | 9447 (8534 to 10359) | 7364 (4150 to 10577) | 58.1 (51.1 to 65.0) | 105.3 (94.1 to 116.4) |  |
| A11_timing_convention | scheduled_start_vl_update_next_month | 50 | 50/50 | 50/50 | 10.0 | median 30 (IQR 30 to 61); mean 127 | median 0 (IQR 0 to 0); mean 0 | 4418 (4163 to 4673) | 2595 (2433 to 2758) | 13285 (12508 to 14063) | 12677 (9847 to 15507) | 83.1 (76.5 to 89.8) | 149.5 (139.4 to 159.7) | Statutory-minimum comparison. ART_START_VL_UPDATE_NEXT_MONTH = True: the first viral-load update of every start processed at a scheduled start month (eligible people after their wait, restarts after re-engagement) is deferred to the next cycle, the convention already applied to starts at diagnosis; analytic period only. |
| A11_timing_convention | scheduled_start_vl_update_next_month_observed_practice | 50 | 50/50 | 50/50 | 9.0 | median 61 (IQR 30 to 91); mean 91 | median 0 (IQR 0 to 0); mean 0 | 2934 (2695 to 3172) | 1784 (1637 to 1931) | 9744 (9009 to 10478) | 11007 (8112 to 13903) | 58.3 (51.7 to 65.0) | 107.0 (97.1 to 116.9) | Observed-practice comparison (A9 empirical_pathways_plus_observed_tail) under the same timing convention; compare with A9 empirical_pathways_plus_observed_tail on the same seeds (the published convention). Optional second setting. |

## A12_age_offset

| analysis | setting | n_seeds | dominant_seeds | NMB_positive_seeds | breakeven_year_median | diag_to_ART_comp_days | diag_to_ART_int_days | Infections averted (95% CI) | Deaths averted (95% CI) | Incremental QALYs, decomposition (95% CI) | Incremental QALYs, direct (95% CI) | Cost saving, billion yen (95% CI) | NMB, billion yen (95% CI) | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A12_age_offset | base (published model, same seeds) | 30 | 30/30 | 30/30 | 12.0 | median 30 (IQR 30 to 61); mean 126 | median 0 (IQR 0 to 0); mean 0 | 3174 (2793 to 3555) | 1845 (1597 to 2093) | 9401 (8031 to 10770) | 6449 (1486 to 11412) | 58.6 (48.6 to 68.5) | 105.6 (89.2 to 121.9) |  |
| A12_age_offset | age_offset_0 | 30 | 30/30 | 30/30 | 12.0 | median 30 (IQR 30 to 61); mean 127 | median 0 (IQR 0 to 0); mean 0 | 4365 (4010 to 4720) | 2085 (1880 to 2290) | 9867 (8972 to 10763) | 12008 (7290 to 16725) | 83.3 (74.9 to 91.7) | 132.6 (120.6 to 144.7) | AGE_MORTALITY_OFFSET = 0.0 (published model 5.0); all other inputs unchanged. |

## A1_factorial

| analysis | setting | n_seeds | dominant_seeds | NMB_positive_seeds | breakeven_year_median | diag_to_ART_comp_days | diag_to_ART_int_days | Infections averted (95% CI) | Deaths averted (95% CI) | Incremental QALYs, decomposition (95% CI) | Incremental QALYs, direct (95% CI) | Cost saving, billion yen (95% CI) | NMB, billion yen (95% CI) | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A1_factorial | base (published model, same seeds) | 100 | 99/100 | 100/100 | 12.0 | median 30 (IQR 30 to 61); mean 126 | median 0 (IQR 0 to 0); mean 0 | 3072 (2893 to 3250) | 1783 (1673 to 1893) | 9228 (8656 to 9800) | 8826 (6602 to 11049) | 56.4 (51.9 to 60.9) | 102.5 (95.4 to 109.7) |  |
| A1_factorial | criteria_removed_one_month | 100 | 66/100 | 74/100 | 6.5 | median 30 (IQR 30 to 61); mean 47 | median 0 (IQR 0 to 0); mean 0 | 596 (448 to 744) | 376 (282 to 471) | 2460 (1940 to 2981) | 3219 (783 to 5656) | 9.7 (5.8 to 13.5) | 22.0 (15.8 to 28.1) | 0 month slot = cell D (full reform); 1 month slot = cell C (criteria removed, one-month delay retained). Contrast C - D = waiting removal given no criteria. |
| A1_factorial | criteria_retained_no_wait | 100 | 73/100 | 80/100 | 9.0 | median 30 (IQR 30 to 61); mean 126 | median 0 (IQR 0 to 0); mean 81 | 760 (593 to 926) | 453 (354 to 552) | 2395 (1854 to 2936) | 3184 (830 to 5539) | 13.3 (9.0 to 17.6) | 25.2 (18.4 to 32.1) | 0 month slot = cell B (criteria retained, waiting interval removed); 1 month slot = cell A (current pathway). Contrast A - B = waiting removal given criteria. |

### A1 factorial contrasts (per-seed paired differences)

| analysis | setting | n_seeds | dominant_seeds | NMB_positive_seeds | breakeven_year_median | Infections averted (95% CI) | Deaths averted (95% CI) | Incremental QALYs, decomposition (95% CI) | Incremental QALYs, direct (95% CI) | Cost saving, billion yen (95% CI) | NMB, billion yen (95% CI) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A1_factorial | Waiting-interval removal, criteria retained (A - B) | 100 | 73/100 | 80/100 |  | 760 (593 to 926) | 453 (354 to 552) | 2395 (1854 to 2936) | 3184 (830 to 5539) | 13.3 (9.0 to 17.6) | 25.2 (18.4 to 32.1) |
| A1_factorial | Waiting-interval removal, criteria removed (C - D) | 100 | 66/100 | 74/100 |  | 596 (448 to 744) | 376 (282 to 471) | 2460 (1940 to 2981) | 3219 (783 to 5656) | 9.7 (5.8 to 13.5) | 22.0 (15.8 to 28.1) |
| A1_factorial | Eligibility-criteria removal, one-month wait retained (A - C) = (A - D) - (C - D) | 100 | 100/100 | 100/100 |  | 2476 (2310 to 2641) | 1407 (1306 to 1507) | 6768 (6208 to 7328) | 5607 (3413 to 7801) | 46.7 (42.4 to 51.1) | 80.6 (73.6 to 87.5) |
| A1_factorial | Eligibility-criteria removal, no wait (B - D) = (A - D) - (A - B) | 100 | 100/100 | 100/100 |  | 2312 (2165 to 2459) | 1330 (1238 to 1423) | 6833 (6342 to 7323) | 5642 (3370 to 7913) | 43.1 (39.5 to 46.8) | 77.3 (71.5 to 83.1) |
| A1_factorial | Full reform (A - D), same seeds | 100 | 99/100 | 100/100 |  | 3072 (2893 to 3250) | 1783 (1673 to 1893) | 9228 (8656 to 9800) | 8826 (6602 to 11049) | 56.4 (51.9 to 60.9) | 102.5 (95.4 to 109.7) |
| A1_factorial | Interaction (A - B) - (C - D) | 100 | 42/100 | 52/100 |  | 164 (-53 to 380) | 76 (-59 to 212) | -65 (-822 to 692) | -35 (-3505 to 3435) | 3.6 (-1.9 to 9.2) | 3.3 (-5.7 to 12.3) |

| contrast | share of full-reform infections averted (%) | share of full-reform NMB (%) |
|---|---|---|
| Waiting-interval removal, criteria retained (A - B) | 24.7 | 24.6 |
| Waiting-interval removal, criteria removed (C - D) | 19.4 | 21.4 |
| Eligibility-criteria removal, one-month wait retained (A - C) = (A - D) - (C - D) | 80.6 | 78.6 |
| Eligibility-criteria removal, no wait (B - D) = (A - D) - (A - B) | 75.3 | 75.4 |

## A2_foi_trajectory

| analysis | setting | n_seeds | dominant_seeds | NMB_positive_seeds | breakeven_year_median | diag_to_ART_comp_days | diag_to_ART_int_days | Infections averted (95% CI) | Deaths averted (95% CI) | Incremental QALYs, decomposition (95% CI) | Incremental QALYs, direct (95% CI) | Cost saving, billion yen (95% CI) | NMB, billion yen (95% CI) | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A2_foi_trajectory | base (published model, same seeds) | 50 | 49/50 | 50/50 | 11.0 | median 30 (IQR 30 to 61); mean 126 | median 0 (IQR 0 to 0); mean 0 | 3129 (2864 to 3393) | 1833 (1662 to 2003) | 9447 (8534 to 10359) | 7364 (4150 to 10577) | 58.1 (51.1 to 65.0) | 105.3 (94.1 to 116.4) |  |
| A2_foi_trajectory | decline | 50 | 50/50 | 50/50 | 11.5 | median 30 (IQR 30 to 61); mean 126 | median 0 (IQR 0 to 0); mean 0 | 2155 (2019 to 2291) | 1299 (1197 to 1402) | 7565 (7009 to 8120) | 9384 (5904 to 12865) | 42.8 (38.6 to 47.0) | 80.6 (74.1 to 87.2) | Multiplier falls linearly from 1.18 (2024) to 1.00 (2043), then constant. |
| A2_foi_trajectory | resurgence | 50 | 50/50 | 50/50 | 12.0 | median 30 (IQR 30 to 61); mean 128 | median 0 (IQR 0 to 0); mean 0 | 4567 (4254 to 4881) | 2555 (2371 to 2738) | 12193 (11295 to 13092) | 12327 (8685 to 15968) | 76.5 (69.1 to 84.0) | 137.5 (126.0 to 149.0) | Multiplier rises linearly from 1.18 (2024) to 1.30 (2033), then constant. |

## A4_time_to_suppression

| analysis | setting | n_seeds | dominant_seeds | NMB_positive_seeds | breakeven_year_median | diag_to_ART_comp_days | diag_to_ART_int_days | Infections averted (95% CI) | Deaths averted (95% CI) | Incremental QALYs, decomposition (95% CI) | Incremental QALYs, direct (95% CI) | Cost saving, billion yen (95% CI) | NMB, billion yen (95% CI) | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A4_time_to_suppression | base (published model, same seeds) | 30 | 30/30 | 30/30 | 12.0 | median 30 (IQR 30 to 61); mean 126 | median 0 (IQR 0 to 0); mean 0 | 3174 (2793 to 3555) | 1845 (1597 to 2093) | 9401 (8031 to 10770) | 6449 (1486 to 11412) | 58.6 (48.6 to 68.5) | 105.6 (89.2 to 121.9) |  |
| A4_time_to_suppression | fast | 30 | 30/30 | 30/30 | 13.5 | median 30 (IQR 30 to 61); mean 126 | median 0 (IQR 0 to 0); mean 0 | 2844 (2597 to 3091) | 1687 (1539 to 1835) | 8992 (8157 to 9826) | 5349 (757 to 9940) | 49.4 (42.5 to 56.3) | 94.4 (83.8 to 104.9) | 1 + 2 x Beta(2,5) months (median about 1.5, maximum 3). |
| A4_time_to_suppression | slow | 30 | 30/30 | 30/30 | 13.0 | median 30 (IQR 30 to 61); mean 128 | median 0 (IQR 0 to 0); mean 0 | 3200 (2825 to 3576) | 1874 (1650 to 2099) | 9936 (8712 to 11160) | 9335 (5048 to 13621) | 58.9 (49.8 to 68.1) | 108.6 (93.6 to 123.7) | 1 + 8 x Beta(2,5) months (median about 3.1, maximum 9). |

## A5_post_diagnosis_behaviour

| analysis | setting | n_seeds | dominant_seeds | NMB_positive_seeds | breakeven_year_median | diag_to_ART_comp_days | diag_to_ART_int_days | Infections averted (95% CI) | Deaths averted (95% CI) | Incremental QALYs, decomposition (95% CI) | Incremental QALYs, direct (95% CI) | Cost saving, billion yen (95% CI) | NMB, billion yen (95% CI) | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A5_post_diagnosis_behaviour | base (published model, same seeds) | 100 | 99/100 | 100/100 | 12.0 | median 30 (IQR 30 to 61); mean 126 | median 0 (IQR 0 to 0); mean 0 | 3072 (2893 to 3250) | 1783 (1673 to 1893) | 9228 (8656 to 9800) | 8826 (6602 to 11049) | 56.4 (51.9 to 60.9) | 102.5 (95.4 to 109.7) |  |
| A5_post_diagnosis_behaviour | factor_0.5 | 100 | 79/100 | 88/100 | 15.0 | median 30 (IQR 30 to 61); mean 121 | median 0 (IQR 0 to 0); mean 0 | 970 (856 to 1085) | 583 (510 to 656) | 3264 (2856 to 3672) | 4257 (2247 to 6267) | 12.7 (9.4 to 15.9) | 29.0 (23.9 to 34.1) | 50% reduction in transmission from diagnosed unsuppressed sources. |
| A5_post_diagnosis_behaviour | factor_0.7 | 100 | 89/100 | 95/100 | 13.0 | median 30 (IQR 30 to 61); mean 124 | median 0 (IQR 0 to 0); mean 0 | 1491 (1358 to 1625) | 897 (814 to 981) | 4655 (4186 to 5124) | 4110 (1939 to 6281) | 23.5 (19.7 to 27.2) | 46.7 (40.9 to 52.5) | 30% reduction in transmission from diagnosed unsuppressed sources. |

## A7_delayed_start

| analysis | setting | n_seeds | dominant_seeds | NMB_positive_seeds | breakeven_year_median | diag_to_ART_comp_days | diag_to_ART_int_days | Infections averted (95% CI) | Deaths averted (95% CI) | Incremental QALYs, decomposition (95% CI) | Incremental QALYs, direct (95% CI) | Cost saving, billion yen (95% CI) | NMB, billion yen (95% CI) | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A7_delayed_start | base (published model, same seeds) | 50 | 49/50 | 50/50 | 11.0 | median 30 (IQR 30 to 61); mean 126 | median 0 (IQR 0 to 0); mean 0 | 3129 (2864 to 3393) | 1833 (1662 to 2003) | 9447 (8534 to 10359) | 7364 (4150 to 10577) | 58.1 (51.1 to 65.0) | 105.3 (94.1 to 116.4) |  |
| A7_delayed_start | start_2027 | 50 | 49/50 | 50/50 | 1.0 | median 30 (IQR 30 to 61); mean 126 | median 0 (IQR 0 to 0); mean 14 | 2516 (2319 to 2713) | 1447 (1328 to 1565) | 6975 (6400 to 7551) | 8357 (5033 to 11681) | 40.5 (35.9 to 45.1) | 75.4 (68.3 to 82.5) | Immediate ART from calendar year 2027 (three-year implementation delay). |
| A7_delayed_start | start_2030 | 50 | 48/50 | 50/50 | 1.0 | median 30 (IQR 30 to 61); mean 126 | median 0 (IQR 0 to 30); mean 32 | 1997 (1850 to 2144) | 1105 (1019 to 1192) | 4786 (4416 to 5156) | 4354 (2049 to 6660) | 27.3 (23.9 to 30.8) | 51.3 (46.3 to 56.3) | Immediate ART from calendar year 2030 (six-year implementation delay). |

### A7 cost of inaction (2024 implementation minus delayed implementation, per seed)

| analysis | setting | n_seeds | dominant_seeds | NMB_positive_seeds | breakeven_year_median | Infections averted (95% CI) | Deaths averted (95% CI) | Incremental QALYs, decomposition (95% CI) | Incremental QALYs, direct (95% CI) | Cost saving, billion yen (95% CI) | NMB, billion yen (95% CI) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A7_delayed_start | cost of inaction: start_2027 versus 2024 | 50 |  |  |  | 612 (392 to 832) | 386 (242 to 531) | 2471 (1664 to 3278) | -993 (-4203 to 2216) | 17.6 (11.5 to 23.6) | 29.9 (20.1 to 39.7) |
| A7_delayed_start | cost of inaction: start_2030 versus 2024 | 50 |  |  |  | 1132 (922 to 1342) | 727 (581 to 873) | 4661 (3804 to 5517) | 3009 (57 to 5962) | 30.7 (24.8 to 36.6) | 54.0 (44.2 to 63.8) |

## A8_status_quo_empirical

| analysis | setting | n_seeds | dominant_seeds | NMB_positive_seeds | breakeven_year_median | diag_to_ART_comp_days | diag_to_ART_int_days | Infections averted (95% CI) | Deaths averted (95% CI) | Incremental QALYs, decomposition (95% CI) | Incremental QALYs, direct (95% CI) | Cost saving, billion yen (95% CI) | NMB, billion yen (95% CI) | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A8_status_quo_empirical | base (published model, same seeds) | 50 | 49/50 | 50/50 | 11.0 | median 30 (IQR 30 to 61); mean 126 | median 0 (IQR 0 to 0); mean 0 | 3129 (2864 to 3393) | 1833 (1662 to 2003) | 9447 (8534 to 10359) | 7364 (4150 to 10577) | 58.1 (51.1 to 65.0) | 105.3 (94.1 to 116.4) |  |
| A8_status_quo_empirical | empirical_pathways | 50 | 50/50 | 50/50 | 12.0 | median 61 (IQR 30 to 91); mean 146 | median 0 (IQR 0 to 0); mean 0 | 3798 (3548 to 4048) | 2209 (2049 to 2369) | 11496 (10655 to 12337) | 8668 (5357 to 11979) | 69.2 (62.7 to 75.7) | 126.7 (116.2 to 137.1) | 1 month slot = empirical status quo: eligibility-to-ART delay gamma(mean 2.3, SD 1.2 months) for the two-test pathway; 65% of individuals with CD4 below 200 (AIDS exemption) start after gamma(0.5, 0.25) months. Monthly steps round delays up to the next cycle, so the modelled medians are 61 to 91 days versus 70 observed. |

## A9_observed_ineligible_wait

| analysis | setting | n_seeds | dominant_seeds | NMB_positive_seeds | breakeven_year_median | diag_to_ART_comp_days | diag_to_ART_int_days | Infections averted (95% CI) | Deaths averted (95% CI) | Incremental QALYs, decomposition (95% CI) | Incremental QALYs, direct (95% CI) | Cost saving, billion yen (95% CI) | NMB, billion yen (95% CI) | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A9_observed_ineligible_wait | base (published model, same seeds) | 100 | 99/100 | 100/100 | 12.0 | median 30 (IQR 30 to 61); mean 126 | median 0 (IQR 0 to 0); mean 0 | 3072 (2893 to 3250) | 1783 (1673 to 1893) | 9228 (8656 to 9800) | 8826 (6602 to 11049) | 56.4 (51.9 to 60.9) | 102.5 (95.4 to 109.7) |  |
| A9_observed_ineligible_wait | criteria_retained_no_wait_plus_observed_tail | 50 | 33/50 | 37/50 | 6.0 | median 30 (IQR 30 to 61); mean 73 | median 0 (IQR 0 to 0); mean 28 | 539 (305 to 772) | 356 (210 to 501) | 2401 (1660 to 3141) | 5400 (2182 to 8617) | 8.3 (1.9 to 14.7) | 20.3 (10.6 to 30.0) | Factorial cell B under the observed ineligible wait: 0 month slot = criteria retained, waiting interval removed; 1 month slot = statutory-minimum comparator; both with the fallback. Contrast = waiting-interval effect under observed assumptions; the eligibility effect is one_month_plus_observed_tail minus this contrast. |
| A9_observed_ineligible_wait | empirical_pathways_plus_observed_tail | 100 | 96/100 | 96/100 | 10.0 | median 61 (IQR 30 to 91); mean 91 | median 0 (IQR 0 to 0); mean 0 | 1839 (1690 to 1989) | 1124 (1024 to 1225) | 6366 (5845 to 6887) | 6390 (4036 to 8745) | 35.6 (31.8 to 39.4) | 67.4 (61.3 to 73.6) | 1 month slot = empirical status quo (A8 pathways) with the ineligible wait capped by the observed fallback distribution: the most realistic current-practice comparator. |
| A9_observed_ineligible_wait | one_month_plus_observed_tail | 50 | 44/50 | 45/50 | 9.0 | median 30 (IQR 30 to 61); mean 73 | median 0 (IQR 0 to 0); mean 0 | 1199 (969 to 1429) | 725 (578 to 871) | 4195 (3462 to 4929) | 3925 (931 to 6919) | 22.5 (16.8 to 28.2) | 43.5 (34.5 to 52.5) | 1 month slot = statutory-minimum comparator (criteria + one-month delay) with the ineligible wait capped by the observed fallback distribution. |

### A9 decomposition under the observed ineligible wait (per-seed paired differences)

| analysis | setting | n_seeds | dominant_seeds | NMB_positive_seeds | breakeven_year_median | Infections averted (95% CI) | Deaths averted (95% CI) | Incremental QALYs, decomposition (95% CI) | Incremental QALYs, direct (95% CI) | Cost saving, billion yen (95% CI) | NMB, billion yen (95% CI) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A9_observed_ineligible_wait | Observed tail: full reform (A' - D) | 50 | 44/50 | 45/50 |  | 1199 (969 to 1429) | 725 (578 to 871) | 4195 (3462 to 4929) | 3925 (931 to 6919) | 22.5 (16.8 to 28.2) | 43.5 (34.5 to 52.5) |
| A9_observed_ineligible_wait | Observed tail: waiting-interval removal, criteria retained (A' - B') | 50 | 33/50 | 37/50 |  | 539 (305 to 772) | 356 (210 to 501) | 2401 (1660 to 3141) | 5400 (2182 to 8617) | 8.3 (1.9 to 14.7) | 20.3 (10.6 to 30.0) |
| A9_observed_ineligible_wait | Observed tail: eligibility-criteria removal (B' - D) = (A' - D) - (A' - B') | 50 | 33/50 | 37/50 |  | 660 (448 to 872) | 369 (236 to 502) | 1795 (1120 to 2469) | -1475 (-4887 to 1937) | 14.2 (8.7 to 19.7) | 23.2 (14.7 to 31.7) |

| contrast | share of full-reform infections averted (%) | share of full-reform NMB (%) |
|---|---|---|
| Observed tail: waiting-interval removal, criteria retained (A' - B') | 44.9 | 46.7 |
| Observed tail: eligibility-criteria removal (B' - D) = (A' - D) - (A' - B') | 55.1 | 53.3 |

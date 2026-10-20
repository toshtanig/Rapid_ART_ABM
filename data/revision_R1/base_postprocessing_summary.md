# Post-processing of the published 200-seed run (whole-population direct estimator)

## Discount-rate sensitivity (same rate for costs and QALYs; year 1 undiscounted)

| comparison | discount_rate | n_seeds | Incremental QALYs, direct (95% CI) | Cost saving, billion yen (95% CI) | NMB, billion yen (95% CI) | dominant_seeds | NMB_positive_seeds |
|---|---|---|---|---|---|---|---|
| 0m vs 1m | 0.0 | 200 | 13527 (11051 to 16002) | 88.4 (83.7 to 93.1) | 156.0 (141.7 to 170.3) | 155/200 | 187/200 |
| 0m vs 1m | 0.02 | 200 | 8389 (6787 to 9992) | 54.4 (51.3 to 57.4) | 96.3 (87.1 to 105.5) | 153/200 | 184/200 |
| 0m vs 1m | 0.03 | 200 | 6706 (5389 to 8022) | 43.2 (40.7 to 45.6) | 76.7 (69.2 to 84.2) | 151/200 | 182/200 |
| 0m vs 1m | 0.05 | 200 | 4410 (3484 to 5336) | 27.8 (26.1 to 29.5) | 49.9 (44.6 to 55.1) | 145/200 | 180/200 |
| 0m vs 2m | 0.0 | 200 | 24214 (21611 to 26817) | 133.6 (128.6 to 138.5) | 254.6 (239.1 to 270.2) | 179/200 | 198/200 |
| 0m vs 2m | 0.02 | 200 | 15171 (13490 to 16853) | 82.5 (79.3 to 85.8) | 158.4 (148.3 to 168.5) | 180/200 | 197/200 |
| 0m vs 2m | 0.03 | 200 | 12196 (10816 to 13577) | 65.7 (63.1 to 68.4) | 126.7 (118.4 to 135.0) | 178/200 | 197/200 |
| 0m vs 2m | 0.05 | 200 | 8123 (7154 to 9093) | 42.6 (40.8 to 44.5) | 83.2 (77.4 to 89.0) | 176/200 | 196/200 |

## 20-year horizon (discounted at 2% per year)

| comparison | horizon_years | n_seeds | Infections averted (95% CI) | Deaths averted (95% CI) | Incremental QALYs, direct (95% CI) | Cost saving, billion yen (95% CI) | NMB, billion yen (95% CI) | dominant_seeds |
|---|---|---|---|---|---|---|---|---|
| 0m vs 1m | 20 | 200 | 1971 (1879 to 2063) | 646 (612 to 681) | 2922 (2068 to 3777) | 18.0 (16.6 to 19.4) | 32.6 (27.9 to 37.3) | 126/200 |
| 0m vs 2m | 20 | 200 | 2962 (2865 to 3059) | 997 (962 to 1032) | 5556 (4683 to 6429) | 28.2 (26.6 to 29.8) | 55.9 (50.8 to 61.1) | 165/200 |

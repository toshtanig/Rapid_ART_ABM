#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Realised timing under the model's monthly conventions (Supplementary Tables S1e and S10a).

Two direct calculations, independent of the population model and of the number of seeds:

1. Time to viral suppression under the implemented viral-load rule (hiv_rapid_art_model.py,
   simulate_viral_load_randomized called with the current viral load): at each monthly update the
   viral load is reduced by LOG_DROP * (t / tau) log10, where t is the months since ART initiation and
   tau the nominal duration drawn at initiation, so the cumulative decline after k updates is
   LOG_DROP * k (k + 1) / (2 tau). Uniform noise and blips are added as in the model. For comparison
   the "described" rule (a linear decline of LOG_DROP over tau from the initial value) is also run.

2. Waiting times realised at monthly resolution: ART starts at the first monthly cycle after the
   drawn gamma delay has elapsed, so the realised wait is the delay rounded up to whole months.

Run as-is in VS Code; settings below. Outputs two CSV files in OUT_DIR.
"""
import os
import numpy as np
import pandas as pd

# =============================================================================
# CONFIG
# =============================================================================
HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(os.path.dirname(HERE), "outputs", "timing_checks")
N_AGENTS = 400_000          # individuals per cell for the viral-load simulation
N_DRAWS = 2_000_000         # draws per waiting distribution
SEED = 1
LOG_DROP = 4.0                                   # total nominal decline, log10 copies/mL
DECAY_PER_MONTH = 0.2 / 365.0 * 30.0             # slow decay after tau (vl_decay_rate_per_month)
NOISE = 0.2                                      # uniform noise, +/- log10
BLIP_PROB, BLIP_SIZE = 0.01, 0.5                 # transient blips
SUPPRESSION_THRESHOLD = 1.5                      # log10 copies/mL (about 30 copies/mL)
MAX_MONTHS = 36
TTS_SETTINGS = {"base (1 + 5 Beta)": 5.0, "shortened (1 + 2 Beta)": 2.0, "lengthened (1 + 8 Beta)": 8.0}   # TTS_RANGE_MONTHS
VL_BANDS = {"3.0-4.0 (CD4 >=500)": (3.0, 4.0), "4.5-5.5 (CD4 200-499)": (4.5, 5.5), "5.5-6.5 (CD4 <200)": (5.5, 6.5)}
WAITS = {   # pathway: (gamma mean, gamma SD) in months
    "statutory minimum, eligible": (1.0, 0.5),
    "two-month comparator, eligible": (2.0, 0.5),
    "observed practice, eligible (two-test pathway)": (2.3, 1.2),
    "observed practice, AIDS exemption": (0.5, 0.25),
    "observed ineligible wait": (14.0, 11.0),
}
DAYS_PER_MONTH = 365.25 / 12.0
# =============================================================================


def simulate_suppression(rng, mode, vl0, tau):
    """Months (monthly updates) until the viral load first falls to the threshold; 99 = not within MAX_MONTHS."""
    n = len(vl0)
    vl = vl0.copy()
    first = np.full(n, 99)
    for m in range(1, MAX_MONTHS + 1):
        t = float(m)
        base = vl if mode == "implemented" else vl0
        before = base - LOG_DROP * np.minimum(t / tau, 1.0)
        after = np.maximum(base - LOG_DROP, 0.0) - DECAY_PER_MONTH * np.maximum(t - tau, 0.0)
        new = np.where(t < tau, before, after)
        new = np.maximum(new, 0.0)
        new = np.maximum(new + rng.uniform(-NOISE, NOISE, n), 0.0)
        new = new + (rng.random(n) < BLIP_PROB) * BLIP_SIZE
        vl = new
        first[(vl <= SUPPRESSION_THRESHOLD) & (first == 99)] = m
    return first


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    rng = np.random.default_rng(SEED)

    rows = []
    for setting, rng_months in TTS_SETTINGS.items():
        for band, (lo, hi) in VL_BANDS.items():
            vl0 = rng.uniform(lo, hi, N_AGENTS)
            tau = 1.0 + rng_months * rng.beta(2.0, 5.0, size=N_AGENTS)
            for mode in ("implemented", "described (linear from initial value)"):
                first = simulate_suppression(rng, mode, vl0, tau)
                sup = first[first < 99]
                rows.append({
                    "setting": setting, "vl_band_at_initiation": band, "rule": mode,
                    "nominal_tau_mean_months": round(float(tau.mean()), 2),
                    "pct_suppressed_by_update_1": round(float((first <= 1).mean() * 100), 1),
                    "pct_suppressed_by_update_2": round(float((first <= 2).mean() * 100), 1),
                    "pct_suppressed_by_update_3": round(float((first <= 3).mean() * 100), 1),
                    "pct_suppressed_by_update_6": round(float((first <= 6).mean() * 100), 1),
                    "pct_not_suppressed_by_12": round(float((first > 12).mean() * 100), 1),
                    "realised_median_cycles": float(np.median(sup)) if len(sup) else np.nan,
                    "realised_mean_cycles": round(float(sup.mean()), 2) if len(sup) else np.nan,
                })
    df1 = pd.DataFrame(rows)
    df1.to_csv(os.path.join(OUT_DIR, "realised_time_to_suppression.csv"), index=False)

    rows = []
    for pathway, (mean, sd) in WAITS.items():
        shape, scale = mean ** 2 / sd ** 2, sd ** 2 / mean
        x = rng.gamma(shape, scale, N_DRAWS)
        c = np.ceil(x)
        rows.append({
            "pathway": pathway, "input_mean_months": mean, "input_sd_months": sd,
            "input_median_months": round(float(np.median(x)), 2),
            "pct_realised_1": round(float((c == 1).mean() * 100), 1), "pct_realised_2": round(float((c == 2).mean() * 100), 1),
            "pct_realised_3": round(float((c == 3).mean() * 100), 1), "pct_realised_4": round(float((c == 4).mean() * 100), 1),
            "pct_realised_5_or_more": round(float((c >= 5).mean() * 100), 1),
            "realised_mean_cycles": round(float(c.mean()), 2), "realised_median_cycles": float(np.median(c)),
            "realised_p25_cycles": float(np.percentile(c, 25)), "realised_p75_cycles": float(np.percentile(c, 75)),
            "realised_mean_days": round(float(c.mean() * DAYS_PER_MONTH)), "realised_median_days": round(float(np.median(c) * DAYS_PER_MONTH)),
        })
    df2 = pd.DataFrame(rows)
    df2.to_csv(os.path.join(OUT_DIR, "realised_waiting_times.csv"), index=False)

    pd.set_option("display.width", 220)
    print(df1.to_string(index=False))
    print()
    print(df2.to_string(index=False))
    print("\nsaved to", OUT_DIR)


if __name__ == "__main__":
    main()

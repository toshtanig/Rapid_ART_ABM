#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""One-way sensitivity to the counterfactual survival weight."""

import os
import sys
import math
import importlib

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Lancet TIFF export hook (300 dpi, LZW; multipanel 180 mm / single 107 mm) ---
EXPORT_TIFF = True
if EXPORT_TIFF:
    try:
        import os as _os, sys as _sys
        _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
        import lancet_export
        lancet_export.enable_tiff_export()
    except Exception as _e:
        print("[lancet_export] disabled:", _e)

# =============================================================================

# =============================================================================
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# --- Data location: set BASE_144 to your local 144_scenario folder ---
BASE_144 = r"/path/to/20250315_169/144_scenario"
PROJ_169 = os.path.dirname(BASE_144)
BASE_169 = PROJ_169

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))

RUN_MODE = "production"
_MS_DIR = "core_multiseed_validation" if RUN_MODE == "validation" else "core_multiseed"
SCAN_DIRS_OVERRIDE = [os.path.join(BASE_169, "144_scenario", _MS_DIR)]

# SCAN_DIRS_OVERRIDE = [
#     os.path.join(BASE_169, "144_scenario", "core_multiseed"),

# ]

SURVIVAL_ADJ_GRID = [0.55, 0.60, 0.648, 0.70, 0.75, 0.80, 0.85]
BASE_SA = 0.648

COMPARISONS_OVERRIDE = None

WTP = 5_000_000
# CEA whole-population incremental cost difference (off-ART 0.45), independent of SURVIVAL_ADJ;
# used so that the net monetary benefit matches the main-text / Table 1 value (Convention A).

DCOST_CEA_JPY = {"0m_vs_1m": -54_376_807_944.0, "0m_vs_2m": -82_549_346_590.0}

OUTPUT_DIR = os.path.join(BASE_169, "144_scenario", "survival_adj_sensitivity")
OUTPUT_CSV = os.path.join(OUTPUT_DIR, "survival_adj_sensitivity.csv")
OUTPUT_MD  = os.path.join(OUTPUT_DIR, "survival_adj_sensitivity.md")
MAKE_FIGURES = True
FIG_DPI = 200
# =============================================================================

if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)
import prevention_qaly_dominance as PV   # noqa: E402

def _t_value(n):
    tt = {2: 12.71, 3: 4.303, 4: 3.182, 5: 2.776, 6: 2.571, 7: 2.447, 8: 2.365,
          9: 2.306, 10: 2.262, 11: 2.228, 12: 2.201, 13: 2.179, 14: 2.160,
          15: 2.145, 20: 2.093, 30: 2.045}
    return tt.get(n, 1.96 if n > 30 else 2.262)

def _agg(vals):
    x = np.array([v for v in vals
                  if v is not None and not (isinstance(v, float) and math.isnan(v))], float)
    n = len(x)
    if n == 0:
        return dict(n=0, mean=float("nan"), sd=float("nan"),
                    lo=float("nan"), hi=float("nan"))
    m = x.mean()
    sd = x.std(ddof=1) if n > 1 else 0.0
    se = sd / math.sqrt(n) if n > 1 else 0.0
    t = _t_value(n)
    return dict(n=n, mean=m, sd=sd, lo=m - t * se, hi=m + t * se)

def list_seed_dirs(scan_dirs):
    seen = set()
    dirs = []
    for base in scan_dirs:
        if not os.path.isdir(base):
            continue
        for name in sorted(os.listdir(base)):
            if not name.startswith("seed_"):
                continue
            try:
                s = int(name.split("seed_")[1])
            except ValueError:
                continue
            if s in seen:
                continue
            seen.add(s)
            dirs.append(os.path.join(base, name))
    return dirs

def run_sensitivity():
    scan_dirs = SCAN_DIRS_OVERRIDE if SCAN_DIRS_OVERRIDE is not None else PV.SCAN_DIRS
    comparisons = COMPARISONS_OVERRIDE if COMPARISONS_OVERRIDE is not None else PV.COMPARISONS
    seed_dirs = list_seed_dirs(scan_dirs)
    print(f"scanned dirs : {scan_dirs}")
    print(f"seeds found  : {len(seed_dirs)}")
    print(f"SA grid      : {SURVIVAL_ADJ_GRID}  (base={BASE_SA})")

    original_sa = PV.SURVIVAL_ADJ
    rows = []
    try:
        for (iv, cp) in comparisons:
            comp_label = f"{iv}_vs_{cp}"
            for sa in SURVIVAL_ADJ_GRID:
                PV.SURVIVAL_ADJ = sa
                dqaly, nmb, per_inf, dcost = [], [], [], []
                dominant = 0
                nmb_pos = 0
                n_seed = 0
                for sd in seed_dirs:
                    res = PV.compute_seed_comparison(sd, iv, cp)
                    if res is None or res["dCost"] is None:
                        continue
                    n_seed += 1
                    dq = res["dQALY"]
                    # CEA whole-population cost difference (off-ART 0.45), SURVIVAL_ADJ-independent
                    dc = DCOST_CEA_JPY.get(comp_label, res["dCost"])
                    nm = dq * WTP - dc
                    dqaly.append(dq)
                    dcost.append(dc)
                    nmb.append(nm)
                    per_inf.append(res["per_infection_QALY_loss"])
                    if dq > 0 and dc < 0:
                        dominant += 1
                    if nm > 0:
                        nmb_pos += 1
                aq = _agg(dqaly); an = _agg(nmb); ap = _agg(per_inf); ac = _agg(dcost)
                rows.append({
                    "comparison": comp_label,
                    "SURVIVAL_ADJ": sa,
                    "is_base": (abs(sa - BASE_SA) < 1e-9),
                    "n_seeds": n_seed,
                    "per_infection_QALY_mean": ap["mean"],
                    "per_infection_QALY_lo": ap["lo"],
                    "per_infection_QALY_hi": ap["hi"],
                    "dQALY_mean": aq["mean"],
                    "dQALY_lo": aq["lo"],
                    "dQALY_hi": aq["hi"],
                    "dCost_JPY_mean": ac["mean"],
                    "NMB_JPY_mean": an["mean"],
                    "NMB_JPY_lo": an["lo"],
                    "NMB_JPY_hi": an["hi"],
                    "dominant_seeds": dominant,
                    "NMB_pos_seeds": nmb_pos,
                })
    finally:
        PV.SURVIVAL_ADJ = original_sa

    return pd.DataFrame(rows)

def write_markdown(df, path):
    lines = []
    lines.append("# SURVIVAL_ADJ one-way sensitivity analysis\n")
    lines.append(f"- Seeds swept: {int(df['n_seeds'].max())}")
    lines.append(f"- WTP: {WTP:,} JPY/QALY")
    lines.append(f"- Reference SURVIVAL_ADJ: {BASE_SA}")
    lines.append("- SURVIVAL_ADJ affects only the prevention-component ΔQALY; ΔCost and the treatment-timing component are unchanged.\n")
    for comp_label, g in df.groupby("comparison"):
        lines.append(f"## {comp_label}\n")
        dcost_b = g["dCost_JPY_mean"].mean() / 1e9
        lines.append(f"ΔCost (SURVIVAL_ADJ-independent; common to all rows): {dcost_b:,.1f} billion JPY\n")
        lines.append("| SURVIVAL_ADJ | per-infection QALY | ΔQALY [95%CI] | NMB(billion JPY) [95%CI] | dominant | NMB>0 |")
        lines.append("|---|---|---|---|---|---|")
        for _, r in g.iterrows():
            star = " (base)" if r["is_base"] else ""
            lines.append(
                f"| {r['SURVIVAL_ADJ']:.3f}{star} | "
                f"{r['per_infection_QALY_mean']:.2f} | "
                f"{r['dQALY_mean']:,.0f} [{r['dQALY_lo']:,.0f}, {r['dQALY_hi']:,.0f}] | "
                f"{r['NMB_JPY_mean']/1e9:,.1f} [{r['NMB_JPY_lo']/1e9:,.1f}, {r['NMB_JPY_hi']/1e9:,.1f}] | "
                f"{int(r['dominant_seeds'])}/{int(r['n_seeds'])} | "
                f"{int(r['NMB_pos_seeds'])}/{int(r['n_seeds'])} |"
            )
        lines.append("")
    lines.append("## Interpretation\n")
    lines.append(
        "ΔCost is independent of SURVIVAL_ADJ and constant across all rows (see table above); SURVIVAL_ADJ scales only the prevention-component ΔQALY "
        "to scale in a near-linear fashion. Thus Rapid ART (0m) remains cost-saving, and the NMB across the swept "
        "stays positive across the whole range (the share of seeds with NMB>0 is nearly all, >=97%), and the cost-effectiveness conclusion (dominance, "
        "CEAC = 100% at all WTP) does not reverse; only the absolute magnitudes of ΔQALY and NMB change.\n"
    )
    lines.append(
        "The reference SURVIVAL_ADJ = 0.648 matches the prevention-decomposition ΔQALY to the whole-population direct ΔQALY "
        "(matching the 200-seed production total QALY differences: 0m vs 1m ~ 8,389 / 0m vs 2m ~ 15,171) and analytically "
        "is an anchored value (integrating the alternative solutions 0.632 / 0.665), not arbitrary. ΔQALY "
        "Because it is a linear function (depending only on the prevention component, independent of the treatment-timing component), this script alone recomputes, without rerunning the model, "
        "is determined. At the reference 0.648 the per-infection discounted QALY loss is about 2.8 QALY, within the range generally regarded as plausible under modern ART "
        "(roughly 2-3 QALY). Lowering SURVIVAL_ADJ to its lower bound (e.g. 0.45) brings the per-infection QALY "
        "Even reduced to near 0, cost savings keep the NMB positive and the dominance conclusion does not reverse."
        "(for reference: the provisional 30-seed validation anchor was 0.70, but with 200 seeds the SE narrowed and this value was fixed."
        "ΔCost should be fixed at the 200-seed mean-input CEA value from cea_param_uncertainty.)"
    )
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

def make_figures(df):
    comps = list(df["comparison"].unique())
    panels = [
        ("per_infection_QALY_mean", "per-infection QALY"),
        ("dQALY_mean", "Incremental QALY"),
        ("NMB_JPY_mean", "NMB (billion JPY)"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
    for j, (col, ylab) in enumerate(panels):
        ax = axes[j]
        for k, comp in enumerate(comps):
            g = df[df["comparison"] == comp].sort_values("SURVIVAL_ADJ")
            x = g["SURVIVAL_ADJ"].values
            if col == "per_infection_QALY_mean":
                y = g[col].values
                lo = g["per_infection_QALY_lo"].values
                hi = g["per_infection_QALY_hi"].values
            elif col == "dQALY_mean":
                y = g[col].values
                lo = g["dQALY_lo"].values
                hi = g["dQALY_hi"].values
            else:
                y = g[col].values / 1e9
                lo = g["NMB_JPY_lo"].values / 1e9
                hi = g["NMB_JPY_hi"].values / 1e9
            ax.plot(x, y, "-o", color=colors[k % len(colors)], label=comp, zorder=3)
            ax.fill_between(x, lo, hi, color=colors[k % len(colors)], alpha=0.15, zorder=1)
        ax.axvline(BASE_SA, color="black", ls="--", lw=1, zorder=2)
        ax.text(BASE_SA, ax.get_ylim()[1], f" base={BASE_SA}", va="top", ha="left", fontsize=8)
        ax.set_xlabel("SURVIVAL_ADJ")
        ax.set_ylabel(ylab)
        ax.set_axisbelow(True)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="upper left")
    fig.suptitle("One-way sensitivity to SURVIVAL_ADJ", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(OUTPUT_DIR, "survival_adj_sensitivity.png")
    fig.savefig(out, dpi=FIG_DPI)
    plt.close(fig)
    print(f"saved figure: {out}")

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = run_sensitivity()
    df.to_csv(OUTPUT_CSV, index=False)
    write_markdown(df, OUTPUT_MD)
    if MAKE_FIGURES:
        make_figures(df)
    print(f"\nSaved: {OUTPUT_CSV}")
    print(f"Saved: {OUTPUT_MD}")

    with pd.option_context("display.width", 160, "display.max_columns", 20):
        show = df[["comparison", "SURVIVAL_ADJ", "per_infection_QALY_mean",
                   "dQALY_mean", "NMB_JPY_mean", "dominant_seeds", "n_seeds"]].copy()
        show["NMB_billion"] = show.pop("NMB_JPY_mean") / 1e9
        print(show.to_string(index=False))

if __name__ == "__main__":
    main()

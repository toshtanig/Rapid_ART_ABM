#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Supplementary Figure S1: national HIV and AIDS notifications in Japan, 1985 to 2024.

Source: AIDS Surveillance Committee, Ministry of Health, Labour and Welfare, annual report 2024
(table 4: annual reports by nationality, transmission route, and diagnosis category).
Cases infected through coagulation-factor products are not included in that table.
Panel a: new reports by diagnosis category (HIV without AIDS at report; AIDS at first report).
Panel b: new reports (HIV plus AIDS) by reported transmission route.
Lancet figure rules: no title inside the artwork, lowercase panel letters, 180 mm width, 300 dpi.
The annual-report Excel tables can be downloaded from https://api-net.jfap.or.jp/status/japan/nenpo.html;
if the file is not present, the figure is drawn from data/revision_R1/figures/FigureS1_national_trends_data.csv.
Run as-is in VS Code (settings below)."""
import os
import numpy as np
import pandas as pd
import openpyxl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# ============================== settings ==============================
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
XLSX = os.path.join(REPO, "data", "external", "2024nenpo1.xlsx")  # annual report 2024 (Excel tables), optional
DATA_CSV = os.path.join(REPO, "data", "revision_R1", "figures", "FigureS1_national_trends_data.csv")
SHEET = "表4"
OUT_DIR = os.path.join(REPO, "outputs", "figures_R1")
OUT_STEM = "FigureS1_national_trends"
YEARS = list(range(1985, 2025))
CALIBRATION_PERIOD = (2004, 2023)
FIG_SIZE_IN = (7.09, 5.2)     # 180 mm wide
FONT_PT = 7
DPI = 300
# validated categorical colours (dataviz reference palette, light surface)
COL_HIV, COL_AIDS = "#2a78d6", "#eb6834"
COL_MSM, COL_HET, COL_OTHER = "#2a78d6", "#eb6834", "#1baf7a"
# ======================================================================

ROUTE_MSM = "同性間の性的接触"
ROUTE_HET = "異性間の性的接触"
ROUTES_OTHER = ["静注薬物使用", "母子感染", "その他", "不明"]


def read_table4():
    ws = openpyxl.load_workbook(XLSX, data_only=True)[SHEET]
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    hdr_i = next(i for i, r in enumerate(rows) if 1985 in r and 2024 in r)
    col = {y: rows[hdr_i].index(y) for y in YEARS}
    diag = nat = None
    recs = []
    for r in rows[hdr_i + 1:]:
        a, b, c = (str(x).strip() if x is not None else "" for x in r[:3])
        if a in ("HIV", "AIDS"):
            diag = a
        if b:
            nat = b
        route = c.replace("　", "").split("*")[0].strip()
        if not route or diag is None or not isinstance(r[col[1985]], (int, float)):
            continue
        recs.append({"diag": diag, "nationality": nat, "route": route,
                     **{y: float(r[col[y]] or 0) for y in YEARS}})
    df = pd.DataFrame(recs)
    tot = df[df["nationality"] == "合計"]
    return tot


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    if os.path.exists(XLSX):
        tot = read_table4()
        def series(diag, routes):
            sub = tot[(tot["diag"].isin(diag)) & (tot["route"].isin(routes))]
            return sub[YEARS].sum(axis=0).to_numpy()
        hiv = series(["HIV"], ["合計"]); aids = series(["AIDS"], ["合計"])
        both = ["HIV", "AIDS"]
        msm = series(both, [ROUTE_MSM]); het = series(both, [ROUTE_HET]); oth = series(both, ROUTES_OTHER)
        unk = series(both, ["不明"])
    else:
        d0 = pd.read_csv(DATA_CSV).set_index("year").loc[YEARS]
        hiv, aids = d0["HIV"].to_numpy(float), d0["AIDS"].to_numpy(float)
        msm, het, oth, unk = (d0[c].to_numpy(float) for c in ("MSM", "heterosexual", "other_or_unknown", "unknown_route"))
    total = hiv + aids
    assert np.allclose(msm + het + oth, total), "route categories do not add up to the total"

    data = pd.DataFrame({"year": YEARS, "HIV": hiv, "AIDS": aids, "total": total,
                         "MSM": msm, "heterosexual": het, "other_or_unknown": oth, "unknown_route": unk})
    data["MSM_share_all"] = data["MSM"] / data["total"]
    data["MSM_share_known_route"] = data["MSM"] / (data["total"] - data["unknown_route"])
    data.to_csv(os.path.join(OUT_DIR, OUT_STEM + "_data.csv"), index=False)

    fs = FONT_PT
    rc = {"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
          "font.size": fs, "axes.labelsize": fs, "xtick.labelsize": fs - 0.5, "ytick.labelsize": fs - 0.5,
          "legend.fontsize": fs, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
          "axes.spines.top": False, "axes.spines.right": False}
    x = np.array(YEARS)
    with plt.rc_context(rc):
        fig, axes = plt.subplots(2, 1, figsize=FIG_SIZE_IN, sharex=True)
        bar = dict(width=0.8, edgecolor="white", linewidth=0.4, zorder=3)

        ax = axes[0]
        ax.axvspan(CALIBRATION_PERIOD[0] - 0.5, CALIBRATION_PERIOD[1] + 0.5, color="#000000", alpha=0.05, lw=0, zorder=0)
        ax.text((CALIBRATION_PERIOD[0] + CALIBRATION_PERIOD[1]) / 2, 1.0, "Model calibration period (2004 to 2023)",
                transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=fs - 0.5, color="#52514e")
        ax.bar(x, hiv, color=COL_HIV, label="HIV (without AIDS when reported)", **bar)
        ax.bar(x, aids, bottom=hiv, color=COL_AIDS, label="AIDS (when first reported)", **bar)
        for y in (int(x[np.argmax(total)]), 2024):
            k = YEARS.index(y)
            ax.text(y, total[k] + 25, f"{int(total[k]):,}", ha="center", va="bottom", fontsize=fs - 0.5, color="#0b0b0b")
        ax.set_ylabel("New reports per year")
        ax.legend(loc="upper left", frameon=False)

        ax = axes[1]
        ax.axvspan(CALIBRATION_PERIOD[0] - 0.5, CALIBRATION_PERIOD[1] + 0.5, color="#000000", alpha=0.05, lw=0, zorder=0)
        ax.bar(x, msm, color=COL_MSM, label="Sexual contact between men", **bar)
        ax.bar(x, het, bottom=msm, color=COL_HET, label="Heterosexual contact", **bar)
        ax.bar(x, oth, bottom=msm + het, color=COL_OTHER, label="Other or unknown route", **bar)
        ax.set_ylabel("New reports per year")
        ax.set_xlabel("Year of report")
        ax.legend(loc="upper left", frameon=False)

        for k, ax in enumerate(axes):
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:,.0f}"))
            ax.grid(True, axis="y", alpha=0.3, linewidth=0.5, zorder=0)
            ax.set_xlim(YEARS[0] - 0.8, YEARS[-1] + 0.8)
            ax.set_xticks(list(range(1985, 2025, 5)) + [2024])
            ax.text(-0.075, 1.02, "ab"[k], transform=ax.transAxes, fontsize=fs + 3, fontweight="bold",
                    va="bottom", ha="left")
        axes[0].set_ylim(0, max(total) * 1.12)
        axes[1].set_ylim(0, max(total) * 1.12)
        fig.tight_layout()
        png = os.path.join(OUT_DIR, OUT_STEM + ".png")
        tif = os.path.join(OUT_DIR, OUT_STEM + ".tiff")
        fig.savefig(png, dpi=DPI)
        fig.savefig(tif, dpi=DPI, pil_kwargs={"compression": "tiff_lzw"})
        plt.close(fig)

    k24 = YEARS.index(2024); kpk = int(np.argmax(total))
    cal = data[(data["year"] >= CALIBRATION_PERIOD[0]) & (data["year"] <= CALIBRATION_PERIOD[1])]
    print(f"saved {png} and .tiff")
    print(f"peak {YEARS[kpk]}: {int(total[kpk]):,}; 2024: {int(total[k24]):,} (HIV {int(hiv[k24])}, AIDS {int(aids[k24])})")
    print(f"MSM share 2024: all {data.loc[k24,'MSM_share_all']:.1%}, known route {data.loc[k24,'MSM_share_known_route']:.1%}")
    print(f"MSM share of known-route reports {CALIBRATION_PERIOD[0]} to {CALIBRATION_PERIOD[1]}: "
          f"min {cal['MSM_share_known_route'].min():.1%} ({int(cal.loc[cal['MSM_share_known_route'].idxmin(),'year'])}), "
          f"max {cal['MSM_share_known_route'].max():.1%}")
    print(f"MSM share of all reports {CALIBRATION_PERIOD[0]} to {CALIBRATION_PERIOD[1]}: min {cal['MSM_share_all'].min():.1%}, max {cal['MSM_share_all'].max():.1%}")


if __name__ == "__main__":
    main()

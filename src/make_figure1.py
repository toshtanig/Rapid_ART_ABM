#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build the schematic of the certification pathway and ART-timing strategies (Figure 1)."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ============================ CONFIG ============================
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# --- Data location: set BASE_144 to your local 144_scenario folder ---
BASE_144 = r"/path/to/project/144_scenario"
PROJ_169 = os.path.dirname(BASE_144)
BASE_169 = BASE_144
OUTPUT_DIR = os.path.join(BASE_169, "figures_tiff")
OUTNAME    = "Figure1"
WIDTH_MM   = 180
DPI        = 300
SAVE_TIFF      = True
SAVE_PNG_PROOF = True
SHOW_PANEL_TITLES = False   # journal rule: no titles inside the artwork (the legend describes the panels)
# ===============================================================

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42

GREY = "#666666"; LINE = "#333333"; FILL = "#FFFFFF"; ACC = "#000000"

def box(ax, x, y, w, h, text, fs=8.0, fill=FILL, edge=GREY, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle="round,pad=0.006,rounding_size=0.012",
                 linewidth=0.9, edgecolor=edge, facecolor=fill, zorder=2))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs,
            zorder=3, fontweight=("bold" if bold else "normal"), color=ACC)

def arrow(ax, x1, y1, x2, y2, style="-|>", color=LINE, lw=1.0, ls="-"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                 mutation_scale=9, lw=lw, color=color, linestyle=ls, zorder=1))

def bracket(ax, x1, x2, ytop, ytick):
    """Square bracket spanning x1..x2 at ytop, with downward ticks to ytick."""
    ax.plot([x1, x1, x2, x2], [ytick, ytop, ytop, ytick],
            color=LINE, lw=0.9, solid_capstyle="butt", zorder=2)

def main():
    mm = 1/25.4
    w_in = WIDTH_MM*mm
    fig = plt.figure(figsize=(w_in, w_in*0.80))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.02], hspace=0.16,
                          left=0.015, right=0.985, top=0.98, bottom=0.02)

    # ----- Panel A: certification pathway (current vs proposed) -----
    axA = fig.add_subplot(gs[0]); axA.set_xlim(0, 1); axA.set_ylim(0, 1); axA.axis("off")
    axA.text(0.0, 0.985, "a", fontsize=11, fontweight="bold", va="top")
    if SHOW_PANEL_TITLES:
        axA.text(0.03, 0.98, "Disability-certification pathway and the ART-initiation strategies",
                 fontsize=8.5, fontweight="bold", va="top")
    axA.text(0.02, 0.905, "CURRENT STANDARD OF CARE  (1–2 month delay)",
             fontsize=7.6, fontweight="bold", color=LINE, va="top")
    # four-week bracket over the two laboratory assessments
    axA.text(0.2925, 0.79, "≥ 4 weeks apart (statutory requirement)",
             ha="center", fontsize=7.6, fontweight="bold", color=ACC)
    bracket(axA, 0.2125, 0.3725, 0.755, 0.72)
    yA = 0.55; hA = 0.17
    box(axA, 0.02, yA, 0.105, hA, "Entry into\ncare", bold=True, edge=ACC)
    box(axA, 0.155, yA, 0.115, hA, "Laboratory\nassessment 1")
    box(axA, 0.315, yA, 0.115, hA, "Laboratory\nassessment 2")
    box(axA, 0.475, yA, 0.115, hA, "Disability\ncertificate")
    box(axA, 0.635, yA, 0.115, hA, "Medical-care\nsubsidy")
    box(axA, 0.795, yA, 0.125, hA, "ART initiation\n(delayed)", bold=True, edge=ACC)
    for x1, x2 in [(0.125, 0.153), (0.27, 0.313), (0.43, 0.473), (0.59, 0.633), (0.75, 0.793)]:
        arrow(axA, x1, yA+hA/2, x2, yA+hA/2)
    axA.text(0.3725, 0.47, "must meet certification criteria*",
             ha="center", fontsize=7.4, style="italic", color=GREY)
    axA.text(0.02, 0.385, "PROPOSED: IMMEDIATE ART  (0 month)",
             fontsize=7.6, fontweight="bold", color=ACC, va="top")
    yA2 = 0.16; hA2 = 0.16
    box(axA, 0.02, yA2, 0.105, hA2, "Entry into\ncare", bold=True, edge=ACC)
    box(axA, 0.30, yA2, 0.30, hA2, "ART initiation at diagnosis\nirrespective of CD4 and HIV-RNA",
        fs=7.4, bold=True, edge=ACC)
    arrow(axA, 0.125, yA2+hA2/2, 0.298, yA2+hA2/2)
    axA.text(0.625, 0.27,
             "Relaxing / reforming the certification pathway:\n"
             "no second assessment, no four-week wait;\n"
             "certificate and subsidy at diagnosis.",
             fontsize=7.0, style="italic", color=GREY, va="top")
    axA.text(0.02, 0.04,
             "*Certification criteria (except AIDS): CD4 ≥500 cells/µL requires HIV-RNA ≥5000 copies/mL;  "
             "HIV-RNA <5000 copies/mL requires CD4 <500 cells/µL.",
             fontsize=6.5, color=GREY, va="top")

    # ----- Panel B: model structure -----
    axB = fig.add_subplot(gs[1]); axB.set_xlim(0, 1); axB.set_ylim(0, 1); axB.axis("off")
    axB.text(0.0, 0.98, "b", fontsize=11, fontweight="bold", va="top")
    if SHOW_PANEL_TITLES:
        axB.text(0.03, 0.975, "Agent-based dynamic transmission and care model",
             fontsize=8.5, fontweight="bold", va="top")
    yB = 0.55; hB = 0.17
    box(axB, 0.02, yB, 0.135, hB, "Susceptible\nMSM + high-risk\nheterosexual", fs=6.6, bold=True, edge=ACC)
    box(axB, 0.185, yB, 0.135, hB, "Force of\ninfection", fs=6.8)
    box(axB, 0.35, yB, 0.125, hB, "Infection\nCD4 / VL", fs=6.8)
    box(axB, 0.505, yB, 0.125, hB, "Testing &\ndiagnosis", fs=6.8)
    box(axB, 0.66, yB, 0.135, hB, "ART initiation\n0 / 1 / 2-month\ndelay", fs=6.6, bold=True, edge=ACC)
    box(axB, 0.825, yB, 0.13, hB, "Viral suppression\n(U = U)", fs=6.8)
    for x1, x2 in [(0.155, 0.183), (0.32, 0.348), (0.475, 0.503), (0.63, 0.658), (0.795, 0.823)]:
        arrow(axB, x1, yB+hB/2, x2, yB+hB/2)
    # treatment as prevention: viral suppression -> force of infection (arc above the row)
    axB.add_patch(FancyArrowPatch((0.89, yB+hB), (0.2525, yB+hB),
                  connectionstyle="arc3,rad=0.15", arrowstyle="-|>",
                  mutation_scale=9, lw=0.9, color=GREY, linestyle=(0, (3, 2)), zorder=1))
    axB.text(0.57, 0.905, "treatment as prevention: viral suppression reduces the force of infection",
             ha="center", fontsize=7.2, style="italic", color=GREY)
    # certification-eligibility gate -> ART initiation (delay arms)
    box(axB, 0.62, 0.16, 0.225, hB,
        "Certification-eligibility check\n(1 / 2-month arms only;\nelse monthly re-evaluation)", fs=7.0)
    arrow(axB, 0.7325, 0.33, 0.7325, yB-0.003, color=GREY, lw=0.8, ls=(0, (3, 2)))
    # background processes
    box(axB, 0.02, 0.16, 0.165, hB, "PrEP: reduces HIV\nacquisition (95% effective)\nMSM / heterosexual", fs=6.4)
    box(axB, 0.205, 0.16, 0.175, hB, "Mortality\nbackground + HIV-related\n(CD4 / AIDS dependent)", fs=6.4)
    box(axB, 0.40, 0.16, 0.185, hB, "Demographic turnover\nnew entrants 2%/yr;\nforeign-born HIV inflow", fs=6.4)
    # PrEP acts on susceptibles
    arrow(axB, 0.1025, 0.33, 0.1025, yB-0.003, color=GREY, lw=0.8, ls=(0, (3, 2)))
    axB.text(0.5, 0.062,
             "Costs (payer perspective) and QALYs accrued over the whole population, discounted 2%/yr.\n"
             "Monthly time step, 20-year burn-in + 40-year horizon.",
             ha="center", va="top", fontsize=6.8, color=ACC, linespacing=1.5)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    if SAVE_TIFF:
        fig.savefig(os.path.join(OUTPUT_DIR, OUTNAME+".tiff"), dpi=DPI, format="tiff",
                    pil_kwargs={"compression": "tiff_lzw"})
    if SAVE_PNG_PROOF:
        fig.savefig(os.path.join(OUTPUT_DIR, OUTNAME+"_proof.png"), dpi=300)
    plt.close(fig)
    print(f"saved {OUTNAME} (width {WIDTH_MM} mm, {DPI} dpi) to {OUTPUT_DIR}")

if __name__ == "__main__":
    main()

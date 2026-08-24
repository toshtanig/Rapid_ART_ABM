#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compose multi-panel figures for the manuscript."""

import os
from PIL import Image, ImageDraw, ImageFont

# ============================ CONFIG ============================
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# --- Data location: set BASE_144 to your local 144_scenario folder ---
BASE_144 = r"/path/to/20250315_169/144_scenario"
PROJ_169 = os.path.dirname(BASE_144)
BASE = BASE_144
PREV_DIR = os.path.join(BASE, "prevention_dominance_figures")
CEA_DIR  = os.path.join(BASE, "cea_param_uncertainty")
OUT_DIR  = os.path.join(BASE, "figures_tiff")          # composites (PNG + TIFF) are written here

TIFF_DPI        = 300
MULTIPANEL_MM   = 180.0     # journal multi-panel width
GAP_PX          = 24        # white gap between stacked panels
LETTER_PX       = 54        # panel-letter font size (relative to panel width ~2200)
# ===============================================================

os.makedirs(OUT_DIR, exist_ok=True)

def _font(size):
    for p in ["/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
              "C:/Windows/Fonts/arialbd.ttf", "C:/Windows/Fonts/Arial.ttf",
              "/Library/Fonts/Arial Bold.ttf"]:
        if os.path.exists(p):
            try:
                return ImageFont.truetype(p, size)
            except Exception:
                pass
    return ImageFont.load_default()

def _load(path):
    im = Image.open(path)
    im.load()
    return im.convert("RGB")

def _scale_w(im, w):
    if im.width == w:
        return im
    h = round(im.height * w / im.width)
    return im.resize((w, h), Image.LANCZOS)

def vstack_with_letters(panels, out_png, out_tiff):
    """panels: list of (image_path, [(letter, x_frac, y_frac), ...]).
    Scales all panels to a common width, stacks vertically, draws panel letters."""
    imgs = [(_load(p), letters) for p, letters in panels]
    W = max(im.width for im, _ in imgs)
    imgs = [(_scale_w(im, W), letters) for im, letters in imgs]
    total_h = sum(im.height for im, _ in imgs) + GAP_PX * (len(imgs) - 1)
    canvas = Image.new("RGB", (W, total_h), "white")
    draw = ImageDraw.Draw(canvas)
    font = _font(LETTER_PX)
    y = 0
    for im, letters in imgs:
        canvas.paste(im, (0, y))
        for letter, xf, yf in letters:
            x = int(xf * W) + 6
            ly = y + int(yf * im.height) + 6
            # white halo for legibility
            for dx, dy in [(-2, 0), (2, 0), (0, -2), (0, 2)]:
                draw.text((x + dx, ly + dy), letter, font=font, fill="white")
            draw.text((x, ly), letter, font=font, fill="black")
        y += im.height + GAP_PX
    canvas.save(out_png, "PNG")
    # publication TIFF: 180 mm wide at 300 dpi
    target_w = int(round(MULTIPANEL_MM / 25.4 * TIFF_DPI))
    tif = _scale_w(canvas, target_w)
    tif.save(out_tiff, "TIFF", compression="tiff_lzw", dpi=(TIFF_DPI, TIFF_DPI))
    print(f"  saved {os.path.basename(out_png)} and {os.path.basename(out_tiff)} "
          f"({canvas.size} -> TIFF {tif.size})")

def single(in_path, out_png, out_tiff, width_mm=107.0):
    im = _load(in_path)
    im.save(out_png, "PNG")
    target_w = int(round(width_mm / 25.4 * TIFF_DPI))
    tif = _scale_w(im, target_w)
    tif.save(out_tiff, "TIFF", compression="tiff_lzw", dpi=(TIFF_DPI, TIFF_DPI))
    print(f"  saved {os.path.basename(out_png)} and {os.path.basename(out_tiff)}")

def main():
    print("Composing main-text figures ...")
    # Figure 3: cost-effectiveness. a = CE plane (top); b, c = QALY / cost-saving bars (bottom).
    vstack_with_letters(
        [(os.path.join(PREV_DIR, "D_ce_plane.jpg"), [("a", 0.0, 0.0)]),
         (os.path.join(PREV_DIR, "D_qaly_cost_bars.jpg"), [("b", 0.0, 0.0), ("c", 0.5, 0.0)])],
        os.path.join(OUT_DIR, "Figure3.png"), os.path.join(OUT_DIR, "Figure3.tiff"))

    # Figure 4: a = tornado (top); b = break-even curve (bottom).
    vstack_with_letters(
        [(os.path.join(CEA_DIR, "fig_tornado.png"), [("a", 0.0, 0.0)]),
         (os.path.join(PREV_DIR, "D_breakeven_curve.jpg"), [("b", 0.0, 0.0)])],
        os.path.join(OUT_DIR, "Figure4.png"), os.path.join(OUT_DIR, "Figure4.tiff"))

    print("Composing supplementary figures ...")
    # Supplementary Figure S6 (PSA): a = CE plane (PSA); b = NMB / cost distribution.
    vstack_with_letters(
        [(os.path.join(CEA_DIR, "fig_ce_plane_psa.png"), [("a", 0.0, 0.0)]),
         (os.path.join(CEA_DIR, "fig_psa_histogram.png"), [("b", 0.0, 0.0)])],
        os.path.join(OUT_DIR, "FigureS6_PSA.png"), os.path.join(OUT_DIR, "FigureS6_PSA.tiff"))

    # Supplementary Figure S4 (CEAC): single panel.
    single(os.path.join(CEA_DIR, "fig_ceac.png"),
           os.path.join(OUT_DIR, "FigureS4_CEAC.png"), os.path.join(OUT_DIR, "FigureS4_CEAC.tiff"))

    print("Done. Composites written to:", OUT_DIR)

if __name__ == "__main__":
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Export tables and figures in the journal submission format."""
import io
import os
import matplotlib
import matplotlib.pyplot as plt

try:
    from PIL import Image
    _HAVE_PIL = True
except Exception:
    _HAVE_PIL = False

WIDTH_DEFAULT_MM = 180
WIDTH_MM = {
    "A_calibration": 180, "B_projection": 180,
    "C_delay_cdf_diag": 107, "C_delay_cdf_eligible": 107, "C_high_cd4_group": 107,
    "D_ce_plane": 107, "D_breakeven_curve": 107, "D_nmb_ci": 107,
    "D_qaly_cost_bars": 107, "D_distributions": 180,
    "fig_tornado": 180, "fig_ceac": 107, "fig_ce_plane_psa": 107,
    "fig_psa_histogram": 107, "fig_breakeven_detailed": 180, "fig_breakeven_summary": 107,
    "convergence": 180, "survival_adj_sensitivity": 180,
    "Figure1": 180,
}

_PATCHED = False

def apply_style():
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

def enable_tiff_export(outdir=None, dpi=300, width_map=None, also_keep_original=True):
    global _PATCHED
    apply_style()
    if _PATCHED:
        return
    if not _HAVE_PIL:
        print("[lancet_export] Pillow (PIL) not available; TIFF output disabled (pip install pillow).")
        return
    wm = dict(WIDTH_MM)
    if width_map:
        wm.update(width_map)
    orig_savefig = matplotlib.figure.Figure.savefig

    def savefig(self, fname, *args, **kwargs):

        result = orig_savefig(self, fname, *args, **kwargs) if also_keep_original else None
        try:
            stem = os.path.splitext(os.path.basename(str(fname)))[0]
            d = outdir or (os.path.dirname(str(fname)) or ".")
            os.makedirs(d, exist_ok=True)
            w0, _h0 = self.get_size_inches()
            target_w_in = wm.get(stem, WIDTH_DEFAULT_MM) / 25.4

            render_dpi = target_w_in * dpi / float(w0)
            buf = io.BytesIO()
            orig_savefig(self, buf, format="tiff", dpi=render_dpi,
                         facecolor=self.get_facecolor())
            buf.seek(0)
            im = Image.open(buf)

            im.save(os.path.join(d, stem + ".tiff"), format="TIFF",
                    dpi=(dpi, dpi), compression="tiff_lzw")
            buf.close()
        except Exception as e:
            print(f"[lancet_export] TIFF export skipped for {fname}: {e}")
        return result

    matplotlib.figure.Figure.savefig = savefig
    _PATCHED = True
    print("[lancet_export] TIFF export enabled (300 dpi, LZW; widths: multipanel 180 mm, single 107 mm; layout unchanged)")

# -*- coding: utf-8 -*-
# R10-R4/R6/R17: merge per-panel TIFFs into one file per figure number,
# flatten RGBA->RGB (white), tag 300 dpi, drop legacy _rev suffixes.
import os
from PIL import Image

SRC = r"results\figures\tiff_300"
OUT_TMP = r"results\figures\tiff_300_final"

FIGS = [
    ("Fig1_flowchart", ["Fig1_flowchart"]),
    ("Fig2_pca_volcano_heatmap", ["Fig2AB_pca_volcano_rev", "Fig2C_heatmap_rev"]),
    ("Fig3_wgcna", ["Fig3_wgcna_rev"]),
    ("Fig4_enrichment", ["Fig4_enrichment_rev"]),
    ("Fig5_candidates", ["Fig5A_venn_rev", "Fig5B_candidates_rev"]),
    ("Fig6_ml_screening", ["Fig6_ml_screening_rev"]),
    ("Fig7_validation", ["Fig7A_validation_roc_rev", "Fig7B_calibration_GSE142025", "Fig7C_dca_GSE142025"]),
    ("Fig8_immune", ["Fig8A_immune_box_rev", "Fig8B_corr_rev"]),
    ("Fig9_scrna_localization", ["Fig9_scrna_pdk_localization"]),
]

def flatten_rgb(im):
    if im.mode == "RGBA":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[3])
        return bg
    return im.convert("RGB")

os.makedirs(OUT_TMP, exist_ok=True)
for out_name, parts in FIGS:
    paths = [os.path.join(SRC, p + ".tif") for p in parts]
    for p in paths:
        assert os.path.exists(p), "missing " + p
    ims = [Image.open(p) for p in paths]
    widths = [im.width for im in ims]
    modes = [im.mode for im in ims]
    target = min(widths)
    if target < 2400:                       # keep >= 8 in @300dpi when possible
        target = min(2400, max(widths))
    scaled = []
    for im in ims:
        if im.width != target:
            h = round(im.height * target / im.width)
            im = im.resize((target, h), Image.LANCZOS)
        scaled.append(flatten_rgb(im))
    gap = 40 if len(scaled) > 1 else 0
    H = sum(im.height for im in scaled) + gap * (len(scaled) - 1)
    canvas = Image.new("RGB", (target, H), (255, 255, 255))
    y = 0
    for im in scaled:
        canvas.paste(im, (0, y))
        y += im.height + gap
    out = os.path.join(OUT_TMP, out_name + ".tif")
    canvas.save(out, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
    print("%-28s <- %-55s %s  src=%s modes=%s" % (out_name, ",".join(parts), canvas.size, widths, modes))

print("--- verify ---")
for f in sorted(os.listdir(OUT_TMP)):
    im = Image.open(os.path.join(OUT_TMP, f))
    print(f, im.size, im.mode, "dpi=", im.info.get("dpi"))
print("DONE", len(os.listdir(OUT_TMP)), "files")

# -*- coding: utf-8 -*-
# R11 assembly: rebuild finals Fig2 (AB+C), Fig6 (flatten), Fig8 (A+B) from fresh PNG sources,
# sync both figure dirs + judge previews.
import os
import shutil
from PIL import Image

SRC = r"results\figures"
FIN = r"results\figures\tiff_300_final"
JUDGE = r"results\figures\judge_png"

def flatten_rgb(im):
    if im.mode == "RGBA":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[3])
        return bg
    return im.convert("RGB")

def stack(names, out_name, gap=40):
    ims = [flatten_rgb(Image.open(os.path.join(SRC, n))) for n in names]
    widths = [im.width for im in ims]
    target = min(widths)
    if target < 2400:
        target = min(2400, max(widths))
    scaled = []
    for im in ims:
        if im.width != target:
            im = im.resize((target, round(im.height * target / im.width)), Image.LANCZOS)
        scaled.append(im)
    H = sum(im.height for im in scaled) + gap * (len(scaled) - 1)
    canvas = Image.new("RGB", (target, H), (255, 255, 255))
    y = 0
    for im in scaled:
        canvas.paste(im, (0, y))
        y += im.height + gap
    out = os.path.join(FIN, out_name)
    canvas.save(out, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
    canvas.resize((1400, round(H * 1400 / target))).save(os.path.join(JUDGE, out_name.replace(".tif", ".png")))
    print(out_name, "<-", names, canvas.size, "src widths", widths)

stack(["Fig2AB_pca_volcano_rev.png", "Fig2C_heatmap_rev.png"], "Fig2_pca_volcano_heatmap.tif")
stack(["Fig6_ml_screening_rev.png"], "Fig6_ml_screening.tif", gap=0)
stack(["Fig8A_immune_box_rev.png", "Fig8B_corr_rev.png"], "Fig8_immune.tif")

sub = [d for d in os.listdir(".") if not d[0].isascii() and os.path.isdir(d)][0]
for name in ["Fig2_pca_volcano_heatmap.tif", "Fig6_ml_screening.tif", "Fig8_immune.tif"]:
    for tgt in [r"results\figures\tiff_300", os.path.join(sub, "02_图件", "tiff_300")]:
        shutil.copy2(os.path.join(FIN, name), os.path.join(tgt, name))
print("ASSEMBLED+SYNCED")

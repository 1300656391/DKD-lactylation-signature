# -*- coding: utf-8 -*-
# R12 assembly: rebuild Fig6 (flatten) + Fig7 (A + new B + C, restamp B/C below titles),
# pixel-verify the two defects are gone, sync dirs + previews.
import os
import shutil
from PIL import Image, ImageDraw, ImageFont

SRC = r"results\figures"
FIN = r"results\figures\tiff_300_final"
JUDGE = r"results\figures\judge_png"
font_path = r"C:\Windows\Fonts\arialbd.ttf"

def flatten_rgb(im):
    if im.mode == "RGBA":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[3])
        return bg
    return im.convert("RGB")

# ---- Fig6: single panel, flatten ----
f6 = flatten_rgb(Image.open(os.path.join(SRC, "Fig6_ml_screening_rev.png")))
p6 = os.path.join(FIN, "Fig6_ml_screening.tif")
f6.save(p6, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
f6.resize((1400, round(f6.height * 1400 / f6.width))).save(os.path.join(JUDGE, "Fig6_ml_screening.png"))
g6 = f6.convert("L"); px6 = g6.load()
edge6 = sum(1 for yy in range(f6.height) for xx in range(f6.width - 6, f6.width) if px6[xx, yy] < 128)
print("Fig6 final:", f6.size, "| rightmost-6px dark:", edge6)

# ---- Fig7: A + new B + C, stamps below title bands ----
A = flatten_rgb(Image.open(os.path.join(SRC, "Fig7A_validation_roc_rev.png")))
B = flatten_rgb(Image.open(os.path.join(SRC, "Fig7B_calibration_GSE142025.png")))
C = flatten_rgb(Image.open(os.path.join(SRC, "Fig7C_dca_GSE142025.png")))
W = A.width; gap = 40; half = (W - gap) // 2
B2 = B.resize((half, round(B.height * half / B.width)), Image.LANCZOS)
C2 = C.resize((half, round(C.height * half / C.width)), Image.LANCZOS)
H = A.height + gap + max(B2.height, C2.height)
fig7 = Image.new("RGB", (W, H), (255, 255, 255))
fig7.paste(A, (0, 0))
yB = A.height + gap
fig7.paste(B2, (0, yB))
fig7.paste(C2, (half + gap, yB))
d = ImageDraw.Draw(fig7)
size7 = max(64, round(W * 0.045))
def band_y(panel):
    band = panel.convert("L").crop((0, 0, min(1200, panel.width), 300))
    bp = band.load()
    rows = [y for y in range(band.height) if any(bp[x, y] < 128 for x in range(0, band.width, 2))]
    return (max(rows) + 14) if rows else 10
yb = band_y(B2); yc = band_y(C2)
d.text((14, yB + yb), "B", fill=(0, 0, 0), font=ImageFont.truetype(font_path, size7))
d.text((half + gap + 14, yB + yc), "C", fill=(0, 0, 0), font=ImageFont.truetype(font_path, size7))
print("B/C y-offsets:", yb, yc)
p7 = os.path.join(FIN, "Fig7_validation.tif")
fig7.save(p7, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
fig7.resize((1400, round(H * 1400 / W))).save(os.path.join(JUDGE, "Fig7_validation.png"))

# ---- sync ----
sub = [d2 for d2 in os.listdir(".") if not d2[0].isascii() and os.path.isdir(d2)][0]
for name in ["Fig6_ml_screening.tif", "Fig7_validation.tif"]:
    for tgt in [r"results\figures\tiff_300", os.path.join(sub, "02_图件", "tiff_300")]:
        shutil.copy2(os.path.join(FIN, name), os.path.join(tgt, name))
print("ASSEMBLED+SYNCED")

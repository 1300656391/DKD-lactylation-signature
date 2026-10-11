# -*- coding: utf-8 -*-
# Rebuild Fig5 final: new Fig5A (177-key venn) + Fig5B, stamp "A" matched to native B (size 50),
# sync dirs + preview.
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

A = flatten_rgb(Image.open(os.path.join(SRC, "Fig5A_venn_rev.png")))
B = flatten_rgb(Image.open(os.path.join(SRC, "Fig5B_candidates_rev.png")))
widths = [A.width, B.width]
target = min(widths)
if target < 2400:
    target = min(2400, max(widths))
def scale_to(im, t):
    if im.width != t:
        return im.resize((t, round(im.height * t / im.width)), Image.LANCZOS)
    return im
A2, B2 = scale_to(A, target), scale_to(B, target)
gap = 40
H = A2.height + gap + B2.height
fig5 = Image.new("RGB", (target, H), (255, 255, 255))
fig5.paste(A2, (0, 0))
fig5.paste(B2, (0, A2.height + gap))
d = ImageDraw.Draw(fig5)
d.text((16, 14), "A", fill=(0, 0, 0), font=ImageFont.truetype(font_path, 50))
out = os.path.join(FIN, "Fig5_candidates.tif")
fig5.save(out, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
fig5.resize((1400, round(H * 1400 / target))).save(os.path.join(JUDGE, "Fig5_candidates.png"))
print("Fig5 rebuilt:", fig5.size, "src widths", widths)

sub = [d for d in os.listdir(".") if not d[0].isascii() and os.path.isdir(d)][0]
for tgt in [r"results\figures\tiff_300", os.path.join(sub, "02_图件", "tiff_300")]:
    shutil.copy2(out, os.path.join(tgt, "Fig5_candidates.tif"))
print("SYNCED")

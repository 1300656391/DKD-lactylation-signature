# -*- coding: utf-8 -*-
# Fig7 relayout v3: rebuild from results\figures PNG sources (panel TIFFs were cleaned),
# [A ROC full width] / [B calibration | C DCA side by side] + bold B/C stamps
# placed automatically below each bottom-panel title band.
import os
from PIL import Image, ImageDraw, ImageFont

SRC = r"results\figures"
OUT = r"results\figures\tiff_300_final\Fig7_validation.tif"

def flatten_rgb(im):
    if im.mode == "RGBA":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[3])
        return bg
    return im.convert("RGB")

A = flatten_rgb(Image.open(os.path.join(SRC, "Fig7A_validation_roc_rev.png")))
B = flatten_rgb(Image.open(os.path.join(SRC, "Fig7B_calibration_GSE142025.png")))
C = flatten_rgb(Image.open(os.path.join(SRC, "Fig7C_dca_GSE142025.png")))
print("src:", A.size, B.size, C.size)

W = A.width
gap = 40
half = (W - gap) // 2
B2 = B.resize((half, round(B.height * half / B.width)), Image.LANCZOS)
C2 = C.resize((half, round(C.height * half / C.width)), Image.LANCZOS)

H = A.height + gap + max(B2.height, C2.height)
canvas = Image.new("RGB", (W, H), (255, 255, 255))
canvas.paste(A, (0, 0))
yB = A.height + gap
canvas.paste(B2, (0, yB))
canvas.paste(C2, (half + gap, yB))

font = ImageFont.truetype(r"C:\Windows\Fonts\arialbd.ttf", max(64, round(W * 0.045)))
draw = ImageDraw.Draw(canvas)

def stamp_below_title(ch, x0, panel_img):
    # find the title text band at the top of the panel; stamp the letter just below it
    band = panel_img.convert("L").crop((0, 0, min(1200, panel_img.width), 300))
    px = band.load()
    rows = [y for y in range(band.height)
            if any(px[x, y] < 128 for x in range(0, band.width, 2))]
    y_local = (max(rows) + 14) if rows else 10
    draw.text((x0 + 14, yB + y_local), ch, fill=(0, 0, 0), font=font)
    return y_local

yb_off = stamp_below_title("B", 0, B2)
yc_off = stamp_below_title("C", half + gap, C2)
print("B stamp y-offset:", yb_off, " C stamp y-offset:", yc_off)

canvas.save(OUT, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
im = Image.open(OUT)
print("Fig7 final:", im.size, im.mode, im.info.get("dpi"))

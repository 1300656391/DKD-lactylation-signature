# -*- coding: utf-8 -*-
# Round-2 fix (recover from bad stamping):
# 1) rebuild Fig5 from PNGs (original merge params), stamp A sized to native B (tight-window measure)
# 2) rebuild Fig7 via py_fig7_relayout.py logic, then unify A size (tight erase + restamp)
import os
from PIL import Image, ImageDraw, ImageFont
import shutil

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

def glyph_bbox(img, box, thr=128):
    """isolated dark-glyph bbox in a TIGHT window; returns None or bbox, plus edge-touch flag"""
    g = img.convert("L").crop(box)
    w, h = g.size
    px = g.load()
    xs, ys = [], []
    for y in range(h):
        for x in range(w):
            if px[x, y] < thr:
                xs.append(x); ys.append(y)
    if not xs:
        return None, True
    bbox = (box[0] + min(xs), box[1] + min(ys), box[0] + max(xs), box[1] + max(ys))
    touches = (min(xs) == 0 or min(ys) == 0 or max(xs) == w - 1 or max(ys) == h - 1)
    return bbox, touches

# ---------- 1) rebuild Fig5 ----------
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
yB = A2.height + gap
# measure native B tag with TIGHT window at bottom-panel top-left corner
b_bbox, b_touch = glyph_bbox(fig5, (0, yB, 140, yB + 150))
print("Fig5 rebuilt", fig5.size, "| tight B bbox:", b_bbox, "| edge-touch:", b_touch)
size_a = 90
if b_bbox and not b_touch:
    bh = b_bbox[3] - b_bbox[1]
    size_a = max(40, int(bh * 1.45))
print("-> A font size:", size_a)
d = ImageDraw.Draw(fig5)
d.text((16, 14), "A", fill=(0, 0, 0), font=ImageFont.truetype(font_path, size_a))
p5 = os.path.join(FIN, "Fig5_candidates.tif")
fig5.save(p5, format="TIFF", dpi=(300, 300), compression="tiff_lzw")

# ---------- 2) rebuild Fig7 (same as py_fig7_relayout v3) ----------
A7 = flatten_rgb(Image.open(os.path.join(SRC, "Fig7A_validation_roc_rev.png")))
B7 = flatten_rgb(Image.open(os.path.join(SRC, "Fig7B_calibration_GSE142025.png")))
C7 = flatten_rgb(Image.open(os.path.join(SRC, "Fig7C_dca_GSE142025.png")))
W7 = A7.width; gap7 = 40; half = (W7 - gap7) // 2
B72 = B7.resize((half, round(B7.height * half / B7.width)), Image.LANCZOS)
C72 = C7.resize((half, round(C7.height * half / C7.width)), Image.LANCZOS)
H7 = A7.height + gap7 + max(B72.height, C72.height)
fig7 = Image.new("RGB", (W7, H7), (255, 255, 255))
fig7.paste(A7, (0, 0))
y7 = A7.height + gap7
fig7.paste(B72, (0, y7))
fig7.paste(C72, (half + gap7, y7))
# measure native A tag on SOURCE A7 (top-left tight window) BEFORE pasting context
a_src_bbox, a_touch = glyph_bbox(A7, (0, 0, 200, 170))
print("Fig7A source tight A bbox:", a_src_bbox, "| edge-touch:", a_touch)
d7 = ImageDraw.Draw(fig7)
# B/C stamps below title bands (same as accepted v3)
def band_stamp(ch, x0, panel):
    band = panel.convert("L").crop((0, 0, min(1200, panel.width), 300))
    px = band.load()
    rows = [y for y in range(band.height) if any(px[x, y] < 128 for x in range(0, band.width, 2))]
    y_local = (max(rows) + 14) if rows else 10
    return y_local
size7 = max(64, round(W7 * 0.045))
yb_off = band_stamp("B", 0, B72); yc_off = band_stamp("C", half + gap7, C72)
d7.text((14, y7 + yb_off), "B", fill=(0, 0, 0), font=ImageFont.truetype(font_path, size7))
d7.text((half + gap7 + 14, y7 + yc_off), "C", fill=(0, 0, 0), font=ImageFont.truetype(font_path, size7))
# unify A: erase native small A (tight bbox from source, panel pasted at origin) then stamp big A
if a_src_bbox and not a_touch:
    x0, y0, x1, y1 = a_src_bbox
    d7.rectangle([x0 - 6, y0 - 6, x1 + 6, y1 + 6], fill=(255, 255, 255))
    print("erased native A at", a_src_bbox)
else:
    print("WARN: native A bbox not isolated; leaving native A, skipping restamp")
    size7 = None
if size7:
    d7.text((14, 10), "A", fill=(0, 0, 0), font=ImageFont.truetype(font_path, size7))
p7 = os.path.join(FIN, "Fig7_validation.tif")
fig7.save(p7, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
print("B/C y-offsets:", yb_off, yc_off, "| A size:", size7)

# previews + sync
for name, im in [("Fig5_candidates", fig5), ("Fig7_validation", fig7)]:
    im.resize((1400, round(im.height * 1400 / im.width))).save(os.path.join(JUDGE, name + ".png"))
sub = [d for d in os.listdir(".") if not d[0].isascii() and os.path.isdir(d)][0]
for name in ["Fig5_candidates.tif", "Fig7_validation.tif"]:
    for tgt in [r"results\figures\tiff_300", os.path.join(sub, "02_图件", "tiff_300")]:
        shutil.copy2(os.path.join(FIN, name), os.path.join(tgt, name))
print("DONE+SYNCED")

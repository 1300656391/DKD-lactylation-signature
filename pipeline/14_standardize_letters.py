# 14_standardize_letters.py ----
# Final lettering pass over the assembled publication TIFFs: remove the
# renderer's native panel tags and redraw them as uniform Arial Bold
# lowercase letters (2.2 mm x-height, same-baseline alignment within rows,
# panel-corner placement), reproducing the submitted figure lettering.
# Fig 3 is first resampled to the submitted canvas (2100 -> 2544 px wide,
# the documented 1.2114x LANCZOS step matching the submission figures).
# Usage: python pipeline/14_standardize_letters.py   (run 13_figures.R and
# assemble_figures.py first; operates on pipeline/figures/*.tif in place).

import os
import sys
import time

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "figures")
FONT_DIRS = [r"C:\Windows\Fonts", "/usr/share/fonts/truetype/dejavu",
             "/usr/share/fonts/truetype/liberation", "/Library/Fonts"]
TARGET_MM = 2.2

# Lettering spec of the submitted figures.
#   scale: canvas resample applied before lettering (Fig 3 only)
#   zones: native-tag scan boxes (x0, x1, y0, y1); an isolated letter-sized
#          ink component inside each box is removed
#   rows : letters sharing a baseline; baseline given in final-canvas px
#   singles: top-anchored single letters (x, top)
#   x    : final x per letter; xh: x-height px at the submitted display width
FIGS = {
    "Fig2_pca_volcano_heatmap.tif": dict(
        xh=40,
        zones={"a": (0, 220, 0, 160), "b": (900, 1400, 0, 160), "c": (0, 220, 962, 1150)},
        rows={"ab": 62}, singles={"c": (34, 984)},
        x={"a": 34, "b": 1084, "c": 34}),
    "Fig3_wgcna.tif": dict(
        xh=38, scale=(2544 / 2100, 1210 / 999),
        zones={},  # renderer draws no panel tags; letters are added only
        rows={}, singles={"a": (21, 21), "b": (1099, 21)},
        x={"a": 21, "b": 1099}),
    "Fig4_candidates.tif": dict(
        xh=35,
        zones={"a": (0, 220, 0, 260), "b": (0, 280, 1870, 2300)},
        rows={}, singles={"a": (19, 19), "b": (19, 1876)},
        x={"a": 19, "b": 19}),
    "Fig5_enrichment.tif": dict(
        xh=32,
        zones={"a": (15, 130, 25, 130), "b": (15, 140, 1065, 1170),
               "c": (15, 130, 2105, 2205)},
        rows={}, singles={"a": (45, 73), "b": (45, 1104), "c": (45, 2168)},
        x={"a": 45, "b": 45, "c": 45}),
    "Fig6_ml_screening.tif": dict(
        xh=32,
        zones={"a": (15, 130, 20, 120), "b": (665, 790, 20, 120),
               "c": (1340, 1480, 20, 120), "d": (15, 130, 820, 920),
               "e": (1290, 1400, 820, 920)},
        rows={"abc": 88, "de": 891}, singles={},
        x={"a": 44, "b": 708, "c": 1382, "d": 44, "e": 1331}),
    "Fig7_validation.tif": dict(
        xh=32,
        zones={"a": (240, 350, 5, 90), "b": (0, 90, 1880, 1980),
               "c": (1040, 1180, 1880, 1980)},
        rows={"bc": 1949}, singles={"a": (17, 17)},
        x={"a": 17, "b": 17, "c": 1087}),
    "Fig8_immune.tif": dict(
        xh=32,
        zones={"a": (0, 200, 0, 150), "b": (0, 280, 1150, 1850)},
        rows={}, singles={"a": (19, 22), "b": (19, 1242)},
        x={"a": 19, "b": 19}),
    "Fig9_scrna_localization.tif": dict(
        xh=32,
        zones={"a": (0, 220, 0, 170), "b": (1050, 1260, 0, 170),
               "c": (0, 160, 780, 940), "d": (950, 1200, 780, 940)},
        rows={"ab": 85, "cd": 832}, singles={},
        x={"a": 41, "b": 1121, "c": 41, "d": 1041}),
}


def find_font(names, size):
    for n in names:
        for d in FONT_DIRS:
            p = os.path.join(d, n)
            if os.path.exists(p):
                try:
                    return ImageFont.truetype(p, size)
                except OSError:
                    pass
    return ImageFont.load_default()


def ascender_ratio():
    """Measured from the font file at render time (never hardcode)."""
    f = find_font(["arialbd.ttf", "Arial_Bold.ttf", "DejaVuSans-Bold.ttf"], 100)
    def probe(ch):
        img = Image.new("L", (300, 300), 255)
        ImageDraw.Draw(img).text((100, 100), ch, font=f, fill=0)
        ys, _ = np.where(np.array(img) < 150)
        return ys.max() - ys.min() + 1
    return probe("b") / probe("a")


def render_letter(ch, xh, asc, fbold):
    want = xh if ch in "ace" else xh * asc
    size = max(10, int(xh / 0.52))
    for _ in range(10):
        img = Image.new("L", (size * 3, size * 3), 255)
        ImageDraw.Draw(img).text((size, size), ch, font=fbold(size), fill=0)
        arr = np.array(img)
        ys, xs = np.where(arr < 150)
        h = ys.max() - ys.min() + 1
        if abs(h - want) <= 1:
            break
        size = max(8, int(size * want / h))
    return arr[ys.min():ys.max() + 1, xs.min():xs.max() + 1]


def clean(gray, wx, wy, gh, gw):
    if wx < 3 or wy < 3 or wy + gh + 2 >= gray.shape[0] or wx + gw + 2 >= gray.shape[1]:
        return False
    return gray[wy - 2:wy + gh + 2, wx - 2:wx + gw + 2].min() >= 245


def remove_tag(gray, rgb, zx0, zx1, zy0, zy1):
    win = gray[zy0:zy1, zx0:zx1]
    dark = win < 150
    lab, n = ndimage.label(dark)
    cands = []
    for i in range(1, n + 1):
        comp = lab == i
        ys, xs = np.where(comp)
        h, w = ys.max() - ys.min() + 1, xs.max() - xs.min() + 1
        if not (14 <= h <= 80 and 8 <= w <= 90) or comp.sum() / (h * w) < 0.08:
            continue
        if (dark & ndimage.binary_dilation(comp, iterations=3) & ~comp).sum() > 0:
            continue
        cands.append((int(xs.min()), int(ys.min()), comp, h))
    if not cands:
        return False
    x0, y0, comp, _ = min(cands, key=lambda t: (t[0], t[1]))
    # erase the whole tag bounding box (2-px margin) plus its anti-aliased
    # fringe; the isolation gate above guarantees no other ink within 3 px
    # of the tag, and the other-ink guard keeps shared-edge pixels intact
    ys0, xs0 = np.where(comp)
    box = np.zeros(win.shape, bool)
    box[max(0, ys0.min() - 2):ys0.max() + 3, max(0, xs0.min() - 2):xs0.max() + 3] = True
    removal = box & (win < 250)
    removal &= ~(ndimage.binary_dilation(dark & ~comp, iterations=1))
    full = np.zeros(gray.shape, bool)
    full[zy0:zy1, zx0:zx1] = removal
    rgb[full] = 255
    gray[full] = 255
    return True


def process(fname, spec, asc, fbold):
    path = os.path.join(SRC, fname)
    im = Image.open(path)
    rgb_im = im.convert("RGB")
    if "scale" in spec:
        sx, sy = spec["scale"]
        rgb_im = rgb_im.resize((round(rgb_im.width * sx), round(rgb_im.height * sy)),
                               Image.LANCZOS)
        dpi = im.info.get("dpi", (300, 300))
    rgb = np.array(rgb_im)
    gray = np.array(rgb_im.convert("L")).copy()
    xh = spec["xh"]
    crops = {t: render_letter(t, xh, asc, fbold) for t in spec["x"]}

    for tag, (zx0, zx1, zy0, zy1) in spec["zones"].items():
        if not remove_tag(gray, rgb, zx0, zx1, zy0, zy1):
            print(f"  {fname}: native tag '{tag}' not found in zone - skipped")

    pos = {}
    for row, baseline in spec["rows"].items():
        for t in row:
            pos[t] = [spec["x"][t], baseline - crops[t].shape[0]]
    for t, (x, top) in spec["singles"].items():
        pos[t] = [x, top]

    for row, baseline in spec["rows"].items():
        shift = 0
        while shift <= 40 and not all(clean(gray, pos[t][0], pos[t][1] + shift,
                                            *crops[t].shape) for t in row):
            shift += 4
        if shift > 40:
            raise SystemExit(f"{fname}: no clean baseline for row {row}")
        for t in row:
            pos[t][1] += shift
    for t in spec["singles"]:
        x, y = pos[t]
        gh, gw = crops[t].shape
        if not clean(gray, x, y, gh, gw):
            for dx in list(range(-24, 25, 4)) + list(range(24, -25, -4)):
                if clean(gray, x + dx, y, gh, gw):
                    pos[t][0] = x + dx
                    break
            else:
                raise SystemExit(f"{fname}: no clean slot for '{t}'")

    for t, (wx, wy) in pos.items():
        gh, gw = crops[t].shape
        rgb[wy:wy + gh, wx:wx + gw] = np.minimum(rgb[wy:wy + gh, wx:wx + gw],
                                                 crops[t][:, :, None])
        gray[wy:wy + gh, wx:wx + gw] = np.minimum(gray[wy:wy + gh, wx:wx + gw], crops[t])
        print(f"  {fname}: {t} -> ({wx},{wy}) baseline={wy + gh}")

    for attempt in range(6):
        try:
            Image.fromarray(rgb).save(path, format="TIFF",
                                      dpi=im.info.get("dpi", (300, 300)),
                                      compression="tiff_lzw")
            break
        except OSError:
            time.sleep(0.8)
    else:
        raise SystemExit(f"{fname}: save failed (file locked?)")
    im.close()


if __name__ == "__main__":
    asc = ascender_ratio()
    print(f"measured ascender/x-height ratio: {asc:.4f}")
    for fname, spec in FIGS.items():
        print(fname)
        process(fname, spec, asc, lambda s: find_font(["arialbd.ttf", "Arial_Bold.ttf",
                                                       "DejaVuSans-Bold.ttf"], s))
    print("LETTERS_DONE")

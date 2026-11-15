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
        xh=80,
        zones={"a": (0, 440, 0, 320), "b": (1800, 2800, 0, 320), "c": (0, 440, 1924, 2300)},
        rows={"ab": 124}, singles={"c": (68, 1968)},
        x={"a": 68, "b": 2168, "c": 68}),
    "Fig3_wgcna.tif": dict(
        xh=76, scale=(5088 / 4200, 2420 / 1998),
        zones={},  # renderer draws no panel tags; letters are added only
        rows={}, singles={"a": (42, 42), "b": (2198, 42)},
        x={"a": 42, "b": 2198}),
    "Fig4_candidates.tif": dict(
        xh=70,
        zones={"a": (0, 440, 0, 520), "b": (0, 560, 3740, 4600)},
        rows={}, singles={"a": (38, 38), "b": (38, 3752)},
        x={"a": 38, "b": 38}),
    "Fig5_enrichment.tif": dict(
        xh=64,
        zones={"a": (30, 260, 50, 260), "b": (30, 280, 2130, 2340),
               "c": (30, 260, 4210, 4410)},
        rows={}, singles={"a": (90, 146), "b": (90, 2208), "c": (90, 4336)},
        x={"a": 90, "b": 90, "c": 90}),
    "Fig6_ml_screening.tif": dict(
        xh=64,
        zones={"a": (30, 260, 40, 240), "b": (1330, 1580, 40, 240),
               "c": (2680, 2960, 40, 240), "d": (30, 260, 1640, 1840),
               "e": (2580, 2800, 1640, 1840)},
        rows={"abc": 176, "de": 1782}, singles={},
        x={"a": 88, "b": 1416, "c": 2764, "d": 88, "e": 2662}),
    "Fig7_validation.tif": dict(
        xh=64,
        zones={"a": (480, 700, 10, 180), "b": (0, 180, 3760, 3960),
               "c": (2080, 2360, 3760, 3960)},
        rows={"bc": 3898}, singles={"a": (34, 34)},
        x={"a": 34, "b": 34, "c": 2174}),
    "Fig8_immune.tif": dict(
        xh=64,
        zones={"a": (0, 400, 0, 300), "b": (0, 560, 2300, 3700)},
        rows={}, singles={"a": (38, 44), "b": (38, 2484)},
        x={"a": 38, "b": 38}),
    "Fig9_scrna_localization.tif": dict(
        xh=64,
        zones={"a": (0, 440, 0, 340), "b": (2100, 2520, 0, 340),
               "c": (0, 320, 1560, 1880), "d": (1900, 2400, 1560, 1880)},
        rows={"ab": 170, "cd": 1664}, singles={},
        x={"a": 82, "b": 2242, "c": 82, "d": 2082}),
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
        if abs(h - want) <= max(1, xh // 32):
            break
        size = max(8, int(size * want / h))
    return arr[ys.min():ys.max() + 1, xs.min():xs.max() + 1]


def clean(gray, wx, wy, gh, gw):
    if wx < 5 or wy < 5 or wy + gh + 4 >= gray.shape[0] or wx + gw + 4 >= gray.shape[1]:
        return False
    return gray[wy - 4:wy + gh + 4, wx - 4:wx + gw + 4].min() >= 245


def remove_tag(gray, rgb, zx0, zx1, zy0, zy1):
    win = gray[zy0:zy1, zx0:zx1]
    dark = win < 150
    lab, n = ndimage.label(dark)
    cands = []
    for i in range(1, n + 1):
        comp = lab == i
        ys, xs = np.where(comp)
        h, w = ys.max() - ys.min() + 1, xs.max() - xs.min() + 1
        if not (28 <= h <= 160 and 16 <= w <= 180) or comp.sum() / (h * w) < 0.08:
            continue
        if (dark & ndimage.binary_dilation(comp, iterations=6) & ~comp).sum() > 0:
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
    box[max(0, ys0.min() - 4):ys0.max() + 5, max(0, xs0.min() - 4):xs0.max() + 5] = True
    removal = box & (win < 250)
    removal &= ~(ndimage.binary_dilation(dark & ~comp, iterations=2))
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
        while shift <= 80 and not all(clean(gray, pos[t][0], pos[t][1] + shift,
                                            *crops[t].shape) for t in row):
            shift += 8
        if shift > 40:
            raise SystemExit(f"{fname}: no clean baseline for row {row}")
        for t in row:
            pos[t][1] += shift
    for t in spec["singles"]:
        x, y = pos[t]
        gh, gw = crops[t].shape
        if not clean(gray, x, y, gh, gw):
            for dx in list(range(-48, 49, 8)) + list(range(48, -49, -8)):
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

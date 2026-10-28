# assemble_figures.py ----
# Assemble the publication TIFF figures (300 dpi, LZW, <=179 mm width) from
# the PNG panels rendered by 13_figures.R.
# Outputs (all in figures/): Fig1_flowchart.tif (600 dpi),
#   Fig2_pca_volcano_heatmap.tif, Fig3_wgcna.tif, Fig4_candidates.tif,
#   Fig5_enrichment.tif, Fig6_ml_screening.tif, Fig7_validation.tif,
#   Fig8_immune.tif, Fig9_scrna_localization.tif

import os
from PIL import Image, ImageDraw, ImageFont

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")

def flatten_rgb(im):
    if im.mode == "RGBA":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[3])
        return bg
    return im.convert("RGB")

def load(name):
    return flatten_rgb(Image.open(os.path.join(SRC, name)))

def save_tif(im, name, dpi=(300, 300)):
    out = os.path.join(SRC, name)
    im.save(out, format="TIFF", dpi=dpi, compression="tiff_lzw")
    w_in = im.size[0] / dpi[0]
    print(f"{name}: {im.size[0]}x{im.size[1]} px, {w_in * 25.4:.0f} mm @ {dpi[0]} dpi")

def stack(ims, target_w=None, gap=40, resize_after=False):
    # target_w = None aligns all panels to the narrowest one; an explicit
    # width (2108 px = 179 mm at 300 dpi) matches the submission figures.
    # resize_after=True stacks at native width first and rescales the whole
    # canvas afterwards (reproduces the submission Fig8 assembly).
    if resize_after:
        canvas = stack(ims, target_w=None, gap=gap)
        if target_w is None or canvas.width == target_w:
            return canvas
        return canvas.resize((target_w, round(canvas.height * target_w / canvas.width)),
                             Image.LANCZOS)
    if target_w is None:
        target_w = min(im.width for im in ims)
    scaled = [im if im.width == target_w
              else im.resize((target_w, round(im.height * target_w / im.width)), Image.LANCZOS)
              for im in ims]
    H = sum(im.height for im in scaled) + gap * (len(scaled) - 1)
    canvas = Image.new("RGB", (target_w, H), (255, 255, 255))
    y = 0
    for im in scaled:
        canvas.paste(im, (0, y))
        y += im.height + gap
    return canvas

# Fig 1: flowchart (source rendered at 600 dpi)
save_tif(load("Fig1_flowchart.png"), "Fig1_flowchart.tif", dpi=(600, 600))

# Fig 2: PCA/volcano (A|B) on top, heatmap (C) below, widths aligned to C
A = load("Fig2AB_pca_volcano_rev.png")
C = load("Fig2C_heatmap_rev.png")
save_tif(stack([A, C], gap=0), "Fig2_pca_volcano_heatmap.tif")

# Fig 3: WGCNA panels
save_tif(load("Fig3_wgcna_rev.png"), "Fig3_wgcna.tif")

# Fig 4: Venn (A) over candidate bar (B)
save_tif(stack([load("Fig4A_venn_rev.png"), load("Fig4B_bar_rev.png")],
               target_w=2108, gap=30), "Fig4_candidates.tif")

# Fig 5: enrichment panels
save_tif(load("Fig5_enrichment_rev.png"), "Fig5_enrichment.tif")

# Fig 6: ML screening grid
save_tif(load("Fig6_ml_screening_rev.png"), "Fig6_ml_screening.tif")

# Fig 7: validation ROC (A) on top, calibration (B) | DCA (C) below.
# The C panel title is drawn on the final canvas so that it matches B's
# rendered size after resizing.
A7 = load("Fig7A_validation_roc_rev.png")
B7 = load("Fig7B_calibration_GSE142025.png")
C7 = load("Fig7C_dca_GSE142025.png")
cdraw = ImageDraw.Draw(C7)
cdraw.rectangle([0, 0, C7.width, 180], fill=(255, 255, 255))
bdraw = ImageDraw.Draw(B7)
bdraw.rectangle([0, 0, 80, 80], fill=(255, 255, 255))
W = A7.width
gap = 40
half = (W - gap) // 2
B2 = B7.resize((half, round(B7.height * half / B7.width)), Image.LANCZOS)
C2 = C7.resize((half, round(C7.height * half / C7.width)), Image.LANCZOS)
H = A7.height + gap + max(B2.height, C2.height)
canvas = Image.new("RGB", (W, H), (255, 255, 255))
canvas.paste(A7, (0, 0))
yB = A7.height + gap
canvas.paste(B2, (0, yB))
canvas.paste(C2, (half + gap, yB))
draw = ImageDraw.Draw(canvas)
FONT_DIRS = [r"C:\Windows\Fonts",
             "/usr/share/fonts/truetype/dejavu",
             "/usr/share/fonts/truetype/liberation",
             "/Library/Fonts"]
def find_font(candidates, size):
    for name in candidates:
        for d in FONT_DIRS:
            p = os.path.join(d, name)
            if os.path.exists(p):
                try:
                    return ImageFont.truetype(p, size)
                except OSError:
                    pass
    return ImageFont.load_default()
c_title_font = find_font(["arialbd.ttf", "Arial_Bold.ttf", "DejaVuSans-Bold.ttf"], 28)
tag_font = find_font(["arial.ttf", "Arial.ttf", "DejaVuSans.ttf"], 36)
c_title_text = "Decision curve analysis \u2013 GSE142025"
bbox = c_title_font.getbbox(c_title_text)
tw = bbox[2] - bbox[0]
tx = half + gap + (half - tw) // 2
draw.text((tx, yB + 15), c_title_text, fill=(0, 0, 0), font=c_title_font)
draw.text((12, yB + 10), "B", fill=(0, 0, 0), font=tag_font)
draw.text((half + gap + 12, yB + 10), "C", fill=(0, 0, 0), font=tag_font)
save_tif(canvas, "Fig7_validation.tif")

# Fig 8: immune boxplots (A) over correlation heatmap (B)
save_tif(stack([load("Fig8A_immune_box_rev.png"), load("Fig8B_corr_rev.png")],
               target_w=2108, resize_after=True), "Fig8_immune.tif")

# Fig 9: snRNA localization
save_tif(load("Fig9_scrna_localization_rev.png"), "Fig9_scrna_localization.tif")

print("\nASSEMBLY_DONE")

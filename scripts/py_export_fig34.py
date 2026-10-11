# -*- coding: utf-8 -*-
# Export freshly regenerated Fig3/Fig4 PNGs to final TIFFs (RGB, 300 dpi)
from PIL import Image

for src, dst in [
    (r"results\figures\Fig3_wgcna_rev.png", r"results\figures\tiff_300_final\Fig3_wgcna.tif"),
    (r"results\figures\Fig4_enrichment_rev.png", r"results\figures\tiff_300_final\Fig4_enrichment.tif"),
]:
    im = Image.open(src)
    if im.mode == "RGBA":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[3])
        im = bg
    else:
        im = im.convert("RGB")
    im.save(dst, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
    print(dst, Image.open(dst).size, Image.open(dst).mode)
print("EXPORT_DONE")

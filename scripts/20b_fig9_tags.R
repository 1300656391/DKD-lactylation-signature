# 20b_fig9_tags.R ----
# R10-R3/R20: Fig9 add panel letters A-D + uppercase UMAP axis labels
# (redrawn from the annotated object; upstream not re-run).
suppressWarnings(suppressMessages({
  library(Seurat); library(ggplot2); library(patchwork)
}))
set.seed(42)

merged <- readRDS("results/revision/GSE131882_seurat_annotated.rds")
Idents(merged) <- "celltype"

pA <- DimPlot(merged, reduction = "umap", label = TRUE, repel = TRUE, pt.size = .3) +
  ggtitle("Cell types (GSE131882 snRNA-seq)") + NoLegend() +
  labs(x = "UMAP1", y = "UMAP2")
pB <- FeaturePlot(merged, "PDK4", reduction = "umap", pt.size = .3, order = TRUE) +
  scale_colour_gradientn(colours = c("grey92", "#FDB863", "#D73027")) +
  ggtitle("PDK4") + NoLegend() +
  labs(x = "UMAP1", y = "UMAP2")
vtheme <- theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5))
pC <- VlnPlot(merged, "PDK4", group.by = "celltype", split.by = "group", split.plot = TRUE, pt.size = 0) & vtheme
pD <- VlnPlot(merged, "PDK2", group.by = "celltype", split.by = "group", split.plot = TRUE, pt.size = 0) & vtheme

fig9 <- (pA + pB + pC + pD + plot_layout(ncol = 2, heights = c(1, 1.25)) +
  plot_annotation(tag_levels = "A")) &
  theme(plot.tag = element_text(size = 26, face = "bold"))

dir.create("results/figures/tiff_300", showWarnings = FALSE)
png("results/figures/Fig9_scrna_pdk_localization.png", type = "windows",
    width = 3900, height = 3600, res = 300)
print(fig9); dev.off()
tiff("results/figures/tiff_300/Fig9_scrna_pdk_localization.tif",
     width = 3900, height = 3600, res = 300, compression = "lzw")
print(fig9); dev.off()
cat("FIG9_DONE tags+UMAP-labels\n")

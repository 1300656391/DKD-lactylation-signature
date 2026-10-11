# 13b_fig2b_fix.R ----
# R10-R11: Fig2B volcano annotation occlusion fix
# Up/Down count annotations aligned to panel edges (hjust), x-range extended,
# repel labels constrained, to avoid occlusion/clipping.
suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
}))
fig_dir <- "results/figures"
png_win <- function(filename, width, height, units = 'in', res = 300, ...)
  grDevices::png(filename, width = width, height = height, units = units, res = res, type = 'windows', ...)
save_fig <- function(p, name, w, h, dpi = 300)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = dpi, device = png_win)

train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
allg <- fread("results/A_all_genes_limma.csv")
deg <- fread("results/A_DEGs.csv")

pca_df <- data.frame(PC1 = prcomp(t(train_expr))$x[, 1],
                     PC2 = prcomp(t(train_expr))$x[, 2], group = train_ph$group)
p1 <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 2.6, alpha = .9) + stat_ellipse(level = .9, linetype = 2, linewidth = .6) +
  scale_color_manual(values = c(Control = "#2166AC", DN = "#B2182B")) +
  theme_bw(base_size = 11) + labs(color = NULL, tag = "A") +
  ggtitle("Training cohort (GSE96804)")

vol <- data.frame(gene = allg$gene, logFC = allg$logFC, neglog10 = -log10(allg$adj.P.Val))
vol$sig <- ifelse(vol$logFC >= .5 & vol$neglog10 > -log10(.05), "Up",
           ifelse(vol$logFC <= -.5 & vol$neglog10 > -log10(.05), "Down", "NS"))
lb <- vol[vol$gene %in% c("PDK4", "PDK2", "ACSF2", "ME1") & vol$sig != "NS", ]
xr <- c(min(vol$logFC) - 1.8, max(vol$logFC) + 1.8)
p2 <- ggplot(vol, aes(logFC, neglog10, color = sig)) +
  geom_point(size = .55, alpha = .45) +
  scale_color_manual(values = c(Up = "#B2182B", Down = "#2166AC", NS = "grey78")) +
  geom_hline(yintercept = -log10(.05), linetype = 3) +
  annotate("text", x = max(vol$logFC) * .55, y = -log10(.05) + 2.5,
           label = "adj.P = 0.05", size = 2.8) +
  geom_vline(xintercept = c(-.5, .5), linetype = 3) +
  geom_text_repel(data = lb, aes(label = gene), size = 2.7, max.overlaps = 20,
                  show.legend = FALSE, xlim = c(xr[1] + .25, xr[2] - .25)) +
  coord_cartesian(xlim = xr) +
  annotate("text", x = xr[2] - .15, y = 14, hjust = 1,
           label = sprintf("Up\nn = %d", sum(vol$sig == "Up")), size = 3, color = "#B2182B") +
  annotate("text", x = xr[1] + .15, y = 14, hjust = 0,
           label = sprintf("Down\nn = %d", sum(vol$sig == "Down")), size = 3, color = "#2166AC") +
  theme_bw(base_size = 11) + labs(tag = "B") +
  ggtitle(sprintf("Volcano (DEGs n = %d)", nrow(deg))) +
  theme(legend.position = "bottom")

save_fig(p1 + p2, "Fig2AB_pca_volcano_rev", 9.6, 4.4)
tiff(file.path(fig_dir, "tiff_300", "Fig2AB_pca_volcano_rev.tif"),
     width = 2880, height = 1320, res = 300, compression = "lzw")
print(p1 + p2); dev.off()
cat("FIG2AB_DONE annotations-repositioned\n")

# 13c_fig4_regen.R ----
# R10 visual-gate fix: tiff_300 vintage mismatch for Fig4 (old 2-panel GO+KEGG render).
# Regenerate the correct 3-panel version (GO BP / Reactome / GSEA-by-NES, tags A/B/C)
# from 13_figures_rev.R code, at 300 dpi.
suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
}))
set.seed(42)
fig_dir <- "results/figures"
png_win <- function(filename, width, height, units = 'in', res = 300, ...)
  grDevices::png(filename, width = width, height = height, units = units, res = res, type = 'windows', ...)
save_fig <- function(p, name, w, h, dpi = 300)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = dpi, device = png_win)

ego <- fread("results/enrich_GO_BP.csv")
er <- fread("results/enrich_Reactome.csv")
gse <- fread("results/gsea_Reactome.csv")
mk_bar <- function(d, n, ttl, tag) {
  if (!hasName(d, "Count")) d[, Count := lengths(strsplit(geneID, "/"))]
  d <- utils::head(d[order(d$p.adjust)], n)
  d$lab <- make.names(substr(d$Description, 1, 52))
  d$lab <- factor(d$lab, levels = rev(d$lab))
  ggplot(d, aes(Count, lab, fill = p.adjust)) + geom_col(width = .7) +
    scale_fill_gradient(low = "#D32F2F", high = "#1976D2", name = "p.adj") +
    theme_bw(base_size = 9.5) + ylab(NULL) + ggtitle(ttl) + labs(tag = tag) +
    theme(plot.title = element_text(size = 10.5))
}
p1 <- mk_bar(ego, 10, "GO BP (top DEGs)", "A")
p2 <- mk_bar(er, 10, "Reactome (top DEGs)", "B")
gse2 <- copy(gse)
gse2$lab <- make.names(paste0(substr(gse2$Description, 1, 42),
                              "  NES=", sprintf("%.2f", gse2$NES)))
topg <- rbind(head(gse2[order(NES), ], 4), head(gse2[order(-NES), ], 4))
topg$lab <- factor(topg$lab, levels = topg$lab[order(topg$NES)])
p3 <- ggplot(topg, aes(NES, lab, fill = p.adjust)) + geom_col(width = .7) +
  scale_fill_gradient(low = "#D32F2F", high = "#1976D2", name = "FDR") +
  theme_bw(base_size = 9) + ylab(NULL) + ggtitle("GSEA (Reactome, top by |NES|)") +
  labs(tag = "C")
save_fig(p1 / p2 / p3 + plot_layout(heights = c(1, 1, .75)), "Fig4_enrichment_rev", 9, 10)
cat("FIG4_DONE 3-panel-ABC-300dpi\n")

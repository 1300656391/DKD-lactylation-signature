# 16d_fig8b_regen.R ----
# R11-F1/S4: re-render Fig8B from the REPAIRED immune_gene_corr.csv (gene symbols restored)
# with row/column names bumped to 9.5pt. Same layout/tag as 16_fig_round3.R.
suppressWarnings(suppressMessages({
  library(ComplexHeatmap); library(circlize); library(grid)
}))
fig_dir <- "results/figures"

immune_corr <- read.csv("results/immune_gene_corr.csv", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(is.character(immune_corr$gene), !immune_corr$gene[[1]] %in% c("1", "1.0"))
cm <- as.matrix(immune_corr[, setdiff(colnames(immune_corr), c("gene", "Normal_mucosa"))])
rownames(cm) <- immune_corr$gene
cat("rows:", nrow(cm), "cols:", ncol(cm), "| first rows:", head(immune_corr$gene, 3), "\n")
col_fun2 <- colorRamp2(c(-.6, 0, .6), c("#2166AC", "white", "#B2182B"))
png(file.path(fig_dir, "Fig8B_corr_rev.png"), type = "windows", width = 2400, height = 1600, res = 250)
draw(Heatmap(cm, name = "Spearman rho", col = col_fun2,
             row_names_gp = gpar(fontsize = 9.5), column_names_gp = gpar(fontsize = 9.5),
             column_title = "Candidate LRGs x immune cells (GSE96804, all samples; 18 cells)"))
grid.text("B", x = unit(2.6, "mm"), y = unit(1, "npc") - unit(2.6, "mm"),
          gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()
cat("FIG8B_DONE symbols+9.5pt\n")

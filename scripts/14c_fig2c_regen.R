# 14c_fig2c_regen.R ----
# R11-S4: re-render Fig2C with row names bumped 6 -> 8.5pt (same filter logic as 16_fig_round3.R).
suppressWarnings(suppressMessages({
  library(data.table); library(ComplexHeatmap); library(circlize); library(grid)
}))
train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
deg <- fread("results/A_DEGs.csv")
pseudogenes_observed <- c("PFN1P2","SLC7A5P2","SMG1P1","CHCHD2P6","CHCHD2P8","NBPF8P")
nc_filter <- function(g) {
  grepl("^(MIR|SNOR|SNORA|SNORD|SCAR|RNU|LINC|RNA|DANCR|VTRNA|NPIPA|NBPF|LOC[0-9]|OTTHUM|ENSG)", g) |
    grepl("-AS[0-9]$", g) | (g %in% pseudogenes_observed)
}
top_pool <- deg[order(deg$adj.P.Val)][1:120]
top50 <- top_pool[!nc_filter(top_pool$gene)][1:50]
mat <- t(scale(t(train_expr[top50$gene, ])))
ha <- HeatmapAnnotation(Group = train_ph$group,
                        col = list(Group = c(Control = "#2166AC", DN = "#B2182B")))
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
png("results/figures/Fig2C_heatmap_rev.png", type = "windows", width = 1800, height = 2400, res = 250)
draw(Heatmap(mat, name = "z-score", col = col_fun, top_annotation = ha,
             show_column_names = FALSE, row_names_gp = gpar(fontsize = 8.5),
             column_title = "Top 50 protein-coding DEGs (GSE96804)"))
grid.text("C", x = unit(2.6, "mm"), y = unit(1, "npc") - unit(2.6, "mm"),
          gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()
cat("FIG2C_DONE 8.5pt rows:", nrow(top50), "\n")

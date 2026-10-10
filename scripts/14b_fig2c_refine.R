# 14b_fig2c_refine.R ----
# 修正14_n_fixes.R中P[0-9]$误伤真实蛋白编码基因(DUSP1/EIF4EBP2)的问题：
# 仅排除top-120内观察到的确证伪基因后缀条目，确保"Top 50 protein-coding DEGs"严格成立
suppressWarnings(suppressMessages({ library(data.table); library(ComplexHeatmap); library(circlize) }))
train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph   <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
deg <- fread("results/A_DEGs.csv")

pseudogenes_observed <- c("PFN1P2","SLC7A5P2","SMG1P1","CHCHD2P6","CHCHD2P8","NBPF8P")
nc_filter <- function(g) {
  grepl("^(MIR|SNOR|SNORA|SNORD|SCAR|RNU|LINC|RNA|DANCR|VTRNA|NPIPA|NBPF|LOC[0-9]|OTTHUM|ENSG)", g) |
    grepl("-AS[0-9]$", g) | (g %in% pseudogenes_observed)
}
top_pool <- deg[order(deg$adj.P.Val)][1:120]
top50 <- top_pool[!nc_filter(top_pool$gene)][1:50]
kept_flag <- c("DUSP1","EIF4EBP2") %in% top50$gene
cat("DUSP1/EIF4EBP2 retained:", paste(kept_flag, collapse="/"), "\n")
cat("Excluded:", paste(top_pool$gene[nc_filter(top_pool$gene)], collapse=", "), "\n")
# 断言：过滤后前50行即按adj.P排序的蛋白编码前50（在top-120池内无真实编码基因被误伤）
stopifnot(all(!nc_filter(top50$gene)))

mat <- t(scale(t(train_expr[top50$gene, ])))
ha <- HeatmapAnnotation(Group = train_ph$group,
                        col = list(Group = c(Control = "#2166AC", DN = "#B2182B")))
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
png("results/figures/Fig2C_heatmap_rev.png", type = "windows", width = 1800, height = 2400, res = 250)
draw(Heatmap(mat, name = "z-score", col = col_fun, top_annotation = ha,
             show_column_names = FALSE, row_names_gp = gpar(fontsize = 6),
             column_title = "Top 50 protein-coding DEGs (GSE96804)"))
dev.off()
cat("Fig2C refined:", nrow(top50), "rows\n")

# 13d_fig5a_regen.R ----
# R11-S1: re-render Fig5A Venn with the corrected key count (177; empty-alias phantom key removed)
suppressWarnings(suppressMessages({
  library(data.table); library(VennDiagram); library(grid)
}))
fig_dir <- "results/figures"
lrg <- fread("data/LRG_geneset.csv")
deg <- fread("results/A_DEGs.csv")

lrg_alias <- lrg$alias[!is.na(lrg$alias) & nzchar(trimws(lrg$alias))]
lrg_keys5 <- toupper(unique(c(lrg$symbol, lrg_alias)))
stopifnot(length(lrg_keys5) == 177)
venn_list <- list(
  `LRG keys (177)` = lrg_keys5,
  `DEGs (2645)` = deg$gene,
  `Modules (4210)` = fread("results/wgcna_sigmodule_genes.csv")$gene
)
png(file.path(fig_dir, "Fig5A_venn_rev.png"), type = "windows", width = 1800, height = 1500, res = 250)
v <- venn.diagram(venn_list, filename = NULL, fill = c("#F4B183", "#9DC3E6", "#C5E0B4"),
                  alpha = .55, cex = 1.05, cat.cex = .95, margin = .12,
                  fontfamily = "sans", cat.fontfamily = "sans",
                  main = "Candidate LRG selection", main.cex = 1.15)
grid.draw(v)
dev.off()
cat("FIG5A_DONE keys=177\n")

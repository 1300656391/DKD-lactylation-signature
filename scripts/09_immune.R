# 09_immune.R ----
# S9: ssGSEA免疫浸润（Bindea细胞类型signature，源自IOBR signature_collection）
#     + 组间比较 + 候选基因与细胞比例相关
# 输出: results/immune_ssgsea.csv + results/immune_tests.csv + results/immune_gene_corr.csv
# 运行: Rscript 09_immune.R

suppressWarnings(suppressMessages({
  library(GSVA)
  library(data.table)
}))

expr <- readRDS("data/GSE96804_expr_symbol.rds")
ph   <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
stopifnot(identical(ph$sample, colnames(expr)))

# 加载signature并筛选Bindea细胞类型
e <- new.env(); load("data/signature_collection.rda", envir = e)
sc <- get("signature_collection", envir = e)
cell_names <- grep("_Bindea_et_al$", names(sc), value = TRUE)
cell_names <- setdiff(cell_names, c("SW480_cancer_cells_Bindea_et_al"))  # 去肿瘤细胞系
gsets <- sc[cell_names]
names(gsets) <- sub("_Bindea_et_al$", "", cell_names)
gsets <- lapply(gsets, function(g) unique(toupper(intersect(g, rownames(expr)))))
gsets <- gsets[vapply(gsets, length, integer(1)) >= 5]
cat("免疫细胞类型数:", length(gsets), "\n")

# ssGSEA（兼容GSVA新旧API）
ss <- tryCatch({
  gsva(as.matrix(expr), gsets, method = "ssgsea",
       kcdf = "Gaussian", verbose = FALSE, ssgsea.norm = TRUE)
}, error = function(e) {
  par <- ssgseaParam(as.matrix(expr), gsets, normalize = TRUE)
  gsva(par, verbose = FALSE)
})
ss <- as.data.frame(t(ss))
ss$sample <- rownames(ss)
fwrite(ss, "results/immune_ssgsea.csv")

# 组间Wilcoxon
grp <- ph$group
tests <- do.call(rbind, lapply(names(gsets), function(cn) {
  p <- wilcox.test(ss[[cn]] ~ factor(grp))$p.value
  data.frame(cell = cn,
             median_DN = median(ss[[cn]][grp == "DN"], na.rm = TRUE),
             median_Ctrl = median(ss[[cn]][grp == "Control"], na.rm = TRUE),
             p = p)
}))
tests$p.adj <- p.adjust(tests$p, method = "BH")
tests <- tests[order(tests$p), ]
fwrite(tests, "results/immune_tests.csv")
cat("显著差异细胞类型(p.adj<0.05):", sum(tests$p.adj < 0.05), "\n")
print(head(tests, 10))

# 候选基因与细胞比例Spearman相关
cand <- fread("results/A_candidates.csv")
cand_expr <- expr[intersect(cand$gene, rownames(expr)), , drop = FALSE]
cor_m <- do.call(rbind, lapply(rownames(cand_expr), function(g) {
  vapply(names(gsets), function(cn)
    suppressWarnings(cor(t(cand_expr[g, ]), ss[[cn]], method = "spearman")),
    numeric(1))
}))
cor_m <- as.data.frame(cor_m)
cor_m$gene <- rownames(cor_m)
fwrite(cor_m, "results/immune_gene_corr.csv")
cat("相关矩阵完成:", nrow(cor_m), "基因 x", length(gsets), "细胞\n")

cat("\nIMMUNE_DONE\n")

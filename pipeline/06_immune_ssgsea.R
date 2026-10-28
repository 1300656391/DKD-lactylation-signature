# 06_immune_ssgsea.R ----
# ssGSEA immune infiltration on the training cohort (Bindea cell-type
# signatures from the IOBR signature collection), group comparisons, and
# candidate-gene x cell-fraction Spearman correlations.
# Outputs: results/immune_ssgsea.csv, results/immune_tests.csv,
#          results/immune_tests_18cells.csv (BH within the 18-cell family),
#          results/immune_gene_corr.csv

suppressWarnings(suppressMessages({
  library(GSVA); library(data.table)
}))
set.seed(42)

expr <- readRDS("data/GSE96804_expr_symbol.rds")
ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
stopifnot(identical(ph$sample, colnames(expr)))

# Load Bindea signatures; drop the SW480 cancer-cell line entry
e <- new.env(); load("data/signature_collection.rda", envir = e)
sc <- get("signature_collection", envir = e)
cell_names <- grep("_Bindea_et_al$", names(sc), value = TRUE)
cell_names <- setdiff(cell_names, "SW480_cancer_cells_Bindea_et_al")
gsets <- sc[cell_names]
names(gsets) <- sub("_Bindea_et_al$", "", cell_names)
gsets <- lapply(gsets, function(g) unique(toupper(intersect(g, rownames(expr)))))
gsets <- gsets[vapply(gsets, length, integer(1)) >= 5]
cat("cell types:", length(gsets), "\n")

# ssGSEA (works with both the old and the new GSVA API)
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

# Group comparison (Wilcoxon DN vs Control, BH across all cell types)
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
cat("significant cell types (all-family BH<0.05):", sum(tests$p.adj < 0.05), "\n")

# 18-cell immune family (tissue-reference entry Normal_mucosa excluded),
# BH recomputed within the 18-test family
it18 <- tests[tests$cell != "Normal_mucosa", ]
it18$p.adj_BH18 <- p.adjust(it18$p, method = "BH")
it18$sig_BH18 <- ifelse(it18$p.adj_BH18 < 0.05, "yes", "no")
fwrite(it18, "results/immune_tests_18cells.csv")
cat("significant cell types (BH18<0.05):", sum(it18$sig_BH18 == "yes"), "\n")

# Candidate gene x cell-fraction Spearman correlations
cand <- fread("results/A_candidates.csv")
cand_expr <- expr[intersect(cand$gene, rownames(expr)), , drop = FALSE]
gene_order <- rownames(cand_expr)
cor_m <- do.call(rbind, lapply(gene_order, function(g) {
  vapply(names(gsets), function(cn)
    suppressWarnings(cor(t(cand_expr[g, ]), ss[[cn]], method = "spearman")),
    numeric(1))
}))
rownames(cor_m) <- gene_order
cor_m <- as.data.frame(cor_m)
cor_m$gene <- rownames(cor_m)
fwrite(cor_m, "results/immune_gene_corr.csv")
stopifnot(is.character(cor_m$gene), all(cor_m$gene %in% rownames(cand_expr)))
cat("correlation matrix:", nrow(cor_m), "genes x", length(gsets), "cells\n")

cat("\nIMMUNE_SSGSEA_DONE\n")

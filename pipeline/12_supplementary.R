# 12_supplementary.R ----
# Supplementary analyses:
#   S20 - candidate-gene x immune-cell Spearman correlations with BH control
#         across the 23 genes x 18 immune cell types (414 tests);
#   S18 - head-to-head discrimination of the PDK2/PDK4 signature versus the
#         PTGS2/NFE2L2 comparator pair under frozen coefficients, across the
#         training context, the eight validation contexts, and GSE273763 under
#         two matrix normalizations.
# Outputs:
#   results/supplementary_summary.csv (pipeline marker),
#   results/supplement/TableS20_fig8b_correlations_BH.csv,
#   results/supplement/TableS18_head_to_head.csv

suppressWarnings(suppressMessages({
  library(data.table); library(org.Hs.eg.db)
}))
set.seed(42)
dir.create("results/supplement", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. S20: gene x cell Spearman + BH (family = 23 genes x 18 cells)
# ============================================================
ssg <- fread("results/immune_ssgsea.csv")
expr <- readRDS("data/GSE96804_expr_symbol.rds")
rho <- fread("results/immune_gene_corr.csv")
genes <- rho$gene
cells <- setdiff(names(rho), c("gene", "Normal_mucosa"))
cat("genes:", length(genes), "| cells (Normal_mucosa excluded):", length(cells), "\n")

P <- Q <- R <- matrix(NA_real_, length(genes), length(cells),
                      dimnames = list(genes, cells))
for (g in genes) for (cc in cells) {
  ct <- suppressWarnings(cor.test(as.numeric(expr[g, ]), ssg[[cc]], method = "spearman"))
  P[g, cc] <- ct$p.value
  R[g, cc] <- ct$estimate
}
Q[] <- p.adjust(as.vector(P), "BH")
cat("BH family size:", length(P), "| BH-significant:", sum(Q < 0.05), "\n")

s20 <- rbindlist(lapply(genes, function(g) data.frame(
  gene = g, cell_type = cells,
  rho = round(R[g, ], 3), p = signif(P[g, ], 3), p_BH = signif(Q[g, ], 3))))
setorder(s20, p_BH)
fwrite(s20, "results/supplement/TableS20_fig8b_correlations_BH.csv")
cat("S20 written:", nrow(s20), "rows\n")

# ============================================================
# 2. S18: head-to-head vs PTGS2/NFE2L2 (frozen coefficients)
# ============================================================
fm <- fread("results/frozen/final_model_coefficients.csv")
bh <- setNames(fm$canonical_deviance, fm$term)
BETA <- c(int = bh[["(Intercept)"]], PDK2 = bh[["PDK2"]], PDK4 = bh[["PDK4"]])
GENES4 <- c("PDK2", "PDK4", "PTGS2", "NFE2L2")

zrow <- function(v) (v - mean(v)) / sd(v)
zmat <- function(Mx) t(apply(Mx, 1, zrow))
# Mann-Whitney AUC in the as-trained direction (higher score = case).
# P values follow the scipy convention: exact test when the smaller group
# has at most 8 samples, otherwise asymptotic with continuity correction.
auc_raw <- function(score, yy) {
  n1 <- sum(yy == 1); n0 <- sum(yy == 0)
  wt <- suppressWarnings(wilcox.test(score[yy == 1], score[yy == 0],
                                     exact = min(n1, n0) <= 8, correct = TRUE))
  c(auc = unname(wt$statistic) / (n1 * n0), p = wt$p.value)
}

# Training context: gene orientation + frozen PTGS2/NFE2L2 logistic
train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
Ztr <- zmat(train_expr[GENES4, , drop = FALSE])
ytr <- as.integer(train_ph$group == "DN")
orient <- vapply(GENES4, function(g) {
  a <- auc_raw(Ztr[g, ], ytr)[["auc"]]
  ifelse(a > 0.5, 1, -1)
}, numeric(1))
gpair <- glm(ytr ~ Ztr["PTGS2", ] + Ztr["NFE2L2", ], family = binomial)
b_pair <- coef(gpair)
cat("orientation:", paste(names(orient), orient, sep = "=", collapse = ", "), "\n")
cat("PTGS2/NFE2L2 frozen logistic: intercept =", round(b_pair[1], 3),
    "PTGS2 =", round(b_pair[2], 3), "NFE2L2 =", round(b_pair[3], 3), "\n")

# Context builder: list(name, M = 4-gene raw matrix, y, orient)
# Default orientation is training-derived; GSE273763 uses the uniform
# training-direction orientation (all four genes downregulated in DN) that
# was pre-specified for that cohort's extension.
ctx <- list()
add_ctx <- function(name, emat, yy, ornt = NULL) {
  hit <- intersect(GENES4, rownames(emat))
  ctx[[name]] <<- list(M = emat[hit, , drop = FALSE], y = yy, orient = ornt)
}
add_ctx("GSE96804_train", train_expr, ytr)
for (acc in c("GSE30528", "GSE30529")) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  p <- read.csv(paste0("data/", acc, "_pheno.csv"), stringsAsFactors = FALSE)
  keep <- p$group %in% c("DKD", "Control")
  add_ctx(paste0(acc, ifelse(acc == "GSE30528", "_glom", "_tub")),
          e[, keep, drop = FALSE], as.integer(p$group[keep] == "DKD"))
}
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])))
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
rownames(e14) <- toupper(rownames(e14))
k_all <- p14$group %in% c("Advanced_DN", "Early_DN", "Control")
add_ctx("GSE142025_allDN", e14[, k_all, drop = FALSE], as.integer(p14$group[k_all] != "Control"))
for (st in c("Early_DN", "Advanced_DN")) {
  k <- p14$group %in% c(st, "Control")
  add_ctx(paste0("GSE142025_", ifelse(st == "Early_DN", "early", "advanced")),
          e14[, k, drop = FALSE], as.integer(p14$group[k] != "Control"))
}
plat <- c(GSE104948_pooled = "GSE104948", `GSE104948_GPL22945` = "GSE104948-GPL22945",
          `GSE104948_GPL24120` = "GSE104948-GPL24120")
for (nm in names(plat)) {
  e <- readRDS(paste0("data/", plat[[nm]], "_expr_symbol.rds"))
  p <- read.csv(paste0("data/", plat[[nm]], "_pheno.csv"), stringsAsFactors = FALSE)
  add_ctx(nm, e[, p$group %in% c("DN", "Control"), drop = FALSE],
          as.integer(p$group[p$group %in% c("DN", "Control")] == "DN"))
}

# GSE273763 under two normalizations (same samples, not independent)
g273_dir <- "data/GSE273763"
M273 <- read.table(gzfile(file.path(g273_dir, "GSE273763_DKD_GRCH38_forGEO_July2026_v2.txt.gz")),
                   sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
colnames(M273) <- gsub("\\s+", "", colnames(M273))
M273 <- as.matrix(M273)
y273 <- as.integer(grepl("DKD", colnames(M273)))
orient_ext <- c(PDK2 = -1, PDK4 = -1, PTGS2 = -1, NFE2L2 = -1)
add_ctx("GSE273763_submitter", M273, y273, orient_ext)
cnt <- read.table(gzfile(file.path(g273_dir, "GSE273763_countsmatrix.txt.gz")),
                  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
colnames(cnt) <- gsub("\\s+", "", colnames(cnt))
sym <- unname(mapIds(org.Hs.eg.db, keys = rownames(cnt), column = "SYMBOL",
                     keytype = "ENSEMBL"))
keep <- !is.na(sym) & sym != ""
cnt <- as.matrix(cnt[keep, ]); sym <- sym[keep]
if (anyDuplicated(sym)) {
  libsum <- rowSums(cnt)
  kp <- tapply(seq_along(sym), sym, function(i) i[which.max(libsum[i])])
  cnt <- cnt[as.integer(kp), ]; sym <- sym[as.integer(kp)]
}
rownames(cnt) <- sym
cpm <- log2(sweep(cnt, 2, colSums(cnt), "/") * 1e6 + 1)
add_ctx("GSE273763_log2cpm", cpm, y273, orient_ext)

TABLE2 <- c(GSE30528_glom = 0.906, GSE30529_tub = 0.650, GSE142025_allDN = 0.992,
            GSE142025_early = 0.889, GSE142025_advanced = 1.000,
            GSE104948_pooled = 1.000, `GSE104948_GPL22945` = 1.000,
            `GSE104948_GPL24120` = 0.600, GSE273763_submitter = 0.805)

rows <- list()
for (nm in names(ctx)) {
  Mx <- ctx[[nm]]$M; yy <- ctx[[nm]]$y
  ov <- if (is.null(ctx[[nm]]$orient)) orient else ctx[[nm]]$orient
  Z <- zmat(Mx)
  feats <- list()
  for (g in rownames(Z)) feats[[g]] <- ov[g] * Z[g, ]
  if (all(c("PTGS2", "NFE2L2") %in% rownames(Z))) {
    feats[["PTGS2+NFE2L2_meanz"]] <- (ov["PTGS2"] * Z["PTGS2", ] +
                                      ov["NFE2L2"] * Z["NFE2L2", ]) / 2
    feats[["PTGS2+NFE2L2_frozen"]] <- b_pair[1] + b_pair[2] * Z["PTGS2", ] +
                                      b_pair[3] * Z["NFE2L2", ]
  }
  if (all(c("PDK2", "PDK4") %in% rownames(Z))) {
    feats[["PDK2/PDK4_signature"]] <- BETA[["int"]] + BETA[["PDK2"]] * Z["PDK2", ] +
                                      BETA[["PDK4"]] * Z["PDK4", ]
  }
  for (fname in names(feats)) {
    a <- auc_raw(as.numeric(feats[[fname]]), yy)
    rows[[length(rows) + 1]] <- data.frame(
      context = nm, n_dn = sum(yy == 1), n_ctrl = sum(yy == 0), feature = fname,
      AUC_as_trained = round(a[["auc"]], 3),
      AUC_oriented = round(max(a[["auc"]], 1 - a[["auc"]]), 3),
      direction_reversed = ifelse(a[["auc"]] < 0.5, "True", "False"),
      p = round(a[["p"]], 6),
      Table2 = if (fname == "PDK2/PDK4_signature" && nm %in% names(TABLE2))
                 TABLE2[[nm]] else NA_real_)
  }
}
s18 <- rbindlist(rows, fill = TRUE)
fwrite(s18, "results/supplement/TableS18_head_to_head.csv")

val <- s18[context != "GSE96804_train" & feature == "PDK2/PDK4_signature"]
sig_mean <- mean(val$AUC_oriented); sig_min <- min(val$AUC_oriented)
mz <- s18[context != "GSE96804_train" & feature == "PTGS2+NFE2L2_meanz"]
cat(sprintf("10-context signature: mean = %.3f, min = %.3f | mean-z pair: min = %.3f\n",
            sig_mean, sig_min, min(mz$AUC_oriented)))
chk <- s18[feature == "PDK2/PDK4_signature" & !is.na(Table2)]
cat("Table 2 cross-check (oriented vs Table 2):",
    all(chk$AUC_oriented == chk$Table2), "\n")

# ============================================================
# 3. Summary marker
# ============================================================
fwrite(data.frame(
  table = c("TableS20_fig8b_correlations_BH", "TableS18_head_to_head"),
  rows = c(nrow(s20), nrow(s18)),
  key_stat = c(sum(Q < 0.05), sprintf("signature mean=%.3f min=%.3f (10 contexts)",
                                      sig_mean, sig_min))),
  "results/supplementary_summary.csv")

cat("\nSUPPLEMENTARY_DONE\n")

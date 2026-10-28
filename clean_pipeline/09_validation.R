# 09_validation.R ----
# Validation of the frozen two-gene signature across cohorts:
#   - AUC (auto-oriented, as pROC/Table 2) with DeLong CI, Youden sensitivity
#     and specificity with exact binomial CIs, for the training cohort, the
#     four microarray/RNA-seq cohorts, and GSE273763 (FFPE RNA-seq);
#   - strict threshold transfer at the training-derived cutoff (S6);
#   - calibration (frozen vs cohort-recalibrated) and DCA in GSE142025 and
#     GSE273763; normalization sensitivities for GSE273763 (log2CPM,
#     DESeq2 median-of-ratios);
#   - GSE294519 tubular DSP scoring context (S15c).
# Outputs:
#   results/validation_summary.csv,
#   results/validation/threshold_transfer.csv, auc_youden_summary.csv,
#     calibration_GSE142025.csv, calibration_GSE273763_dual.csv,
#     dca_GSE142025.csv, dca_GSE273763.csv, gse273763_normalization_sensitivity.csv,
#     gse294519_scoring.csv,
#   results/supplement/TableS5_missing_genes_by_cohort.csv,
#     TableS6_training_cutoff_transfer.csv,
#     TableS15c_gse294519_dsp.csv, TableS16_gse273763_calibration_dual.csv

suppressWarnings(suppressMessages({
  library(data.table); library(pROC); library(rms); library(org.Hs.eg.db)
}))
set.seed(42)
dir.create("results/validation", showWarnings = FALSE, recursive = TRUE)
dir.create("results/supplement", showWarnings = FALSE, recursive = TRUE)

fm <- fread("results/frozen/final_model_coefficients.csv")
bh <- setNames(fm$canonical_deviance, fm$term)
hub <- readLines("results/hub_genes.txt")
BETA <- c(int = bh[["(Intercept)"]], PDK2 = bh[["PDK2"]], PDK4 = bh[["PDK4"]])

zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
zrow <- function(v) (v - mean(v)) / sd(v)
sc <- function(e) {
  hg <- intersect(hub, rownames(e))
  Xv <- t(zrows(e[hg, , drop = FALSE]))
  as.vector(Xv %*% bh[hg] + bh[["(Intercept)"]])
}
bci <- function(k, n) binom.test(k, n)$conf.int

# ============================================================
# 1. Cohort contexts
# ============================================================
train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
y_tr <- as.integer(train_ph$group == "DN")

ctx <- list()
ctx$GSE96804_train <- list(y = y_tr, s = sc(train_expr))
for (acc in c("GSE30528", "GSE30529")) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  p <- read.csv(paste0("data/", acc, "_pheno.csv"), stringsAsFactors = FALSE)
  keep <- p$group %in% c("DKD", "Control")
  ctx[[acc]] <- list(y = as.integer(p$group[keep] == "DKD"), s = sc(e[, keep, drop = FALSE]))
}
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])))
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
rownames(e14) <- toupper(rownames(e14))
k3 <- p14$group %in% c("Advanced_DN", "Early_DN", "Control")
k4 <- p14$group %in% c("Early_DN", "Control")
ctx$GSE142025_all   <- list(y = as.integer(p14$group[k3] != "Control"), s = sc(e14[, k3, drop = FALSE]))
ctx$GSE142025_early <- list(y = as.integer(p14$group[k4] == "Early_DN"), s = sc(e14[, k4, drop = FALSE]))
e10 <- readRDS("data/GSE104948_expr_symbol.rds")
p10 <- read.csv("data/GSE104948_pheno.csv", stringsAsFactors = FALSE)
k5 <- p10$group %in% c("DN", "Control")
ctx$GSE104948 <- list(y = as.integer(p10$group[k5] == "DN"), s = sc(e10[, k5, drop = FALSE]))

# GSE273763 (FFPE RNA-seq): submitter-normalized matrix is the primary scheme;
# counts-based log2CPM and DESeq2 median-of-ratios are sensitivities
g273_dir <- "data/GSE273763"
M <- read.table(gzfile(file.path(g273_dir, "GSE273763_DKD_GRCH38_forGEO_July2026_v2.txt.gz")),
                sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
colnames(M) <- gsub("\\s+", "", colnames(M))
M <- as.matrix(M)
storage.mode(M) <- "numeric"
stopifnot(ncol(M) == 63)
y273 <- as.integer(grepl("DKD", colnames(M)))
stopifnot(sum(y273 == 1) == 50, sum(y273 == 0) == 13)
s273 <- BETA[["int"]] + BETA[["PDK2"]] * zrow(as.numeric(M["PDK2", ])) +
        BETA[["PDK4"]] * zrow(as.numeric(M["PDK4", ]))
ctx$GSE273763 <- list(y = y273, s = s273)

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
norm_cpm <- function(m) log2(t(t(m) / colSums(m)) * 1e6 + 1)
norm_mor <- function(m) {
  # DESeq2 estimateSizeFactors median-of-ratios: pseudo-reference per gene =
  # exp(mean of log counts) over genes positive in all samples; per-sample
  # size factor = median count/reference ratio.
  g <- m[rowSums(m > 0) == ncol(m), , drop = FALSE]
  ref <- exp(rowMeans(log(g)))
  size <- apply(sweep(g, 1, ref, "/"), 2, median)
  log2(sweep(m, 2, size, "/") + 1)
}
s273_cpm <- BETA[["int"]] + BETA[["PDK2"]] * zrow(norm_cpm(cnt)["PDK2", ]) +
            BETA[["PDK4"]] * zrow(norm_cpm(cnt)["PDK4", ])
s273_mor <- BETA[["int"]] + BETA[["PDK2"]] * zrow(norm_mor(cnt)["PDK2", ]) +
            BETA[["PDK4"]] * zrow(norm_mor(cnt)["PDK4", ])

# ============================================================
# 2. AUC + Youden table (Table 2 / Table 3 basis)
# ============================================================
auc_row <- function(key, label) {
  r <- roc(ctx[[key]]$y, ctx[[key]]$s, quiet = TRUE)
  ci <- ci.auc(r)
  cut <- coords(r, x = "best", ret = "threshold", best.method = "youden")
  cutv <- as.numeric(cut[1])
  yy <- ctx[[key]]$y; ss <- ctx[[key]]$s
  pred <- as.integer(ss >= cutv)
  tp <- sum(pred == 1 & yy == 1); fn <- sum(pred == 0 & yy == 1)
  fp <- sum(pred == 1 & yy == 0); tn <- sum(pred == 0 & yy == 0)
  data.frame(cohort = label, n = length(yy), n_case = sum(yy), n_ctrl = sum(!as.logical(yy)),
             AUC = round(as.numeric(auc(r)), 3),
             CI = sprintf("%.3f-%.3f", ci[1], ci[3]),
             sens = round(tp / (tp + fn), 3),
             sens_lo = round(bci(tp, tp + fn)[1], 3),
             sens_hi = round(bci(tp, tp + fn)[2], 3),
             spec = round(tn / (tn + fp), 3),
             spec_lo = round(bci(tn, tn + fp)[1], 3),
             spec_hi = round(bci(tn, tn + fp)[2], 3))
}
tab <- rbind(
  auc_row("GSE96804_train", "GSE96804_train"),
  auc_row("GSE30528", "GSE30528_glom"),
  auc_row("GSE30529", "GSE30529_tub"),
  auc_row("GSE142025_all", "GSE142025_allDN"),
  auc_row("GSE142025_early", "GSE142025_early"),
  auc_row("GSE104948", "GSE104948_pooled"),
  auc_row("GSE273763", "GSE273763"))
fwrite(tab, "results/validation/auc_youden_summary.csv")
print(tab)

# ============================================================
# 3. Threshold transfer at the training-derived cutoff (S6)
# ============================================================
roc_tr <- roc(y_tr, ctx$GSE96804_train$s, quiet = TRUE)
cut0 <- as.numeric(coords(roc_tr, x = "best", ret = "threshold", best.method = "youden"))
cat("training Youden cutoff =", round(cut0, 4), "\n")
tt_row <- function(key, label) {
  yy <- ctx[[key]]$y; ss <- ctx[[key]]$s
  pred <- as.integer(ss >= cut0)
  tp <- sum(pred == 1 & yy == 1); fn <- sum(pred == 0 & yy == 1)
  fp <- sum(pred == 1 & yy == 0); tn <- sum(pred == 0 & yy == 0)
  data.frame(cohort = label, n_case = sum(yy), n_ctrl = sum(!as.logical(yy)),
             sensitivity = tp / (tp + fn),
             specificity = tn / (tn + fp),
             PPV = tp / (tp + fp), NPV = tn / (tn + fn))
}
s6 <- rbind(tt_row("GSE96804_train", "GSE96804_train"),
            tt_row("GSE30528", "GSE30528"), tt_row("GSE30529", "GSE30529"),
            tt_row("GSE142025_all", "GSE142025_all"), tt_row("GSE142025_early", "GSE142025_early"),
            tt_row("GSE104948", "GSE104948"), tt_row("GSE273763", "GSE273763"))
fwrite(s6, "results/validation/threshold_transfer.csv")
file.copy("results/validation/threshold_transfer.csv",
          "results/supplement/TableS6_training_cutoff_transfer.csv", overwrite = TRUE)
print(s6)

# ============================================================
# 4. Calibration and DCA
# ============================================================
calib_stats <- function(p, yy) {
  vp <- val.prob(p, yy, pl = FALSE)
  c(intercept = unname(vp["Intercept"]), slope = unname(vp["Slope"]),
    MAE = unname(vp["Eavg"]), Emax = unname(vp["Emax"]),
    Spiegelhalter_z = unname(vp["S:z"]), Brier = mean((p - yy)^2))
}
nb_curve <- function(p, yy) {
  do.call(rbind, lapply(seq(0.05, 0.95, by = 0.05), function(th) {
    nb <- function(pp) mean((pp >= th) * yy - (pp >= th) * (1 - yy) * th / (1 - th))
    data.frame(threshold = th, model = round(nb(p), 4),
               treat_all = round(mean(yy) - mean(1 - yy) * th / (1 - th), 4),
               treat_none = 0)
  }))
}

# GSE142025 (all DN vs control): recalibrated score
y3 <- ctx$GSE142025_all$y
g1 <- glm(y3 ~ ctx$GSE142025_all$s, family = binomial)
p_reca <- as.numeric(predict(g1, type = "response"))
p_frozen_142 <- plogis(ctx$GSE142025_all$s)
cal142 <- rbind(
  data.frame(probability = "frozen", t(round(calib_stats(p_frozen_142, y3), 3))),
  data.frame(probability = "recalibrated", t(round(calib_stats(p_reca, y3), 3))))
fwrite(cal142, "results/validation/calibration_GSE142025.csv")
dca142 <- nb_curve(p_reca, y3)
fwrite(dca142, "results/validation/dca_GSE142025.csv")

# GSE273763: dual report (frozen vs recalibrated)
p_frozen_273 <- plogis(s273)
reca273 <- glm(y273 ~ s273, family = binomial)
p_reca273 <- as.numeric(predict(reca273, type = "response"))
cal273 <- rbind(
  data.frame(probability = "frozen", t(round(calib_stats(p_frozen_273, y273), 3))),
  data.frame(probability = "recalibrated", t(round(calib_stats(p_reca273, y273), 3))))
fwrite(cal273, "results/validation/calibration_GSE273763_dual.csv")
file.copy("results/validation/calibration_GSE273763_dual.csv",
          "results/supplement/TableS16_gse273763_calibration_dual.csv", overwrite = TRUE)
dca273 <- nb_curve(p_reca273, y273)
setnames(dca273, "model", "recalibrated")
dca273$frozen <- nb_curve(p_frozen_273, y273)$model
setcolorder(dca273, c("threshold", "recalibrated", "frozen", "treat_all", "treat_none"))
fwrite(dca273, "results/validation/dca_GSE273763.csv")

# ============================================================
# 5. GSE273763 normalization sensitivities + single-gene rows
# ============================================================
as_trained_auc <- function(score, yy) as.numeric(auc(roc(yy, score, quiet = TRUE, direction = "<")))
# Single genes are oriented by their training-cohort direction (all four
# genes are downregulated in DN in training)
GENES4 <- c("PDK2", "PDK4", "PTGS2", "NFE2L2")
Ztr <- t(apply(train_expr[GENES4, , drop = FALSE], 1, zrow))
orient <- vapply(GENES4, function(g) {
  a <- as_trained_auc(Ztr[g, ], y_tr)
  ifelse(a > 0.5, 1, -1)
}, numeric(1))
feat_auc <- function(Mx, tag) {
  zn <- function(g) zrow(as.numeric(Mx[g, ]))
  rows <- list()
  for (g in GENES4) {
    if (!g %in% rownames(Mx)) next
    a <- as_trained_auc(orient[[g]] * zn(g), y273)
    rows[[length(rows) + 1]] <- data.frame(norm = tag, feature = g,
                                           AUC_as_trained = round(a, 3),
                                           AUC_oriented = round(max(a, 1 - a), 3))
  }
  sigv <- BETA[["int"]] + BETA[["PDK2"]] * zn("PDK2") + BETA[["PDK4"]] * zn("PDK4")
  a <- as_trained_auc(sigv, y273)
  rows[[length(rows) + 1]] <- data.frame(norm = tag, feature = "PDK2/PDK4_signature",
                                         AUC_as_trained = round(a, 3),
                                         AUC_oriented = round(max(a, 1 - a), 3))
  do.call(rbind, rows)
}
norm_sens <- rbind(feat_auc(M, "submitter"),
                   feat_auc(norm_cpm(cnt), "log2cpm"),
                   feat_auc(norm_mor(cnt), "medianRatio"))
fwrite(norm_sens, "results/validation/gse273763_normalization_sensitivity.csv")
print(norm_sens)

# ============================================================
# 6. GSE294519 tubular DSP context (S15c)
# ============================================================
# The PTGS2/NFE2L2 comparator is a frozen logistic fit on the training
# z-scores (single-gene orientation as in section 5).
b_pair <- coef(glm(y_tr ~ Ztr["PTGS2", ] + Ztr["NFE2L2", ], family = binomial))

d294 <- read.table(gzfile("data/GSE294519/GSE294519_Normalization_data_allsamples.txt.gz"),
                   sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
d294 <- as.matrix(d294)
rn <- toupper(trimws(as.character(rownames(d294))))
if (anyDuplicated(rn)) {
  mu <- rowMeans(d294, na.rm = TRUE)
  keep <- tapply(seq_along(rn), rn, function(i) i[which.max(mu[i])])
  d294 <- d294[as.integer(keep), , drop = FALSE]; rn <- rn[as.integer(keep)]
}
rownames(d294) <- rn
y294 <- as.integer(startsWith(colnames(d294), "DKD"))
z294 <- function(g) if (g %in% rownames(d294)) zrow(as.numeric(d294[g, ])) else NULL
feats294 <- list()
for (g in GENES4) {
  v <- z294(g)
  if (!is.null(v)) feats294[[g]] <- orient[[g]] * v
}
if (!is.null(z294("PTGS2")) && !is.null(z294("NFE2L2"))) {
  feats294[["PTGS2+NFE2L2_meanz"]] <- (feats294[["PTGS2"]] + feats294[["NFE2L2"]]) / 2
  feats294[["PTGS2+NFE2L2_frozen"]] <- b_pair[1] + b_pair[2] * z294("PTGS2") +
                                       b_pair[3] * z294("NFE2L2")
}
if (!is.null(z294("PDK2")) && !is.null(z294("PDK4")))
  feats294[["PDK2/PDK4_signature"]] <- BETA[["int"]] + BETA[["PDK2"]] * z294("PDK2") +
                                       BETA[["PDK4"]] * z294("PDK4")
s294 <- do.call(rbind, lapply(names(feats294), function(fname) {
  v <- as.numeric(feats294[[fname]])
  a <- as_trained_auc(v, y294)
  # asymptotic rank-sum test with continuity correction (scipy convention,
  # as in the published table)
  wt <- suppressWarnings(wilcox.test(v[y294 == 1], v[y294 == 0],
                                     exact = FALSE, correct = TRUE))
  data.frame(cohort = "GSE294519", n_dn = sum(y294 == 1), n_ctrl = sum(y294 == 0),
             feature = fname, AUC_as_trained = round(a, 3),
             AUC_oriented = round(max(a, 1 - a), 3), p = round(wt$p.value, 3))
}))
fwrite(s294, "results/validation/gse294519_scoring.csv")
file.copy("results/validation/gse294519_scoring.csv",
          "results/supplement/TableS15c_gse294519_dsp.csv", overwrite = TRUE)
i <- match("PDK2/PDK4_signature", s294$feature)
cat("GSE294519 signature: as-trained AUC =", s294$AUC_as_trained[i],
    " p =", s294$p[i], "\n")

# ============================================================
# 7. Validation summary + missing-genes table (S5)
# ============================================================
# S5: hub-gene availability per cohort ("None" when all five are present)
s5_rows <- list(
  list("GSE96804_train", length(y_tr), train_expr),
  list("GSE30528", length(ctx$GSE30528$y), readRDS("data/GSE30528_expr_symbol.rds")),
  list("GSE30529", length(ctx$GSE30529$y), readRDS("data/GSE30529_expr_symbol.rds")),
  list("GSE142025_DN_all", sum(k3), e14),
  list("GSE142025_Early", sum(k4), e14),
  list("GSE142025_Advanced", sum(p14$group %in% c("Advanced_DN", "Control")), e14),
  list("GSE104948", length(ctx$GSE104948$y), e10),
  list("GSE273763", length(y273), M))
s5 <- do.call(rbind, lapply(s5_rows, function(r) {
  miss <- setdiff(hub, rownames(r[[3]]))
  data.frame(cohort = r[[1]], n = r[[2]],
             missing_genes = if (length(miss) == 0) "None" else paste(miss, collapse = ","))
}))
fwrite(s5, "results/supplement/TableS5_missing_genes_by_cohort.csv")
print(s5)

fwrite(tab[, c("cohort", "n", "AUC", "CI", "sens", "spec")],
       "results/validation_summary.csv")

cat("\nVALIDATION_DONE\n")

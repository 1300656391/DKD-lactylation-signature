# 12_revision_analyses.R ----
# 按review_round1.md执行的补充分析
# M3: hub基因完整limma统计 | M7: GSE104948分平台AUC | m6/M4: 固定cutoff敏感度/特异度
# M6: GSE142025校准曲线+DCA（非完美分离队列）| M9: 免疫剔除Normal_mucosa
# 输出: results/revision/*.csv + figures/Fig7B/C/D新版
suppressWarnings(suppressMessages({
  library(data.table); library(limma); library(glmnet); library(pROC)
  library(rms); library(dcurves); library(ggplot2)
}))
set.seed(42)
dir.create("results/revision", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# ---------- 复现两基因模型（与10/11脚本同种子同流程） ----------
train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
hub <- readLines("results/hub_genes.txt")
zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
Xh <- data.frame(t(zrows(train_expr[hub, , drop = FALSE])), check.names = FALSE)
y <- as.integer(train_ph$group == "DN")
cv2 <- cv.glmnet(as.matrix(Xh), y, family = "binomial", alpha = 1, nfolds = 10)
bh <- as.matrix(coef(cv2, s = "lambda.1se"))
score_tr <- as.vector(as.matrix(Xh) %*% bh[-1, , drop = FALSE] + bh["(Intercept)", 1])
roc_tr <- roc(y, score_tr, quiet = TRUE)
cat("惩罚模型系数:\n"); print(bh)

# ---------- M3: hub基因完整统计 ----------
allg <- fread("results/A_all_genes_limma.csv")
hub_tab <- allg[gene %in% hub, .(gene, logFC, adj.P.Val, AveExpr)]
hub_tab <- hub_tab[match(hub, gene)]
fwrite(hub_tab, "results/revision/M3_hub_full_stats.csv")
cat("\nM3 hub统计:\n"); print(hub_tab)

# ---------- 通用评分函数 ----------
cohort_data <- function(e, p, case_groups, ctrl_groups) {
  keep <- p$group %in% c(case_groups, ctrl_groups)
  list(e = e[, keep, drop = FALSE], y = as.integer(p$group[keep] %in% case_groups))
}
score_cohort <- function(e, yy) {
  hg <- intersect(hub, rownames(e))
  Xv <- data.frame(t(zrows(e[hg, , drop = FALSE])), check.names = FALSE)
  s <- rep(bh["(Intercept)", 1], ncol(e))
  for (g in hg) s <- s + bh[g, 1] * Xv[[g]]
  list(score = s, y = yy, missing = setdiff(hub, hg))
}

# ---------- m6/M4: 固定cutoff敏感度/特异度 ----------
coord <- coords(roc_tr, x = "best", ret = c("threshold", "sensitivity", "specificity"),
                best.method = "youden")
cut0 <- as.numeric(coord["threshold"]); cat("\n训练集Youden cutoff:", cut0, "\n")
perf_at_cut <- function(score, yy, label) {
  pred <- as.integer(score >= cut0)
  tp <- sum(pred == 1 & yy == 1); fp <- sum(pred == 1 & yy == 0)
  fn <- sum(pred == 0 & yy == 1); tn <- sum(pred == 0 & yy == 0)
  data.frame(cohort = label, n_case = tp + fn, n_ctrl = fp + tn,
             sensitivity = tp / (tp + fn), specificity = tn / (tn + fp),
             PPV = tp / (tp + fp), NPV = tn / (tn + fn))
}
perf_list <- list(perf_at_cut(score_tr, y, "GSE96804_train"))

# 各队列数据
vsets <- list()
for (acc in c("GSE30528", "GSE30529")) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  p <- read.csv(paste0("data/", acc, "_pheno.csv"), stringsAsFactors = FALSE)
  vsets[[acc]] <- cohort_data(e, p, "DKD", "Control")
}
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])),
                  stringsAsFactors = FALSE)
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
rownames(e14) <- toupper(rownames(e14))
e10 <- readRDS("data/GSE104948_expr_symbol.rds")
p10 <- read.csv("data/GSE104948_pheno.csv", stringsAsFactors = FALSE)

sc1 <- score_cohort(vsets$GSE30528$e, vsets$GSE30528$y)
sc2 <- score_cohort(vsets$GSE30529$e, vsets$GSE30529$y)
cd3 <- cohort_data(e14, p14, c("Advanced_DN", "Early_DN"), "Control"); sc3 <- score_cohort(cd3$e, cd3$y)
cd4 <- cohort_data(e14, p14, "Early_DN", "Control"); sc4 <- score_cohort(cd4$e, cd4$y)
cd5 <- cohort_data(e10, p10, "DN", "Control"); sc5 <- score_cohort(cd5$e, cd5$y)

perf_list <- c(perf_list,
  list(perf_at_cut(sc1$score, sc1$y, "GSE30528"),
       perf_at_cut(sc2$score, sc2$y, "GSE30529"),
       perf_at_cut(sc3$score, sc3$y, "GSE142025_all"),
       perf_at_cut(sc4$score, sc4$y, "GSE142025_early"),
       perf_at_cut(sc5$score, sc5$y, "GSE104948")))
perf_tab <- rbindlist(perf_list)
fwrite(perf_tab, "results/revision/M6_fixed_cutoff_performance.csv")
cat("\n固定cutoff性能:\n"); print(perf_tab)

# ---------- M7: GSE104948分平台AUC ----------
clean_v <- function(x) trimws(gsub('^"|"$', "", x))
split_line <- function(l) {
  if (grepl("=", l, fixed = TRUE)) s <- sub("^[^=]*=", "", l)
  else s <- sub("^[^\t]*\t", "", l)
  clean_v(strsplit(s, "\t", fixed = TRUE)[[1]])
}
plat_auc <- lapply(c("GSE104948-GPL22945", "GSE104948-GPL24120"), function(acc) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  ml <- readLines(paste0("data/", acc, "_meta.txt"), warn = FALSE)
  ch <- ml[grepl("^!Sample_characteristics_ch1", ml)]
  vals <- lapply(ch, split_line)
  n <- ncol(e); maxlen <- max(vapply(vals, length, integer(1)))
  mat <- do.call(rbind, lapply(vals, function(v) { length(v) <- maxlen; v }))
  diagv <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    h <- mat[, i][grepl("diagnosis:", mat[, i])]
    if (length(h)) diagv[i] <- sub("^diagnosis:\\s*", "", h[1])
  }
  keep <- diagv %in% c("Diabetic Nephropathy", "Tumor Nephrectomy")
  e2 <- e[, keep, drop = FALSE]
  yy <- as.integer(diagv[keep] == "Diabetic Nephropathy")
  sc <- score_cohort(e2, yy)
  r <- roc(yy, sc$score, quiet = TRUE)
  data.frame(platform = acc, n = sum(keep), n_case = sum(yy), n_ctrl = sum(!yy),
             AUC = as.numeric(auc(r)), lower = ci.auc(r)[1], upper = ci.auc(r)[3])
})
plat_tab <- rbindlist(plat_auc)
fwrite(plat_tab, "results/revision/M7_platform_stratified_AUC.csv")
cat("\nM7分平台AUC:\n"); print(plat_tab)

# ---------- M6: GSE142025校准+DCA ----------
png("results/figures/Fig7B_calibration_GSE142025.png", type = "windows",
    width = 1400, height = 1400, res = 240)
val.prob(sc3$score, sc3$y, rule = "youden", riskdist = FALSE,
         main = "Calibration — GSE142025 (all DN vs control)", cex.lab = 1.1)
dev.off()

dca_df <- data.frame(y = sc3$y, Signature = fitted(glm(sc3$y ~ sc3$score, family = binomial)))
png("results/figures/Fig7C_dca_GSE142025.png", type = "windows",
    width = 1500, height = 1300, res = 240)
tryCatch(print(plot(dca(y ~ Signature, data = dca_df), smooth = FALSE)),
         error = function(e) cat("DCA失败:", conditionMessage(e), "\n"))
dev.off()

# ---------- M9: 免疫剔除Normal_mucosa ----------
it <- fread("results/immune_tests.csv")
it2 <- it[cell != "Normal_mucosa"]
it2$sig <- it2$p.adj < 0.05
fwrite(it2, "results/revision/M9_immune_tests_no_mucosa.csv")
cat("\nM9剔除Normal_mucosa后: 细胞类型", length(it2), " 显著", sum(it2$sig), "\n")
cat("显著细胞:\n"); print(it2[sig == TRUE, .(cell, median_DN, median_Ctrl, p.adj)])

cat("\nREVISION_ANALYSES_DONE\n")

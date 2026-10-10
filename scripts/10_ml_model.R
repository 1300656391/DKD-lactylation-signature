# 10_ml_model.R ----
# S10: LASSO + SVM-RFE + XGBoost 三算法特征筛选 → hub基因
# S11: logistic诊断模型 + nomogram + 校准 + DCA + 多队列ROC验证
# S12: SHAP可解释性
# 输出: results/hub_genes.txt + results/AUC_summary.csv + results/AUC_single_gene.csv
#       + results/figures/Fig6_ml_draft.png + Fig7_model_draft.png + Fig8_shap_draft.png
# 运行: Rscript 10_ml_model.R

suppressWarnings(suppressMessages({
  library(data.table); library(glmnet); library(e1071); library(xgboost)
  library(pROC); library(rms); library(dcurves); library(shapviz); library(ggplot2)
}))
set.seed(42)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# ============ 数据准备 ============
train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph   <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
cand <- fread("results/A_candidates.csv")$gene
genes <- intersect(cand, rownames(train_expr))
cat("候选基因进入ML:", length(genes), ":", paste(genes, collapse = ", "), "\n")

zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
Xdf <- data.frame(t(zrows(train_expr[genes, , drop = FALSE])), check.names = FALSE)
y <- as.integer(train_ph$group == "DN")

# ============ S10-1 LASSO ============
cvfit <- cv.glmnet(as.matrix(Xdf), y, family = "binomial", alpha = 1,
                   nfolds = 10, type.measure = "deviance")
cf <- as.matrix(coef(cvfit, s = "lambda.1se"))
lasso_genes <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")
cat("\nLASSO入选(", length(lasso_genes), "):", paste(lasso_genes, collapse = ", "), "\n")

# ============ S10-2 SVM-RFE ============
svm_rfe_rank <- function(X, y) {
  feats <- colnames(X); ranked <- character(0)
  while (length(feats) > 1) {
    m <- svm(X[, feats, drop = FALSE], y, kernel = "linear", scale = FALSE)
    w <- (t(m$coefs) %*% m$SV)^2
    wr <- as.vector(w); names(wr) <- feats
    rem <- names(which.min(wr))
    ranked <- c(rem, ranked)
    feats <- setdiff(feats, rem)
  }
  c(ranked, feats)  # 从最不重要到最重要
}
yf <- factor(y, levels = c(0, 1), labels = c("Ctrl", "DN"))
rk <- svm_rfe_rank(as.matrix(Xdf), yf)
sizes <- 2:min(12, ncol(Xdf))
folds <- sample(rep(1:10, length.out = nrow(Xdf)))
acc <- vapply(sizes, function(k) {
  a <- vapply(1:10, function(f) {
    tr <- folds != f
    rk2 <- svm_rfe_rank(as.matrix(Xdf[tr, , drop = FALSE]), yf[tr])
    sel <- tail(rk2, k)
    m <- svm(as.matrix(Xdf[tr, sel, drop = FALSE]), yf[tr], kernel = "linear")
    mean(predict(m, as.matrix(Xdf[!tr, sel, drop = FALSE])) == yf[!tr])
  }, numeric(1))
  mean(a)
}, numeric(1))
best_k <- sizes[which.max(acc)]
svm_genes <- tail(rk, best_k)
cat("SVM-RFE入选(", best_k, "):", paste(rev(svm_genes), collapse = ", "), "\n")

# ============ S10-3 XGBoost ============
dtrain <- xgb.DMatrix(as.matrix(Xdf), label = y)
# 分层5折（避免单一类别折导致AUC=NaN）
mk_folds <- function(y, k = 5) {
  f0 <- sample(rep(1:k, length.out = sum(y == 0)))
  f1 <- sample(rep(1:k, length.out = sum(y == 1)))
  lapply(1:k, function(i) c(which(y == 0)[f0 == i], which(y == 1)[f1 == i]))
}
folds5 <- mk_folds(y, 5)
xpar <- list(objective = "binary:logistic", eval_metric = "auc",
             max_depth = 3, eta = 0.1)
cvb <- xgb.cv(data = dtrain, params = xpar, nrounds = 300, folds = folds5,
              early_stopping_rounds = 20, verbose = 0)
bst <- xgboost(x = as.matrix(Xdf), y = factor(y, levels = c(0, 1)),
               max_depth = 3, learning_rate = 0.1, objective = "binary:logistic",
               nrounds = max(cvb$best_iteration, 10), nthread = 2)
imp <- xgb.importance(model = bst)
xgb_genes <- head(imp$Feature, 10)
cat("XGBoost top:", paste(xgb_genes, collapse = ", "), "\n")

# ============ 交集 → hub ============
hub <- Reduce(intersect, list(lasso_genes, svm_genes, xgb_genes))
decision <- "三算法交集"
if (length(hub) < 2) {
  cnt <- table(c(lasso_genes, svm_genes, xgb_genes))
  hub <- names(cnt)[cnt >= 2]
  decision <- "至少两算法入选"
}
cat("\n决策规则:", decision, "\nHUB基因(", length(hub), "):", paste(hub, collapse = ", "), "\n")
writeLines(hub, "results/hub_genes.txt")

# ============ S11 模型（惩罚logistic，避免完美分离） ============
hub <- intersect(hub, colnames(Xdf))
Xh <- as.matrix(Xdf[, hub, drop = FALSE])
cv2 <- cv.glmnet(Xh, y, family = "binomial", alpha = 1, nfolds = 10,
                 type.measure = "auc")
bh <- as.matrix(coef(cv2, s = "lambda.1se"))       # 行名: (Intercept), hub...
keep_g <- rownames(bh)[bh[, 1] != 0]
cat("惩罚模型保留基因:", paste(setdiff(keep_g, "(Intercept)"), collapse = ", "), "\n")
score_tr <- as.vector(Xh %*% bh[-1, , drop = FALSE] + bh["(Intercept)", 1])
roc_tr <- roc(y, score_tr, quiet = TRUE)
# 诚实10折CV-AUC（分层重复glmnet）
set.seed(42)
folds_cv <- sample(rep(1:10, length.out = length(y)))
cv_auc <- vapply(1:10, function(f) {
  tr <- folds_cv != f
  cvg <- cv.glmnet(Xh[tr, , drop = FALSE], y[tr], family = "binomial",
                   alpha = 1, nfolds = 10)
  bb <- as.matrix(coef(cvg, s = "lambda.1se"))
  s_te <- as.vector(Xh[!tr, , drop = FALSE] %*% bb[-1, , drop = FALSE] + bb[1, 1])
  as.numeric(auc(roc(y[!tr], s_te, quiet = TRUE)))
}, numeric(1))
cat("训练集AUC(表观):", round(auc(roc_tr), 3),
    " 10折CV-AUC(均值):", round(mean(cv_auc), 3),
    "±", round(sd(cv_auc), 3), "\n")

# 单基因ROC（训练集）
sg_auc <- t(sapply(hub, function(g) {
  r <- roc(y, Xdf[[g]], quiet = TRUE)
  c(AUC = as.numeric(auc(r)),
    lower = ci.auc(r)[1], upper = ci.auc(r)[3])
}))
sg_auc <- as.data.frame(sg_auc); sg_auc$gene <- rownames(sg_auc)

# ============ 验证函数（惩罚模型系数） ============
validate_cohort <- function(expr, pheno, case_groups, ctrl_groups, label) {
  keep <- pheno$group %in% c(case_groups, ctrl_groups)
  e <- expr[, keep, drop = FALSE]
  ph <- pheno[keep, ]
  hg <- intersect(hub, rownames(e))
  Xv <- data.frame(t(zrows(e[hg, , drop = FALSE])), check.names = FALSE)
  score <- rep(bh["(Intercept)", 1], ncol(e))
  for (g in hg) score <- score + bh[g, 1] * Xv[[g]]
  yy <- as.integer(ph$group %in% case_groups)
  if (length(unique(yy)) < 2) return(NULL)
  r <- roc(yy, score, quiet = TRUE)
  cat(sprintf("%-28s n=%d  AUC=%.3f [%.3f-%.3f]\n", label, length(yy),
              auc(r), ci.auc(r)[1], ci.auc(r)[3]))
  data.frame(cohort = label, n = length(yy), AUC = as.numeric(auc(r)),
             lower = ci.auc(r)[1], upper = ci.auc(r)[3],
             missing_genes = paste(setdiff(hub, hg), collapse = ","))
}

auc_tab <- data.frame(cohort = "GSE96804_train", n = length(y),
                      AUC = as.numeric(auc(roc_tr)),
                      lower = ci.auc(roc_tr)[1], upper = ci.auc(roc_tr)[3],
                      missing_genes = "")

# GSE30528 / GSE30529
for (acc in c("GSE30528", "GSE30529")) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  p <- read.csv(paste0("data/", acc, "_pheno.csv"), stringsAsFactors = FALSE)
  auc_tab <- rbind(auc_tab, validate_cohort(e, p, "DKD", "Control", acc))
}

# GSE142025（重建pheno：GSM↔title对齐）
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
ttl <- cl(ml[grepl("^!Sample_title", ml)][1])
gsm <- cl(ml[grepl("^!Sample_geo_accession", ml)][1])
grp <- cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])
grp <- sub("^group:\\s*", "", grp)
p14 <- data.frame(sample = ttl, group = grp, gsm = gsm, stringsAsFactors = FALSE)
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
stopifnot(!any(is.na(match(colnames(e14), p14$gsm))))
rownames(e14) <- toupper(rownames(e14))
auc_tab <- rbind(auc_tab,
  validate_cohort(e14, p14, c("Advanced_DN", "Early_DN"), "Control", "GSE142025_DN_all"),
  validate_cohort(e14, p14, "Early_DN", "Control", "GSE142025_Early"),
  validate_cohort(e14, p14, "Advanced_DN", "Control", "GSE142025_Advanced"))

# GSE104948（小样本，仅报告）
e10 <- readRDS("data/GSE104948_expr_symbol.rds")
p10 <- read.csv("data/GSE104948_pheno.csv", stringsAsFactors = FALSE)
auc_tab <- rbind(auc_tab, validate_cohort(e10, p10, "DN", "Control", "GSE104948"))

fwrite(auc_tab, "results/AUC_summary.csv")
fwrite(sg_auc, "results/AUC_single_gene.csv")

# ============ nomogram + 校准 + DCA ============
dfm <- data.frame(y = y, Xdf[, hub, drop = FALSE])
lrm_fit <- tryCatch({
  dd <- datadist(dfm); options(datadist = "dd")
  lrm(y ~ ., data = dfm)
}, error = function(e) { cat("lrm拟合失败(可能分离):", conditionMessage(e), "\n"); NULL })
if (!is.null(lrm_fit)) {
  nom <- nomogram(lrm_fit, fun = plogis, fun.at = c(.1, .25, .5, .75, .9),
                  lp = FALSE)
  png("results/figures/Fig7_nom_draft.png", width = 900, height = 600, res = 120)
  plot(nom); dev.off()
  cal <- tryCatch(calibrate(lrm_fit, method = "boot", B = 1000),
                  error = function(e) NULL)
  if (!is.null(cal)) {
    png("results/figures/Fig7_calib_draft.png", width = 700, height = 700, res = 120)
    plot(cal, xlab = "Predicted probability", ylab = "Observed", subtitles = FALSE)
    dev.off()
  }
} else {
  # 回退：以score为唯一预测子
  lrm2 <- lrm(y ~ score_tr)
  png("results/figures/Fig7_nom_draft.png", width = 900, height = 500, res = 120)
  plot(nomogram(lrm2, fun = plogis, fun.at = c(.1, .3, .5, .7, .9))); dev.off()
  cal2 <- tryCatch(calibrate(lrm2, method = "boot", B = 1000), error = function(e) NULL)
  if (!is.null(cal2)) {
    png("results/figures/Fig7_calib_draft.png", width = 700, height = 700, res = 120)
    plot(cal2, xlab = "Predicted probability", ylab = "Observed", subtitles = FALSE)
    dev.off()
  }
}

# DCA（基于score的预测概率）
dfd <- data.frame(y = y,
                  Signature = fitted(glm(y ~ score_tr, family = binomial)))
png("results/figures/Fig7_dca_draft.png", width = 800, height = 650, res = 120)
tryCatch(print(plot(dca(y ~ Signature, data = dfd), smooth = FALSE)),
         error = function(e) cat("DCA绘制失败:", conditionMessage(e), "\n"))
dev.off()

# ============ S12 SHAP ============
shp <- tryCatch({
  m <- as.matrix(predict(bst, newdata = as.matrix(Xdf), predcontrib = TRUE))
  colnames(m) <- c(colnames(Xdf), "bias")
  m
}, error = function(e) NULL)
if (!is.null(shp) && ncol(shp) > 1) {
  p_sh <- sv_importance(shapviz(shp[, -ncol(shp), drop = FALSE]))
  ggsave("results/figures/Fig8_shap_draft.png", p_sh, width = 8, height = 6, dpi = 200)
} else {
  cat("SHAP输出维度异常，跳过（不影响主线）\n")
}

cat("\nML_MODEL_DONE\n")

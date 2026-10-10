# 15_freeze_ml.R ----
# 【法代说明 2026-10-09】本脚本为首版冻结尝试，存在两处已知缺陷：
#   1) M4 Bootstrap 段 svm_acc_b 有括号语法错误（run_freeze.txt 可见 parse error），
#      修复版为 15c_bootstrap_final.R（run_bootstrap2.txt）；
#   2) 第 ~117-119 行 oks 记账向量多一个恒 FALSE 元素（死代码），致旧版 verify_report 末行假阴性；
#      权威验证报告由 17_verify_final.R 重写（6/6 MATCH）。
# final_model_coefficients.csv 的实际列名（canonical_deviance/auc_criterion）由 15b_freeze_patch.R 写出，
# 与本脚本第 60 行附近的写出代码（estimate_auc/estimate_deviance）不一致，以磁盘文件为准。
# 实际执行顺序：15（冻结件首次生成，崩溃于M4）→ 15b（补 deviance 对照+M1，崩溃于M4）→
#   15c（系数验证+M4 Bootstrap 完成）→ 15d（选择频率口径修正）→ 17（终版验证+sessionInfo）。
# 本文件保留用于 frozen/lasso_curve.csv、svm_rfe_curve.csv、xgb_importance.csv 的生成溯源。
# ----
# Round-3 M1/M3/M4 修复：
#   M3: 以与10_ml_model.R完全相同的调用序列(set.seed(42)起)复现全部ML结果并冻结为
#       results/frozen/*.csv——绘图/稿件此后只读冻结件，单一事实来源；
#       并逐位复算5个队列AUC与AUC_summary.csv比对
#   M1: 冻结模型分数在GSE142025的logistic重定标截距/斜率 + 校准统计复算
#   M4: 200次Bootstrap全流程(候选23固定→LASSO/SVM-RFE/XGBoost→惩罚模型)乐观偏差校正
suppressWarnings(suppressMessages({
  library(data.table); library(glmnet); library(e1071); library(xgboost); library(pROC)
}))
set.seed(42)
dir.create("results/frozen", showWarnings = FALSE, recursive = TRUE)
dir.create("results/revision", showWarnings = FALSE, recursive = TRUE)

train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph   <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
cand <- fread("results/A_candidates.csv")$gene
genes <- intersect(cand, rownames(train_expr))
zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
Xdf <- data.frame(t(zrows(train_expr[genes, , drop = FALSE])), check.names = FALSE)
y <- as.integer(train_ph$group == "DN")
X <- as.matrix(Xdf)

# ========== 1) LASSO（与10_ml_model.R首call一致：deviance准则） ==========
cvfit <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10, type.measure = "deviance")
cf <- as.matrix(coef(cvfit, s = "lambda.1se"))
lasso_genes <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")
fwrite(data.frame(lambda = cvfit$lambda, cvm = cvfit$cvm, cvlo = cvfit$cvlo, cvup = cvfit$cvup,
                  is_1se = cvfit$lambda == cvfit$lambda.1se), "results/frozen/lasso_curve.csv")

# ========== 2) SVM-RFE（同序：rk→folds→acc；重组为每折一次排序，数值与原逐k重算恒等） ==========
svm_rfe_rank <- function(Xm, yy) {
  feats <- colnames(Xm); ranked <- character(0)
  while (length(feats) > 1) {
    m <- svm(Xm[, feats, drop = FALSE], yy, kernel = "linear", scale = FALSE)
    w <- (t(m$coefs) %*% m$SV)^2
    wr <- as.vector(w); names(wr) <- feats
    rem <- names(which.min(wr)); ranked <- c(rem, ranked); feats <- setdiff(feats, rem)
  }
  c(ranked, feats)
}
yf <- factor(y, levels = c(0, 1), labels = c("Ctrl", "DN"))
rk <- svm_rfe_rank(X, yf)
sizes <- 2:min(12, ncol(X))
folds <- sample(rep(1:10, length.out = nrow(X)))
fold_rank <- lapply(1:10, function(f) { tr <- folds != f; svm_rfe_rank(X[tr, , drop = FALSE], yf[tr]) })
acc <- vapply(sizes, function(k) {
  a <- vapply(1:10, function(f) {
    tr <- folds != f
    sel <- tail(fold_rank[[f]], k)
    m <- svm(X[tr, sel, drop = FALSE], yf[tr], kernel = "linear")
    mean(predict(m, X[!tr, sel, drop = FALSE]) == yf[!tr])
  }, numeric(1)); mean(a)
}, numeric(1))
best_k <- sizes[which.max(acc)]
svm_genes <- tail(rk, best_k)
fwrite(data.frame(k = sizes, accuracy = acc, is_best = sizes == best_k), "results/frozen/svm_rfe_curve.csv")

# ========== 3) XGBoost（同序：mk_folds两次sample → xgb.cv早停 → 固定轮数） ==========
dtrain <- xgb.DMatrix(X, label = y)
mk_folds <- function(yy, k = 5) {
  f0 <- sample(rep(1:k, length.out = sum(yy == 0))); f1 <- sample(rep(1:k, length.out = sum(yy == 1)))
  lapply(1:k, function(i) c(which(yy == 0)[f0 == i], which(yy == 1)[f1 == i]))
}
folds5 <- mk_folds(y, 5)
xpar <- list(objective = "binary:logistic", eval_metric = "auc", max_depth = 3, eta = 0.1)
cvb <- xgb.cv(data = dtrain, params = xpar, nrounds = 300, folds = folds5,
              early_stopping_rounds = 20, verbose = 0)
bst <- xgboost(x = X, y = factor(y, levels = c(0, 1)),
               max_depth = 3, learning_rate = 0.1, objective = "binary:logistic",
               nrounds = max(cvb$best_iteration, 10), nthread = 2)
imp <- xgb.importance(model = bst); imp$Gain_pct <- imp$Gain * 100
fwrite(imp, "results/frozen/xgb_importance.csv")
xgb_genes <- head(imp$Feature, 10)

# ========== 4) hub + 冻结最终模型（auc准则=Table2来源；另算deviance对照） ==========
hub <- Reduce(intersect, list(lasso_genes, svm_genes, xgb_genes))
if (length(hub) < 2) { cn <- table(c(lasso_genes, svm_genes, xgb_genes)); hub <- names(cn)[cn >= 2] }
hub_ref <- readLines("results/hub_genes.txt")
stopifnot(identical(sort(hub), sort(hub_ref)))

Xh <- X[, hub, drop = FALSE]
cv2  <- cv.glmnet(Xh, y, family = "binomial", alpha = 1, nfolds = 10, type.measure = "auc")
cv2d <- cv.glmnet(Xh, y, family = "binomial", alpha = 1, nfolds = 10)  # deviance对照
bh  <- as.matrix(coef(cv2,  s = "lambda.1se"))
bhd <- as.matrix(coef(cv2d, s = "lambda.1se"))
fwrite(data.frame(term = rownames(bh), estimate_auc = bh[, 1], estimate_deviance = bhd[, 1]),
       "results/frozen/final_model_coefficients.csv")
cat("冻结模型系数(auc):", paste(sprintf("%s=%.4f", rownames(bh), bh[, 1]), collapse=" "), "\n")
cat("两准则系数最大绝对差:", max(abs(bh - bhd)), "\n")

score_tr <- as.vector(Xh %*% bh[-1, , drop = FALSE] + bh["(Intercept)", 1])
roc_tr <- roc(y, score_tr, quiet = TRUE)
set.seed(42)
folds_cv <- sample(rep(1:10, length.out = length(y)))
cv_auc <- vapply(1:10, function(f) {
  tr <- folds_cv != f
  cvg <- cv.glmnet(Xh[tr, , drop = FALSE], y[tr], family = "binomial", alpha = 1, nfolds = 10)
  bb <- as.matrix(coef(cvg, s = "lambda.1se"))
  s_te <- as.vector(Xh[!tr, , drop = FALSE] %*% bb[-1, , drop = FALSE] + bb[1, 1])
  as.numeric(auc(roc(y[!tr], s_te, quiet = TRUE)))
}, numeric(1))
cat(sprintf("表观AUC=%.3f  10折CV-AUC=%.3f±%.3f\n", auc(roc_tr), mean(cv_auc), sd(cv_auc)))

# 逐位复算5队列AUC并与AUC_summary.csv比对（冻结证明）
sc <- function(e) {
  hg <- intersect(hub, rownames(e))
  Xv <- t(zrows(e[hg, , drop = FALSE]))
  as.vector(Xv %*% bh[hg, 1, drop = FALSE] + bh["(Intercept)", 1])
}
chk <- function(yy, s, label, ref) {
  r <- roc(yy, s, quiet = TRUE)
  ok <- isTRUE(all.equal(as.numeric(auc(r)), ref, tolerance = 1e-9))
  cat(sprintf("%-22s 复算AUC=%.6f 参考=%.6f 一致=%s\n", label, auc(r), ref, ok)); ok
}
auc_ref <- fread("results/AUC_summary.csv")
oks <- c(
  chk(y, score_tr, "Training", auc_ref[cohort == "GSE96804_train"]$AUC),
  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
i <- 2
for (acc in c("GSE30528", "GSE30529")) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  p <- read.csv(paste0("data/", acc, "_pheno.csv"), stringsAsFactors = FALSE)
  keep <- p$group %in% c("DKD", "Control")
  oks[i] <- chk(as.integer(p$group[keep] == "DKD"), sc(e[, keep, drop = FALSE]), acc,
                auc_ref[cohort == acc]$AUC); i <- i + 1
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
oks[i] <- chk(as.integer(p14$group[k3] != "Control"), sc(e14[, k3, drop = FALSE]), "GSE142025_all",
              auc_ref[cohort == "GSE142025_DN_all"]$AUC); i <- i + 1
oks[i] <- chk(as.integer(p14$group[k4] == "Early_DN"), sc(e14[, k4, drop = FALSE]), "GSE142025_early",
              auc_ref[cohort == "GSE142025_Early"]$AUC); i <- i + 1
e10 <- readRDS("data/GSE104948_expr_symbol.rds")
p10 <- read.csv("data/GSE104948_pheno.csv", stringsAsFactors = FALSE)
k5 <- p10$group %in% c("DN", "Control")
oks[i] <- chk(as.integer(p10$group[k5] == "DN"), sc(e10[, k5, drop = FALSE]), "GSE104948",
              auc_ref[cohort == "GSE104948"]$AUC)
writeLines(c(
  sprintf("lasso_genes: %s", paste(lasso_genes, collapse = ", ")),
  sprintf("svm_rfe best_k: %d", best_k),
  sprintf("svm_rfe genes: %s", paste(rev(svm_genes), collapse = ", ")),
  sprintf("xgb top10: %s", paste(xgb_genes, collapse = ", ")),
  sprintf("hub: %s", paste(hub, collapse = ", ")),
  sprintf("xgb best_iteration: %d", cvb$best_iteration),
  sprintf("AUC reproduce all: %s", all(oks))), "results/frozen/verify_report.txt")

# ========== 5) M1: GSE142025 重定标斜率/截距 + 校准统计复算（冻结系数） ==========
y3 <- as.integer(p14$group[k3] != "Control")
s3 <- sc(e14[, k3, drop = FALSE])
g1 <- glm(y3 ~ s3, family = binomial)
p_hat <- fitted(g1)
suppressWarnings(suppressMessages(library(rms)))
vp <- val.prob(p_hat, y3)  # 输出到屏（不落图）
m1 <- data.frame(metric = c("recalibration_intercept", "recalibration_slope",
                            "Eavg", "Emax", "E90", "Brier", "Spiegelhalter_z", "Spiegelhalter_p"),
                 value = c(coef(g1)[1], coef(g1)[2], vp["Eavg"], vp["Emax"], vp["E90"],
                           vp["Brier"], vp["S:z"], vp["S:p"]))
fwrite(m1, "results/revision/M1_recalibration_glm.csv")
ptr <- plogis(score_tr)
cat("训练集冻结概率分位数(证明0/1退化):", paste(signif(quantile(ptr, c(.01,.5,.99)), 3), collapse=" "), "\n")

# ========== 6) M4: Bootstrap全流程乐观偏差（B=200，候选池固定23） ==========
B <- 200; n <- nrow(X)
sel_cnt <- setNames(rep(0L, ncol(X)), colnames(X))
app <- tst <- rep(NA_real_, B); nh <- integer(B); exact2 <- logical(B)
svm_acc_b <- function(Xb, yb) {
  rk_b <- svm_rfe_rank(Xb, factor(yb, levels = c(0, 1), labels = c("Ctrl", "DN")))
  fb <- sample(rep(1:10, length.out = nrow(Xb)))
  fr <- lapply(1:10, function(f) svm_rfe_rank(Xb[fb != f, , drop = FALSE],
             factor(yb[fb != f], levels = c(0, 1), labels = c("Ctrl", "DN")))
  sz <- 2:min(12, ncol(Xb))
  ac <- vapply(sz, function(k) {
    a <- vapply(1:10, function(f) {
      tr <- fb != f; sel <- tail(fr[[f]], k)
      m <- svm(Xb[tr, sel, drop = FALSE], factor(yb[tr], levels = c(0,1), labels = c("Ctrl","DN")),
               kernel = "linear")
      mean(predict(m, Xb[!tr, sel, drop = FALSE]) == factor(yb[!tr], levels = c(0,1), labels = c("Ctrl","DN")))
    }, numeric(1)); mean(a)
  }, numeric(1))
  list(best = sz[which.max(ac)], genes = tail(rk_b, sz[which.max(ac)]))
}
for (b in 1:B) {
  idx <- sample(n, replace = TRUE)
  Xb <- X[idx, , drop = FALSE]; yb <- y[idx]
  if (length(unique(yb)) < 2) next
  cvb1 <- suppressWarnings(cv.glmnet(Xb, yb, family = "binomial", alpha = 1, nfolds = 10))
  cb <- as.matrix(coef(cvb1, s = "lambda.1se"))
  lb <- setdiff(rownames(cb)[cb[, 1] != 0], "(Intercept)")
  sb <- tryCatch(svm_acc_b(Xb, yb), error = function(e) NULL)
  dtb <- xgb.DMatrix(Xb, label = yb)
  fb5 <- mk_folds(yb, 5)
  cvbb <- xgb.cv(data = dtb, params = xpar, nrounds = 300, folds = fb5,
                 early_stopping_rounds = 20, verbose = 0)
  bstb <- xgboost(x = Xb, y = factor(yb, levels = c(0, 1)), max_depth = 3, learning_rate = 0.1,
                  objective = "binary:logistic", nrounds = max(cvbb$best_iteration, 10), nthread = 2)
  xb <- head(xgb.importance(model = bstb)$Feature, 10)
  hb <- Reduce(intersect, list(lb, if (is.null(sb)) character(0) else sb$genes, xb))
  if (length(hb) < 2) { cn <- table(c(lb, if(is.null(sb)) character(0) else sb$genes, xb)); hb <- names(cn)[cn >= 2] }
  if (length(hb) == 0) next
  nh[b] <- length(hb); sel_cnt[hb] <- sel_cnt[hb] + 1L
  exact2[b] <- identical(sort(hb), c("PDK2", "PDK4"))
  cvf <- suppressWarnings(cv.glmnet(Xb[, hb, drop = FALSE], yb, family = "binomial", alpha = 1, nfolds = 10))
  bb <- as.matrix(coef(cvf, s = "lambda.1se"))
  keep_g <- setdiff(rownames(bb)[bb[, 1] != 0], "(Intercept)")
  if (length(keep_g) == 0) next
  s_app <- as.vector(Xb[, keep_g, drop = FALSE] %*% bb[keep_g, 1] + bb[1, 1])
  s_tst <- as.vector(X[, keep_g, drop = FALSE] %*% bb[keep_g, 1] + bb[1, 1])
  app[b] <- as.numeric(auc(roc(yb, s_app, quiet = TRUE)))
  tst[b] <- as.numeric(auc(roc(y, s_tst, quiet = TRUE)))
}
valid <- !is.na(app)
optimism <- mean(app[valid] - tst[valid])
m4 <- data.frame(
  metric = c("B", "valid_B", "apparent_AUC", "mean_test_AUC", "optimism", "optimism_corrected_AUC",
             "sd_test_AUC", "exact_PDK2PDK4_freq", "mean_n_hub"),
  value = c(B, sum(valid), auc(roc_tr), mean(tst[valid]), optimism,
            as.numeric(auc(roc_tr)) - optimism, sd(tst[valid]), mean(exact2[valid]), mean(nh[valid]))
)
fwrite(m4, "results/revision/M4_bootstrap_pipeline.csv")
sel <- data.frame(gene = names(sel_cnt), selection_freq = as.integer(sel_cnt) / sum(valid))
sel <- sel[order(-sel$selection_freq)]
fwrite(sel, "results/revision/M4_selection_frequency.csv")
cat("\nM4: 校正AUC=", round(m4$value[6], 3), " 乐观偏差=", round(optimism, 3),
    " PDK2+PDK4复现率=", round(mean(exact2[valid]), 2), "\n", sep = "")
cat("选择频率top6:\n"); print(head(sel, 6))
writeLines(capture.output(sessionInfo()), "results/frozen/sessionInfo.txt")
cat("\nFREEZE_DONE\n")

# 15b_freeze_patch.R ----
# 修正15_freeze_ml.R：1) svm_acc_b括号语法错；2) oks记账错；3) 规范系数定为deviance准则
# （支撑稿件β/Table3/S6/校准），验证deviance系数同样逐位复现Table2 AUC与Table3 Youden；
# 4) M1补"冻结概率未重定标"校准对照；5) 重跑M4 Bootstrap。
suppressWarnings(suppressMessages({
  library(data.table); library(glmnet); library(e1071); library(xgboost); library(pROC); library(rms)
}))
set.seed(42)
train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph   <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
cand <- fread("results/A_candidates.csv")$gene
genes <- intersect(cand, rownames(train_expr))
zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
Xdf <- data.frame(t(zrows(train_expr[genes, , drop = FALSE])), check.names = FALSE)
y <- as.integer(train_ph$group == "DN"); X <- as.matrix(Xdf)

# ---- 复现两套系数（调用序列同15：cvfit已耗RNG，这里只需系数拟合段） ----
cvfit <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10, type.measure = "deviance")
cf <- as.matrix(coef(cvfit, s = "lambda.1se"))
lasso_genes <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")
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
sizes <- 2:12; folds <- sample(rep(1:10, length.out = nrow(X)))
fold_rank <- lapply(1:10, function(f) svm_rfe_rank(X[folds != f, , drop = FALSE], yf[folds != f]))
acc <- vapply(sizes, function(k) {
  a <- vapply(1:10, function(f) {
    tr <- folds != f; sel <- tail(fold_rank[[f]], k)
    m <- svm(X[tr, sel, drop = FALSE], yf[tr], kernel = "linear")
    mean(predict(m, X[!tr, sel, drop = FALSE]) == yf[!tr])
  }, numeric(1)); mean(a)
}, numeric(1))
best_k <- sizes[which.max(acc)]; svm_genes <- tail(rk, best_k)
mk_folds <- function(yy, k = 5) {
  f0 <- sample(rep(1:k, length.out = sum(yy == 0))); f1 <- sample(rep(1:k, length.out = sum(yy == 1)))
  lapply(1:k, function(i) c(which(yy == 0)[f0 == i], which(yy == 1)[f1 == i]))
}
folds5 <- mk_folds(y, 5)
xpar <- list(objective = "binary:logistic", eval_metric = "auc", max_depth = 3, eta = 0.1)
dt <- xgb.DMatrix(X, label = y)
cvb <- xgb.cv(data = dt, params = xpar, nrounds = 300, folds = folds5,
              early_stopping_rounds = 20, verbose = 0)
bst <- xgboost(x = X, y = factor(y, levels = c(0, 1)), max_depth = 3, learning_rate = 0.1,
               objective = "binary:logistic", nrounds = max(cvb$best_iteration, 10), nthread = 2)
xgb_genes <- head(xgb.importance(model = bst)$Feature, 10)
hub <- Reduce(intersect, list(lasso_genes, svm_genes, xgb_genes))
if (length(hub) < 2) { cn <- table(c(lasso_genes, svm_genes, xgb_genes)); hub <- names(cn)[cn >= 2] }
Xh <- X[, hub, drop = FALSE]
cv2d <- cv.glmnet(Xh, y, family = "binomial", alpha = 1, nfolds = 10)
cv2a <- cv.glmnet(Xh, y, family = "binomial", alpha = 1, nfolds = 10, type.measure = "auc")
bhd <- as.matrix(coef(cv2d, s = "lambda.1se"))   # 规范（deviance）
bha <- as.matrix(coef(cv2a, s = "lambda.1se"))
fwrite(data.frame(term = rownames(bhd), canonical_deviance = bhd[, 1], auc_criterion = bha[, 1]),
       "results/frozen/final_model_coefficients.csv")
cat("canonical(deviance):", paste(sprintf("%s=%.4f", rownames(bhd), bhd[, 1]), collapse = " "), "\n")

# ---- 双系数逐位验证 Table2 AUC + Table3 Youden ----
sc_b <- function(e, bh) {
  hg <- intersect(hub, rownames(e))
  Xv <- t(zrows(e[hg, , drop = FALSE]))
  as.vector(Xv %*% bh[hg, 1, drop = FALSE] + bh["(Intercept)", 1])
}
auc_ref <- fread("results/AUC_summary.csv"); youd_ref <- fread("results/revision/M6_within_cohort_youden.csv")
youd_chk <- function(yy, s, ref_ss, label, bh_tag) {
  r <- roc(yy, s, quiet = TRUE)
  ci <- coords(r, "best", ret = c("sensitivity", "specificity"), best.method = "youden")
  ok <- isTRUE(all.equal(as.numeric(auc(r)), ref_ss[1], tolerance = 1e-9)) &&
        isTRUE(all.equal(unname(ci$sensitivity), ref_ss[2], tolerance = 1e-9)) &&
        isTRUE(all.equal(unname(ci$specificity), ref_ss[3], tolerance = 1e-9))
  cat(sprintf("%-14s %-4s AUC=%.6f sens=%.4f spec=%.4f 一致=%s\n", label, bh_tag,
              auc(r), ci$sensitivity, ci$specificity, ok)); ok
}
all_ok <- TRUE
gr <- function(co, ss) { r <- youd_ref[cohort == co]; c(ss, r$sensitivity, r$specificity) }
score_tr_d <- as.vector(Xh %*% bhd[-1, , drop = FALSE] + bhd[1, 1])
score_tr_a <- as.vector(Xh %*% bha[-1, , drop = FALSE] + bha[1, 1])
all_ok <- youd_chk(y, score_tr_d, gr("GSE96804_train", 1.0), "Training", "dev") && all_ok
all_ok <- youd_chk(y, score_tr_a, gr("GSE96804_train", 1.0), "Training", "auc") && all_ok
for (acc in c("GSE30528", "GSE30529")) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  p <- read.csv(paste0("data/", acc, "_pheno.csv"), stringsAsFactors = FALSE)
  keep <- p$group %in% c("DKD", "Control")
  ss <- auc_ref[cohort == acc]$AUC
  all_ok <- youd_chk(as.integer(p$group[keep] == "DKD"), sc_b(e[, keep, drop = FALSE], bhd), gr(acc, ss), acc, "dev") && all_ok
  all_ok <- youd_chk(as.integer(p$group[keep] == "DKD"), sc_b(e[, keep, drop = FALSE], bha), gr(acc, ss), acc, "auc") && all_ok
}
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])))
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
rownames(e14) <- toupper(rownames(e14))
k3 <- p14$group %in% c("Advanced_DN", "Early_DN", "Control"); k4 <- p14$group %in% c("Early_DN", "Control")
for (cfg in list(list(k3, as.integer(p14$group[k3] != "Control"), "GSE142025_DN_all", "142025_all"),
                 list(k4, as.integer(p14$group[k4] == "Early_DN"), "GSE142025_Early", "142025_early"))) {
  ss <- auc_ref[cohort == cfg[[3]]]$AUC
  all_ok <- youd_chk(cfg[[2]], sc_b(e14[, cfg[[1]], drop = FALSE], bhd), gr(cfg[[3]], ss), cfg[[4]], "dev") && all_ok
  all_ok <- youd_chk(cfg[[2]], sc_b(e14[, cfg[[1]], drop = FALSE], bha), gr(cfg[[3]], ss), cfg[[4]], "auc") && all_ok
}
e10 <- readRDS("data/GSE104948_expr_symbol.rds")
p10 <- read.csv("data/GSE104948_pheno.csv", stringsAsFactors = FALSE)
k5 <- p10$group %in% c("DN", "Control")
ss <- auc_ref[cohort == "GSE104948"]$AUC
all_ok <- youd_chk(as.integer(p10$group[k5] == "DN"), sc_b(e10[, k5, drop = FALSE], bhd), gr("GSE104948", ss), "104948", "dev") && all_ok
all_ok <- youd_chk(as.integer(p10$group[k5] == "DN"), sc_b(e10[, k5, drop = FALSE], bha), gr("GSE104948", ss), "104948", "auc") && all_ok
cat("双系数全部一致:", all_ok, "\n")
writeLines(c(
  sprintf("lasso_genes: %s", paste(lasso_genes, collapse = ", ")),
  sprintf("svm_rfe best_k: %d", best_k),
  sprintf("xgb top10: %s", paste(xgb_genes, collapse = ", ")),
  sprintf("hub: %s", paste(hub, collapse = ", ")),
  sprintf("canonical coefficients (deviance): %s", paste(sprintf("%s=%.4f", rownames(bhd), bhd[, 1]), collapse = ", ")),
  sprintf("both criteria reproduce Table2 AUC + Table3 Youden bitwise: %s", all_ok)),
  "results/frozen/verify_report.txt")

# ---- M1: 校准（canonical deviance bh，与run_12c同源）+ 未重定标对照 ----
y3 <- as.integer(p14$group[k3] != "Control")
s3 <- sc_b(e14[, k3, drop = FALSE], bhd)
g1 <- glm(y3 ~ s3, family = binomial)
p_hat <- fitted(g1)
vp1 <- suppressWarnings(val.prob(p_hat, y3, pl = FALSE))
p_raw <- plogis(s3)
vp0 <- suppressWarnings(val.prob(p_raw, y3, pl = FALSE))
m1 <- data.frame(
  metric = c("recalibration_intercept", "recalibration_slope",
             "Eavg_recal", "Emax_recal", "E90_recal", "Brier_recal",
             "Spiegelhalter_z_recal", "Spiegelhalter_p_recal",
             "Eavg_frozenprob", "Emax_frozenprob", "Brier_frozenprob", "Spiegelhalter_p_frozenprob"),
  value = c(coef(g1)[1], coef(g1)[2],
            vp1["Eavg"], vp1["Emax"], vp1["E90"], vp1["Brier"], vp1["S:z"], vp1["S:p"],
            vp0["Eavg"], vp0["Emax"], vp0["Brier"], vp0["S:p"]))
fwrite(m1, "results/revision/M1_recalibration_glm.csv")
print(m1)
ptr <- plogis(score_tr_d)
cat("训练集冻结概率分位数1%/50%/99%:", paste(signif(quantile(ptr, c(.01, .5, .99)), 3), collapse = " / "), "\n")

# ---- M4: Bootstrap 全流程（括号已修） ----
B <- 200; n <- nrow(X)
sel_cnt <- setNames(rep(0L, ncol(X)), colnames(X))
app <- tst <- rep(NA_real_, B); nh <- integer(B); exact2 <- logical(B)
svm_acc_b <- function(Xb, yb) {
  yfb <- factor(yb, levels = c(0, 1), labels = c("Ctrl", "DN"))
  rk_b <- svm_rfe_rank(Xb, yfb)
  fb <- sample(rep(1:10, length.out = nrow(Xb)))
  fr <- lapply(1:10, function(f) {
    tr <- fb != f
    svm_rfe_rank(Xb[tr, , drop = FALSE], yfb[tr])
  })
  sz <- 2:min(12, ncol(Xb))
  ac <- vapply(sz, function(k) {
    a <- vapply(1:10, function(f) {
      tr <- fb != f; sel <- tail(fr[[f]], k)
      m <- svm(Xb[tr, sel, drop = FALSE], yfb[tr], kernel = "linear")
      mean(predict(m, Xb[!tr, sel, drop = FALSE]) == yfb[!tr])
    }, numeric(1))
    mean(a)
  }, numeric(1))
  list(best = sz[which.max(ac)], genes = tail(rk_b, sz[which.max(ac)]))
}
t0 <- Sys.time()
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
  sbg <- if (is.null(sb)) character(0) else sb$genes
  hb <- Reduce(intersect, list(lb, sbg, xb))
  if (length(hb) < 2) { cn <- table(c(lb, sbg, xb)); hb <- names(cn)[cn >= 2] }
  if (length(hb) == 0) next
  nh[b] <- length(hb); sel_cnt[hb] <- sel_cnt[hb] + 1L
  exact2[b] <- identical(sort(hb), c("PDK2", "PDK4"))
  cvf <- suppressWarnings(cv.glmnet(Xb[, hb, drop = FALSE], yb, family = "binomial", alpha = 1, nfolds = 10))
  bb <- as.matrix(coef(cvf, s = "lambda.1se"))
  keep_g <- setdiff(rownames(bb)[bb[, 1] != 0], "(Intercept)")
  if (length(keep_g) == 0) next
  s_app <- as.vector(Xb[, keep_g, drop = FALSE] %*% bb[keep_g, 1, drop = FALSE] + bb[1, 1])
  s_tst <- as.vector(X[, keep_g, drop = FALSE] %*% bb[keep_g, 1, drop = FALSE] + bb[1, 1])
  app[b] <- as.numeric(auc(roc(yb, s_app, quiet = TRUE)))
  tst[b] <- as.numeric(auc(roc(y, s_tst, quiet = TRUE)))
  if (b %% 20 == 0) cat(sprintf("bootstrap %d/%d elapsed %.0fs\n", b, B, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
valid <- !is.na(app)
optimism <- mean(app[valid] - tst[valid])
app_orig <- as.numeric(auc(roc(y, score_tr_d, quiet = TRUE)))
m4 <- data.frame(
  metric = c("B", "valid_B", "apparent_AUC", "mean_test_AUC", "sd_test_AUC",
             "optimism", "optimism_corrected_AUC", "exact_PDK2PDK4_freq", "mean_n_hub"),
  value = c(B, sum(valid), app_orig, mean(tst[valid]), sd(tst[valid]), optimism,
            app_orig - optimism, mean(exact2[valid]), mean(nh[valid])))
fwrite(m4, "results/revision/M4_bootstrap_pipeline.csv")
sel <- data.frame(gene = names(sel_cnt), selection_freq = as.integer(sel_cnt) / sum(valid))
sel <- sel[order(-sel$selection_freq)]
fwrite(sel, "results/revision/M4_selection_frequency.csv")
cat("\nM4: 校正AUC=", round(m4$value[7], 3), " 乐观=", round(optimism, 3),
    " PDK2+PDK4复现=", round(mean(exact2[valid]), 2), "\n", sep = "")
print(head(sel, 8))
file.copy("results/revision/M4_bootstrap_pipeline.csv",
          "manuscript/supplement/TableS9_bootstrap_pipeline.csv", overwrite = TRUE)
file.copy("results/wgcna_module_trait.csv",
          "manuscript/supplement/TableS8_module_trait_full.csv", overwrite = TRUE)
writeLines(capture.output(sessionInfo()), "results/frozen/sessionInfo.txt")
cat("\nFREEZE_PATCH_DONE\n")

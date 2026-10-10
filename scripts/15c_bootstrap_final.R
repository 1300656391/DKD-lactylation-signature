# 15c_bootstrap_final.R ----
# 修复：1) Youden参照改为硬编码（稿件Table3+AUC_summary定值）；2) M4单基因hub时跳过建模
suppressWarnings(suppressMessages({
  library(data.table); library(glmnet); library(e1071); library(xgboost); library(pROC)
}))
set.seed(42)
train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph   <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
cand <- fread("results/A_candidates.csv")$gene
genes <- intersect(cand, rownames(train_expr))
zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
Xdf <- data.frame(t(zrows(train_expr[genes, , drop = FALSE])), check.names = FALSE)
y <- as.integer(train_ph$group == "DN"); X <- as.matrix(Xdf)

# ---- 与15b同序复现canonical系数 ----
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
hub <- c("ACSF2", "HIST1H4B", "ME1", "PDK2", "PDK4")
bhd <- as.matrix(coef(cv.glmnet(X[, hub, drop = FALSE], y, family = "binomial", alpha = 1, nfolds = 10),
                      s = "lambda.1se"))
cat("canonical:", paste(sprintf("%s=%.4f", rownames(bhd), bhd[, 1]), collapse = " "), "\n")

# ---- 双系数参照验证（硬编码稿件定值） ----
sc_b <- function(e) {
  hg <- intersect(hub, rownames(e))
  Xv <- t(zrows(e[hg, , drop = FALSE]))
  as.vector(Xv %*% bhd[hg, 1, drop = FALSE] + bhd["(Intercept)", 1])
}
REF <- list(
  Training      = list(yy = y, ss = c(1.000, 1.000, 1.000)),
  GSE30528      = list(ss = c(0.905982905982906, 0.777777777777778, 1)),
  GSE30529      = list(ss = c(0.65, 0.6, 1)),
  GSE142025_all = list(ss = c(0.991769547325103, 0.962962962962963, 1)),
  GSE142025_early = list(ss = c(0.888888888888889, 1, 0.777777777777778)),
  GSE104948     = list(ss = c(1, 1, 1)))
youd_chk <- function(yy, s, ref, label) {
  r <- roc(yy, s, quiet = TRUE)
  ci <- coords(r, "best", ret = c("sensitivity", "specificity"), best.method = "youden")
  ok <- isTRUE(all.equal(round(as.numeric(auc(r)), 9), round(ref[1], 9))) &&
        isTRUE(all.equal(round(unname(ci$sensitivity), 6), round(ref[2], 6))) &&
        isTRUE(all.equal(round(unname(ci$specificity), 6), round(ref[3], 6)))
  cat(sprintf("%-15s AUC=%.6f sens=%.4f spec=%.4f ok=%s\n", label, auc(r), ci$sensitivity, ci$specificity, ok)); ok
}
all_ok <- youd_chk(y, as.vector(X[, hub] %*% bhd[-1, , drop = FALSE] + bhd[1, 1]), REF$Training$ss, "Training")
for (acc in c("GSE30528", "GSE30529")) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  p <- read.csv(paste0("data/", acc, "_pheno.csv"), stringsAsFactors = FALSE)
  keep <- p$group %in% c("DKD", "Control")
  all_ok <- youd_chk(as.integer(p$group[keep] == "DKD"), sc_b(e[, keep, drop = FALSE]), REF[[acc]]$ss, acc) && all_ok
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
all_ok <- youd_chk(as.integer(p14$group[k3] != "Control"), sc_b(e14[, k3, drop = FALSE]),
                   REF$GSE142025_all$ss, "GSE142025_all") && all_ok
all_ok <- youd_chk(as.integer(p14$group[k4] == "Early_DN"), sc_b(e14[, k4, drop = FALSE]),
                   REF$GSE142025_early$ss, "GSE142025_early") && all_ok
e10 <- readRDS("data/GSE104948_expr_symbol.rds")
p10 <- read.csv("data/GSE104948_pheno.csv", stringsAsFactors = FALSE)
k5 <- p10$group %in% c("DN", "Control")
all_ok <- youd_chk(as.integer(p10$group[k5] == "DN"), sc_b(e10[, k5, drop = FALSE]), REF$GSE104948$ss, "GSE104948") && all_ok
cat("canonical系数复现Table2+Table3全部一致:", all_ok, "\n")

# ---- M4 Bootstrap（单基因hub跳过建模） ----
mk_folds <- function(yy, k = 5) {
  f0 <- sample(rep(1:k, length.out = sum(yy == 0))); f1 <- sample(rep(1:k, length.out = sum(yy == 1)))
  lapply(1:k, function(i) c(which(yy == 0)[f0 == i], which(yy == 1)[f1 == i]))
}
xpar <- list(objective = "binary:logistic", eval_metric = "auc", max_depth = 3, eta = 0.1)
svm_acc_b <- function(Xb, yb) {
  yfb <- factor(yb, levels = c(0, 1), labels = c("Ctrl", "DN"))
  rk_b <- svm_rfe_rank(Xb, yfb)
  fb <- sample(rep(1:10, length.out = nrow(Xb)))
  fr <- lapply(1:10, function(f) { tr <- fb != f; svm_rfe_rank(Xb[tr, , drop = FALSE], yfb[tr]) })
  sz <- 2:min(12, ncol(Xb))
  ac <- vapply(sz, function(k) {
    a <- vapply(1:10, function(f) {
      tr <- fb != f; sel <- tail(fr[[f]], k)
      m <- svm(Xb[tr, sel, drop = FALSE], yfb[tr], kernel = "linear")
      mean(predict(m, Xb[!tr, sel, drop = FALSE]) == yfb[!tr])
    }, numeric(1)); mean(a)
  }, numeric(1))
  list(best = sz[which.max(ac)], genes = tail(rk_b, sz[which.max(ac)]))
}
B <- 200; n <- nrow(X)
sel_cnt <- setNames(rep(0L, ncol(X)), colnames(X))
app <- tst <- rep(NA_real_, B); nh <- integer(B); exact2 <- logical(B)
t0 <- Sys.time()
for (b in 1:B) {
  idx <- sample(n, replace = TRUE)
  Xb <- X[idx, , drop = FALSE]; yb <- y[idx]
  if (length(unique(yb)) < 2) next
  cvb1 <- suppressWarnings(cv.glmnet(Xb, yb, family = "binomial", alpha = 1, nfolds = 10))
  cb <- as.matrix(coef(cvb1, s = "lambda.1se"))
  lb <- setdiff(rownames(cb)[cb[, 1] != 0], "(Intercept)")
  sb <- tryCatch(svm_acc_b(Xb, yb), error = function(e) NULL)
  fb5 <- mk_folds(yb, 5)
  cvbb <- xgb.cv(data = xgb.DMatrix(Xb, label = yb), params = xpar, nrounds = 300, folds = fb5,
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
  if (length(hb) < 2) next
  cvf <- suppressWarnings(cv.glmnet(Xb[, hb, drop = FALSE], yb, family = "binomial", alpha = 1, nfolds = 10))
  bb <- as.matrix(coef(cvf, s = "lambda.1se"))
  keep_g <- setdiff(rownames(bb)[bb[, 1] != 0], "(Intercept)")
  if (length(keep_g) == 0) next
  s_app <- as.vector(Xb[, keep_g, drop = FALSE] %*% bb[keep_g, 1, drop = FALSE] + bb[1, 1])
  s_tst <- as.vector(X[, keep_g, drop = FALSE] %*% bb[keep_g, 1, drop = FALSE] + bb[1, 1])
  app[b] <- as.numeric(auc(roc(yb, s_app, quiet = TRUE)))
  tst[b] <- as.numeric(auc(roc(y, s_tst, quiet = TRUE)))
  if (b %% 25 == 0) cat(sprintf("bootstrap %d/%d elapsed %.0fs\n", b, B,
                                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
valid <- !is.na(app)
optimism <- mean(app[valid] - tst[valid])
app_orig <- 1.0
m4 <- data.frame(
  metric = c("B", "valid_B", "apparent_AUC", "mean_test_AUC", "sd_test_AUC",
             "optimism", "optimism_corrected_AUC", "exact_PDK2PDK4_freq", "mean_n_hub"),
  value = c(B, sum(valid), app_orig, mean(tst[valid]), sd(tst[valid]), optimism,
            app_orig - optimism, mean(exact2[valid]), mean(nh[valid])))
fwrite(m4, "results/revision/M4_bootstrap_pipeline.csv")
sel <- data.frame(gene = names(sel_cnt), selection_freq = as.integer(sel_cnt) / sum(valid))
sel <- sel[order(sel$selection_freq, decreasing = TRUE), , drop = FALSE]
fwrite(sel, "results/revision/M4_selection_frequency.csv")
cat("\nM4 summary:\n"); print(m4)
cat("Selection frequency top8:\n"); print(head(sel, 8))
file.copy("results/revision/M4_bootstrap_pipeline.csv",
          "manuscript/supplement/TableS9_bootstrap_pipeline.csv", overwrite = TRUE)
file.copy("results/wgcna_module_trait.csv",
          "manuscript/supplement/TableS8_module_trait_full.csv", overwrite = TRUE)
cat("\nBOOTSTRAP_DONE\n")

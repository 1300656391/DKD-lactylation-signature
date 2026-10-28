# 07_ml_screening.R ----
# Three-algorithm feature screening (LASSO + SVM-RFE + XGBoost) on the 23
# candidate genes, consensus hub-gene selection, penalized diagnostic model,
# and initial validation across the training and four microarray/RNA-seq
# cohorts (GSE273763 and GSE294519 validation lives in 09_validation.R).
# Outputs: results/hub_genes.txt, results/AUC_summary.csv,
#          results/AUC_single_gene.csv, results/ml_screening.txt

suppressWarnings(suppressMessages({
  library(data.table); library(glmnet); library(e1071); library(xgboost); library(pROC)
}))
set.seed(42)

train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
cand <- fread("results/A_candidates.csv")$gene
genes <- intersect(cand, rownames(train_expr))
cat("candidates entering ML:", length(genes), ":", paste(genes, collapse = ", "), "\n")

zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
Xdf <- data.frame(t(zrows(train_expr[genes, , drop = FALSE])), check.names = FALSE)
y <- as.integer(train_ph$group == "DN")
X <- as.matrix(Xdf)

# ---------- LASSO (10-fold CV, binomial deviance) ----------
cvfit <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10,
                   type.measure = "deviance")
cf <- as.matrix(coef(cvfit, s = "lambda.1se"))
lasso_genes <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")
cat("LASSO selected (", length(lasso_genes), "):", paste(lasso_genes, collapse = ", "), "\n")

# ---------- SVM-RFE (linear kernel, 10-fold CV accuracy over k=2..12) ----------
svm_rfe_rank <- function(Xm, yy) {
  feats <- colnames(Xm); ranked <- character(0)
  while (length(feats) > 1) {
    m <- svm(Xm[, feats, drop = FALSE], yy, kernel = "linear", scale = FALSE)
    w <- (t(m$coefs) %*% m$SV)^2
    wr <- as.vector(w); names(wr) <- feats
    rem <- names(which.min(wr)); ranked <- c(rem, ranked); feats <- setdiff(feats, rem)
  }
  c(ranked, feats)  # least important first
}
yf <- factor(y, levels = c(0, 1), labels = c("Ctrl", "DN"))
rk <- svm_rfe_rank(X, yf)
sizes <- 2:min(12, ncol(X))
folds <- sample(rep(1:10, length.out = nrow(X)))
acc <- vapply(sizes, function(k) {
  a <- vapply(1:10, function(f) {
    tr <- folds != f
    rk2 <- svm_rfe_rank(X[tr, , drop = FALSE], yf[tr])
    sel <- tail(rk2, k)
    m <- svm(X[tr, sel, drop = FALSE], yf[tr], kernel = "linear")
    mean(predict(m, X[!tr, sel, drop = FALSE]) == yf[!tr])
  }, numeric(1))
  mean(a)
}, numeric(1))
best_k <- sizes[which.max(acc)]
svm_genes <- tail(rk, best_k)
cat("SVM-RFE selected (k =", best_k, "):", paste(rev(svm_genes), collapse = ", "), "\n")

# ---------- XGBoost (stratified 5-fold early stopping, top-10 gain) ----------
dtrain <- xgb.DMatrix(X, label = y)
mk_folds <- function(yy, k = 5) {
  f0 <- sample(rep(1:k, length.out = sum(yy == 0)))
  f1 <- sample(rep(1:k, length.out = sum(yy == 1)))
  lapply(1:k, function(i) c(which(yy == 0)[f0 == i], which(yy == 1)[f1 == i]))
}
folds5 <- mk_folds(y, 5)
xpar <- list(objective = "binary:logistic", eval_metric = "auc",
             max_depth = 3, eta = 0.1)
cvb <- xgb.cv(data = dtrain, params = xpar, nrounds = 300, folds = folds5,
              early_stopping_rounds = 20, verbose = 0)
bst <- xgboost(x = X, y = factor(y, levels = c(0, 1)),
               max_depth = 3, learning_rate = 0.1, objective = "binary:logistic",
               nrounds = max(cvb$best_iteration, 10), nthread = 2)
imp <- xgb.importance(model = bst)
xgb_genes <- head(imp$Feature, 10)
cat("XGBoost top:", paste(xgb_genes, collapse = ", "), "\n")

# ---------- Consensus -> hub genes ----------
hub <- Reduce(intersect, list(lasso_genes, svm_genes, xgb_genes))
decision <- "three-algorithm intersection"
if (length(hub) < 2) {
  cnt <- table(c(lasso_genes, svm_genes, xgb_genes))
  hub <- names(cnt)[cnt >= 2]
  decision <- "selected by at least two algorithms"
}
cat("decision rule:", decision, "\nhub genes (", length(hub), "):",
    paste(hub, collapse = ", "), "\n")
writeLines(hub, "results/hub_genes.txt")

# ---------- Penalized logistic model on hub genes (AUC criterion) ----------
hub <- intersect(hub, colnames(Xdf))
Xh <- as.matrix(Xdf[, hub, drop = FALSE])
cv2 <- cv.glmnet(Xh, y, family = "binomial", alpha = 1, nfolds = 10,
                 type.measure = "auc")
bh <- as.matrix(coef(cv2, s = "lambda.1se"))
score_tr <- as.vector(Xh %*% bh[-1, , drop = FALSE] + bh["(Intercept)", 1])
roc_tr <- roc(y, score_tr, quiet = TRUE)

# Honest 10-fold CV AUC (refit glmnet inside every fold)
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
cat("training AUC (apparent):", round(auc(roc_tr), 3),
    "| 10-fold CV AUC:", round(mean(cv_auc), 3), "+/-", round(sd(cv_auc), 3), "\n")

# Single-gene ROC (training)
sg_auc <- t(sapply(hub, function(g) {
  r <- roc(y, Xdf[[g]], quiet = TRUE)
  c(AUC = as.numeric(auc(r)), lower = ci.auc(r)[1], upper = ci.auc(r)[3])
}))
sg_auc <- as.data.frame(sg_auc); sg_auc$gene <- rownames(sg_auc)
fwrite(sg_auc, "results/AUC_single_gene.csv")

# ---------- Validation across cohorts ----------
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

for (acc in c("GSE30528", "GSE30529")) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  p <- read.csv(paste0("data/", acc, "_pheno.csv"), stringsAsFactors = FALSE)
  auc_tab <- rbind(auc_tab, validate_cohort(e, p, "DKD", "Control", acc))
}

# GSE142025: align GSM column names via the meta file
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(sample = cl(ml[grepl("^!Sample_title", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])),
                  gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  stringsAsFactors = FALSE)
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
stopifnot(!any(is.na(match(colnames(e14), p14$gsm))))
rownames(e14) <- toupper(rownames(e14))
auc_tab <- rbind(auc_tab,
  validate_cohort(e14, p14, c("Advanced_DN", "Early_DN"), "Control", "GSE142025_DN_all"),
  validate_cohort(e14, p14, "Early_DN", "Control", "GSE142025_Early"),
  validate_cohort(e14, p14, "Advanced_DN", "Control", "GSE142025_Advanced"))

e10 <- readRDS("data/GSE104948_expr_symbol.rds")
p10 <- read.csv("data/GSE104948_pheno.csv", stringsAsFactors = FALSE)
auc_tab <- rbind(auc_tab, validate_cohort(e10, p10, "DN", "Control", "GSE104948"))

fwrite(auc_tab, "results/AUC_summary.csv")

# Run log for the record
writeLines(c(
  sprintf("lasso_genes: %s", paste(lasso_genes, collapse = ", ")),
  sprintf("svm_rfe best_k: %d", best_k),
  sprintf("svm_rfe genes: %s", paste(rev(svm_genes), collapse = ", ")),
  sprintf("xgb top10: %s", paste(xgb_genes, collapse = ", ")),
  sprintf("hub: %s", paste(hub, collapse = ", ")),
  sprintf("xgb best_iteration: %d", cvb$best_iteration),
  sprintf("training AUC apparent: %.3f  10-fold CV: %.3f +/- %.3f",
          as.numeric(auc(roc_tr)), mean(cv_auc), sd(cv_auc)),
  sprintf("decision rule: %s", decision)), "results/ml_screening.txt")

cat("\nML_SCREENING_DONE\n")

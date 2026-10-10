# 12b_revision_patch.R ----
# 修复: val.prob参数 + 队列内Youden敏感度/特异度表
suppressWarnings(suppressMessages({
  library(data.table); library(glmnet); library(pROC); library(rms); library(dcurves)
}))
set.seed(42)
dir.create("results/revision", showWarnings = FALSE, recursive = TRUE)

train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
hub <- readLines("results/hub_genes.txt")
zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
Xh <- data.frame(t(zrows(train_expr[hub, , drop = FALSE])), check.names = FALSE)
y <- as.integer(train_ph$group == "DN")
cv2 <- cv.glmnet(as.matrix(Xh), y, family = "binomial", alpha = 1, nfolds = 10)
bh <- as.matrix(coef(cv2, s = "lambda.1se"))
score_tr <- as.vector(as.matrix(Xh) %*% bh[-1, , drop = FALSE] + bh["(Intercept)", 1])

score_cohort <- function(e, yy) {
  hg <- intersect(hub, rownames(e))
  Xv <- data.frame(t(zrows(e[hg, , drop = FALSE])), check.names = FALSE)
  s <- rep(bh["(Intercept)", 1], ncol(e))
  for (g in hg) s <- s + bh[g, 1] * Xv[[g]]
  s
}
cohort_data <- function(e, p, case_groups, ctrl_groups) {
  keep <- p$group %in% c(case_groups, ctrl_groups)
  list(e = e[, keep, drop = FALSE], y = as.integer(p$group[keep] %in% case_groups))
}

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

sets <- list(
  list(label = "GSE96804_train", score = score_tr, yy = y),
  list(label = "GSE30528", score = score_cohort(vsets$GSE30528$e, vsets$GSE30528$y), yy = vsets$GSE30528$y),
  list(label = "GSE30529", score = score_cohort(vsets$GSE30529$e, vsets$GSE30529$y), yy = vsets$GSE30529$y),
  list(label = "GSE142025_all", score = {
    cd <- cohort_data(e14, p14, c("Advanced_DN", "Early_DN"), "Control")
    score_cohort(cd$e, cd$y) }, yy = {
    cd <- cohort_data(e14, p14, c("Advanced_DN", "Early_DN"), "Control"); cd$y }),
  list(label = "GSE142025_early", score = {
    cd <- cohort_data(e14, p14, "Early_DN", "Control")
    score_cohort(cd$e, cd$y) }, yy = {
    cd <- cohort_data(e14, p14, "Early_DN", "Control"); cd$y }),
  list(label = "GSE104948", score = {
    cd <- cohort_data(e10, p10, "DN", "Control"); score_cohort(cd$e, cd$y) }, yy = {
    cd <- cohort_data(e10, p10, "DN", "Control"); cd$y })
)

# 队列内Youden: 敏感度/特异度
tab <- rbindlist(lapply(sets, function(s) {
  r <- roc(s$yy, s$score, quiet = TRUE)
  co <- coords(r, x = "best", ret = c("threshold", "sensitivity", "specificity"),
               best.method = "youden")
  data.frame(cohort = s$label, AUC = as.numeric(auc(r)),
             cutoff = as.numeric(co["threshold"]),
             sensitivity = as.numeric(co["sensitivity"]),
             specificity = as.numeric(co["specificity"]))
}))
fwrite(tab, "results/revision/M6_within_cohort_youden.csv")
cat("队列内Youden性能:\n"); print(tab)

# 校准曲线（GSE142025_all）—— 简化参数
sc3 <- sets[[4]]$score; y3 <- sets[[4]]$yy
png("results/figures/Fig7B_calibration_GSE142025.png", type = "windows",
    width = 1400, height = 1400, res = 240)
val.prob(sc3, y3)
title("Calibration - GSE142025 (all DN vs control)")
dev.off()

dca_df <- data.frame(y = y3, Signature = fitted(glm(y3 ~ sc3, family = binomial)))
png("results/figures/Fig7C_dca_GSE142025.png", type = "windows",
    width = 1500, height = 1300, res = 240)
tryCatch(print(plot(dca(y ~ Signature, data = dca_df), smooth = FALSE)),
         error = function(e) cat("DCA失败:", conditionMessage(e), "\n"))
dev.off()

cat("\nPATCH_DONE\n")

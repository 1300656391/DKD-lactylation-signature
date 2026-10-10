# 12c_calibration_fix.R ----
# 校准曲线：score→概率后传入val.prob
suppressWarnings(suppressMessages({ library(pROC); library(rms); library(dcurves) }))
set.seed(42)
source_lines <- NULL
# 复用12b的数据构造（简化重算）
library(data.table); library(glmnet)
train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
hub <- readLines("results/hub_genes.txt")
zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
Xh <- data.frame(t(zrows(train_expr[hub, , drop = FALSE])), check.names = FALSE)
y <- as.integer(train_ph$group == "DN")
cv2 <- cv.glmnet(as.matrix(Xh), y, family = "binomial", alpha = 1, nfolds = 10)
bh <- as.matrix(coef(cv2, s = "lambda.1se"))
sc <- function(e, yy) {
  hg <- intersect(hub, rownames(e))
  Xv <- data.frame(t(zrows(e[hg, , drop = FALSE])), check.names = FALSE)
  s <- rep(bh["(Intercept)", 1], ncol(e))
  for (g in hg) s <- s + bh[g, 1] * Xv[[g]]
  s
}
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])),
                  stringsAsFactors = FALSE)
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
rownames(e14) <- toupper(rownames(e14))
keep <- p14$group %in% c("Advanced_DN", "Early_DN", "Control")
e14k <- e14[, keep]; y3 <- as.integer(p14$group[keep] != "Control")
s3 <- sc(e14k, y3)

p_hat <- fitted(glm(y3 ~ s3, family = binomial))
png("results/figures/Fig7B_calibration_GSE142025.png", type = "windows",
    width = 1400, height = 1400, res = 240)
val.prob(p_hat, y3)
title("Calibration - GSE142025 (DN vs control)")
dev.off()

dca_df <- data.frame(y = y3, Signature = p_hat)
png("results/figures/Fig7C_dca_GSE142025.png", type = "windows",
    width = 1500, height = 1300, res = 240)
tryCatch(print(plot(dca(y ~ Signature, data = dca_df), smooth = FALSE)),
         error = function(e) cat("DCA失败:", conditionMessage(e), "\n"))
dev.off()
cat("CALIBRATION_DONE\n")

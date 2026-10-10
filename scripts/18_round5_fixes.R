# 18_round5_fixes.R ----
# Round-5 修复分析：
#   M2: 各队列 PDK2-only / PDK4-only / 两基因签名 AUC 并列 + DeLong(签名 vs PDK4)
#       同时输出 PDK4/PDK2 中位数方向（供 M1 方向冲突讨论：区室解释的直接证据）
#   M4: Table3 sens/spec 精确二项 95% CI（Youden 截断为样本内最优，题注说明）
#   m3: immune_gene_corr 各细胞平均|r|（验证"myeloid and T-helper"措辞）
suppressWarnings(suppressMessages({ library(data.table); library(pROC) }))
hub <- readLines("results/hub_genes.txt")
fm <- fread("results/frozen/final_model_coefficients.csv")
bh <- setNames(fm$canonical_deviance, fm$term)
train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph   <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
sig  <- function(e) { hg <- intersect(hub, rownames(e)); Xv <- t(zrows(e[hg, , drop = FALSE]))
                      as.vector(Xv %*% bh[hg] + bh[["(Intercept)"]]) }
one  <- function(e, g) as.numeric(zrows(as.matrix(e[g, , drop = FALSE]))[1, ])
med  <- function(e, g, yy) { m <- as.matrix(e[g, , drop = FALSE])
  c(median(m[1, yy == 1]), median(m[1, yy == 0])) }

row_out <- list(); i <- 0
addrow <- function(label, e, yy, n_case, n_ctrl) {
  s <- sig(e); p4 <- one(e, "PDK4"); p2 <- one(e, "PDK2")
  m4v <- med(e, "PDK4", yy); m2v <- med(e, "PDK2", yy)
  p4o <- p4 * ifelse(m4v[1] < m4v[2], -1, 1)   # 方向对齐：病例侧高分，DeLong 有效
  r_s <- roc(yy, s, quiet = TRUE); r4 <- roc(yy, p4o, quiet = TRUE)
  r2 <- roc(yy, p2, direction = "auto", quiet = TRUE)
  dl <- tryCatch(roc.test(r_s, r4, method = "delong")$p.value, error = function(e) NA)
  i <<- i + 1
  row_out[[i]] <<- data.frame(
    cohort = label, n_case = n_case, n_ctrl = n_ctrl,
    AUC_signature = round(as.numeric(auc(r_s)), 3),
    AUC_PDK4 = round(as.numeric(auc(r4)), 3), PDK4_dir = ifelse(m4v[1] < m4v[2], "down", "up"),
    AUC_PDK2 = round(as.numeric(auc(r2)), 3), PDK2_dir = ifelse(m2v[1] < m2v[2], "down", "up"),
    DeLong_sig_vs_PDK4 = ifelse(is.na(dl), NA, signif(dl, 3)),
    PDK4_median_case = round(m4v[1], 2), PDK4_median_ctrl = round(m4v[2], 2),
    PDK2_median_case = round(m2v[1], 2), PDK2_median_ctrl = round(m2v[2], 2))
}
y <- as.integer(train_ph$group == "DN")
addrow("GSE96804_train", train_expr, y, sum(y), sum(1 - y))
for (acc in c("GSE30528", "GSE30529")) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  p <- read.csv(paste0("data/", acc, "_pheno.csv"), stringsAsFactors = FALSE)
  keep <- p$group %in% c("DKD", "Control")
  yy <- as.integer(p$group[keep] == "DKD")
  addrow(acc, e[, keep, drop = FALSE], yy, sum(yy), sum(1 - yy))
}
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])))
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
rownames(e14) <- toupper(rownames(e14))
k3 <- p14$group %in% c("Advanced_DN", "Early_DN", "Control"
)
k4 <- p14$group %in% c("Early_DN", "Control"); k5 <- p14$group %in% c("Advanced_DN", "Control")
addrow("GSE142025_all", e14[, k3, drop = FALSE], as.integer(p14$group[k3] != "Control"), sum(p14$group[k3] != "Control"), sum(p14$group[k3] == "Control"))
addrow("GSE142025_early", e14[, k4, drop = FALSE], as.integer(p14$group[k4] == "Early_DN"), sum(p14$group[k4] == "Early_DN"), sum(p14$group[k4] == "Control"))
addrow("GSE142025_advanced", e14[, k5, drop = FALSE], as.integer(p14$group[k5] == "Advanced_DN"), sum(p14$group[k5] == "Advanced_DN"), sum(p14$group[k5] == "Control"))
e10 <- readRDS("data/GSE104948_expr_symbol.rds")
p10 <- read.csv("data/GSE104948_pheno.csv", stringsAsFactors = FALSE)
k6 <- p10$group %in% c("DN", "Control")
addrow("GSE104948", e10[, k6, drop = FALSE], as.integer(p10$group[k6] == "DN"), sum(p10$group[k6] == "DN"), sum(p10$group[k6] == "Control"))
m2 <- rbindlist(row_out)
fwrite(m2, "results/revision/M2_single_vs_signature.csv")
file.copy("results/revision/M2_single_vs_signature.csv",
          "manuscript/supplement/TableS10_single_vs_signature.csv", overwrite = TRUE)
print(m2)

# M4: Table3 精确二项 CI
youd <- fread("results/revision/M6_within_cohort_youden.csv")
ci_row <- function(tab) {
  y <- as.integer(train_ph$group == "DN"); s <- sig(train_expr)
  r <- roc(y, s, quiet = TRUE); ci <- coords(r, "best", ret = c("sensitivity", "specificity"), best.method = "youden")
  data.frame(sens_lo = binom.test(round(ci$sensitivity * sum(y)), sum(y))$conf.int[1],
             sens_hi = binom.test(round(ci$sensitivity * sum(y)), sum(y))$conf.int[2],
             spec_lo = binom.test(round(ci$specificity * sum(1 - y)), sum(1 - y))$conf.int[1],
             spec_hi = binom.test(round(ci$specificity * sum(1 - y)), sum(1 - y))$conf.int[2])
}
out <- list()
for (j in seq_len(nrow(youd))) {
  lab <- youd$cohort[j]
  if (lab == "GSE96804_train") { yy <- as.integer(train_ph$group == "DN"); e <- train_expr; nc <- sum(yy); nk <- sum(1 - yy) }
  else if (lab %in% c("GSE30528", "GSE30529")) {
    e <- readRDS(paste0("data/", lab, "_expr_symbol.rds"))
    p <- read.csv(paste0("data/", lab, "_pheno.csv"), stringsAsFactors = FALSE)
    keep <- p$group %in% c("DKD", "Control"); e <- e[, keep, drop = FALSE]
    yy <- as.integer(p$group[keep] == "DKD"); nc <- sum(yy); nk <- sum(1 - yy)
  } else if (lab == "GSE142025_all") { yy <- as.integer(p14$group[k3] != "Control"); e <- e14[, k3, drop = FALSE]; nc <- sum(yy); nk <- sum(1 - yy) }
  else if (lab == "GSE142025_early") { yy <- as.integer(p14$group[k4] == "Early_DN"); e <- e14[, k4, drop = FALSE]; nc <- sum(yy); nk <- sum(1 - yy) }
  else { yy <- as.integer(p10$group[k6] == "DN"); e <- e10[, k6, drop = FALSE]; nc <- sum(yy); nk <- sum(1 - yy) }
  r <- roc(yy, sig(e), quiet = TRUE)
  ci <- coords(r, "best", ret = c("sensitivity", "specificity"), best.method = "youden")
  ns <- round(ci$sensitivity * nc); nn <- round((1 - ci$specificity) * nk)
  bs <- binom.test(ns, nc); bn <- binom.test(round(ci$specificity * nk), nk)
  out[[j]] <- data.frame(cohort = lab, sens = youd$sensitivity[j],
    sens_CI = sprintf("%.2f-%.2f", bs$conf.int[1], bs$conf.int[2]),
    spec = youd$specificity[j],
    spec_CI = sprintf("%.2f-%.2f", bn$conf.int[1], bn$conf.int[2]))
}
m4 <- rbindlist(out)
fwrite(m4, "results/revision/M4b_youden_ci.csv")
print(m4)

# m3: 各细胞平均|r|
ic <- read.csv("results/immune_gene_corr.csv", stringsAsFactors = FALSE, check.names = FALSE)
cm <- as.matrix(ic[, setdiff(colnames(ic), c("gene", "Normal_mucosa"))])
ar <- sort(round(colMeans(abs(cm)), 3), decreasing = TRUE)
cat("\nMean |rho| by cell:\n"); print(head(ar, 6))
cat("\nROUND5_FIXES_DONE\n")

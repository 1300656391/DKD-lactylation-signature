# 19c_train_z_sensitivity.R ----
# Sensitivity companion to the Fig7A scoring convention:
# re-scores the four RNA validation cohorts with TRAINING-SET PDK2/PDK4-panel
# means and SDs for the z step (instead of the within-cohort z used by
# 13_figures.R / Fig7A), keeping the frozen canonical coefficients.
# Expected: cross-platform AUC degradation -- reported as-is; demonstrates
# that within-cohort z is a scoring convention, not a single-sample assay.
# Run inside the project root (same prerequisites as 13_figures.R:
# data/GSE96804_expr_symbol.rds, data/GSE30528-*, GSE142025_*,
# GSE104948_*, results/hub_genes.txt, results/frozen/final_model_coefficients.csv).
# Output: results/revision/M5c_train_z_sensitivity.csv
suppressWarnings(suppressMessages({ library(data.table); library(pROC) }))

hub <- readLines("results/hub_genes.txt")
fm <- fread("results/frozen/final_model_coefficients.csv")
bhf <- setNames(fm$canonical_deviance, fm$term)

zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))

train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
hg0 <- hub[hub %in% rownames(train_expr)]
mu_tr <- apply(train_expr[hg0, , drop = FALSE], 1, mean)
sd_tr <- apply(train_expr[hg0, , drop = FALSE], 1, sd)

score_within <- function(e) {
  hg <- intersect(hub, rownames(e))
  Xv <- t(zrows(e[hg, , drop = FALSE]))
  as.vector(Xv %*% bhf[hg] + bhf[["(Intercept)"]])
}
score_trainz <- function(e) {
  hg <- intersect(hub, rownames(e))
  m <- e[hg, , drop = FALSE]
  Z <- (m - mu_tr[hg]) / sd_tr[hg]
  as.vector(t(Z) %*% bhf[hg] + bhf[["(Intercept)"]])
}

res <- list(); i <- 0
add_row <- function(cohort, y, s_w, s_t, note) {
  i <<- i + 1
  res[[i]] <<- data.frame(
    cohort = cohort,
    n_case = sum(y == 1), n_ctrl = sum(y == 0),
    AUC_within_cohort_z = round(as.numeric(auc(roc(y, s_w, quiet = TRUE))), 3),
    AUC_training_z = round(as.numeric(auc(roc(y, s_t, quiet = TRUE))), 3),
    delta = round(as.numeric(auc(roc(y, s_t, quiet = TRUE))) -
                  as.numeric(auc(roc(y, s_w, quiet = TRUE))), 3),
    note = note, stringsAsFactors = FALSE)
}

# Training reference: within-z and training-z are the same vectors here --
# kept as an internal consistency check.
train_ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
y_tr <- as.integer(train_ph$group == "DN")
add_row("Training_GSE96804", y_tr, score_within(train_expr), score_trainz(train_expr),
        "training cohort (reference; z conventions coincide)")

for (acc in c("GSE30528", "GSE30529")) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  p <- read.csv(paste0("data/", acc, "_pheno.csv"), stringsAsFactors = FALSE)
  keep <- p$group %in% c("DKD", "Control")
  y <- as.integer(p$group[keep] == "DKD")
  lbl <- ifelse(acc == "GSE30528", "GSE30528_glomeruli", "GSE30529_tubules")
  add_row(lbl, y, score_within(e[, keep, drop = FALSE]),
          score_trainz(e[, keep, drop = FALSE]), "within-cohort z vs training-z")
}
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])))
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
rownames(e14) <- toupper(rownames(e14))
k3 <- p14$group %in% c("Advanced_DN", "Early_DN", "Control")
y14 <- as.integer(p14$group[k3] != "Control")
add_row("GSE142025_allDN", y14, score_within(e14[, k3, drop = FALSE]),
        score_trainz(e14[, k3, drop = FALSE]), "all DN vs controls")
e10 <- readRDS("data/GSE104948_expr_symbol.rds")
p10 <- read.csv("data/GSE104948_pheno.csv", stringsAsFactors = FALSE)
y10 <- as.integer(p10$group == "DN")
add_row("GSE104948_pooled", y10, score_within(e10), score_trainz(e10),
        "DN vs LD/HT/TN controls, both platforms")

out <- rbindlist(res)
fwrite(out, "results/revision/M5c_train_z_sensitivity.csv")
print(out)
# sanity: training row must have identical AUCs under both conventions
stopifnot(abs(out$AUC_within_cohort_z[1] - out$AUC_training_z[1]) < 1e-6)
cat("\nTRAINZ_DONE\n")

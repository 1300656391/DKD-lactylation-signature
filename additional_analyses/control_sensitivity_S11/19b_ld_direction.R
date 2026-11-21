# 19b_ld_direction.R ----
# Companion to 19_ld_sensitivity.R (Table S11):
# adds as-trained, observed-oriented, and observed-direction AUCs for PDK4 and
# PDK2 in the GSE104948 living-donor / hypertensive / nephrectomy comparisons,
# all platforms, plus a direction_reversed flag (observed != training direction;
# training direction is DOWN in DN for both genes).
# C-2 (2026-11-21): the former AUC_PDK4_as_trained column was computed with
# pROC's default direction="auto" and therefore carried the ORIENTED value;
# it is now split into AUC_PDK4_as_trained (fixed training direction) and
# AUC_PDK4_oriented (observed-direction), with PDK4_dir renamed to
# PDK4_dir_observed and a direction_reversed column appended.
# Run inside the project root (prerequisites: data/GSE104948-*_expr_symbol.rds,
# results/hub_genes.txt, results/frozen/final_model_coefficients.csv).
# Outputs: results/revision/M5b_gse104948_direction.csv
#          updates 投稿包_FIG/03_补充材料/ESM_13.csv when that file exists
#          (matched by comparison; pooled row keeps blank single-gene cells).
suppressWarnings(suppressMessages({ library(data.table); library(pROC) }))

hub <- readLines("results/hub_genes.txt")
fm <- fread("results/frozen/final_model_coefficients.csv")
bh <- setNames(fm$canonical_deviance, fm$term)

zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
p_of <- function(e, g) as.numeric(zrows(as.matrix(e[g, , drop = FALSE]))[1, ])
# AUC with cases higher than controls (deterministic; no pROC auto-flip)
auc_hi <- function(score, y) as.numeric(auc(roc(y, score, quiet = TRUE, direction = ">")))

grp_from_title <- function(p) {
  ml <- readLines(paste0("data/GSE104948-", p, "_meta.txt"), warn = FALSE)
  gsm <- trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", ml[grepl("^!Sample_geo_accession", ml)][1]), "\t", fixed = TRUE)[[1]]))
  ttl <- trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", ml[grepl("^!Sample_title", ml)][1]), "\t", fixed = TRUE)[[1]]))
  data.frame(gsm = gsm, grp = sub("^[A-Z][0-9]+-Glom-", "", sub("[0-9]+$", "", ttl)), stringsAsFactors = FALSE)
}
load_expr <- function(p) readRDS(paste0("data/GSE104948-", p, "_expr_symbol.rds"))

res <- list(); i <- 0
add_dir <- function(label, e, case_idx, ctrl_idx, note) {
  es <- as.matrix(e[, c(case_idx, ctrl_idx), drop = FALSE])
  y <- c(rep(1, length(case_idx)), rep(0, length(ctrl_idx)))
  # calibers (empirically verified: pROC direction=">" returns the AUC of the
  # "cases lower than controls" ordering, which IS the as-trained caliber --
  # training direction is PDK2/PDK4 down in DN):
  a4_at <- auc_hi(p_of(es, "PDK4"), y)
  a2_at <- auc_hi(p_of(es, "PDK2"), y)
  a4_or <- max(a4_at, 1 - a4_at)   # observed-direction oriented
  a2_or <- max(a2_at, 1 - a2_at)   # = 1 - a4_at / 1 - a2_at by construction
  dir4 <- ifelse(median(es["PDK4", y == 1]) < median(es["PDK4", y == 0]), "down", "up")
  dir2 <- ifelse(median(es["PDK2", y == 1]) < median(es["PDK2", y == 0]), "down", "up")
  rev4 <- dir4 != "down"
  stopifnot(abs(a4_or - max(a4_at, 1 - a4_at)) < 1e-9,
            abs(a2_or - max(a2_at, 1 - a2_at)) < 1e-9)  # oriented = max of the two calibers
  i <<- i + 1
  res[[i]] <<- data.frame(
    comparison = label, n_case = length(case_idx), n_ctrl = length(ctrl_idx),
    AUC_PDK4_as_trained = round(a4_at, 3),
    AUC_PDK4_oriented   = round(a4_or, 3),
    PDK4_dir_observed   = dir4,
    PDK4_median_case = round(median(es["PDK4", y == 1]), 3),
    PDK4_median_ctrl = round(median(es["PDK4", y == 0]), 3),
    AUC_PDK2_as_trained = round(a2_at, 3),
    AUC_PDK2_oriented   = round(a2_or, 3),
    PDK2_dir_observed   = dir2,
    PDK2_median_case = round(median(es["PDK2", y == 1]), 3),
    PDK2_median_ctrl = round(median(es["PDK2", y == 0]), 3),
    direction_reversed = ifelse(rev4, "True", "False"),
    note = note, stringsAsFactors = FALSE)
}

for (p in c("GPL22945", "GPL24120")) {
  g <- grp_from_title(p); e <- load_expr(p)
  dn <- which(g$grp == "DN"); ld <- which(g$grp == "LD_"); ht <- which(g$grp == "HT"); tn <- which(g$grp == "TN")
  cat(p, ": DN", length(dn), "LD", length(ld), "HT", length(ht), "TN", length(tn), "\n")
  if (length(ld) > 0) add_dir(paste0("DN_vs_LD_", p), e, dn, ld, "living-donor control (direction disclosure)")
  if (length(ht) > 0) add_dir(paste0("DN_vs_HT_", p), e, dn, ht, "hypertensive control (direction disclosure)")
  if (length(tn) > 0) add_dir(paste0("DN_vs_TN_", p), e, dn, tn, "tumor-nephrectomy control (cross-check vs Table S18)")
}

out <- rbindlist(res)
if (!dir.exists("results/revision")) dir.create("results/revision", recursive = TRUE)
fwrite(out, "results/revision/M5b_gse104948_direction.csv")
print(out)

# ---- update Table S11 (ESM_13.csv) when staged: split the per-gene AUC
# columns and add direction_reversed; matched by comparison; pooled row keeps
# blank single-gene cells. NOTE: staged through ASCII intermediate paths
# (results/revision/ESM_13_in.csv -> ESM_13_updated.csv) because the project's
# CJK package path is not representable in the R session's locale on this
# machine; the caller copies the staged file back over the package file. ----
esm13_in <- "results/revision/ESM_13_in.csv"
esm13_out <- "results/revision/ESM_13_updated.csv"
if (file.exists(esm13_in)) {
  e13 <- read.csv(esm13_in, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot("AUC_PDK4_as_trained" %in% names(e13), "PDK4_dir" %in% names(e13),
            "PDK2_dir" %in% names(e13))
  m <- match(e13$comparison, out$comparison)
  e13$AUC_PDK4_as_trained[!is.na(m)] <- out$AUC_PDK4_as_trained[m[!is.na(m)]]
  e13$AUC_PDK4_oriented <- ifelse(is.na(m), NA, out$AUC_PDK4_oriented[m])
  names(e13)[names(e13) == "PDK4_dir"] <- "PDK4_dir_observed"
  e13$PDK4_dir_observed <- ifelse(is.na(m), NA, out$PDK4_dir_observed[m])
  e13$PDK4_median_case <- ifelse(is.na(m), NA, out$PDK4_median_case[m])
  e13$PDK4_median_ctrl <- ifelse(is.na(m), NA, out$PDK4_median_ctrl[m])
  e13$AUC_PDK2_as_trained[!is.na(m)] <- out$AUC_PDK2_as_trained[m[!is.na(m)]]
  e13$AUC_PDK2_oriented <- ifelse(is.na(m), NA, out$AUC_PDK2_oriented[m])
  names(e13)[names(e13) == "PDK2_dir"] <- "PDK2_dir_observed"
  e13$PDK2_dir_observed <- ifelse(is.na(m), NA, out$PDK2_dir_observed[m])
  e13$PDK2_median_case <- ifelse(is.na(m), NA, out$PDK2_median_case[m])
  e13$PDK2_median_ctrl <- ifelse(is.na(m), NA, out$PDK2_median_ctrl[m])
  e13$direction_reversed <- ifelse(is.na(m), "", out$direction_reversed[m])
  # column order: as_trained, oriented, dir_observed, medians per gene; flag after PDK2 medians
  want <- c("comparison", "n_case", "n_ctrl", "AUC_signature", "CI_lower", "CI_upper", "AUC_PDK4",
            "AUC_PDK4_as_trained", "AUC_PDK4_oriented", "PDK4_dir_observed",
            "PDK4_median_case", "PDK4_median_ctrl",
            "AUC_PDK2_as_trained", "AUC_PDK2_oriented", "PDK2_dir_observed",
            "PDK2_median_case", "PDK2_median_ctrl", "direction_reversed",
            "youden_sens", "youden_spec", "note")
  stopifnot(setequal(names(e13), want))
  e13 <- e13[, want]
  write.csv(e13, esm13_out, row.names = FALSE, na = "", quote = TRUE)
  cat("ESM_13 updated:", esm13_out, "\n")
  print(e13[, c("comparison", "AUC_PDK4_as_trained", "AUC_PDK4_oriented", "PDK4_dir_observed",
                "AUC_PDK2_as_trained", "AUC_PDK2_oriented", "PDK2_dir_observed", "direction_reversed")])
} else {
  cat("ESM_13 staging copy not found at", esm13_in, "- skipped Table S11 update\n")
}
cat("\nDIRECTION_DONE\n")

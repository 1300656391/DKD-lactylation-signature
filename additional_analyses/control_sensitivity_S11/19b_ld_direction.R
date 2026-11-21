# 19b_ld_direction.R ----
# Companion to 19_ld_sensitivity.R (Table S11):
# adds AS-TRAINED direction and medians for PDK4 and PDK2 in the GSE104948
# living-donor / hypertensive / nephrectomy comparisons, all platforms.
# The original S11 reported only a self-oriented (direction-invariant) PDK4 AUC;
# this script discloses the raw direction. Run inside the pipeline workdir
# (same prerequisites as 19_ld_sensitivity.R: data/GSE104948-*_expr_symbol.rds,
# results/hub_genes.txt, results/frozen/final_model_coefficients.csv).
# Output: results/revision/M5b_gse104948_direction.csv
suppressWarnings(suppressMessages({ library(data.table); library(pROC) }))

hub <- readLines("results/hub_genes.txt")
fm <- fread("results/frozen/final_model_coefficients.csv")
bh <- setNames(fm$canonical_deviance, fm$term)

zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
p_of <- function(e, g) as.numeric(zrows(as.matrix(e[g, , drop = FALSE]))[1, ])

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
  a4 <- as.numeric(auc(roc(y, p_of(es, "PDK4"), quiet = TRUE)))  # as-trained: >0.5 = PDK4 higher in DN
  a2 <- as.numeric(auc(roc(y, p_of(es, "PDK2"), quiet = TRUE)))
  i <<- i + 1
  res[[i]] <<- data.frame(
    comparison = label, n_case = length(case_idx), n_ctrl = length(ctrl_idx),
    AUC_PDK4_as_trained = round(a4, 3),
    PDK4_dir   = ifelse(median(es["PDK4", y == 1]) < median(es["PDK4", y == 0]), "down", "up"),
    PDK4_median_case = round(median(es["PDK4", y == 1]), 3),
    PDK4_median_ctrl = round(median(es["PDK4", y == 0]), 3),
    AUC_PDK2_as_trained = round(a2, 3),
    PDK2_dir   = ifelse(median(es["PDK2", y == 1]) < median(es["PDK2", y == 0]), "down", "up"),
    PDK2_median_case = round(median(es["PDK2", y == 1]), 3),
    PDK2_median_ctrl = round(median(es["PDK2", y == 0]), 3),
    note = note)
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
fwrite(out, "results/revision/M5b_gse104948_direction.csv")
print(out)
cat("\nDIRECTION_DONE\n")
# Decision rule: DN_vs_LD rows with PDK4_dir == "down" (and PDK2 same direction as main cohorts)
#   -> v3 wording may be narrowed to "preparation-associated" and cite updated Table S11;
#   PDK4_dir == "up" vs living donors -> keep "preparation- and control-source-associated"
#   (the shipped v3 abstract/Discussion wording is safe under either outcome).

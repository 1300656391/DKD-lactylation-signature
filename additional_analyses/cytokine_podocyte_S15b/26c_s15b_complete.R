# 26c_s15b_complete.R ----
# Round-2 fixes: complete the S15b evidence base.
#   - cytokine panel: replication-cohort within-DKD rows for IL6/IL1B (TNF absent) with all readouts
#   - p values for ALL readouts (score/PDK2/PDK4; podocyte PDK4/score) in every context
# Output: results/revision/M9c_s15b_complete.csv (tidy, pure CSV)
suppressWarnings(suppressMessages({ library(data.table) }))
set.seed(42)

hub <- readLines("results/hub_genes.txt")
fm <- fread("results/frozen/final_model_coefficients.csv")
bh <- setNames(fm$canonical_deviance, fm$term)
zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
sig_of <- function(e) {
  hg <- intersect(hub, rownames(e))
  Xv <- t(zrows(e[hg, , drop = FALSE]))
  as.vector(Xv %*% bh[hg] + bh[["(Intercept)"]])
}
rho <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))
rhop <- function(x, y) suppressWarnings(tryCatch(cor.test(x, y, method = "spearman")$p.value, error = function(e) NA))

out <- list(); i <- 0
addrow <- function(...) { i <<- i + 1; out[[i]] <<- data.frame(...) }
cyto_genes <- c("IL6", "TNF", "IL1B", "NLRP3")
podo_genes <- c("NPHS1", "NPHS2", "PODXL", "SYNPO")

run_cyto <- function(cohort, ctx, e, is_dn = NULL, note = "") {
  sc <- sig_of(e)
  idx <- if (is.null(is_dn)) seq_len(ncol(e)) else which(is_dn)
  for (g in cyto_genes) {
    if (!g %in% rownames(e)) {
      addrow(analysis = "cytokine_panel", cohort = cohort, context = ctx, feature = g,
             readout = "score", value = NA, p = NA, n = length(idx), note = "gene absent from processed expression matrix")
      next
    }
    x <- as.numeric(e[g, idx])
    sc_i <- if (is.null(is_dn)) sc else sc[idx]
    for (rd in c("score", "PDK2", "PDK4")) {
      y <- if (rd == "score") sc_i else as.numeric(e[rd, idx])
      addrow(analysis = "cytokine_panel", cohort = cohort, context = ctx, feature = g,
             readout = rd, value = round(rho(y, x), 3), p = signif(rhop(y, x), 3),
             n = length(idx), note = note)
    }
  }
}

run_podo <- function(cohort, ctx, e, is_dn = NULL, note = "") {
  avail <- podo_genes[podo_genes %in% rownames(e)]
  idx <- if (is.null(is_dn)) seq_len(ncol(e)) else which(is_dn)
  pm_all <- colMeans(zrows(e[avail, , drop = FALSE]))
  pm <- pm_all[idx]
  sc <- sig_of(e)
  sc_i <- if (is.null(is_dn)) sc else sc[idx]
  for (rd in c("score", "PDK4")) {
    y <- if (rd == "score") sc_i else as.numeric(e[rd, idx])
    nt <- note
    if (rd == "score") nt <- paste0(note, ifelse(note == "", "", "; "),
                                    "disease polarity dominates; not localization evidence")
    if (rd == "PDK4" && is.null(is_dn)) nt <- paste0(note, ifelse(note == "", "", "; "),
                                    "partially reflects disease polarity; DN-only subset shown separately")
    addrow(analysis = "podocyte_module", cohort = cohort, context = ctx,
           feature = paste(avail, collapse = ";"), readout = rd,
           value = round(rho(y, pm), 3), p = signif(rhop(y, pm), 3),
           n = length(idx), note = nt)
  }
}

# training
te <- readRDS("data/GSE96804_expr_symbol.rds")
ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
dn <- ph$group == "DN"
run_cyto("GSE96804", "all-sample", as.matrix(te))
run_cyto("GSE96804", "DN-only (subset-z)", as.matrix(te[, dn]), note = "subset-internal z")
run_podo("GSE96804", "all-sample", as.matrix(te))
run_podo("GSE96804", "DN-only", as.matrix(te[, dn]), note = "polarity-minimal read-out")

# replication
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^\"|\"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])),
                  stringsAsFactors = FALSE)
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
rownames(e14) <- toupper(rownames(e14))
keep <- p14$group %in% c("Advanced_DN", "Early_DN", "Control")
e14k <- as.matrix(e14[, keep])
is_dn14 <- p14$group[keep] != "Control"
run_cyto("GSE142025", "all-sample", e14k)
run_cyto("GSE142025", "within-DKD (full-cohort scoring)", e14k, is_dn14,
         note = "score from full-cohort standardization; subset used for correlation")
run_podo("GSE142025", "all-sample", e14k)

# glomerular validation cohorts
for (acc in c("GSE30528", "GSE30529")) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  run_podo(acc, "all-sample", as.matrix(e))
}
for (p in c("GPL22945", "GPL24120")) {
  ep <- readRDS(paste0("data/GSE104948-", p, "_expr_symbol.rds"))
  mlp <- readLines(paste0("data/GSE104948-", p, "_meta.txt"), warn = FALSE)
  ttlp <- trimws(gsub('^\"|\"$', "", strsplit(sub("^[^\t]*\t", "", mlp[grepl("^!Sample_title", mlp)][1]), "\t", fixed = TRUE)[[1]]))
  gp <- sub("^[A-Z][0-9]+-Glom-", "", sub("[0-9]+$", "", ttlp))
  idx <- which(gp %in% c("DN", "TN", "HT", "LD_"))
  run_podo(paste0("GSE104948_", p), "all-sample; all annotated tissue states", as.matrix(ep[, idx]),
           note = "DN/TN/HT/LD pooled (as in 2.6 sensitivity analyses)")
}

outdf <- rbindlist(out, fill = TRUE)
fwrite(outdf, "results/revision/M9c_s15b_complete.csv")
print(outdf[analysis == "cytokine_panel" & context != "all-sample" & cohort == "GSE142025"])
print(outdf[analysis == "podocyte_module"])
cat("\nS15B_COMPLETE_DONE\n")

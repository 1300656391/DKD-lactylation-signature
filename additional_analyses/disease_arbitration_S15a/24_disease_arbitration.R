# 24_disease_arbitration.R ----
# Science-fresh-review arbitration analyses (existing data only, frozen conventions):
#   A2-i   HT vs living-donor glomeruli within GPL24120 (generic-injury alternative)
#   A2-ii  Signature score vs fibrosis/tubular proxies (training GSE96804 + GSE142025)
#   A5     Dependent-correlation test for PDK2/PDK4-macrophage asymmetry in GSE142025
#          (bootstrap, 10k, seed 42) + control-subgroup bootstrap CI
#   (rec)  Partial Spearman PDK2-macrophage in GSE142025 (replication of training analysis)
suppressWarnings(suppressMessages({
  library(data.table); library(pROC); library(GSVA)
}))
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
mk_proxies <- function(e) {
  g <- function(x) if (x %in% rownames(e)) as.numeric(e[x, ]) else NA
  tub <- (scale(g("LRP2")) + scale(g("CUBN"))) / 2
  fib <- (scale(g("COL1A1")) + scale(g("COL3A1"))) / 2
  list(tub = tub, fib = fib)
}

out <- list(); i <- 0
addrow <- function(...) { i <<- i + 1; out[[i]] <<- data.frame(...) }

# ---------- A2-i: HT vs LD within GPL24120 ----------
grp_from_title <- function(p) {
  ml <- readLines(paste0("data/GSE104948-", p, "_meta.txt"), warn = FALSE)
  gsm <- trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", ml[grepl("^!Sample_geo_accession", ml)][1]), "\t", fixed = TRUE)[[1]]))
  ttl <- trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", ml[grepl("^!Sample_title", ml)][1]), "\t", fixed = TRUE)[[1]]))
  data.frame(gsm = gsm, grp = sub("^[A-Z][0-9]+-Glom-", "", sub("[0-9]+$", "", ttl)), stringsAsFactors = FALSE)
}
g <- grp_from_title("GPL24120"); e <- readRDS("data/GSE104948-GPL24120_expr_symbol.rds")
stopifnot(all(colnames(e) == g$gsm))
ht <- which(g$grp == "HT"); ld <- which(g$grp == "LD_"); dn <- which(g$grp == "DN"); tn <- which(g$grp == "TN")
cat("GPL24120: DN", length(dn), "HT", length(ht), "LD", length(ld), "TN", length(tn), "\n")
es <- as.matrix(e[, c(ht, ld), drop = FALSE])
y <- c(rep(1, length(ht)), rep(0, length(ld)))
r <- roc(y, sig_of(es), quiet = TRUE)
cat(sprintf("HT_vs_LD_GPL24120: AUC=%.3f CI=(%.3f-%.3f) n=%d vs %d\n",
            as.numeric(auc(r)), ci.auc(r)[1], ci.auc(r)[3], length(ht), length(ld)))
addrow(analysis = "HT_vs_LD_GPL24120", estimate = round(as.numeric(auc(r)), 3),
       lo = round(ci.auc(r)[1], 3), hi = round(ci.auc(r)[3], 3),
       n1 = length(ht), n2 = length(ld), note = "frozen signature, subset-internal z")
# descriptive medians per group (z within full platform)
sc_all <- sig_of(as.matrix(e))
for (gp in c("DN", "HT", "LD_", "TN")) {
  idx <- which(g$grp == gp)
  if (length(idx)) addrow(analysis = paste0("median_score_", gp, "_GPL24120"),
                          estimate = round(median(sc_all[idx]), 3), lo = NA, hi = NA,
                          n1 = length(idx), n2 = NA, note = "descriptive")
}
# PDK4 single gene for context (direction-oriented)
p4 <- as.numeric(zrows(as.matrix(e["PDK4", c(ht, ld), drop = FALSE]))[1, ])
r4 <- roc(y, p4, quiet = TRUE)
addrow(analysis = "HT_vs_LD_GPL24120_PDK4", estimate = round(as.numeric(auc(r4)), 3),
       lo = round(ci.auc(r4)[1], 3), hi = round(ci.auc(r4)[3], 3),
       n1 = length(ht), n2 = length(ld), note = "PDK4 single gene")

# ---------- A2-ii: signature vs fibrosis/tubular proxies ----------
# training
te <- readRDS("data/GSE96804_expr_symbol.rds")
sc_tr <- sig_of(as.matrix(te))
pr <- mk_proxies(te)
for (v in c("tub", "fib")) {
  rho <- suppressWarnings(cor(sc_tr, pr[[v]], method = "spearman"))
  pv <- suppressWarnings(cor.test(sc_tr, pr[[v]], method = "spearman")$p.value)
  addrow(analysis = paste0("score_vs_", v, "_GSE96804"), estimate = round(rho, 3),
         lo = NA, hi = NA, n1 = ncol(te), n2 = NA, note = sprintf("P=%.4g", pv))
  cat(sprintf("training score~%s: rho=%.3f P=%.3g\n", v, rho, pv))
}
# within-DN only (severity axis beyond case status)
ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
is_dn_tr <- ph$group == "DN"
for (v in c("tub", "fib")) {
  rho <- suppressWarnings(cor(sc_tr[is_dn_tr], pr[[v]][is_dn_tr], method = "spearman"))
  pv <- suppressWarnings(cor.test(sc_tr[is_dn_tr], pr[[v]][is_dn_tr], method = "spearman")$p.value)
  addrow(analysis = paste0("score_vs_", v, "_GSE96804_DNonly"), estimate = round(rho, 3),
         lo = NA, hi = NA, n1 = sum(is_dn_tr), n2 = NA, note = sprintf("P=%.4g", pv))
  cat(sprintf("training DN-only score~%s: rho=%.3f P=%.3g (n=%d)\n", v, rho, pv, sum(is_dn_tr)))
}
# GSE142025
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])),
                  stringsAsFactors = FALSE)
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
rownames(e14) <- toupper(rownames(e14))
keep <- p14$group %in% c("Advanced_DN", "Early_DN", "Control")
e14k <- as.matrix(e14[, keep])
sc14 <- sig_of(e14k)
pr14 <- mk_proxies(e14k)
for (v in c("tub", "fib")) {
  if (!any(is.na(pr14[[v]]))) {
    rho <- suppressWarnings(cor(sc14, pr14[[v]], method = "spearman"))
    pv <- suppressWarnings(cor.test(sc14, pr14[[v]], method = "spearman")$p.value)
    addrow(analysis = paste0("score_vs_", v, "_GSE142025"), estimate = round(rho, 3),
           lo = NA, hi = NA, n1 = ncol(e14k), n2 = NA, note = sprintf("P=%.4g", pv))
    cat(sprintf("GSE142025 score~%s: rho=%.3f P=%.3g\n", v, rho, pv))
  }
}

# ---------- macrophage scores for GSE142025 (as in 21) ----------
ev <- new.env(); load("data/signature_collection.rda", envir = ev)
scoll <- get("signature_collection", envir = ev)
cn <- grep("_Bindea_et_al$", names(scoll), value = TRUE)
cn <- setdiff(cn, c("SW480_cancer_cells_Bindea_et_al"))
gsets <- scoll[cn]; names(gsets) <- sub("_Bindea_et_al$", "", cn)
bindea_cells <- c("B_cells","T_cells","T_helper_cells","Tcm","Tem","Th1_cells","Th2_cells","TFH",
  "CD8_T_cells","Cytotoxic_cells","NK_cells","NK_CD56dim_cells","DC","iDC","Eosinophils",
  "Macrophages","Mast_cells","Neutrophils")
sig18 <- gsets[intersect(bindea_cells, names(gsets))]
sig18 <- lapply(sig18, function(x) unique(toupper(x)))
par <- ssgseaParam(exprData = e14k, geneSets = sig18, normalize = TRUE)
ss <- gsva(par, verbose = FALSE)
mac <- as.numeric(ss["Macrophages", ])
pdk2 <- as.numeric(e14k["PDK2", ]); pdk4 <- as.numeric(e14k["PDK4", ])
grp14 <- p14$group[keep]

r12 <- suppressWarnings(cor(pdk2, mac, method = "spearman"))
r13 <- suppressWarnings(cor(pdk4, mac, method = "spearman"))
r23 <- suppressWarnings(cor(pdk2, pdk4, method = "spearman"))
cat(sprintf("\nGSE142025: rho(PDK2,mac)=%.3f rho(PDK4,mac)=%.3f rho(PDK2,PDK4)=%.3f\n", r12, r13, r23))

# bootstrap difference (10k)
n <- length(mac)
B <- 10000
diffs <- replicate(B, {
  idx <- sample.int(n, n, replace = TRUE)
  d2 <- suppressWarnings(cor(pdk2[idx], mac[idx], method = "spearman"))
  d3 <- suppressWarnings(cor(pdk4[idx], mac[idx], method = "spearman"))
  if (is.na(d2) || is.na(d3)) NA else d2 - d3
})
diffs <- diffs[!is.na(diffs)]
pboot <- 2 * min(mean(diffs <= 0), mean(diffs >= 0))
addrow(analysis = "PDK2mac_minus_PDK4mac_GSE142025", estimate = round(r12 - r13, 3),
       lo = round(quantile(diffs, .025), 3), hi = round(quantile(diffs, .975), 3),
       n1 = n, n2 = B, note = sprintf("bootstrap two-sided P=%.4g", pboot))
cat(sprintf("bootstrap diff: est=%.3f CI=(%.3f-%.3f) P=%.3g\n",
            r12 - r13, quantile(diffs, .025), quantile(diffs, .975), pboot))

# Steiger z on ranks (dependent overlapping correlations) as cross-check
rr <- data.frame(apply(cbind(pdk2, pdk4, mac), 2, rank))
z12 <- atanh(cor(rr$pdk2, rr$mac)); z13 <- atanh(cor(rr$pdk4, rr$mac)); r23r <- cor(rr$pdk2, rr$pdk4)
detR <- 1 - cor(rr$pdk2, rr$mac)^2 - cor(rr$pdk4, rr$mac)^2 - r23r^2 + 2 * cor(rr$pdk2, rr$mac) * cor(rr$pdk4, rr$mac) * r23r
se <- sqrt((n - 3) / ((n - 1) * (1 + r23r)) * detR) / (n - 3)  # Steiger 1980 SE on z-scale ~
zstat <- (z12 - z13) / sqrt(2 * (1 - r23r) / (n - 3))
psteiger <- 2 * pnorm(-abs(zstat))
addrow(analysis = "PDK2mac_vs_PDK4mac_Steiger_GSE142025", estimate = round(zstat, 2),
       lo = NA, hi = NA, n1 = n, n2 = NA, note = sprintf("P=%.4g", psteiger))
cat(sprintf("Steiger z=%.2f P=%.3g\n", zstat, psteiger))

# control-subgroup PDK2-mac rho + bootstrap CI (n=9)
ctl <- grp14 == "Control"
rctl <- suppressWarnings(cor(pdk2[ctl], mac[ctl], method = "spearman"))
bc <- replicate(B, {
  idx <- sample.int(sum(ctl), sum(ctl), replace = TRUE)
  suppressWarnings(cor(pdk2[ctl][idx], mac[ctl][idx], method = "spearman"))
})
bc <- bc[!is.na(bc)]
addrow(analysis = "PDK2mac_controlsub_GSE142025", estimate = round(rctl, 3),
       lo = round(quantile(bc, .025), 3), hi = round(quantile(bc, .975), 3),
       n1 = sum(ctl), n2 = B, note = "bootstrap CI")
cat(sprintf("control-subgroup rho=%.3f CI=(%.3f-%.3f)\n", rctl, quantile(bc, .025), quantile(bc, .975)))

# ---------- partial Spearman PDK2-mac in GSE142025 ----------
df14 <- data.frame(PDK2 = pdk2, Macrophage = mac, tub = pr14$tub, fib = pr14$fib)
rk <- rank
partial_spearman <- function(x, y, ...) {
  d2 <- data.frame(x = rk(x), y = rk(y), ...)
  r1 <- residuals(lm(x ~ ., d2[, setdiff(names(d2), "y"), drop = FALSE]))
  r2 <- residuals(lm(y ~ ., d2[, setdiff(names(d2), "x"), drop = FALSE]))
  cor(r1, r2, method = "spearman")
}
if (!any(is.na(df14$tub))) {
  rho_raw <- cor(df14$PDK2, df14$Macrophage, method = "spearman")
  rho_adj <- partial_spearman(df14$PDK2, df14$Macrophage, tub = rk(df14$tub), fib = rk(df14$fib))
  addrow(analysis = "PDK2mac_partial_GSE142025", estimate = round(rho_adj, 3), lo = NA, hi = NA,
         n1 = n, n2 = NA, note = sprintf("raw=%.3f (training: -0.585 -> -0.305)", rho_raw))
  cat(sprintf("GSE142025 partial: raw=%.3f -> adjusted=%.3f\n", rho_raw, rho_adj))
}

outdf <- rbindlist(out)
fwrite(outdf, "results/revision/M8_disease_arbitration.csv")
print(outdf)
cat("\nARBITRATION_DONE\n")

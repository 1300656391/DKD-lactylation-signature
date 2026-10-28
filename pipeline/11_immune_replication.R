# 11_immune_replication.R ----
# Replication of the immune-microenvironment analyses:
#   A. 18-cell ssGSEA in the independent RNA-seq cohort GSE142025 (DN vs
#      control, BH within 18 tests) and PDK2/PDK4 x myeloid correlations;
#   B. immune-compartment proportions in the annotated snRNA-seq object;
#   C. partial Spearman analysis in the training cohort (PDK-macrophage
#      association controlled for tubular-content and fibrosis proxies) and
#      macrophage state-module correlations;
#   D. MCP-counter cross-method sensitivity (training + replication).
# Outputs:
#   results/replication_summary.csv (pipeline marker),
#   results/replication/gse142025_immune.csv (S13a), gse142025_pdk_myeloid.csv (S13b),
#     scrna_immune_proportion.csv (S13c), partial_spearman_macrostate.csv (S14),
#   results/supplement/TableS21_mcp_counter.csv

suppressWarnings(suppressMessages({
  library(data.table); library(GSVA); library(Seurat); library(MCPcounter)
}))
set.seed(42)
dir.create("results/replication", showWarnings = FALSE, recursive = TRUE)
dir.create("results/supplement", showWarnings = FALSE, recursive = TRUE)

# ---------- Shared: Bindea 18-cell signature ----------
e <- new.env(); load("data/signature_collection.rda", envir = e)
sc <- get("signature_collection", envir = e)
cn <- grep("_Bindea_et_al$", names(sc), value = TRUE)
cn <- setdiff(cn, "SW480_cancer_cells_Bindea_et_al")
gsets <- sc[cn]; names(gsets) <- sub("_Bindea_et_al$", "", cn)
bindea_cells <- c("B_cells","T_cells","T_helper_cells","Tcm","Tem","Th1_cells","Th2_cells","TFH",
                  "CD8_T_cells","Cytotoxic_cells","NK_cells","NK_CD56dim_cells","DC","iDC",
                  "Eosinophils","Macrophages","Mast_cells","Neutrophils")
sig18 <- gsets[intersect(bindea_cells, names(gsets))]
sig18 <- lapply(sig18, function(g) unique(toupper(g)))
stopifnot(length(sig18) == 18)

# ---------- A. GSE142025 replication ----------
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])))
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
rownames(e14) <- toupper(rownames(e14))
keep <- p14$group %in% c("Advanced_DN", "Early_DN", "Control")
e14k <- as.matrix(e14[, keep])
grp <- p14$group[keep]
is_dn <- grp != "Control"
stage <- ifelse(grp == "Early_DN", "Early", ifelse(grp == "Advanced_DN", "Advanced", "Control"))

par18 <- ssgseaParam(exprData = e14k, geneSets = sig18, normalize = TRUE)
ss <- as.data.frame(t(gsva(par18, verbose = FALSE)))
ss$group <- grp; ss$stage <- stage

res <- do.call(rbind, lapply(bindea_cells, function(cc) {
  x <- ss[[cc]]
  data.frame(cell = cc, median_DN = median(x[is_dn]), median_Ctrl = median(x[!is_dn]),
             p = wilcox.test(x[is_dn], x[!is_dn])$p.value)
}))
res$p.adj_BH18 <- p.adjust(res$p, "BH")
res$direction <- ifelse(res$median_DN > res$median_Ctrl, "up", "down")
fwrite(res, "results/replication/gse142025_immune.csv")
file.copy("results/replication/gse142025_immune.csv",
          "results/supplement/TableS13a_gse142025_immune.csv", overwrite = TRUE)
cat("GSE142025 DN vs Control (18 cells):", sum(res$p.adj_BH18 < 0.05), "significant\n")

mye <- c("Macrophages", "Neutrophils", "DC", "iDC", "Eosinophils")
cor_res <- do.call(rbind, lapply(c("PDK2", "PDK4"), function(g) {
  ge <- as.numeric(e14k[g, ])
  do.call(rbind, lapply(mye, function(cc) {
    data.frame(gene = g, cell = cc,
               rho = round(suppressWarnings(cor(ge, ss[[cc]], method = "spearman")), 3))
  }))
}))

# Within-group and per-stage PDK x macrophage correlations (GSE142025):
# all samples / DN only / controls only / early / advanced
mac <- ss$Macrophages
pdk2 <- as.numeric(e14k["PDK2", ]); pdk4 <- as.numeric(e14k["PDK4", ])
subsets <- list(
  list(scope = "All (DN+Ctrl n=36)", idx = seq_along(mac)),
  list(scope = "DN only (n=27)",      idx = which(grp != "Control")),
  list(scope = "Control only (n=9)",  idx = which(grp == "Control")),
  list(scope = "Early DN (n=6)",      idx = which(grp == "Early_DN")),
  list(scope = "Advanced DN (n=21)",  idx = which(grp == "Advanced_DN")))
within <- do.call(rbind, lapply(subsets, function(s) {
  i <- s$idx
  data.frame(gene_or_subset = "PDK2+PDK4", cell_or_scope = s$scope,
             PDK2_rho = round(suppressWarnings(cor(pdk2[i], mac[i], method = "spearman")), 3),
             PDK4_rho = round(suppressWarnings(cor(pdk4[i], mac[i], method = "spearman")), 3),
             note = "within-group Spearman")
}))

# S13b: all-sample myeloid correlations + within-group macrophage rows
s13b <- rbind(
  data.frame(gene_or_subset = cor_res$gene[cor_res$gene == "PDK2"],
             cell_or_scope = cor_res$cell[cor_res$gene == "PDK2"],
             PDK2_rho = cor_res$rho[cor_res$gene == "PDK2"], PDK4_rho = NA_real_,
             note = "all-sample Spearman"),
  data.frame(gene_or_subset = cor_res$gene[cor_res$gene == "PDK4"],
             cell_or_scope = cor_res$cell[cor_res$gene == "PDK4"],
             PDK2_rho = NA_real_,
             PDK4_rho = cor_res$rho[cor_res$gene == "PDK4"],
             note = "all-sample Spearman"),
  within)
fwrite(s13b, "results/replication/gse142025_pdk_myeloid.csv")
file.copy("results/replication/gse142025_pdk_myeloid.csv",
          "results/supplement/TableS13b_pdk_myeloid.csv", overwrite = TRUE)
fwrite(within, "results/replication/within_group_macro_corr.csv")

# Per-stage macrophage medians
stg <- data.frame(group = c("Control", "Early DN", "Advanced DN"),
                  n = c(sum(grp == "Control"), sum(grp == "Early_DN"),
                        sum(grp == "Advanced_DN")),
                  median_macrophage = c(round(median(mac[grp == "Control"]), 3),
                                        round(median(mac[grp == "Early_DN"]), 3),
                                        round(median(mac[grp == "Advanced_DN"]), 3)))
fwrite(stg, "results/replication/perstage_macrophage.csv")

w <- suppressWarnings(wilcox.test(ss$Macrophages[stage == "Early"], ss$Macrophages[stage == "Advanced"]))
cat("Macrophages Early vs Advanced: p =", signif(w$p.value, 3), "\n")

# ---------- B. snRNA immune-compartment proportions + per-stage medians (S13c) ----------
merged <- readRDS("results/snrna/seurat_annotated.rds")
ct <- as.character(merged$celltype); gp <- as.character(merged$group)
dn <- sum(gp == "DN"); ctl <- sum(gp == "Control")
props <- do.call(rbind, lapply(c("Macro", "T/NK", "B"), function(t) {
  data.frame(type_or_group = t,
             pct_DN_or_median = round(100 * sum(ct == t & gp == "DN") / dn, 2),
             pct_Ctrl_or_n = round(100 * sum(ct == t & gp == "Control") / ctl, 2),
             source = "GSE131882 snRNA immune proportions (pct of total nuclei)")
}))
s13c <- rbind(props, data.frame(
  type_or_group = stg$group, pct_DN_or_median = stg$median_macrophage,
  pct_Ctrl_or_n = stg$n,
  source = "GSE142025 per-stage macrophage medians (median, n)"))
fwrite(s13c, "results/replication/scrna_immune_proportion.csv")
file.copy("results/replication/scrna_immune_proportion.csv",
          "results/supplement/TableS13c_scrna_immune_proportion.csv", overwrite = TRUE)
print(s13c)

# ---------- C. Partial Spearman + macrophage state modules (training) ----------
train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
immune_ss <- as.data.frame(fread("results/immune_ssgsea.csv"))
rownames(immune_ss) <- immune_ss$sample; immune_ss$sample <- NULL
ov <- intersect(colnames(train_expr), rownames(immune_ss))
gz <- function(g) as.numeric(train_expr[g, ov])
zsc <- function(v) (v - mean(v)) / sd(v)
df <- data.frame(
  PDK2 = gz("PDK2"), PDK4 = gz("PDK4"),
  Macrophage = immune_ss[ov, "Macrophages"],
  tub  = (zsc(gz("LRP2")) + zsc(gz("CUBN"))) / 2,
  fib  = (zsc(gz("COL1A1")) + zsc(gz("COL3A1"))) / 2,
  homeo = (zsc(gz("ARG1")) + zsc(gz("MRC1")) + zsc(gz("CD163"))) / 3,
  LAM = (zsc(gz("TREM2")) + zsc(gz("APOE")) + zsc(gz("LPL"))) / 3)

rk <- rank
partial_spearman <- function(x, y, ...) {
  df2 <- data.frame(x = rk(x), y = rk(y), ...)
  r1 <- residuals(lm(x ~ ., df2[, setdiff(names(df2), "y"), drop = FALSE]))
  r2 <- residuals(lm(y ~ ., df2[, setdiff(names(df2), "x"), drop = FALSE]))
  list(rho = suppressWarnings(cor(r1, r2, method = "spearman")),
       p = suppressWarnings(cor.test(r1, r2, method = "spearman")$p.value))
}
spear <- function(x, y) list(rho = suppressWarnings(cor(x, y, method = "spearman")),
                             p = suppressWarnings(cor.test(x, y, method = "spearman")$p.value))

ps <- list(
  PDK2_Macro_raw = spear(df$PDK2, df$Macrophage),
  PDK2_Macro_ctrl_tub = partial_spearman(df$PDK2, df$Macrophage, tub = rk(df$tub)),
  PDK2_Macro_ctrl_fib = partial_spearman(df$PDK2, df$Macrophage, fib = rk(df$fib)),
  PDK2_Macro_ctrl_both = partial_spearman(df$PDK2, df$Macrophage, tub = rk(df$tub), fib = rk(df$fib)),
  PDK4_Macro_raw = spear(df$PDK4, df$Macrophage),
  PDK4_Macro_ctrl_both = partial_spearman(df$PDK4, df$Macrophage, tub = rk(df$tub), fib = rk(df$fib)),
  PDK2_homeo = spear(df$PDK2, df$homeo), PDK2_LAM = spear(df$PDK2, df$LAM),
  PDK4_homeo = spear(df$PDK4, df$homeo), PDK4_LAM = spear(df$PDK4, df$LAM),
  Macro_homeo = spear(df$Macrophage, df$homeo), Macro_LAM = spear(df$Macrophage, df$LAM))
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 4))
s14 <- data.frame(analysis = names(ps),
                  spearman_rho = round(vapply(ps, `[[`, 0, "rho"), 3),
                  p_value = fmt_p(vapply(ps, `[[`, 0, "p")))
fwrite(s14, "results/replication/partial_spearman_macrostate.csv")
file.copy("results/replication/partial_spearman_macrostate.csv",
          "results/supplement/TableS14_partial_spearman.csv", overwrite = TRUE)
print(s14)

# ---------- D. MCP-counter cross-method sensitivity ----------
# The package's estimate() wrapper errors on this R version; the exported
# appendSignatures workhorse reproduces it (verified equivalent).
mcp_genes <- read.table("data/mcp_genes.tsv", sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE, colClasses = "character",
                        check.names = FALSE)
for (j in seq_len(ncol(mcp_genes))) mcp_genes[[j]] <- trimws(gsub('^"|"$', '', mcp_genes[[j]]))
mcp_estimate <- function(expr) {
  f <- mcp_genes[mcp_genes$`HUGO symbols` %in% rownames(expr), ]
  features <- split(f$`HUGO symbols`, f$`Cell population`)
  t(appendSignatures(expr, features))
}
run_cohort <- function(tag, expr, is_case) {
  mcp <- mcp_estimate(expr)
  rows <- lapply(rownames(mcp), function(p) {
    x <- mcp[p, ]
    wt <- suppressWarnings(wilcox.test(x[is_case], x[!is_case]))
    data.frame(cohort = tag, population = p,
               median_DN = round(median(x[is_case]), 3),
               median_ctrl = round(median(x[!is_case]), 3),
               W = unname(wt$statistic), p = wt$p.value,
               rho_PDK2 = suppressWarnings(cor(x, as.numeric(expr["PDK2", ]), method = "spearman")),
               p_PDK2 = suppressWarnings(cor.test(x, as.numeric(expr["PDK2", ]), method = "spearman")$p.value),
               rho_PDK4 = suppressWarnings(cor(x, as.numeric(expr["PDK4", ]), method = "spearman")),
               p_PDK4 = suppressWarnings(cor.test(x, as.numeric(expr["PDK4", ]), method = "spearman")$p.value),
               rho_PDK2_withinDN = suppressWarnings(cor(x[is_case], as.numeric(expr["PDK2", is_case]), method = "spearman")),
               p_PDK2_withinDN = suppressWarnings(cor.test(x[is_case], as.numeric(expr["PDK2", is_case]), method = "spearman")$p.value))
  })
  d <- rbindlist(rows)
  d[, p_BH := p.adjust(p, "BH")]
  d[, p_PDK2_BH := p.adjust(p_PDK2, "BH")]
  d[, p_PDK4_BH := p.adjust(p_PDK4, "BH")]
  d
}
train_mcp <- run_cohort("GSE96804", train_expr, train_ph$group == "DN")
rep_mcp <- run_cohort("GSE142025", e14k, is_dn)
mcp_out <- rbind(train_mcp, rep_mcp)
setcolorder(mcp_out, c("cohort", "population", "median_DN", "median_ctrl", "W", "p",
                       "rho_PDK2", "p_PDK2", "rho_PDK4", "p_PDK4",
                       "p_BH", "p_PDK2_BH", "p_PDK4_BH",
                       "rho_PDK2_withinDN", "p_PDK2_withinDN"))
setorder(mcp_out, cohort, p)
fwrite(mcp_out, "results/supplement/TableS21_mcp_counter.csv")

# ---------- Summary marker ----------
fwrite(data.frame(
  analysis = c("gse142025_immune_18cells", "gse142025_pdk_myeloid", "scrna_immune_proportion",
               "partial_spearman_macrostate", "mcp_counter"),
  significant_or_key = c(sum(res$p.adj_BH18 < 0.05),
                         nrow(s13b), nrow(s13c), nrow(s14), nrow(mcp_out)),
  output = c("results/replication/gse142025_immune.csv",
             "results/replication/gse142025_pdk_myeloid.csv",
             "results/replication/scrna_immune_proportion.csv",
             "results/replication/partial_spearman_macrostate.csv",
             "results/supplement/TableS21_mcp_counter.csv")),
  "results/replication_summary.csv")

cat("\nIMMUNE_REPLICATION_DONE\n")

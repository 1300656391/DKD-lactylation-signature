# 16e_fig6_fig7b_regen.R ----
# R12-R2/R3: Fig6 re-render at 12.5in (10.5in caused coord_equal panel overflow:
# clipped E title/legend) and remove the contradictory k-annotation;
# Fig7B re-render with a short title (long title clipped at canvas right edge).
suppressWarnings(suppressMessages({
  library(data.table); library(pROC); library(ggplot2); library(patchwork); library(rms)
}))
set.seed(42)
fig_dir <- "results/figures"
png_win <- function(filename, width, height, units = 'in', res = 300, ...)
  grDevices::png(filename, width = width, height = height, units = units, res = res, type = 'windows', ...)
save_fig <- function(p, name, w, h, dpi = 300)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = dpi, device = png_win)

train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph   <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
allg <- fread("results/A_all_genes_limma.csv")
hub <- readLines("results/hub_genes.txt")
zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
cand <- fread("results/A_candidates.csv")$gene
genes <- intersect(cand, rownames(train_expr))
Xdf <- data.frame(t(zrows(train_expr[genes, , drop = FALSE])), check.names = FALSE)
y <- as.integer(train_ph$group == "DN")

# ---- Fig6 (12.5in, no k-annotation) ----
fl <- fread("results/frozen/lasso_curve.csv")
p1 <- ggplot(fl, aes(lambda, cvm)) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), color = "grey70") +
  geom_point(color = "#B2182B", size = .9) + scale_x_log10() +
  geom_vline(xintercept = fl$lambda[fl$is_1se], linetype = 2) +
  theme_bw(base_size = 11) + ylab("Binomial deviance") + xlab("lambda") +
  ggtitle("LASSO CV (lambda.1se)") + labs(tag = "A")
fs <- fread("results/frozen/svm_rfe_curve.csv")
p2 <- ggplot(fs, aes(k, accuracy)) +
  geom_line(color = "#2166AC") + geom_point(size = 2) +
  geom_vline(xintercept = fs$k[fs$is_best], linetype = 2) +
  theme_bw(base_size = 11) + xlab("Number of features") + ylab("CV accuracy") +
  ggtitle(sprintf("SVM-RFE (optimum k = %d)", fs$k[fs$is_best])) + labs(tag = "B") +
  ylim(min(fs$accuracy) - .02, 1)
fx <- fread("results/frozen/xgb_importance.csv")
p3 <- ggplot(head(fx, 10), aes(reorder(Feature, Gain_pct), Gain_pct)) +
  geom_col(fill = "#4CAF50", width = .7) + coord_flip() +
  theme_bw(base_size = 11) + xlab(NULL) + ylab("Gain (%)") +
  ggtitle("XGBoost importance") + labs(tag = "C")
hubexp <- do.call(rbind, lapply(hub, function(g)
  data.frame(gene = g, expr = as.numeric(train_expr[g, ]), group = train_ph$group)))
allgk <- allg[gene %in% hub]
pv <- setNames(allgk$adj.P.Val, allgk$gene)
star <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "ns")))
hubexp$gene <- factor(hubexp$gene, levels = hub)
ann <- data.frame(gene = factor(hub, levels = hub), lab = star(pv[hub]),
                  expr = max(hubexp$expr) * 1.06, group = "DN")
p4 <- ggplot(hubexp, aes(gene, expr, fill = group)) +
  geom_boxplot(outlier.size = .5, linewidth = .3) +
  geom_text(data = ann, aes(x = gene, y = expr, label = lab), size = 3.4, inherit.aes = FALSE) +
  scale_fill_manual(values = c(Control = "#2166AC", DN = "#B2182B")) +
  theme_bw(base_size = 11) + xlab(NULL) + ylab("Expression (log2)") +
  ggtitle("Hub gene expression (training)") + labs(tag = "D") +
  theme(legend.position = "bottom") +
  annotate("text", x = .5, y = max(hubexp$expr) * 1.12, hjust = 0,
           label = "ns / * / ** / *** vs adj.P 0.05/0.01/0.001", size = 2.8)
rocs <- lapply(hub, function(g) roc(y, Xdf[[g]], quiet = TRUE))
names(rocs) <- hub
p5 <- ggroc(rocs, linewidth = .7) +
  geom_abline(slope = 1, intercept = 1, linetype = 3) +
  scale_color_brewer(palette = "Set1") + theme_bw(base_size = 11) +
  coord_equal() + ggtitle("Single-gene ROC (training)") + labs(color = NULL, tag = "E") +
  theme(legend.position = "bottom", legend.box.margin = margin(0, 30, 0, 0, "pt"))
save_fig((p1 | p2 | p3) / (p4 | p5), "Fig6_ml_screening_rev", 12.5, 7.0)
cat("FIG6_DONE 12.5in-noannot\n")

# ---- Fig7B (short title, balanced) ----
fm <- fread("results/frozen/final_model_coefficients.csv")
bh <- setNames(fm$canonical_deviance, fm$term)
sc <- function(e) {
  hg <- intersect(hub, rownames(e))
  Xv <- t(zrows(e[hg, , drop = FALSE]))
  as.vector(Xv %*% bh[hg] + bh[["(Intercept)"]])
}
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])))
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
k3 <- p14$group %in% c("Advanced_DN", "Early_DN", "Control")
y3 <- as.integer(p14$group[k3] != "Control")
s3 <- sc(e14[, k3, drop = FALSE])
g1 <- glm(y3 ~ s3, family = binomial)
p_hat <- fitted(g1)
png(file.path(fig_dir, "Fig7B_calibration_GSE142025.png"), type = "windows", width = 1400, height = 1400, res = 240)
vp <- val.prob(p_hat, y3)
title("Calibration - GSE142025 (recalibrated)")
dev.off()
cat("FIG7B_DONE shorttitle; Eavg", vp["Eavg"], "Emax", vp["Emax"], "Brier", vp["Brier"], "S:p", vp["S:p"], "\n")

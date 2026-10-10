# 16_fig_round3.R ----
# Round-3 图件修复：
#   Fig1: 3处文字（LRDEG术语/Nomogram/19→18细胞）
#   Fig3: 悬空括注改为指向Table S8
#   Fig6: 全部面板只读 results/frozen/*.csv（不再重算ML）；E面板加宽右边距
#   Fig7A: 图例加AUC数值；Fig7B/7C: canonical系数重生成；7C刻度显式0.2步进
#   Fig8A: 星号改用BH18口径（TableS7）；Fig8B: 去Normal_mucosa列+改图题(all samples)
#   Fig2C: 面板字母C移至左上（与A/B统一）
suppressWarnings(suppressMessages({
  library(data.table); library(glmnet); library(e1071); library(pROC)
  library(ggplot2); library(ggrepel); library(ComplexHeatmap); library(circlize)
  library(patchwork); library(dcurves); library(rms); library(grid)
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

# ========== Fig1 流程图（3处文字修正） ==========
png(file.path(fig_dir, "Fig1_flowchart.png"), type = "windows", width = 2400, height = 1500, res = 220)
par(mar = c(1, 1, 2, 1)); plot.new(); plot.window(c(0, 100), c(0, 62), xaxs = "i", yaxs = "i")
box_at <- function(x, y, w, h, label, fill = "white", cex = 0.95) {
  symbols(x, y, rectangles = matrix(c(w, h), 1), inches = FALSE, add = TRUE,
          fg = "grey25", bg = fill, lwd = 1.4)
  text(x, y, label, cex = cex)
}
arr <- function(x1, y1, x2, y2) arrows(x1, y1, x2, y2, lwd = 1.6, length = 0.12)
title("Study workflow", cex.main = 1.5, adj = 0)
box_at(16, 55, 26, 7, "GEO cohorts\nTrain: GSE96804 (DN 41 / Ctrl 20)\nValidation: 4 cohorts (n=97)", "#DCEEFB")
box_at(50, 55, 24, 7, "Lactylation-related genes\nGO + Reactome + literature\n(171 genes)", "#DCEEFB")
box_at(84, 55, 24, 7, "DEGs (limma)\n2,645 genes\n+ WGCNA modules", "#DCEEFB")
arr(29, 55, 38, 55); arr(62, 55, 72, 55)
box_at(50, 41, 44, 7, "Candidate LRGs: 23  (LRG \u2229 (DEG \u222a disease modules))", "#FDF2D0")
arr(84, 51.5, 60, 44.8); arr(50, 51.5, 50, 44.7)
box_at(22, 30, 34, 8, "3 ML algorithms\nLASSO | SVM-RFE | XGBoost\nconsensus", "#E8F5E9")
box_at(70, 30, 40, 8, "Signature: PDK2 + PDK4 (penalized logistic)\nMulti-cohort ROC | Calibration | DCA", "#E8F5E9")
arr(42, 37, 36, 34); arr(58, 37, 64, 34)
box_at(22, 15, 34, 8, "ROC validation\nGSE30528 / GSE30529 / GSE142025\nincl. early-stage DN | GSE104948", "#FCE4EC")
box_at(70, 15, 40, 8, "Immune infiltration (ssGSEA, 18 cell types)\nGene-cell correlation", "#FCE4EC")
arr(22, 26, 22, 19); arr(70, 26, 70, 19)
box_at(46, 4, 50, 6, "Lactylation-related diagnostic signature for DKD", "#D1C4E9", 1.05)
arr(32, 11, 42, 7); arr(66, 11, 58, 7)
dev.off()
cat("Fig1 done\n")

# ========== Fig3（标题括注指向Table S8） ==========
wgc_mt <- fread("results/wgcna_module_trait.csv")
wgc_p <- wgc_mt[module != "grey"]
png(file.path(fig_dir, "Fig3_wgcna_rev.png"), type = "windows", width = 2400, height = 1100, res = 240)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2) + .1, oma = c(0, 2, 0, 0))
suppressWarnings(suppressMessages(library(WGCNA)))
mads_f <- apply(train_expr, 1, mad)
dat_f <- t(train_expr[order(-mads_f)[1:5000], ])
sft_f <- pickSoftThreshold(dat_f, powerVector = c(1:10, seq(12, 20, 2)),
                           networkType = "signed", verbose = 0, corFnc = "cor")
fi <- sft_f$fitIndices
plot(fi$Power, -sign(fi$slope) * fi$SFT.R.sq, xlab = "Soft threshold (power)",
     ylab = "Scale-free fit R2", main = "WGCNA soft threshold", pch = 16, col = "#B2182B",
     cex.lab = 1, cex.axis = .9)
abline(h = .85, col = "red", lty = 2)
text(fi$Power, -sign(fi$slope) * fi$SFT.R.sq, labels = fi$Power, cex = .65, pos = 4)
par(mar = c(5, 9, 4, 2) + .1)
mt2 <- matrix(wgc_p$V1, ncol = 1, dimnames = list(wgc_p$module, "DN"))
lab <- matrix(paste0(signif(wgc_p$V1, 2), "\n(", signif(wgc_p$p, 1), ")"),
              ncol = 1, dimnames = list(wgc_p$module, "DN"))
labeledHeatmap(Matrix = mt2, yLabels = wgc_p$module, ySymbols = wgc_p$module,
               colorLabels = FALSE, colors = blueWhiteRed(50),
               setStdMargins = FALSE, cex.text = .5, xLabels = "DN",
               xLabelsPosition = "bottom", textMatrix = lab,
               main = "Module-trait correlations (grey in Table S8)")
dev.off()
cat("Fig3 done\n")

# ========== Fig6（只读冻结件） ==========
fl <- fread("results/frozen/lasso_curve.csv")
p1 <- ggplot(fl, aes(lambda, cvm)) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), color = "grey70") +
  geom_point(color = "#B2182B", size = .9) + scale_x_log10() +
  geom_vline(xintercept = fl$lambda[fl$is_1se], linetype = 2) +
  theme_bw(base_size = 10) + ylab("Binomial deviance") + xlab("lambda") +
  ggtitle("LASSO CV (lambda.1se)") + labs(tag = "A")
fs <- fread("results/frozen/svm_rfe_curve.csv")
p2 <- ggplot(fs, aes(k, accuracy)) +
  geom_line(color = "#2166AC") + geom_point(size = 2) +
  geom_vline(xintercept = fs$k[fs$is_best], linetype = 2) +
  theme_bw(base_size = 10) + xlab("Number of features") + ylab("CV accuracy") +
  ggtitle(sprintf("SVM-RFE (optimum k = %d)", fs$k[fs$is_best])) + labs(tag = "B") +
  annotate("text", x = fs$k[fs$is_best], y = min(fs$accuracy) + .01, label = sprintf("k = %d", fs$k[fs$is_best]),
           vjust = -0.2, hjust = -0.15, size = 3) +
  ylim(min(fs$accuracy) - .02, 1)
fx <- fread("results/frozen/xgb_importance.csv")
fx$Gain_pct <- fx$Gain_pct
p3 <- ggplot(head(fx, 10), aes(reorder(Feature, Gain_pct), Gain_pct)) +
  geom_col(fill = "#4CAF50", width = .7) + coord_flip() +
  theme_bw(base_size = 10) + xlab(NULL) + ylab("Gain (%)") +
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
  geom_text(data = ann, aes(x = gene, y = expr, label = lab), size = 3.2, inherit.aes = FALSE) +
  scale_fill_manual(values = c(Control = "#2166AC", DN = "#B2182B")) +
  theme_bw(base_size = 10) + xlab(NULL) + ylab("Expression (log2)") +
  ggtitle("Hub gene expression (training)") + labs(tag = "D") +
  theme(legend.position = "bottom") +
  annotate("text", x = .5, y = max(hubexp$expr) * 1.12, hjust = 0,
           label = "ns / * / ** / *** vs adj.P 0.05/0.01/0.001", size = 2.6)
rocs <- lapply(hub, function(g) roc(y, Xdf[[g]], quiet = TRUE))
names(rocs) <- hub
p5 <- ggroc(rocs, linewidth = .7) +
  geom_abline(slope = 1, intercept = 1, linetype = 3) +
  scale_color_brewer(palette = "Set1") + theme_bw(base_size = 10) +
  coord_equal() + ggtitle("Single-gene ROC (training)") + labs(color = NULL, tag = "E") +
  theme(legend.position = "bottom", legend.box.margin = margin(0, 30, 0, 0, "pt"))
save_fig((p1 | p2 | p3) / (p4 | p5), "Fig6_ml_screening_rev", 14.2, 7.8)
cat("Fig6 done\n")

# ========== Fig7A（canonical系数 + 图例AUC） ==========
fm <- fread("results/frozen/final_model_coefficients.csv")
bh <- setNames(fm$canonical_deviance, fm$term)
sc <- function(e) {
  hg <- intersect(hub, rownames(e))
  Xv <- t(zrows(e[hg, , drop = FALSE]))
  as.vector(Xv %*% bh[hg] + bh[["(Intercept)"]])
}
Xh <- data.frame(t(zrows(train_expr[hub, , drop = FALSE])), check.names = FALSE)
roc_list <- list(Training = roc(y, as.vector(as.matrix(Xh) %*% bh[hub] + bh[["(Intercept)"]]), quiet = TRUE))
for (acc in c("GSE30528", "GSE30529")) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  p <- read.csv(paste0("data/", acc, "_pheno.csv"), stringsAsFactors = FALSE)
  keep <- p$group %in% c("DKD", "Control")
  roc_list[[acc]] <- roc(as.integer(p$group[keep] == "DKD"), sc(e[, keep, drop = FALSE]), quiet = TRUE)
}
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])))
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
rownames(e14) <- toupper(rownames(e14))
k3 <- p14$group %in% c("Advanced_DN", "Early_DN", "Control")
k4 <- p14$group %in% c("Early_DN", "Control")
roc_list[["GSE142025 (all DN)"]] <- roc(as.integer(p14$group[k3] != "Control"), sc(e14[, k3, drop = FALSE]), quiet = TRUE)
roc_list[["GSE142025 (early)"]] <- roc(as.integer(p14$group[k4] == "Early_DN"), sc(e14[, k4, drop = FALSE]), quiet = TRUE)
e10 <- readRDS("data/GSE104948_expr_symbol.rds")
p10 <- read.csv("data/GSE104948_pheno.csv", stringsAsFactors = FALSE)
k5 <- p10$group %in% c("DN", "Control")
roc_list[["GSE104948"]] <- roc(as.integer(p10$group[k5] == "DN"), sc(e10[, k5, drop = FALSE]), quiet = TRUE)
nm <- names(roc_list)
names(roc_list) <- sprintf("%s (AUC %s)", nm, formatC(sapply(roc_list, function(r) as.numeric(auc(r))), format = "f", digits = 3))
p7 <- ggroc(roc_list, linewidth = .75) + geom_abline(slope = 1, intercept = 1, linetype = 3) +
  coord_equal() + theme_bw(base_size = 10.5) +
  scale_color_brewer(palette = "Dark2") + labs(color = NULL, tag = "A") +
  ggtitle("PDK2+PDK4 signature ROC across cohorts") +
  theme(legend.position = c(.38, .16),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.key.width = unit(22, "pt"), legend.text = element_text(size = 8.2))
save_fig(p7, "Fig7A_validation_roc_rev", 7.6, 6.8)
cat("Fig7A done; AUCs:", paste(sprintf("%.3f", sapply(roc_list, function(r) as.numeric(auc(r)))), collapse = ", "), "\n")

# ========== Fig7B/7C（canonical s3） ==========
y3 <- as.integer(p14$group[k3] != "Control")
s3 <- sc(e14[, k3, drop = FALSE])
g1 <- glm(y3 ~ s3, family = binomial)
p_hat <- fitted(g1)
png(file.path(fig_dir, "Fig7B_calibration_GSE142025.png"), type = "windows", width = 1400, height = 1400, res = 240)
vp <- val.prob(p_hat, y3)
title("Calibration - GSE142025 (DN vs control, recalibrated score)")
dev.off()
cat("Fig7B stats: Eavg", vp["Eavg"], "Emax", vp["Emax"], "Brier", vp["Brier"], "S:p", vp["S:p"], "\n")
dca_df <- data.frame(y = y3, Signature = p_hat)
png(file.path(fig_dir, "Fig7C_dca_GSE142025.png"), type = "windows", width = 1500, height = 1300, res = 240)
pd <- tryCatch({
  pl <- plot(dca(y ~ Signature, data = dca_df), smooth = FALSE)
  pl + scale_x_continuous(breaks = seq(0, 1, .2)) + ggtitle("Decision curve analysis - GSE142025")
}, error = function(e) { cat("DCA gg adjust fail:", conditionMessage(e), "\n"); plot(dca(y ~ Signature, data = dca_df, smooth = FALSE)) })
print(pd)
dev.off()
cat("Fig7C done\n")

# ========== Fig8A（BH18星号） ==========
immune_ss <- fread("results/immune_ssgsea.csv")
it7 <- fread("manuscript/supplement/TableS7_immune_tests_18cells.csv")
it2 <- it7[p.adj_BH18 < 0.05]
sig_cells <- it2[order(p.adj_BH18)]$cell
iss <- as.data.frame(immune_ss); rownames(iss) <- iss$sample; iss$sample <- NULL
grp <- train_ph$group
dfi <- do.call(rbind, lapply(sig_cells, function(cn)
  data.frame(cell = cn, score = as.numeric(iss[[cn]]), group = grp)))
dfi$cell <- factor(dfi$cell, levels = sig_cells)
starvec <- setNames(ifelse(it2$p.adj_BH18 < .001, "***", ifelse(it2$p.adj_BH18 < .01, "**", "*")), it2$cell)
ann8 <- data.frame(cell = factor(sig_cells, levels = sig_cells), lab = starvec[sig_cells],
                   score = max(dfi$score, na.rm = TRUE) * 1.08, group = "Control")
p8a <- ggplot(dfi, aes(cell, score, fill = group)) +
  geom_boxplot(outlier.size = .4, linewidth = .3) +
  geom_text(data = ann8, aes(x = cell, y = score, label = lab), size = 2.9, inherit.aes = FALSE) +
  scale_fill_manual(values = c(Control = "#2166AC", DN = "#B2182B")) +
  theme_bw(base_size = 10) + xlab(NULL) + ylab("ssGSEA score") +
  ggtitle(sprintf("Immune infiltration (significant cells, n = %d of 18; BH across 18 tests)", length(sig_cells))) +
  labs(tag = "A") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
save_fig(p8a, "Fig8A_immune_box_rev", 9.2, 5.2)
cat("Fig8A done; stars:", paste(sig_cells, starvec[sig_cells], sep = "=", collapse = " "), "\n")

# ========== Fig8B（去Normal_mucosa + all samples） ==========
immune_corr <- read.csv("results/immune_gene_corr.csv", stringsAsFactors = FALSE, check.names = FALSE)
cm <- as.matrix(immune_corr[, setdiff(colnames(immune_corr), c("gene", "Normal_mucosa"))])
rownames(cm) <- immune_corr$gene
col_fun2 <- colorRamp2(c(-.6, 0, .6), c("#2166AC", "white", "#B2182B"))
png(file.path(fig_dir, "Fig8B_corr_rev.png"), type = "windows", width = 2400, height = 1600, res = 250)
draw(Heatmap(cm, name = "Spearman rho", col = col_fun2,
             row_names_gp = gpar(fontsize = 7.5), column_names_gp = gpar(fontsize = 7.5),
             column_title = "Candidate LRGs x immune cells (GSE96804, all samples; 18 cells)"))
grid.text("B", x = unit(2.6, "mm"), y = unit(1, "npc") - unit(2.6, "mm"),
          gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()
cat("Fig8B done; columns:", ncol(cm), "\n")

# ========== Fig2C（字母C移左上；沿用14b精确过滤） ==========
deg <- fread("results/A_DEGs.csv")
pseudogenes_observed <- c("PFN1P2","SLC7A5P2","SMG1P1","CHCHD2P6","CHCHD2P8","NBPF8P")
nc_filter <- function(g) {
  grepl("^(MIR|SNOR|SNORA|SNORD|SCAR|RNU|LINC|RNA|DANCR|VTRNA|NPIPA|NBPF|LOC[0-9]|OTTHUM|ENSG)", g) |
    grepl("-AS[0-9]$", g) | (g %in% pseudogenes_observed)
}
top_pool <- deg[order(deg$adj.P.Val)][1:120]
top50 <- top_pool[!nc_filter(top_pool$gene)][1:50]
mat <- t(scale(t(train_expr[top50$gene, ])))
ha <- HeatmapAnnotation(Group = train_ph$group,
                        col = list(Group = c(Control = "#2166AC", DN = "#B2182B")))
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
png("results/figures/Fig2C_heatmap_rev.png", type = "windows", width = 1800, height = 2400, res = 250)
draw(Heatmap(mat, name = "z-score", col = col_fun, top_annotation = ha,
             show_column_names = FALSE, row_names_gp = gpar(fontsize = 6),
             column_title = "Top 50 protein-coding DEGs (GSE96804)"))
grid.text("C", x = unit(2.6, "mm"), y = unit(1, "npc") - unit(2.6, "mm"),
          gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()
cat("Fig2C done:", nrow(top50), "rows\n")
cat("\nFIG_ROUND3_DONE\n")

# 13_figures.R ----
# Render the nine main figures (300 dpi PNG panels) from the frozen results
# of steps 3-12. TIFF assembly is done afterwards by assemble_figures.py.
# Outputs (all in figures/): Fig1_flowchart.png, Fig2AB_pca_volcano_rev.png,
#   Fig2C_heatmap_rev.png, Fig3_wgcna_rev.png, Fig4A_venn_rev.png,
#   Fig4B_bar_rev.png, Fig5_enrichment_rev.png, Fig6_ml_screening_rev.png,
#   Fig7A_validation_roc_rev.png, Fig7B_calibration_GSE142025.png,
#   Fig7C_dca_GSE142025.png, Fig8A_immune_box_rev.png, Fig8B_corr_rev.png,
#   Fig9_scrna_localization_rev.png

suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
  library(WGCNA); library(ComplexHeatmap); library(circlize)
  library(VennDiagram); library(pROC); library(rms); library(viridis)
  library(Seurat); library(grid)
}))
set.seed(42)
fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
png_dev <- function(filename, width, height, units = "in", res = 300, ...)
  grDevices::png(filename, width = width, height = height, units = units,
                 res = res, ...)

train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
allg <- fread("results/A_all_genes_limma.csv")
deg <- fread("results/A_DEGs.csv")
hub <- readLines("results/hub_genes.txt")
cand <- fread("results/A_candidates.csv")
zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
genes <- intersect(cand$gene, rownames(train_expr))
Xdf <- data.frame(t(zrows(train_expr[genes, , drop = FALSE])), check.names = FALSE)
y <- as.integer(train_ph$group == "DN")

# ============================================================
# Fig 1: study flowchart
# ============================================================
# Canvas 4216 x 2720 @600 dpi; the y-range extends to -2 so the conclusion
# box (bottom edge at y = -1.5) renders with its full border at unchanged
# scale - all boxes, arrows, and labels keep their original coordinates.
png(file.path(fig_dir, "Fig1_flowchart.png"),
    width = 4216, height = 2720, res = 600, pointsize = 10)
par(mar = c(1, 1, 1, 1)); plot.new(); plot.window(c(0, 100), c(-2, 62.5), xaxs = "i", yaxs = "i")
box_at <- function(x, y, w, h, label, fill = "white", cex = 0.85, ty = 0) {
  symbols(x, y, rectangles = matrix(c(w, h), 1), inches = FALSE, add = TRUE,
          fg = "grey25", bg = fill, lwd = 1.6)
  text(x, y + ty, label, cex = cex)
}
arr <- function(x1, y1, x2, y2) arrows(x1, y1, x2, y2, lwd = 1.7, length = 0.06, angle = 25)

box_at(21, 57, 34, 7.5, "GEO cohorts\nTrain: GSE96804 (DN 41/Control 20)\nExternal cohorts: 5 (n=160)", "#DCEEFB")
box_at(76, 57, 30, 7.5, "Lactylation-related genes\nGO + Reactome + literature\n(171 genes)", "#DCEEFB", ty = 0.06)
box_at(21, 44, 34, 7, "Differential expression & modules\nDEGs: 2,645 (limma)\nWGCNA disease modules", "#DCEEFB")
box_at(76, 44, 30, 7, "Candidate LRGs: 23\nLRG \u2229 (DEG \u222a modules)", "#FDF2D0")
box_at(21, 29, 34, 8, "3 ML algorithms\nLASSO | SVM-RFE | XGBoost\nconsensus (\u22652) \u2192 5 hub genes", "#C8E6C9")
box_at(76, 29, 30, 8, "Two-gene signature\nPDK2 + PDK4\n(penalized logistic regression)", "#C8E6C9")
box_at(15.5, 13.5, 30, 12.6, paste0(
  "Validation (frozen signature)\n",
  "GSE30528 | GSE30529 (compartment)\n",
  "GSE142025 (incl. early DN)\n",
  "GSE104948 | GSE273763 (FFPE)\n",
  "ROC: all cohorts\n",
  "Calibration/DCA: GSE142025+GSE273763"), "#FCE4EC", cex = 0.70, ty = -0.14)
box_at(50.5, 14, 24, 6, "snRNA-seq localization\n(GSE131882; 29,048 nuclei)", "#FCE4EC")
box_at(84, 14, 25, 8, "Immune infiltration\nssGSEA (18 cell types)\nGSE142025 replication\nGSE273763 (FFPE): Suppl. S17", "#FCE4EC", cex = 0.75)
box_at(48.5, 1.5, 55, 6, "Lactylation-related diagnostic signature for DKD", "#D1C4E9", 0.9)

arr(21, 53.2, 21, 47.6);        text(22.3, 50.4, "training", cex = 0.58, adj = c(0, -0.2))
arr(76, 53.2, 76, 47.6)
arr(38, 44, 61, 44)
arr(66.5, 40.5, 30.5, 33.9)
arr(38, 29, 61, 29)
arr(63, 25, 18.6, 20.6)
arr(72, 25, 50.5, 17.2)
arr(84, 25, 84, 18.2)
arr(16, 7.2, 38.7, 4.9)
arr(50.5, 11, 46.5, 4.8)
arr(84, 10, 58.3, 4.9)
dev.off()
cat("Fig1 done\n")

# ============================================================
# Fig 2A/B: PCA + volcano (training)
# ============================================================
pca_res <- prcomp(t(train_expr))
var_pct <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], group = train_ph$group)
p1 <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 2.6, alpha = .9) + stat_ellipse(level = .9, linetype = 2, linewidth = .6) +
  scale_color_manual(values = c(Control = "#2166AC", DN = "#B2182B")) +
  theme_bw(base_size = 12) +
  labs(color = NULL, tag = "A",
       x = sprintf("PC1 (%.1f%%)", var_pct[1]), y = sprintf("PC2 (%.1f%%)", var_pct[2])) +
  ggtitle("Training cohort (GSE96804)") +
  theme(legend.position = "bottom")

vol <- data.frame(gene = allg$gene, logFC = allg$logFC, neglog10 = -log10(allg$adj.P.Val))
vol$sig <- ifelse(vol$logFC >= .5 & vol$neglog10 > -log10(.05), "Up",
           ifelse(vol$logFC <= -.5 & vol$neglog10 > -log10(.05), "Down", "NS"))
lb <- vol[vol$gene %in% c("PDK4", "PDK2", "ACSF2") & vol$sig != "NS", ]
p2 <- ggplot(vol, aes(logFC, neglog10, color = sig)) +
  geom_point(size = .55, alpha = .45) +
  scale_color_manual(values = c(Up = "#B2182B", Down = "#2166AC", NS = "grey78")) +
  geom_hline(yintercept = -log10(.05), linetype = 3) +
  geom_vline(xintercept = 0.5, linetype = 3, color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = -0.5, linetype = 3, color = "grey40", linewidth = 0.5) +
  annotate("label", x = 1.5, y = 2.8, label = "adjusted P = 0.05",
           size = 3.5, hjust = 1, color = "black", fill = "white", alpha = 0.85,
           label.size = 0, label.padding = unit(0.12, "lines")) +
  geom_label_repel(data = lb, aes(label = gene), size = 3.0, max.overlaps = 30,
                   show.legend = FALSE, seed = 123, force = 5,
                   box.padding = 0.6, min.segment.length = 0,
                   fill = "white", alpha = 0.85, label.size = 0.15,
                   max.time = 2, max.iter = 10000) +
  annotate("label", x = 3, y = 14, label = sprintf("Up\nn = %d", sum(vol$sig == "Up")),
           size = 3.2, color = "#B2182B", hjust = 0, fill = "white", alpha = 0.85,
           label.size = 0, label.padding = unit(0.12, "lines")) +
  annotate("label", x = -3, y = 14, label = sprintf("Down\nn = %d", sum(vol$sig == "Down")),
           size = 3.2, color = "#2166AC", hjust = 1, fill = "white", alpha = 0.85,
           label.size = 0, label.padding = unit(0.12, "lines")) +
  theme_bw(base_size = 12) + labs(tag = "B",
       x = "log2FC", y = expression(-log[10](adjusted~P))) +
  ggtitle(sprintf("Volcano (DEGs n = %d)", nrow(deg))) +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = -4:4, limits = c(-4.5, 4.5))
ggsave(file.path(fig_dir, "Fig2AB_pca_volcano_rev.png"), p1 | p2,
       width = 9.6, height = 4.4, dpi = 300, device = png_dev)
cat("Fig2AB done\n")

# ============================================================
# Fig 2C: top-50 protein-coding DEG heatmap
# ============================================================
top_pool <- deg[order(deg$adj.P.Val)][1:120]
nc_prefixes <- "^(MIR|SNOR|RNU|LINC|RNA|NPIPA|NBPF|LOC|OTTHUM|VTRNA|DANCR|PFN1P|SLC7A5P|SMG1P|NEBL|CHCHD2P)"
antisense <- grepl("-AS[0-9]$", top_pool$gene)
top50 <- top_pool[!grepl(nc_prefixes, top_pool$gene) & !antisense][1:50]
stopifnot(all(c("DUSP1", "EIF4EBP2", "PDK2", "PDK4", "MXD1") %in% top50$gene))
mat <- t(scale(t(train_expr[top50$gene, ])))
ha <- HeatmapAnnotation(Group = train_ph$group,
                        col = list(Group = c(Control = "#2166AC", DN = "#B2182B")))
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
png(file.path(fig_dir, "Fig2C_heatmap_rev.png"),
    width = 1800, height = 2400, res = 250)
draw(Heatmap(mat, name = "z-score", col = col_fun, top_annotation = ha,
             show_column_names = FALSE, row_names_gp = gpar(fontsize = 6),
             column_title = "Top 50 protein-coding DEGs (GSE96804)"))
grid.text("C", x = unit(2.2, "mm"), y = unit(1, "npc") - unit(2.2, "mm"),
          gp = gpar(fontsize = 10, fontface = "bold"))
dev.off()
cat("Fig2C done\n")

# ============================================================
# Fig 3: WGCNA soft threshold + module-trait (ggplot wide tiles)
# ============================================================
wgc_mt <- fread("results/wgcna_module_trait.csv")
wgc_p <- wgc_mt[module != "grey"]
mads_f <- apply(train_expr, 1, mad)
dat_f <- t(train_expr[order(-mads_f)[1:5000], ])
sft_f <- pickSoftThreshold(dat_f, powerVector = c(1:10, seq(12, 20, 2)),
                           networkType = "signed", verbose = 0, corFnc = "cor")
fi <- sft_f$fitIndices
fi$sft <- -sign(fi$slope) * fi$SFT.R.sq
fp1 <- ggplot(fi, aes(Power, sft)) +
  geom_hline(yintercept = .85, linetype = 2, color = "red") +
  geom_point(color = "#B2182B", shape = 16, size = 1.4) +
  geom_text(aes(label = Power), cex = 2.2, vjust = -0.9) +
  theme_bw(base_size = 9) +
  xlab("Soft threshold (power)") + ylab("Scale-free fit R2") +
  ggtitle("WGCNA soft threshold") +
  theme(plot.title = element_text(size = 9))
mt <- copy(wgc_p)
setnames(mt, c("corr_r", "p"), c("r", "pval"))
mt[, lab := paste0(signif(r, 2), "\n(", signif(pval, 1), ")")]
mt[, module := factor(module, levels = rev(module))]
rng <- range(mt$r)
fp2 <- ggplot(mt, aes(x = "DN", y = module, fill = r)) +
  geom_tile(width = 0.98, height = 0.98, color = "white") +
  geom_text(aes(label = lab), size = 2.4, lineheight = 0.85) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(rng[1], rng[2]), name = "Pearson r") +
  theme_bw(base_size = 9) + xlab(NULL) + ylab(NULL) +
  ggtitle("Module-trait correlations") +
  theme(axis.text.x = element_text(size = 8), panel.grid = element_blank(),
        plot.title = element_text(size = 9), legend.key.height = unit(0.8, "cm"))
grDevices::png(file.path(fig_dir, "Fig3_wgcna_rev.png"),
               width = 7.0, height = 3.33, units = "in", res = 300)
print(fp1 + fp2 + plot_layout(widths = c(1, 1.15)))
grDevices::dev.off()
cat("Fig3 done\n")

# ============================================================
# Fig 4A: candidate-source Venn; Fig 4B: candidate log2FC bar
# ============================================================
lrg <- fread("data/LRG_geneset.csv")
lrg_alias <- lrg$alias[!is.na(lrg$alias) & nzchar(trimws(lrg$alias))]
lrg_entries <- toupper(unique(c(lrg$symbol, lrg_alias)))
deg_genes <- deg$gene
mod_genes <- fread("results/wgcna_sigmodule_genes.csv")$gene
stopifnot(length(lrg_entries) == 177, length(deg_genes) == 2645, length(mod_genes) == 4210)
venn_list <- list(
  `LRG entries (177)` = lrg_entries,
  `DEGs (2645)` = deg_genes,
  `Modules (4210)` = mod_genes)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
png(file.path(fig_dir, "Fig4A_venn_rev.png"),
    width = 1500, height = 1300, res = 250)
v <- venn.diagram(venn_list, filename = NULL,
                  fill = c("#F4B183", "#9DC3E6", "#C5E0B4"), alpha = .55,
                  cex = 1.2, cat.cex = .9, margin = .10,
                  fontfamily = "sans", cat.fontfamily = "sans")
grid.draw(v)
grid.text("A", x = unit(2, "mm"), y = unit(1, "npc") - unit(2, "mm"),
          gp = gpar(fontsize = 10, fontface = "bold"))
dev.off()

cd <- cand[!is.na(cand$logFC)]
cd$lab <- ifelse(cd$adj.P.Val < 0.001, "***", ifelse(cd$adj.P.Val < 0.01, "**",
          ifelse(cd$adj.P.Val < 0.05, "*", "")))
cd <- cd[order(cd$logFC)]
p_bar <- ggplot(cd, aes(reorder(gene, logFC), logFC, fill = LRG_sources)) +
  geom_col(width = .72, color = "grey25", linewidth = .18) +
  geom_text(aes(label = lab, y = logFC + sign(logFC) * 0.18), size = 5.5, na.rm = TRUE) +
  coord_flip() + theme_bw(base_size = 11) +
  labs(x = NULL, y = "log2FC (DN versus Control)", fill = "LRG source", tag = "B",
       title = "Differentially expressed candidate LRGs (n = 11)") +
  scale_y_continuous(breaks = -3:1, limits = c(-3.5, 1.2)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))
ggsave(file.path(fig_dir, "Fig4B_bar_rev.png"), p_bar,
       width = 7.0, height = 5.5, dpi = 300, device = png_dev)
cat("Fig4 done\n")

# ============================================================
# Fig 5: enrichment panels (cividis gradient)
# ============================================================
ego <- fread("results/enrich_GO_BP.csv")
er <- fread("results/enrich_Reactome.csv")
gse <- fread("results/gsea_Reactome.csv")
mk_bar <- function(d, n, ttl, tag) {
  if (!hasName(d, "Count")) d[, Count := lengths(strsplit(geneID, "/"))]
  d <- utils::head(d[order(d$p.adjust)], n)
  d$lab <- make.names(substr(d$Description, 1, 52))
  d$lab <- factor(d$lab, levels = rev(d$lab))
  d$logFDR <- log10(d$p.adjust)
  ggplot(d, aes(Count, lab, fill = logFDR)) + geom_col(width = .7) +
    scale_fill_viridis_c(option = "cividis",
                         name = expression(log[10](adjusted~P))) +
    theme_bw(base_size = 10) + ylab(NULL) + ggtitle(ttl) + labs(tag = tag) +
    theme(plot.title = element_text(size = 10.5))
}
fp1 <- mk_bar(ego, 10, "GO BP (top DEGs)", "A")
fp2 <- mk_bar(er, 10, "Reactome (top DEGs)", "B")
gse2 <- copy(gse)
gse2$lab <- make.names(paste0(substr(gse2$Description, 1, 48), "  NES=", sprintf("%.2f", gse2$NES)))
topg <- rbind(head(gse2[order(NES), ], 4), head(gse2[order(-NES), ], 4))
topg$lab <- factor(topg$lab, levels = topg$lab[order(topg$NES)])
fp3 <- ggplot(topg, aes(NES, lab, fill = log10(p.adjust))) + geom_col(width = .7) +
  scale_fill_viridis_c(option = "cividis", name = expression(log[10](FDR))) +
  theme_bw(base_size = 9) + ylab(NULL) +
  ggtitle("GSEA (Reactome, top 4 up- and top 4 down-regulated)") + labs(tag = "C")
grDevices::png(file.path(fig_dir, "Fig5_enrichment_rev.png"),
               width = 7.0, height = 10, units = "in", res = 300)
print(fp1 / fp2 / fp3 + plot_layout(heights = c(1, 1, .8)))
grDevices::dev.off()
cat("Fig5 done\n")

# ============================================================
# Fig 6: ML screening (A-E; E with deterministic label stack)
# ============================================================
fl <- fread("results/frozen/lasso_curve.csv")
lambda_range <- range(fl$lambda)
brks <- c(1e-5, 1e-4, 1e-3, 0.01, 0.1, 1)
brks <- brks[brks >= lambda_range[1] & brks <= lambda_range[2]]
if (length(brks) < 2) brks <- pretty(log10(fl$lambda))
mp1 <- ggplot(fl, aes(lambda, cvm)) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), color = "grey70") +
  geom_point(color = "#D55E00", size = .9) + scale_x_log10(breaks = brks) +
  geom_vline(xintercept = fl$lambda[fl$is_1se], linetype = 2) +
  theme_bw(base_size = 9) + ylab("Binomial deviance") + xlab("lambda (log scale)") +
  ggtitle("LASSO CV (lambda.1se)") + labs(tag = "A") +
  theme(plot.title = element_text(size = 8.5))
fs <- fread("results/frozen/svm_rfe_curve.csv")
mp2 <- ggplot(fs, aes(k, accuracy)) +
  geom_line(color = "#0072B2") + geom_point(size = 1.5) +
  geom_vline(xintercept = fs$k[fs$is_best], linetype = 2) +
  theme_bw(base_size = 9) + xlab("Number of features") + ylab("CV accuracy") +
  ggtitle(sprintf("SVM-RFE (optimum k = %d)", fs$k[fs$is_best])) + labs(tag = "B") +
  theme(plot.title = element_text(size = 8.5))
fx <- fread("results/frozen/xgb_importance.csv")
mp3 <- ggplot(head(fx, 10), aes(reorder(Feature, Gain_pct), Gain_pct)) +
  geom_col(fill = "#E69F00", width = .7) + coord_flip() +
  theme_bw(base_size = 9) + xlab(NULL) + ylab("Gain (%)") +
  ggtitle("XGBoost gain (top 10)") + labs(tag = "C") +
  theme(plot.title = element_text(size = 8.5))
hubexp <- do.call(rbind, lapply(hub, function(g) {
  v <- as.numeric(train_expr[g, ])
  data.frame(gene = g, expr = v, group = train_ph$group)
}))
allgk <- allg[gene %in% hub]
pv <- setNames(allgk$adj.P.Val, allgk$gene)
star <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "ns")))
hubexp$gene <- factor(hubexp$gene, levels = hub)
ann <- data.frame(gene = factor(hub, levels = hub), lab = star(pv[hub]),
                  expr = max(hubexp$expr) * 1.06, group = "DN")
mp4 <- ggplot(hubexp, aes(gene, expr, fill = group)) +
  geom_boxplot(outlier.size = .4, linewidth = .25) +
  geom_text(data = ann, aes(x = gene, y = expr, label = lab), size = 2.8, inherit.aes = FALSE) +
  scale_fill_manual(values = c(Control = "#0072B2", DN = "#D55E00")) +
  theme_bw(base_size = 9) + xlab(NULL) + ylab("Expression (log2)") +
  ggtitle("Hub gene expression (training)") + labs(tag = "D") +
  theme(legend.position = "bottom", plot.title = element_text(size = 8.5))
rocs <- lapply(hub, function(g) roc(y, Xdf[[g]], quiet = TRUE))
names(rocs) <- hub
oi <- c(PDK4 = "#000000", PDK2 = "#0072B2", ACSF2 = "#E69F00",
        HIST1H4B = "#56B4E9", ME1 = "#CC79A7")
lt <- c(PDK4 = "solid", PDK2 = "solid", ACSF2 = "longdash",
        HIST1H4B = "dashed", ME1 = "dotdash")
present <- intersect(names(oi), hub)
roc_df <- do.call(rbind, lapply(present, function(g) {
  r <- rocs[[g]]
  data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities,
             gene = g, auc = as.numeric(auc(r)))
}))
best_pts <- do.call(rbind, lapply(present, function(g) {
  coords(rocs[[g]], "best", ret = c("specificity", "sensitivity"), best.method = "youden")
}))
best_df <- data.frame(gene = present, fpr = 1 - best_pts$specificity,
                      tpr = best_pts$sensitivity)
best_df$auc <- sapply(present, function(g) as.numeric(auc(rocs[[g]])))
lab_ord <- best_df[order(-best_df$auc), ]
rows_y <- c(0.36, 0.27, 0.18, 0.10, 0.03)[seq_len(nrow(lab_ord))]
x0 <- 0.42
lab_df <- data.frame(gene = lab_ord$gene, auc = lab_ord$auc,
                     lx = x0, ly = rows_y, sx = lab_ord$fpr, sy = lab_ord$tpr)
mp5 <- ggplot(roc_df, aes(fpr, tpr, color = gene, linetype = gene)) +
  geom_line(linewidth = .6) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey50") +
  scale_color_manual(values = oi[present], guide = "none") +
  scale_linetype_manual(values = lt[present], guide = "none") +
  geom_segment(data = lab_df, aes(x = sx, y = sy, xend = lx - 0.012, yend = ly),
               inherit.aes = FALSE, color = "grey35", linewidth = 0.25, alpha = 0.55) +
  geom_text(data = lab_df, aes(x = lx, y = ly, label = sprintf("%s  AUC=%.3f", gene, auc)),
            hjust = 0, vjust = 0.5, size = 2.35, show.legend = FALSE,
            inherit.aes = FALSE, color = oi[lab_df$gene]) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  coord_equal() +
  theme_bw(base_size = 9) + xlab("1 \u2212 Specificity") + ylab("Sensitivity") +
  ggtitle("Single-gene ROC (training)") + labs(tag = "E") +
  theme(plot.title = element_text(size = 8.5))
grDevices::png(file.path(fig_dir, "Fig6_ml_screening_rev.png"),
               width = 7.0, height = 6.0, units = "in", res = 300)
print((mp1 | mp2 | mp3) / (mp4 | mp5))
grDevices::dev.off()
cat("Fig6 done\n")

# ============================================================
# Fig 7A: validation ROC across cohorts (frozen canonical coefficients)
# ============================================================
fm <- fread("results/frozen/final_model_coefficients.csv")
bhf <- setNames(fm$canonical_deviance, fm$term)
scf <- function(e) {
  hg <- intersect(hub, rownames(e))
  Xv <- t(zrows(e[hg, , drop = FALSE]))
  as.vector(Xv %*% bhf[hg] + bhf[["(Intercept)"]])
}
roc_list <- list(Training = roc(y, scf(train_expr), quiet = TRUE))
for (acc in c("GSE30528", "GSE30529")) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  p <- read.csv(paste0("data/", acc, "_pheno.csv"), stringsAsFactors = FALSE)
  keep <- p$group %in% c("DKD", "Control")
  lbl <- ifelse(acc == "GSE30528", "GSE30528 (glomeruli)", "GSE30529 (tubules)")
  roc_list[[lbl]] <- roc(as.integer(p$group[keep] == "DKD"),
                         scf(e[, keep, drop = FALSE]), quiet = TRUE)
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
roc_list[["GSE142025 (all DN)"]] <- roc(as.integer(p14$group[k3] != "Control"),
                                        scf(e14[, k3, drop = FALSE]), quiet = TRUE)
roc_list[["GSE142025 (early)"]] <- roc(as.integer(p14$group[k4] == "Early_DN"),
                                       scf(e14[, k4, drop = FALSE]), quiet = TRUE)
e10 <- readRDS("data/GSE104948_expr_symbol.rds")
p10 <- read.csv("data/GSE104948_pheno.csv", stringsAsFactors = FALSE)
roc_list[["GSE104948"]] <- roc(as.integer(p10$group == "DN"), scf(e10), quiet = TRUE)
oi8 <- c("#000000", "#E69F00", "#56B4E9", "#0072B2", "#D55E00", "#F0E442", "#CC79A7", "#999999")
lt8 <- c("solid", "solid", "solid", "longdash", "longdash", "dashed", "dashed", "dotdash")
nm <- names(roc_list)
names(roc_list) <- sprintf("%s (AUC %.3f)", nm, sapply(roc_list, function(r) as.numeric(auc(r))))
cols <- setNames(rep(oi8, length.out = length(roc_list)), names(roc_list))
lts <- setNames(rep(lt8, length.out = length(roc_list)), names(roc_list))
p7 <- ggroc(roc_list, linewidth = .65) +
  geom_abline(slope = 1, intercept = 1, linetype = 3) +
  scale_color_manual(values = cols) + scale_linetype_manual(values = lts) +
  guides(linetype = "none") +
  theme_bw(base_size = 9) + coord_equal() +
  ggtitle("PDK2+PDK4 signature ROC across cohorts") +
  labs(color = NULL, tag = "A") +
  theme(legend.position = "bottom", legend.text = element_text(size = 6.5))
grDevices::png(file.path(fig_dir, "Fig7A_validation_roc_rev.png"),
               width = 7.0, height = 6.2, units = "in", res = 300)
print(p7)
grDevices::dev.off()
cat("Fig7A done\n")

# ============================================================
# Fig 7B: calibration GSE142025 (recalibrated score)
# ============================================================
y3 <- as.integer(p14$group[k3] != "Control")
s3 <- scf(e14[, k3, drop = FALSE])
g1 <- glm(y3 ~ s3, family = binomial)
p_hat <- fitted(g1)
png(file.path(fig_dir, "Fig7B_calibration_GSE142025.png"),
    width = 1400, height = 1400, res = 240)
vp <- val.prob(p_hat, y3)
title("Calibration - GSE142025 (recalibrated)")
dev.off()
cat("Fig7B done; Eavg", vp["Eavg"], "Emax", vp["Emax"], "Brier", vp["Brier"], "\n")

# ============================================================
# Fig 7C: DCA GSE142025 (manual net benefit, Okabe-Ito)
# ============================================================
nb <- function(thr, pred, outcome) {
  tp <- sum(pred >= thr & outcome == 1); fp <- sum(pred >= thr & outcome == 0)
  n <- length(outcome)
  tp / n - fp / n * (thr / (1 - thr))
}
thresholds <- seq(0.05, 0.95, by = 0.05)
nb_sig <- sapply(thresholds, nb, pred = p_hat, outcome = y3)
nb_all <- sapply(thresholds, function(t) { prev <- mean(y3); prev - (1 - prev) * t / (1 - t) })
df_dca <- rbind(
  data.frame(thr = thresholds, nb = nb_sig, curve = "Signature"),
  data.frame(thr = thresholds, nb = nb_all, curve = "Treat all"),
  data.frame(thr = thresholds, nb = 0, curve = "Treat none"))
p7c <- ggplot(df_dca, aes(thr, nb, color = curve, linetype = curve)) +
  geom_line(linewidth = .8) +
  scale_color_manual(values = c(Signature = "#0072B2", `Treat all` = "#D55E00",
                                `Treat none` = "#000000")) +
  scale_linetype_manual(values = c(Signature = "solid", `Treat all` = "longdash",
                                   `Treat none` = "dashed")) +
  theme_bw(base_size = 11) +
  xlab("Threshold probability") + ylab("Net benefit") +
  ggtitle("Decision curve analysis - GSE142025") + labs(tag = "C") +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(fig_dir, "Fig7C_dca_GSE142025.png"), p7c,
       width = 7.0, height = 5.5, dpi = 300, device = png_dev)
cat("Fig7C done\n")

# ============================================================
# Fig 8A: immune cell boxplots (BH18 stars); Fig 8B: correlation heatmap
# ============================================================
it18 <- fread("results/immune_tests_18cells.csv")
it2 <- it18[p.adj_BH18 < 0.05]
sig_cells <- it2[order(p.adj_BH18)]$cell
immune_ss <- as.data.frame(fread("results/immune_ssgsea.csv"))
rownames(immune_ss) <- immune_ss$sample; immune_ss$sample <- NULL
grp <- train_ph$group
dfi <- do.call(rbind, lapply(sig_cells, function(cn)
  data.frame(cell = cn, score = as.numeric(immune_ss[[cn]]), group = grp)))
dfi$cell <- factor(dfi$cell, levels = sig_cells)
starvec <- setNames(ifelse(it2$p.adj_BH18 < .001, "***",
                    ifelse(it2$p.adj_BH18 < .01, "**", "*")), it2$cell)
ann8 <- data.frame(cell = factor(sig_cells, levels = sig_cells), lab = starvec[sig_cells],
                   score = max(dfi$score, na.rm = TRUE) * 1.08, group = "Control")
p8a <- ggplot(dfi, aes(cell, score, fill = group)) +
  geom_boxplot(outlier.size = .4, linewidth = .3) +
  geom_text(data = ann8, aes(x = cell, y = score, label = lab), size = 2.9, inherit.aes = FALSE) +
  scale_fill_manual(values = c(Control = "#2166AC", DN = "#B2182B")) +
  theme_bw(base_size = 10) + xlab(NULL) + ylab("ssGSEA score") +
  ggtitle(sprintf("Immune infiltration (significant cells, n = %d of 18; BH across 18 tests)",
                  length(sig_cells))) +
  labs(tag = "A") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
ggsave(file.path(fig_dir, "Fig8A_immune_box_rev.png"), p8a, width = 9.2, height = 5.2,
       dpi = 300, device = png_dev)

immune_corr <- read.csv("results/immune_gene_corr.csv", stringsAsFactors = FALSE,
                        check.names = FALSE)
stopifnot(is.character(immune_corr$gene))
cm <- as.matrix(immune_corr[, setdiff(colnames(immune_corr), c("gene", "Normal_mucosa"))])
rownames(cm) <- immune_corr$gene
col_fun2 <- colorRamp2(c(-.6, 0, .6), c("#2166AC", "white", "#B2182B"))
png(file.path(fig_dir, "Fig8B_corr_rev.png"),
    width = 2400, height = 1600, res = 250)
draw(Heatmap(cm, name = "Spearman rho", col = col_fun2,
             row_names_gp = gpar(fontsize = 9.5), column_names_gp = gpar(fontsize = 9.5),
             column_title = "Candidate LRGs x immune cells (GSE96804, all samples; 18 cells)"))
grid.text("B", x = unit(2.6, "mm"), y = unit(1, "npc") - unit(2.6, "mm"),
          gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()
cat("Fig8 done; BH18-significant cells:", length(sig_cells), "\n")

# ============================================================
# Fig 9: snRNA localization (UMAP + detection-rate bars)
# ============================================================
s12 <- fread("results/snrna_pdk_celltype.csv")
pdk4_rates <- s12[gene == "PDK4", .(type = celltype, rate = as.numeric(pct_exp))]
pdk2_rates <- s12[gene == "PDK2", .(type = celltype, rate = as.numeric(pct_exp))]
ord <- pdk4_rates[order(-rate)]$type
pdk4_rates[, type := factor(type, levels = ord)]
pdk2_rates[, type := factor(type, levels = ord)]
merged <- readRDS("results/snrna/seurat_annotated.rds")
umap_emb <- Embeddings(merged, "umap")
expr_data <- GetAssayData(merged, assay = "RNA", layer = "data")
pdk4_vec <- as.numeric(expr_data["PDK4", ])
ct <- as.character(merged$celltype)
umap_df <- data.frame(UMAP1 = umap_emb[, 1], UMAP2 = umap_emb[, 2],
                      celltype = ct, PDK4_expr = pdk4_vec)
np1 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = celltype)) +
  geom_point(size = 0.2, alpha = 0.6) +
  scale_color_manual(values = c(
    "#000000", "#E69F00", "#56B4E9", "#0072B2", "#D55E00",
    "#F0E442", "#CC79A7", "#009E73", "#999999", "#8B4513",
    "#FF6600", "#32CD32", "#8A2BE2", "#20B2AA", "#DC143C"), name = NULL) +
  theme_bw(base_size = 8) + ggtitle("Cell types") + labs(tag = "A") +
  theme(legend.position = "right", legend.text = element_text(size = 5.5))
np2 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = PDK4_expr)) +
  geom_point(size = 0.2, alpha = 0.7) +
  scale_color_gradient(low = "grey92", high = "#E41A1C", name = "PDK4") +
  theme_bw(base_size = 8) + ggtitle("PDK4 expression") + labs(tag = "B") +
  theme(legend.position = "right", legend.key.width = unit(6, "pt"))
np3 <- ggplot(pdk4_rates, aes(x = type, y = rate)) +
  geom_col(fill = "#E41A1C", width = 0.7) + theme_bw(base_size = 8) +
  xlab(NULL) + ylab("Detection rate (%)") + ggtitle("PDK4 detection rate") +
  labs(tag = "C") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  ylim(0, 75)
np4 <- ggplot(pdk2_rates, aes(x = type, y = rate)) +
  geom_col(fill = "#377EB8", width = 0.7) + theme_bw(base_size = 8) +
  xlab(NULL) + ylab("Detection rate (%)") + ggtitle("PDK2 detection rate") +
  labs(tag = "D") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  ylim(0, 10)
ggsave(file.path(fig_dir, "Fig9_scrna_localization_rev.png"), (np1 | np2) / (np3 | np4),
       width = 7.0, height = 5.5, dpi = 300, device = png_dev)
cat("Fig9 done\n")

cat("\nFIGURES_DONE\n")

# 13_figures_rev.R ----
# 按review_round1.md修订全部图件
suppressWarnings(suppressMessages({
  library(data.table); library(glmnet); library(pROC); library(ggplot2)
  library(ggrepel); library(ComplexHeatmap); library(circlize)
  library(VennDiagram); library(patchwork); library(grid)
}))
set.seed(42)
fig_dir <- "results/figures"
png_win <- function(filename, width, height, units = 'in', res = 300, ...)
  grDevices::png(filename, width = width, height = height, units = units, res = res, type = 'windows', ...)
save_fig <- function(p, name, w, h, dpi = 300)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = dpi, device = png_win)

train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
allg <- fread("results/A_all_genes_limma.csv")
deg <- fread("results/A_DEGs.csv")
lrg <- fread("data/LRG_geneset.csv")
immune_corr <- read.csv("results/immune_gene_corr.csv", stringsAsFactors = FALSE, check.names = FALSE)

# ---------- Fig2 AB: 面板字母 + 阈值线标注 + 数量 ----------
pca_df <- data.frame(PC1 = prcomp(t(train_expr))$x[, 1],
                     PC2 = prcomp(t(train_expr))$x[, 2], group = train_ph$group)
p1 <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 2.6, alpha = .9) + stat_ellipse(level = .9, linetype = 2, linewidth = .6) +
  scale_color_manual(values = c(Control = "#2166AC", DN = "#B2182B")) +
  theme_bw(base_size = 11) + labs(color = NULL, tag = "A") +
  ggtitle("Training cohort (GSE96804)")
vol <- data.frame(gene = allg$gene, logFC = allg$logFC, neglog10 = -log10(allg$adj.P.Val))
vol$sig <- ifelse(vol$logFC >= .5 & vol$neglog10 > -log10(.05), "Up",
           ifelse(vol$logFC <= -.5 & vol$neglog10 > -log10(.05), "Down", "NS"))
lb <- vol[vol$gene %in% c("PDK4", "PDK2", "ACSF2", "ME1") & vol$sig != "NS", ]
p2 <- ggplot(vol, aes(logFC, neglog10, color = sig)) +
  geom_point(size = .55, alpha = .45) +
  scale_color_manual(values = c(Up = "#B2182B", Down = "#2166AC", NS = "grey78")) +
  geom_hline(yintercept = -log10(.05), linetype = 3) +
  annotate("text", x = max(vol$logFC) * .55, y = -log10(.05) + 2.5,
           label = "adj.P = 0.05", size = 2.8) +
  geom_vline(xintercept = c(-.5, .5), linetype = 3) +
  geom_text_repel(data = lb, aes(label = gene), size = 2.7, max.overlaps = 20, show.legend = FALSE) +
  annotate("text", x = 4.2, y = 14, label = sprintf("Up\nn = %d", sum(vol$sig == "Up")),
           size = 3, color = "#B2182B") +
  annotate("text", x = -4.2, y = 14, label = sprintf("Down\nn = %d", sum(vol$sig == "Down")),
           size = 3, color = "#2166AC") +
  theme_bw(base_size = 11) + labs(tag = "B") +
  ggtitle(sprintf("Volcano (DEGs n = %d)", nrow(deg))) +
  theme(legend.position = "bottom")
save_fig(p1 + p2, "Fig2AB_pca_volcano_rev", 9.6, 4.4)

# ---------- Fig2C: 过滤非编码 + 面板字母 ----------
top50 <- deg[order(deg$adj.P.Val)][1:80]
top50 <- top50[!grepl("^(MIR|SNOR|RNU|LINC|RNA)", top50$gene)][1:50]
mat <- t(scale(t(train_expr[top50$gene, ])))
ha <- HeatmapAnnotation(Group = train_ph$group,
                        col = list(Group = c(Control = "#2166AC", DN = "#B2182B")))
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
png(file.path(fig_dir, "Fig2C_heatmap_rev.png"), type = "windows", width = 1800, height = 2400, res = 250)
draw(Heatmap(mat, name = "z-score", col = col_fun, top_annotation = ha,
             show_column_names = FALSE, row_names_gp = gpar(fontsize = 6),
             column_title = "Top 50 protein-coding DEGs (GSE96804)"))
grid.text("C", x = unit(1, "npc") - unit(2.6, "mm"), y = unit(2.6, "mm"),
          gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

# ---------- Fig3: 边距修复 + grey排除说明 ----------
wgc_mt <- fread("results/wgcna_module_trait.csv")
wgc_p <- wgc_mt[module != "grey"]  # 展示排除grey
png(file.path(fig_dir, "Fig3_wgcna_rev.png"), type = "windows", width = 2400, height = 1100, res = 240)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2) + .1, oma = c(0, 2, 0, 0))
suppressWarnings(suppressMessages({
  library(WGCNA)
}))
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
               main = "Module-trait correlations (grey shown in Suppl.)")
dev.off()

# ---------- Fig4: 面板字母 + GSEA NES/FDR ----------
ego <- fread("results/enrich_GO_BP.csv")
er <- fread("results/enrich_Reactome.csv")
gse <- fread("results/gsea_Reactome.csv")
mk_bar <- function(d, n, ttl, tag) {
  if (!hasName(d, "Count")) d[, Count := lengths(strsplit(geneID, "/"))]
  d <- utils::head(d[order(d$p.adjust)], n)
  d$lab <- make.names(substr(d$Description, 1, 52))
  d$lab <- factor(d$lab, levels = rev(d$lab))
  ggplot(d, aes(Count, lab, fill = p.adjust)) + geom_col(width = .7) +
    scale_fill_gradient(low = "#D32F2F", high = "#1976D2", name = "p.adj") +
    theme_bw(base_size = 9.5) + ylab(NULL) + ggtitle(ttl) + labs(tag = tag) +
    theme(plot.title = element_text(size = 10.5))
}
p1 <- mk_bar(ego, 10, "GO BP (top DEGs)", "A")
p2 <- mk_bar(er, 10, "Reactome (top DEGs)", "B")
gse2 <- copy(gse)
gse2$lab <- make.names(paste0(substr(gse2$Description, 1, 42),
                              "  NES=", sprintf("%.2f", gse2$NES)))
topg <- rbind(head(gse2[order(NES), ], 4), head(gse2[order(-NES), ], 4))
topg$lab <- factor(topg$lab, levels = topg$lab[order(topg$NES)])
p3 <- ggplot(topg, aes(NES, lab, fill = p.adjust)) + geom_col(width = .7) +
  scale_fill_gradient(low = "#D32F2F", high = "#1976D2", name = "FDR") +
  theme_bw(base_size = 9) + ylab(NULL) + ggtitle("GSEA (Reactome, top by |NES|)") +
  labs(tag = "C")
save_fig(p1 / p2 / p3 + plot_layout(heights = c(1, 1, .75)), "Fig4_enrichment_rev", 9, 10)

# ---------- Fig5A: Venn（177键口径+标签修复；空别名不再计入键集） ----------
lrg_alias <- lrg$alias[!is.na(lrg$alias) & nzchar(trimws(lrg$alias))]
lrg_keys5 <- toupper(unique(c(lrg$symbol, lrg_alias)))
venn_list <- list(
  `LRG keys (177)` = lrg_keys5,
  `DEGs (2645)` = deg$gene,
  `Modules (4210)` = fread("results/wgcna_sigmodule_genes.csv")$gene
)
png(file.path(fig_dir, "Fig5A_venn_rev.png"), type = "windows", width = 1800, height = 1500, res = 250)
v <- venn.diagram(venn_list, filename = NULL, fill = c("#F4B183", "#9DC3E6", "#C5E0B4"),
                  alpha = .55, cex = 1.05, cat.cex = .95, margin = .12,
                  fontfamily = "sans", cat.fontfamily = "sans",
                  main = "Candidate LRG selection", main.cex = 1.15)
grid.draw(v)
dev.off()

# ---------- Fig5B: 标题修正 ----------
cand <- fread("results/A_candidates.csv")
cd <- cand[!is.na(cand$logFC)]
cd$lab <- ifelse(cd$adj.P.Val < 0.001, "***", ifelse(cd$adj.P.Val < 0.01, "**",
         ifelse(cd$adj.P.Val < 0.05, "*", "")))
p5 <- ggplot(cd, aes(reorder(gene, logFC), logFC, fill = LRG_sources)) +
  geom_col(width = .72, color = "grey25", linewidth = .18) +
  geom_text(aes(label = lab, y = logFC + .14 * sign(logFC)), size = 3) +
  coord_flip() + theme_bw(base_size = 10) +
  labs(x = NULL, y = "log2FC (DN vs Control)", fill = "LRG source",
       tag = "B", title = "Differentially expressed candidate LRGs (n = 11)") +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))
save_fig(p5, "Fig5B_candidates_rev", 8, 6.2)

# ---------- Fig6: 字母+星标 ----------
genes <- intersect(fread("results/A_candidates.csv")$gene, rownames(train_expr))
zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
Xdf <- data.frame(t(zrows(train_expr[genes, , drop = FALSE])), check.names = FALSE)
y <- as.integer(train_ph$group == "DN")
library(e1071); library(xgboost)
cvfit <- cv.glmnet(as.matrix(Xdf), y, family = "binomial", alpha = 1, nfolds = 10,
                   type.measure = "deviance")
p1 <- ggplot(data.frame(lambda = cvfit$lambda, dev = cvfit$cvm,
                        lo = cvfit$cvlo, up = cvfit$cvup),
             aes(lambda, dev)) + geom_errorbar(aes(ymin = lo, ymax = up), color = "grey70") +
  geom_point(color = "#B2182B", size = .9) + scale_x_log10() +
  geom_vline(xintercept = cvfit$lambda.1se, linetype = 2) +
  theme_bw(base_size = 10) + ylab("Binomial deviance") +
  ggtitle("LASSO CV (lambda.1se)") + labs(tag = "A")
sizes <- 2:12
folds <- sample(rep(1:10, length.out = nrow(Xdf)))
svm_rfe_rank <- function(X, yy) {
  feats <- colnames(X); ranked <- character(0)
  while (length(feats) > 1) {
    m <- svm(X[, feats, drop = FALSE], yy, kernel = "linear", scale = FALSE)
    w <- (t(m$coefs) %*% m$SV)^2
    wr <- as.vector(w); names(wr) <- feats
    rem <- names(which.min(wr)); ranked <- c(rem, ranked); feats <- setdiff(feats, rem)
  }
  c(ranked, feats)
}
acc <- vapply(sizes, function(k) {
  a <- vapply(1:10, function(f) {
    tr <- folds != f
    rk2 <- svm_rfe_rank(as.matrix(Xdf[tr, , drop = FALSE]), factor(y[tr]))
    sel <- tail(rk2, k)
    m <- svm(as.matrix(Xdf[tr, sel, drop = FALSE]), factor(y[tr]), kernel = "linear")
    mean(as.integer(as.character(predict(m, as.matrix(Xdf[!tr, sel, drop = FALSE])))) == y[!tr])
  }, numeric(1)); mean(a)
}, numeric(1))
p2 <- ggplot(data.frame(k = sizes, acc = acc), aes(k, acc)) +
  geom_line(color = "#2166AC") + geom_point(size = 2) +
  geom_vline(xintercept = sizes[which.max(acc)], linetype = 2) +
  theme_bw(base_size = 10) + xlab("Number of features") + ylab("CV accuracy") +
  ggtitle("SVM-RFE") + labs(tag = "B") +
  ylim(min(acc) - .02, 1)
dtrain <- xgb.DMatrix(as.matrix(Xdf), label = y)
bst <- xgboost(x = as.matrix(Xdf), y = factor(y, levels = c(0, 1)),
               max_depth = 3, learning_rate = 0.1, objective = "binary:logistic",
               nrounds = 50, nthread = 2)
imp <- xgb.importance(model = bst); imp$Gain <- imp$Gain * 100
p3 <- ggplot(head(imp, 10), aes(reorder(Feature, Gain), Gain)) +
  geom_col(fill = "#4CAF50", width = .7) + coord_flip() +
  theme_bw(base_size = 10) + xlab(NULL) + ylab("Gain (%)") +
  ggtitle("XGBoost importance") + labs(tag = "C")
hub <- readLines("results/hub_genes.txt")
hubexp <- do.call(rbind, lapply(hub, function(g)
  data.frame(gene = g, expr = as.numeric(train_expr[g, ]), group = train_ph$group)))
allgk <- allg[gene %in% hub]
pv <- setNames(allgk$adj.P.Val, allgk$gene)
star <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "ns")))
hubexp$gene <- factor(hubexp$gene, levels = hub)
ann <- data.frame(gene = factor(hub, levels = hub),
                  lab = star(pv[hub]), expr = max(hubexp$expr) * 1.06, group = "DN")
p4 <- ggplot(hubexp, aes(gene, expr, fill = group)) +
  geom_boxplot(outlier.size = .5, linewidth = .3) +
  geom_text(data = ann, aes(x = gene, y = expr, label = lab), size = 3.2, inherit.aes = FALSE) +
  scale_fill_manual(values = c(Control = "#2166AC", DN = "#B2182B")) +
  theme_bw(base_size = 10) + xlab(NULL) + ylab("Expression (log2)") +
  ggtitle("Hub genes (training; vs adj.P") + labs(tag = "D") +
  theme(legend.position = "bottom") +
  annotate("text", x = .5, y = max(hubexp$expr) * 1.12, hjust = 0,
           label = "ns / * / ** / *** vs adj.P 0.05/0.01/0.001", size = 2.6)
p4 <- p4 + ggtitle("Hub gene expression (training)")
rocs <- lapply(hub, function(g) roc(y, Xdf[[g]], quiet = TRUE))
names(rocs) <- hub
p5 <- ggroc(rocs, linewidth = .7) +
  geom_abline(slope = 1, intercept = 1, linetype = 3) +
  scale_color_brewer(palette = "Set1") + theme_bw(base_size = 10) +
  coord_equal() + ggtitle("Single-gene ROC (training)") + labs(color = NULL, tag = "E") +
  theme(legend.position = "bottom", legend.box.margin = margin(0, 0, 0, 0, "pt"))
save_fig((p1 | p2 | p3) / (p4 | p5), "Fig6_ml_screening_rev", 13.5, 7.8)

# ---------- Fig7A: 多队列ROC（字母A） ----------
Xh <- data.frame(t(zrows(train_expr[hub, , drop = FALSE])), check.names = FALSE)
cv2 <- cv.glmnet(as.matrix(Xh), y, family = "binomial", alpha = 1, nfolds = 10)
bh <- as.matrix(coef(cv2, s = "lambda.1se"))
sc <- function(e, yy) {
  hg <- intersect(hub, rownames(e))
  Xv <- data.frame(t(zrows(e[hg, , drop = FALSE])), check.names = FALSE)
  s <- rep(bh["(Intercept)", 1], ncol(e))
  for (g in hg) s <- s + bh[g, 1] * Xv[[g]]
  s
}
roc_list <- list(Training = roc(y, as.vector(as.matrix(Xh) %*% bh[-1, , drop = FALSE] +
                                              bh["(Intercept)", 1]), quiet = TRUE))
for (acc in c("GSE30528", "GSE30529")) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  p <- read.csv(paste0("data/", acc, "_pheno.csv"), stringsAsFactors = FALSE)
  keep <- p$group %in% c("DKD", "Control")
  roc_list[[acc]] <- roc(as.integer(p$group[keep] == "DKD"),
                         sc(e[, keep, drop = FALSE]), quiet = TRUE)
}
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])),
                  stringsAsFactors = FALSE)
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
rownames(e14) <- toupper(rownames(e14))
k3 <- p14$group %in% c("Advanced_DN", "Early_DN", "Control")
roc_list[["GSE142025 (all DN)"]] <- roc(as.integer(p14$group[k3] != "Control"),
                                        sc(e14[, k3, drop = FALSE]), quiet = TRUE)
k4 <- p14$group %in% c("Early_DN", "Control")
roc_list[["GSE142025 (early)"]] <- roc(as.integer(p14$group[k4] == "Early_DN"),
                                       sc(e14[, k4, drop = FALSE]), quiet = TRUE)
e10 <- readRDS("data/GSE104948_expr_symbol.rds")
p10 <- read.csv("data/GSE104948_pheno.csv", stringsAsFactors = FALSE)
k5 <- p10$group %in% c("DN", "Control")
roc_list[["GSE104948"]] <- roc(as.integer(p10$group[k5] == "DN"),
                               sc(e10[, k5, drop = FALSE]), quiet = TRUE)
p7 <- ggroc(roc_list, linewidth = .75) + geom_abline(slope = 1, intercept = 1, linetype = 3) +
  coord_equal() + theme_bw(base_size = 11) +
  scale_color_brewer(palette = "Dark2") + labs(color = NULL, tag = "A") +
  ggtitle("PDK2+PDK4 signature ROC across cohorts") +
  theme(legend.position = c(.42, .18),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.key.width = unit(22, "pt"))
save_fig(p7, "Fig7A_validation_roc_rev", 7.2, 6.6)

# ---------- Fig8: 剔除Normal_mucosa + 星标 + 修复基因名 ----------
immune_ss <- fread("results/immune_ssgsea.csv")
it <- fread("results/immune_tests.csv")
it2 <- it[cell != "Normal_mucosa" & p.adj < 0.05]
sig_cells <- it2[order(p.adj)]$cell
iss <- as.data.frame(immune_ss); rownames(iss) <- iss$sample; iss$sample <- NULL
grp <- train_ph$group
dfi <- do.call(rbind, lapply(sig_cells, function(cn)
  data.frame(cell = cn, score = as.numeric(iss[[cn]]), group = grp)))
dfi$cell <- factor(dfi$cell, levels = sig_cells)
starvec <- setNames(ifelse(it2$p.adj < .001, "***", ifelse(it2$p.adj < .01, "**", "*")), it2$cell)
ann8 <- data.frame(cell = factor(sig_cells, levels = sig_cells),
                   lab = starvec[sig_cells],
                   score = max(dfi$score, na.rm = TRUE) * 1.08, group = "Control")
p1 <- ggplot(dfi, aes(cell, score, fill = group)) +
  geom_boxplot(outlier.size = .4, linewidth = .3) +
  geom_text(data = ann8, aes(x = cell, y = score, label = lab), size = 2.9, inherit.aes = FALSE) +
  scale_fill_manual(values = c(Control = "#2166AC", DN = "#B2182B")) +
  theme_bw(base_size = 10) + xlab(NULL) + ylab("ssGSEA score") +
  ggtitle(sprintf("Immune infiltration (significant cells, n = %d of 18)", length(sig_cells))) +
  labs(tag = "A") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
save_fig(p1, "Fig8A_immune_box_rev", 9, 5.2)

cm <- as.matrix(immune_corr[, setdiff(colnames(immune_corr), "gene")])
rownames(cm) <- immune_corr$gene
col_fun2 <- colorRamp2(c(-.6, 0, .6), c("#2166AC", "white", "#B2182B"))
png(file.path(fig_dir, "Fig8B_corr_rev.png"), type = "windows", width = 2400, height = 1600, res = 250)
draw(Heatmap(cm, name = "Spearman rho", col = col_fun2,
             row_names_gp = gpar(fontsize = 7.5), column_names_gp = gpar(fontsize = 7.5),
             column_title = "Candidate LRGs x immune cells (DN, GSE96804)"))
grid.text("B", x = unit(1, "npc") - unit(2.6, "mm"), y = unit(2.6, "mm"),
          gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

cat("\nFIGURES_REV_DONE\n")

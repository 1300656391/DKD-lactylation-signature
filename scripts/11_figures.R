# 11_figures.R ----
# 论文A统一出图（出版级，300dpi PNG + PDF）
# Fig1 流程图 / Fig2 QC+DEG / Fig3 WGCNA / Fig4 富集 / Fig5 候选基因
# Fig6 ML筛选+hub / Fig7 模型+验证ROC / Fig8 免疫浸润
# 运行: Rscript 11_figures.R

suppressWarnings(suppressMessages({
  library(data.table); library(limma); library(glmnet); library(e1071)
  library(xgboost); library(pROC); library(rms); library(dcurves)
  library(ggplot2); library(ggrepel); library(ComplexHeatmap); library(circlize)
  library(VennDiagram); library(patchwork); library(scales)
}))
set.seed(42)
fig_dir <- "results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
png_win <- function(filename, width, height, units = 'in', res = 300, ...) grDevices::png(filename, width = width, height = height, units = units, res = res, type = 'windows', ...)
save_fig <- function(p, name, w, h, dpi = 300) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = dpi, device = png_win)
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
}

# ===== 数据重载（确定性重算signature） =====
train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
allg <- fread("results/A_all_genes_limma.csv")
deg <- fread("results/A_DEGs.csv")
cand <- fread("results/A_candidates.csv")
lrg <- fread("data/LRG_geneset.csv")
wgc_mt <- fread("results/wgcna_module_trait.csv")
immune_ss <- fread("results/immune_ssgsea.csv")
immune_tests <- fread("results/immune_tests.csv")
immune_corr <- fread("results/immune_gene_corr.csv")

genes <- intersect(cand$gene, rownames(train_expr))
zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
Xdf <- data.frame(t(zrows(train_expr[genes, , drop = FALSE])), check.names = FALSE)
y <- as.integer(train_ph$group == "DN")

# 复现S10筛选（同种子）
cvfit <- cv.glmnet(as.matrix(Xdf), y, family = "binomial", alpha = 1,
                   nfolds = 10, type.measure = "deviance")
cf <- as.matrix(coef(cvfit, s = "lambda.1se"))
lasso_genes <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")
yf <- factor(y, levels = c(0, 1), labels = c("Ctrl", "DN"))
svm_rfe_rank <- function(X, yy) {
  feats <- colnames(X); ranked <- character(0)
  while (length(feats) > 1) {
    m <- svm(X[, feats, drop = FALSE], yy, kernel = "linear", scale = FALSE)
    w <- (t(m$coefs) %*% m$SV)^2
    wr <- as.vector(w); names(wr) <- feats
    rem <- names(which.min(wr)); ranked <- c(rem, ranked)
    feats <- setdiff(feats, rem)
  }
  c(ranked, feats)
}
rk <- svm_rfe_rank(as.matrix(Xdf), yf)
sizes <- 2:min(12, ncol(Xdf)); folds <- sample(rep(1:10, length.out = nrow(Xdf)))
acc <- vapply(sizes, function(k) {
  a <- vapply(1:10, function(f) {
    tr <- folds != f
    rk2 <- svm_rfe_rank(as.matrix(Xdf[tr, , drop = FALSE]), yf[tr])
    sel <- tail(rk2, k)
    m <- svm(as.matrix(Xdf[tr, sel, drop = FALSE]), yf[tr], kernel = "linear")
    mean(predict(m, as.matrix(Xdf[!tr, sel, drop = FALSE])) == yf[!tr])
  }, numeric(1)); mean(a)
}, numeric(1))
best_k <- sizes[which.max(acc)]; svm_genes <- tail(rk, best_k)
dtrain <- xgb.DMatrix(as.matrix(Xdf), label = y)
mk_folds <- function(yy, k = 5) {
  f0 <- sample(rep(1:k, length.out = sum(yy == 0))); f1 <- sample(rep(1:k, length.out = sum(yy == 1)))
  lapply(1:k, function(i) c(which(yy == 0)[f0 == i], which(yy == 1)[f1 == i]))
}
xpar <- list(objective = "binary:logistic", eval_metric = "auc", max_depth = 3, eta = 0.1)
cvb <- xgb.cv(data = dtrain, params = xpar, nrounds = 300, folds = mk_folds(y),
              early_stopping_rounds = 20, verbose = 0)
bst <- xgboost(x = as.matrix(Xdf), y = factor(y, levels = c(0, 1)),
               max_depth = 3, learning_rate = 0.1, objective = "binary:logistic",
               nrounds = max(cvb$best_iteration, 10), nthread = 2)
imp <- xgb.importance(model = bst)
xgb_genes <- head(imp$Feature, 10)
hub <- Reduce(intersect, list(lasso_genes, svm_genes, xgb_genes))
if (length(hub) < 2) {
  cnt <- table(c(lasso_genes, svm_genes, xgb_genes)); hub <- names(cnt)[cnt >= 2]
}
hub <- intersect(hub, colnames(Xdf))
cat("复现hub:", paste(hub, collapse = ", "), "\n")
stopifnot(identical(sort(hub), sort(readLines("results/hub_genes.txt"))))

# 惩罚模型系数
Xh <- as.matrix(Xdf[, hub, drop = FALSE])
cv2 <- cv.glmnet(Xh, y, family = "binomial", alpha = 1, nfolds = 10)
bh <- as.matrix(coef(cv2, s = "lambda.1se"))
score_tr <- as.vector(Xh %*% bh[-1, , drop = FALSE] + bh["(Intercept)", 1])
roc_tr <- roc(y, score_tr, quiet = TRUE)

# ===== Fig1 流程图 =====
pdf(NULL); dev.off()  # 清理设备
png(file.path(fig_dir, "Fig1_flowchart.png"), type = "windows", width = 2400, height = 1500, res = 220)
par(mar = c(1, 1, 2, 1)); plot.new(); plot.window(c(0, 100), c(0, 62), xaxs = "i", yaxs = "i")
box_at <- function(x, y, w, h, label, fill = "white", cex = 0.95) {
  roundrect <- function(x, y, w, h, col) {
    symbols(x, y, rectangles = matrix(c(w, h), 1), inches = FALSE, add = TRUE,
            fg = "grey25", bg = col, lwd = 1.4)
  }
  roundrect(x, y, w, h, fill)
  text(x, y, label, cex = cex)
}
arr <- function(x1, y1, x2, y2) arrows(x1, y1, x2, y2, lwd = 1.6, length = 0.12)
title("Study workflow", cex.main = 1.5, adj = 0)
box_at(16, 55, 26, 7, "GEO cohorts\nTrain: GSE96804 (DN 41 / Ctrl 20)\nValidation: 4 cohorts (n=97)", "#DCEEFB")
box_at(50, 55, 24, 7, "Lactylation-related genes\nGO + Reactome + literature\n(171 genes)", "#DCEEFB")
box_at(84, 55, 24, 7, "DEGs (limma)\n2,645 genes\n+ WGCNA modules", "#DCEEFB")
arr(29, 55, 38, 55); arr(62, 55, 72, 55)
box_at(50, 41, 44, 7, "Candidate LRDEGs: 23  (LRG \u2229 (DEG \u222a disease modules))", "#FDF2D0")
arr(84, 51.5, 60, 44.8); arr(50, 51.5, 50, 44.7)
box_at(22, 30, 34, 8, "3 ML algorithms\nLASSO | SVM-RFE | XGBoost\nconsensus", "#E8F5E9")
box_at(70, 30, 40, 8, "Signature: PDK2 + PDK4 (penalized logistic)\nNomogram | Calibration | DCA", "#E8F5E9")
arr(42, 37, 36, 34); arr(58, 37, 64, 34)
box_at(22, 15, 34, 8, "ROC validation\nGSE30528 / GSE30529 / GSE142025\nincl. early-stage DN | GSE104948", "#FCE4EC")
box_at(70, 15, 40, 8, "Immune infiltration (ssGSEA, 19 cell types)\nGene-cell correlation", "#FCE4EC")
arr(22, 26, 22, 19); arr(70, 26, 70, 19)
box_at(46, 4, 50, 6, "Lactylation-related diagnostic signature for DKD", "#D1C4E9", 1.05)
arr(32, 11, 42, 7); arr(66, 11, 58, 7)
dev.off()

# ===== Fig2 QC + DEG =====
pca_df <- data.frame(PC1 = prcomp(t(train_expr))$x[, 1],
                     PC2 = prcomp(t(train_expr))$x[, 2], group = train_ph$group)
p1 <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 2.6, alpha = .9) + stat_ellipse(level = .9, linetype = 2, linewidth = .6) +
  scale_color_manual(values = c(Control = "#2166AC", DN = "#B2182B")) +
  theme_bw(base_size = 11) + labs(color = NULL, title = "Training cohort (GSE96804)")
vol <- data.frame(gene = allg$gene, logFC = allg$logFC, neglog10 = -log10(allg$adj.P.Val))
vol$sig <- ifelse(vol$logFC >= .5 & vol$neglog10 > -log10(.05), "Up",
           ifelse(vol$logFC <= -.5 & vol$neglog10 > -log10(.05), "Down", "NS"))
lb <- vol[vol$gene %in% c("PDK4", "PDK2", "ME1", "ACSF2", "HIST1H4B") & vol$sig != "NS", ]
p2 <- ggplot(vol, aes(logFC, neglog10, color = sig)) +
  geom_point(size = .55, alpha = .45) +
  scale_color_manual(values = c(Up = "#B2182B", Down = "#2166AC", NS = "grey78")) +
  geom_hline(yintercept = -log10(.05), linetype = 3) + geom_vline(xintercept = c(-.5, .5), linetype = 3) +
  geom_text_repel(data = lb, aes(label = gene), size = 2.7, max.overlaps = 20, show.legend = FALSE) +
  theme_bw(base_size = 11) + labs(title = sprintf("Volcano (DEGs n=%d)", nrow(deg))) +
  theme(legend.position = "bottom")
save_fig(p1 + p2, "Fig2AB_pca_volcano", 9, 4.4)

# 热图 top50
top50 <- deg[order(deg$adj.P.Val)][1:50]
mat <- train_expr[top50$gene, ]
mat <- t(scale(t(mat)))
ha <- HeatmapAnnotation(Group = train_ph$group,
                        col = list(Group = c(Control = "#2166AC", DN = "#B2182B")))
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
png(file.path(fig_dir, "Fig2C_heatmap.png"), type = "windows", width = 1800, height = 2400, res = 250)
draw(Heatmap(mat, name = "z-score", col = col_fun, top_annotation = ha,
             show_column_names = FALSE, row_names_gp = gpar(fontsize = 6),
             column_title = "Top 50 DEGs (GSE96804)"))
dev.off()

# ===== Fig3 WGCNA =====
suppressWarnings(suppressMessages(library(WGCNA)))
mads_f <- apply(train_expr, 1, mad)
dat_f <- t(train_expr[order(-mads_f)[1:5000], ])
sft_f <- pickSoftThreshold(dat_f, powerVector = c(1:10, seq(12, 20, 2)),
                           networkType = "signed", verbose = 0, corFnc = "cor")
png(file.path(fig_dir, "Fig3_wgcna.png"), type = "windows", width = 2400, height = 1000, res = 240)
par(mfrow = c(1, 2))
fi <- sft_f$fitIndices
plot(fi$Power, -sign(fi$slope) * fi$SFT.R.sq, xlab = "Soft threshold (power)",
     ylab = "Scale-free topology fit R2",
     main = "WGCNA soft threshold", pch = 16, col = "#B2182B")
abline(h = .85, col = "red", lty = 2)
text(fi$Power, -sign(fi$slope) * fi$SFT.R.sq, labels = fi$Power, cex = .7, pos = 4)
mt2 <- matrix(wgc_mt$V1, ncol = 1, dimnames = list(wgc_mt$module, "DN"))
lab <- matrix(paste0(signif(wgc_mt$V1, 2), "\n(", signif(wgc_mt$p, 1), ")"),
              ncol = 1, dimnames = list(wgc_mt$module, "DN"))
labeledHeatmap(Matrix = mt2, yLabels = wgc_mt$module, ySymbols = wgc_mt$module,
               colorLabels = FALSE, colors = blueWhiteRed(50),
               setStdMargins = FALSE, cex.text = .55, xLabels = "DN",
               xLabelsPosition = "bottom", textMatrix = lab,
               main = "Module-trait correlations")
dev.off()

# ===== Fig4 富集 =====
ego <- fread("results/enrich_GO_BP.csv")
er <- fread("results/enrich_Reactome.csv")
gse <- fread("results/gsea_Reactome.csv")
mk_enrich_bar <- function(d, n = 10, ttl) {
  d <- utils::head(d[order(d$p.adjust)], n)
  d$lab <- make.names(substr(d$Description, 1, 58))
  d$lab <- factor(d$lab, levels = rev(d$lab))
  ggplot(d, aes(Count, lab, fill = p.adjust)) + geom_col(width = .7) +
    scale_fill_gradient(low = "#D32F2F", high = "#1976D2", name = "p.adjust") +
    theme_bw(base_size = 9.5) + ylab(NULL) + ggtitle(ttl) +
    theme(plot.title = element_text(size = 10.5))
}
p1 <- mk_enrich_bar(ego[, Count := NA], 10, "GO BP")
if (!hasName(ego, "Count")) ego[, Count := NA]
# Count列不存在时用geneID计数
if (all(is.na(ego$Count))) {
  ego$Count <- lengths(strsplit(ego$geneID, "/"))
  p1 <- mk_enrich_bar(ego, 10, "GO BP")
}
if (!hasName(er, "Count")) er[, Count := lengths(strsplit(er$geneID, "/"))]
p2 <- mk_enrich_bar(er, 10, "Reactome pathway")
gse2 <- copy(gse)
gse2$lab <- make.names(substr(gse2$Description, 1, 50))
topg <- rbind(head(gse2[order(NES), ], 5), head(gse2[order(-NES), ], 5))
topg$lab <- factor(topg$lab, levels = topg$lab[order(topg$NES)])
p3 <- ggplot(topg, aes(NES, lab, fill = p.adjust)) + geom_col(width = .7) +
  scale_fill_gradient(low = "#D32F2F", high = "#1976D2") +
  theme_bw(base_size = 9) + ylab(NULL) + ggtitle("GSEA (Reactome)")
save_fig(p1 / p2 / p3 + plot_layout(heights = c(1, 1, .8)), "Fig4_enrichment", 9, 10)

# ===== Fig5 候选基因 =====
venn_list <- list(
  `Lactylation-related (171)` = unique(c(lrg$symbol, na.omit(lrg$alias))),
  DEGs = deg$gene,
  `WGCNA modules` = fread("results/wgcna_sigmodule_genes.csv")$gene
)
png(file.path(fig_dir, "Fig5A_venn.png"), type = "windows", width = 1600, height = 1400, res = 250)
v <- venn.diagram(venn_list, filename = NULL, fill = c("#F4B183", "#9DC3E6", "#C5E0B4"),
                  alpha = .55, cex = 1.1, cat.cex = 1, margin = .06, fontfamily = "sans",
                  main = "Candidate gene selection", main.cex = 1.2)
grid::grid.draw(v)
dev.off()

# 候选基因表达点图（DN vs Ctrl logFC + significance）
cd <- cand[!is.na(cand$logFC)]
cd$lab <- ifelse(cd$adj.P.Val < 0.001, "**", ifelse(cd$adj.P.Val < 0.05, "*", ""))
p5 <- ggplot(cd, aes(reorder(gene, logFC), logFC, fill = LRG_sources)) +
  geom_col(width = .72, color = "grey25", linewidth = .18) +
  geom_text(aes(label = lab, y = logFC + .12 * sign(logFC)), size = 3) +
  coord_flip() + theme_bw(base_size = 10) +
  labs(x = NULL, y = "log2FC (DN vs Control)",
       fill = "LRG source", title = "Candidate LRDEGs (n=23)") +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))
save_fig(p5, "Fig5B_candidates", 8, 6.4)

# ===== Fig6 ML筛选 =====
p1 <- ggplot(data.frame(lambda = cvfit$lambda, dev = cvfit$cvm,
                        lo = cvfit$cvlo, up = cvfit$cvup),
             aes(lambda, dev)) + geom_errorbar(aes(ymin = lo, ymax = up), color = "grey70") +
  geom_point(color = "#B2182B", size = .9) + scale_x_log10() +
  geom_vline(xintercept = cvfit$lambda.1se, linetype = 2) +
  theme_bw(base_size = 10) + ylab("Binomial deviance") + ggtitle("LASSO cross-validation")
p2 <- ggplot(data.frame(k = sizes, acc = acc), aes(k, acc)) +
  geom_line(color = "#2166AC") + geom_point(size = 2) +
  geom_vline(xintercept = best_k, linetype = 2) +
  annotate("text", x = best_k, y = min(acc) + .02, label = paste("best k =", best_k), size = 3) +
  theme_bw(base_size = 10) + xlab("Number of features") + ylab("CV accuracy (10-fold)") +
  ggtitle("SVM-RFE")
imp$Gain <- imp$Gain * 100
p3 <- ggplot(head(imp, 10), aes(reorder(Feature, Gain), Gain)) +
  geom_col(fill = "#4CAF50", width = .7) + coord_flip() +
  theme_bw(base_size = 10) + xlab(NULL) + ylab("Gain (%)") + ggtitle("XGBoost importance")
hubexp <- do.call(rbind, lapply(hub, function(g)
  data.frame(gene = g, expr = as.numeric(train_expr[g, ]), group = train_ph$group)))
p4 <- ggplot(hubexp, aes(gene, expr, fill = group)) +
  geom_boxplot(outlier.size = .5, linewidth = .3) +
  scale_fill_manual(values = c(Control = "#2166AC", DN = "#B2182B")) +
  theme_bw(base_size = 10) + xlab(NULL) + ylab("Expression (log2)") +
  ggtitle("Hub genes (training)") + theme(legend.position = "bottom")
rocs <- lapply(hub, function(g) roc(y, Xdf[[g]], quiet = TRUE))
p5 <- ggroc(setNames(rocs, hub), linewidth = .7) +
  geom_abline(slope = 1, intercept = 1, linetype = 3) +
  scale_color_brewer(palette = "Set1") + theme_bw(base_size = 10) +
  coord_equal() + ggtitle("Single-gene ROC (training)") + labs(color = NULL) +
  theme(legend.position = "bottom")
save_fig((p1 | p2 | p3) / (p4 | p5), "Fig6_ml_screening", 13, 7.6)

# ===== Fig7 模型 + 验证 =====
dfm <- data.frame(y = y, Xdf[, hub, drop = FALSE])
dd <- datadist(dfm); options(datadist = "dd")
lrm_fit <- tryCatch(lrm(y ~ ., data = dfm), error = function(e) NULL)
png(file.path(fig_dir, "Fig7A_nomogram.png"), type = "windows", width = 2000, height = 1300, res = 240)
if (!is.null(lrm_fit)) {
  plot(nomogram(lrm_fit, fun = plogis, fun.at = c(.1, .25, .5, .75, .9), lp = FALSE),
       cex.axis = .65, cex.var = .8)
} else {
  plot(nomogram(lrm(y ~ score_tr), fun = plogis, fun.at = c(.1, .3, .5, .7, .9)))
}
dev.off()
png(file.path(fig_dir, "Fig7B_calibration.png"), type = "windows", width = 1400, height = 1400, res = 240)
cal <- tryCatch(calibrate(lrm_fit %||% lrm(y ~ score_tr), method = "boot", B = 1000),
                error = function(e) NULL)
if (!is.null(cal)) plot(cal, xlab = "Predicted probability", ylab = "Observed", subtitles = FALSE)
dev.off()
dfd <- data.frame(y = y, Signature = fitted(glm(y ~ score_tr, family = binomial)))
png(file.path(fig_dir, "Fig7C_dca.png"), type = "windows", width = 1500, height = 1300, res = 240)
tryCatch(print(plot(dca(y ~ Signature, data = dfd), smooth = FALSE)), error = function(e) NULL)
dev.off()

# 验证ROC（重算各队列score）
valid_sets <- list()
for (acc in c("GSE30528", "GSE30529")) {
  e <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  p <- read.csv(paste0("data/", acc, "_pheno.csv"), stringsAsFactors = FALSE)
  valid_sets[[acc]] <- list(e = e, p = p, case = "DKD", ctrl = "Control")
}
e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(
  gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])),
  stringsAsFactors = FALSE)
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
rownames(e14) <- toupper(rownames(e14))
e10 <- readRDS("data/GSE104948_expr_symbol.rds")
p10 <- read.csv("data/GSE104948_pheno.csv", stringsAsFactors = FALSE)

cohort_score <- function(e, p, case_groups, ctrl_groups) {
  keep <- p$group %in% c(case_groups, ctrl_groups)
  e <- e[, keep, drop = FALSE]; p2 <- p[keep, ]
  hg <- intersect(hub, rownames(e))
  Xv <- data.frame(t(zrows(e[hg, , drop = FALSE])), check.names = FALSE)
  s <- rep(bh["(Intercept)", 1], ncol(e))
  for (g in hg) s <- s + bh[g, 1] * Xv[[g]]
  list(score = s, y = as.integer(p2$group %in% case_groups))
}
roc_list <- list(Training = roc_tr)
s1 <- cohort_score(valid_sets$GSE30528$e, valid_sets$GSE30528$p, "DKD", "Control")
roc_list[["GSE30528 (glomeruli)"]] <- roc(s1$y, s1$score, quiet = TRUE)
s2 <- cohort_score(valid_sets$GSE30529$e, valid_sets$GSE30529$p, "DKD", "Control")
roc_list[["GSE30529 (tubules)"]] <- roc(s2$y, s2$score, quiet = TRUE)
s3 <- cohort_score(e14, p14, c("Advanced_DN", "Early_DN"), "Control")
roc_list[["GSE142025 (all DN)"]] <- roc(s3$y, s3$score, quiet = TRUE)
s4 <- cohort_score(e14, p14, "Early_DN", "Control")
roc_list[["GSE142025 (early DN)"]] <- roc(s4$y, s4$score, quiet = TRUE)
s5 <- cohort_score(e10, p10, "DN", "Control")
roc_list[["GSE104948"]] <- roc(s5$y, s5$score, quiet = TRUE)
p7 <- ggroc(roc_list, linewidth = .75) + geom_abline(slope = 1, intercept = 1, linetype = 3) +
  coord_equal() + theme_bw(base_size = 11) +
  scale_color_brewer(palette = "Dark2") + labs(color = NULL) +
  ggtitle("Signature ROC across cohorts") + theme(legend.position = c(.65, .22),
                                                  legend.background = element_rect(fill = "white", color = "grey80"))
save_fig(p7, "Fig7D_validation_roc", 7, 6.4)

# ===== Fig8 免疫 =====
iss <- as.data.frame(immune_ss)
rownames(iss) <- iss$sample; iss$sample <- NULL
grp <- train_ph$group
sig_cells <- immune_tests[p.adj < 0.05]$cell
dfi <- do.call(rbind, lapply(sig_cells, function(cn)
  data.frame(cell = cn, score = as.numeric(iss[[cn]]), group = grp)))
dfi$cell <- factor(dfi$cell, levels = sig_cells)
p1 <- ggplot(dfi, aes(cell, score, fill = group)) +
  geom_boxplot(outlier.size = .4, linewidth = .3) +
  scale_fill_manual(values = c(Control = "#2166AC", DN = "#B2182B")) +
  theme_bw(base_size = 10) + xlab(NULL) + ylab("ssGSEA score") +
  ggtitle(sprintf("Immune infiltration (%d significant cell types)", length(sig_cells))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
cm <- as.matrix(immune_corr[, -1]); rownames(cm) <- immune_corr$gene
col_fun2 <- colorRamp2(c(-.6, 0, .6), c("#2166AC", "white", "#B2182B"))
png(file.path(fig_dir, "Fig8B_corr.png"), type = "windows", width = 2200, height = 1500, res = 250)
draw(Heatmap(cm, name = "Spearman \u03c1", col = col_fun2,
             row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
             column_title = "Candidate genes \u00d7 immune cells"))
dev.off()
save_fig(p1, "Fig8A_immune_box", 9, 5)

cat("\nFIGURES_DONE\n")


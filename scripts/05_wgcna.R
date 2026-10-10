# 05_wgcna.R ----
# S5: WGCNA 共表达网络与疾病相关模块（训练集 GSE96804）
# 输入: data/GSE96804_expr_symbol.rds + data/GSE96804_pheno.csv
# 输出: results/wgcna_modules.csv（基因-模块-模块隶属度）
#       results/wgcna_module_trait.csv（模块-性状相关）
#       results/figures/Fig3_wgcna_draft.png
# 运行: Rscript 05_wgcna.R

suppressWarnings(suppressMessages({
  library(WGCNA)
  library(ggplot2)
  library(data.table)
}))
options(stringsAsFactors = FALSE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

expr <- readRDS("data/GSE96804_expr_symbol.rds")
ph   <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
stopifnot(identical(ph$sample, colnames(expr)))

# 方差过滤：MAD top 5000
mads <- apply(expr, 1, mad)
dat <- t(expr[order(-mads)[1:min(5000, nrow(expr))], ])
cat("WGCNA输入:", nrow(dat), "样本 x", ncol(dat), "基因\n")

# 软阈值
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(dat, powerVector = powers, networkType = "signed",
                         verbose = 0, corFnc = "cor")
sp <- sft$powerEstimate
if (is.na(sp)) sp <- 6
cat("选定软阈值 power =", sp, "\n")

# 建网
net <- blockwiseModules(
  dat, power = sp, networkType = "signed", TOMType = "signed",
  minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  maxBlockSize = 6000, verbose = 0)
moduleColors <- labels2colors(net$colors)
cat("模块数:", length(unique(moduleColors)), "\n")

# 模块-性状相关
trait <- as.numeric(ph$group == "DN")
MEs <- orderMEs(moduleEigengenes(dat, moduleColors)$eigengenes)
mod_trait <- as.data.frame(cor(MEs, trait, use = "p"))
mod_trait$p <- corPvalueStudent(mod_trait[, 1], nrow(dat))
mod_trait$module <- sub("^ME", "", rownames(mod_trait))
fwrite(mod_trait, "results/wgcna_module_trait.csv")
print(mod_trait[order(-abs(mod_trait[, 1])), ])

# 基因-模块-隶属度表
kme <- signedKME(dat, MEs)
out <- data.frame(gene = colnames(dat), module = moduleColors)
for (m in colnames(kme)) out[[m]] <- kme[, m]
fwrite(out, "results/wgcna_modules.csv")

# 疾病显著模块的基因清单
sig_mods <- mod_trait$module[abs(mod_trait[, 1]) > 0.3 & mod_trait$p < 0.05]
sig_mods <- setdiff(sig_mods, "grey")
cat("\n疾病相关模块(|r|>0.3且p<0.05):", paste(sig_mods, collapse = ", "), "\n")
sig_genes <- out$gene[out$module %in% sig_mods]
fwrite(data.frame(gene = sig_genes), "results/wgcna_sigmodule_genes.csv")
cat("疾病相关模块基因数:", length(sig_genes), "\n")

# 草稿图：软阈值曲线+模块-性状热图
png("results/figures/Fig3_wgcna_draft.png", width = 1400, height = 600, res = 130)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold", ylab = "Scale-free R2", type = "n",
     main = "Soft threshold")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = .8)
abline(h = 0.85, col = "red")
labeledHeatmap(Matrix = mod_trait[, 1, drop = FALSE],
               yLabels = mod_trait$module,
               ySymbols = mod_trait$module,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = paste0(signif(mod_trait[, 1], 2), "\n(", signif(mod_trait$p, 1), ")"),
               setStdMargins = FALSE, cex.text = 0.6,
               main = "Module-trait (DN)")
dev.off()

cat("\nWGCNA_DONE\n")

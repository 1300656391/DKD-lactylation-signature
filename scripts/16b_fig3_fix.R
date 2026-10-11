# 16b_fig3_fix.R ----
# R10 visual-gate fix: Fig3 right-panel title truncated at right edge.
# Drop the in-figure parenthetical (legend already points to Table S8), re-render at 300 dpi.
suppressWarnings(suppressMessages({
  library(data.table); library(WGCNA)
}))
fig_dir <- "results/figures"
train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
wgc_mt <- fread("results/wgcna_module_trait.csv")
wgc_p <- wgc_mt[module != "grey"]

png(file.path(fig_dir, "Fig3_wgcna_rev.png"), type = "windows",
    width = 3000, height = 1375, res = 300)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2) + .1, oma = c(0, 2, 0, 0))
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
               main = "Module-trait correlations")
dev.off()
cat("FIG3_DONE title-fixed-300dpi\n")

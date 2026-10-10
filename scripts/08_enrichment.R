# 08_enrichment.R ----
# S8: GO/KEGG富集(DEGs) + GSEA(全基因排序)
# 输出: results/enrich_*.csv + results/gsea_*.csv + figures/Fig4_enrich_draft.png
# 运行: Rscript 08_enrichment.R

suppressWarnings(suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(enrichplot)
  library(ggplot2)
  library(data.table)
}))

deg <- fread("results/A_DEGs.csv")
allg <- fread("results/A_all_genes_limma.csv")
bg_genes <- allg$gene

eg <- bitr(deg$gene, fromType = "SYMBOL", toType = "ENTREZID",
           OrgDb = org.Hs.eg.db)
cat("映射到Entrez:", nrow(eg), "/", nrow(deg), "\n")
bg_eg <- bitr(bg_genes, fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)

# GO BP（本地）
ego <- enrichGO(gene = eg$ENTREZID, universe = bg_eg$ENTREZID,
                OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1,
                readable = TRUE)
fwrite(as.data.frame(ego), "results/enrich_GO_BP.csv")
cat("GO BP显著条目:", nrow(ego), "\n")

# Reactome通路（本地reactome.db，替代被墙的KEGG在线库）
er <- enrichPathway(gene = eg$ENTREZID, universe = bg_eg$ENTREZID,
                    organism = "human", pAdjustMethod = "BH",
                    pvalueCutoff = 0.05, readable = TRUE)
fwrite(as.data.frame(er), "results/enrich_Reactome.csv")
cat("Reactome显著条目:", nrow(er), "\n")

# GSEA（Reactome，按limma t值排序）
allg_t <- allg[!is.na(allg$t), ]
gl <- allg_t$t; names(gl) <- allg_t$gene
gl <- sort(gl, decreasing = TRUE)
gl <- na.omit(gl[!duplicated(names(gl))])
gl_eg <- bitr(names(gl), fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)
gl2 <- gl[gl_eg$SYMBOL]; names(gl2) <- gl_eg$ENTREZID
gl2 <- sort(gl2, decreasing = TRUE)
gser <- gsePathway(geneList = gl2, organism = "human", minGSSize = 10,
                   maxGSSize = 500, pvalueCutoff = 0.1, verbose = FALSE)
fwrite(as.data.frame(gser), "results/gsea_Reactome.csv")
cat("GSEA Reactome显著条目:", nrow(as.data.frame(gser)), "\n")

# 草稿图：GO/KEGG top条目
mk_bar <- function(obj, title, n = 10) {
  d <- as.data.frame(obj)
  if (nrow(d) == 0) return(ggplot() + annotate("text", x=1, y=1, label="no terms") + theme_void())
  d <- head(d, n); d$Description <- factor(d$Description, levels = rev(d$Description))
  ggplot(d, aes(Count, Description, fill = p.adjust)) +
    geom_col() + scale_fill_gradient(low = "red", high = "blue") +
    theme_bw(base_size = 10) + ggtitle(title)
}
p1 <- mk_bar(ego, "GO BP")
p2 <- mk_bar(er, "Reactome")
png("results/figures/Fig4_enrich_draft.png", width = 1500, height = 900, res = 130)
print(patchwork::wrap_plots(p1, p2, ncol = 1))
dev.off()

cat("\nENRICHMENT_DONE\n")

# 05_enrichment.R ----
# GO BP / Reactome over-representation and Reactome GSEA on the training DEGs,
# plus per-WGCNA-module GO BP enrichment (supplementary table S19).
# Outputs: results/enrich_GO_BP.csv, results/enrich_Reactome.csv,
#          results/gsea_Reactome.csv,
#          results/supplement/TableS19_module_enrichment.csv

suppressWarnings(suppressMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(ReactomePA)
  library(data.table)
}))
set.seed(42)
dir.create("results/supplement", showWarnings = FALSE, recursive = TRUE)

deg <- fread("results/A_DEGs.csv")
allg <- fread("results/A_all_genes_limma.csv")

eg <- bitr(deg$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("DEGs mapped to Entrez:", nrow(eg), "/", nrow(deg), "\n")
bg_eg <- bitr(allg$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO BP over-representation
ego <- enrichGO(gene = eg$ENTREZID, universe = bg_eg$ENTREZID,
                OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1,
                readable = TRUE)
fwrite(as.data.frame(ego), "results/enrich_GO_BP.csv")
cat("GO BP significant terms:", nrow(ego), "\n")

# Reactome pathways
er <- enrichPathway(gene = eg$ENTREZID, universe = bg_eg$ENTREZID,
                    organism = "human", pAdjustMethod = "BH",
                    pvalueCutoff = 0.05, readable = TRUE)
fwrite(as.data.frame(er), "results/enrich_Reactome.csv")
cat("Reactome significant terms:", nrow(er), "\n")

# GSEA on the full ranked gene list (limma t statistic)
allg_t <- allg[!is.na(allg$t), ]
gl <- allg_t$t; names(gl) <- allg_t$gene
gl <- sort(gl, decreasing = TRUE)
gl <- na.omit(gl[!duplicated(names(gl))])
gl_eg <- bitr(names(gl), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gl2 <- gl[gl_eg$SYMBOL]; names(gl2) <- gl_eg$ENTREZID
gl2 <- sort(gl2, decreasing = TRUE)
gser <- gsePathway(geneList = gl2, organism = "human", minGSSize = 10,
                   maxGSSize = 500, pvalueCutoff = 0.1, verbose = FALSE)
fwrite(as.data.frame(gser), "results/gsea_Reactome.csv")
cat("GSEA Reactome significant terms:", nrow(as.data.frame(gser)), "\n")

# ---- Per-module GO BP enrichment (S19), grey module excluded ----
mods <- fread("results/wgcna_modules.csv")[, .(gene, module)]
mods <- mods[module != "grey"]
cat("modules for S19:", paste(sort(unique(mods$module)), collapse = ", "), "\n")

map_ids <- function(genes) {
  suppressWarnings(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID)
}

all_rows <- list()
for (m in sort(unique(mods$module))) {
  genes <- mods[module == m]$gene
  ez <- map_ids(genes)
  if (length(ez) < 10) { cat(m, ": too few mapped genes, skipped\n"); next }
  egm <- suppressWarnings(enrichGO(gene = ez, OrgDb = org.Hs.eg.db,
                                   keyType = "ENTREZID", ont = "BP",
                                   pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.1, readable = TRUE))
  if (nrow(as.data.frame(egm)) == 0) { cat(m, ": no significant terms\n"); next }
  d <- as.data.frame(egm)
  d$module <- m
  d$module_n_mapped <- length(ez)
  all_rows[[length(all_rows) + 1]] <- d
  cat(m, ":", nrow(d), "terms; top:", d$Description[1], "\n")
}

out <- rbindlist(all_rows, fill = TRUE)
setcolorder(out, c("module", "module_n_mapped"))
setorder(out, module, p.adjust)
out <- out[, .(module, module_n_mapped, ID, Description, GeneRatio, BgRatio,
               pvalue, p.adjust, qvalue, Count, geneID)]
fwrite(out, "results/supplement/TableS19_module_enrichment.csv")
cat("S19 written:", nrow(out), "rows,", length(unique(out$module)), "modules\n")

cat("\nENRICHMENT_DONE\n")

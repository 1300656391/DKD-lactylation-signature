# 20_scrna_localization.R ----
# Round-6 A1-(i)：GSE131882（早期DN snRNA-seq，3对照+3糖尿病，Wilson PNAS 2019）
# PDK2/PDK4 细胞类型定位 + 疾病对照比较
# 数据格式：双重gzip的rds；x$umicount$inex$all = ENSG行×细胞列的稀疏计数（intron+exon）
suppressWarnings(suppressMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(data.table)
  library(org.Hs.eg.db); library(patchwork)
}))
set.seed(42)

files <- c(
  CTRL_s1 = "data/scrna/GSM3823939_control.s1.dgecounts.rds.gz",
  CTRL_s2 = "data/scrna/GSM3823940_control.s2.dgecounts.rds.gz",
  CTRL_s3 = "data/scrna/GSM3823941_control.s3.dgecounts.rds.gz",
  DN_s1   = "data/scrna/GSM3823942_diabetes.s1.dgecounts.rds.gz",
  DN_s2   = "data/scrna/GSM3823943_diabetes.s2.dgecounts.rds.gz",
  DN_s3   = "data/scrna/GSM3823944_diabetes.s3.dgecounts.rds.gz")

# ENSG→SYMBOL 映射（mapIds，ENSEMBL键型）
sym_of <- function(genes) unname(mapIds(org.Hs.eg.db, keys = genes,
                                         column = "SYMBOL", keytype = "ENSEMBL"))

objs <- list()
for (nm in names(files)) {
  con <- gzcon(gzfile(files[[nm]], "rb"))
  x <- readRDS(con); close(con)
  cm <- x$umicount$inex$all
  rm(x); gc()
  cat(nm, "raw:", nrow(cm), "genes x", ncol(cm), "nuclei\n")
  sym <- sym_of(rownames(cm))
  keep <- !is.na(sym) & sym != ""
  cm <- cm[keep, ]
  sym <- sym[keep]
  # 同名基因UMI求和折叠：t(G) %*% counts（稀疏指示矩阵按组求和）
  u <- unique(sym); gi <- match(sym, u)
  cn <- colnames(cm)
  G <- Matrix::sparseMatrix(i = seq_along(sym), j = gi, x = 1,
                            dims = c(length(sym), length(u)))
  prod <- as(t(G) %*% cm, "dgCMatrix")
  rn <- make.unique(as.character(u))
  cn2 <- make.unique(as.character(cn))
  dimnames(prod) <- list(rn, cn2)
  storage.mode(prod@x) <- "numeric"
  cm <- prod
  cat(nm, "after symbol collapse:", nrow(cm), "genes x", ncol(cm), "nuclei; PDK4 present:",
      "PDK4" %in% rownames(cm), "\n")
  s <- CreateSeuratObject(counts = cm, project = nm, min.cells = 3, min.features = 200)
  s$group <- if (grepl("^CTRL", nm)) "Control" else "DN"
  s$sample <- nm
  s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
  s <- subset(s, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
  objs[[nm]] <- s
  rm(cm, s); gc()
  cat(nm, "QC-passed:", ncol(objs[[nm]]), "nuclei\n")
}
stopifnot(all(c("PDK2", "PDK4") %in% rownames(objs[[1]])))

# ---- 合并+标准流程（带断点缓存）----
cache <- "results/revision/GSE131882_merged_preannot.rds"
if (file.exists(cache)) {
  merged <- readRDS(cache)
  cat("从缓存载入 merged 对象\n")
} else {
  merged <- merge(objs[[1]], objs[-1], add.cell.ids = names(objs))
merged <- JoinLayers(merged)
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
merged <- FindVariableFeatures(merged, nfeatures = 2000, verbose = FALSE)
merged <- ScaleData(merged, vars.to.regress = "nCount_RNA", verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- FindNeighbors(merged, dims = 1:20, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.6, verbose = FALSE)
  merged <- RunUMAP(merged, dims = 1:20, verbose = FALSE)
  cat("总核数:", ncol(merged), "| 聚类数:", nlevels(merged$seurat_clusters), "\n")
  saveRDS(merged, cache)
  cat("缓存已保存:", cache, "\n")
}

# ---- 标记基因注释 ----
markers <- c(Podo = "NPHS1,NPHS2,WT1,PLCE1", PT = "LRP2,CUBN,VCAM1,SLC34A1",
             `TAL/LOH` = "UMOD,SLC12A1", DCT = "SLC12A3,TRPV5,CALB1",
             `CD-PC` = "AQP2,AQP3", `CD-IC` = "ATP6V1B1,SLC26A4",
             Endo = "PECAM1,VWF,KDR", Mes = "PDGFRB,DES,NOTCH3",
             Fibro = "COL1A1,COL3A1,DCN", MyoFibro = "ACTA2,TAGLN,MYH11",
             Macro = "CD68,CSF1R,LST1", `T/NK` = "CD3D,NKG7,CD3E", B = "MS4A1,CD79A",
             PEC = "CLDN1,ITGA8,PAX8")
genes_all <- unique(unlist(strsplit(unlist(markers), ",")))
pdk <- c("PDK2", "PDK4")
keep_genes <- intersect(unique(c(genes_all, pdk)), rownames(merged))
sub <- as.matrix(GetAssayData(merged, layer = "data")[keep_genes, ])
cl <- as.character(Idents(merged))
cl_means <- t(sapply(split(seq_along(cl), cl), function(ii) rowMeans(sub[, ii, drop = FALSE])))
# 按标签聚合标记基因得分（cluster × label），再取argmax
label_score <- sapply(names(markers), function(lb) {
  gs <- intersect(unlist(strsplit(markers[[lb]], ",")), colnames(cl_means))
  if (length(gs) == 0) return(rep(NA_real_, nrow(cl_means)))
  rowMeans(cl_means[, gs, drop = FALSE])
})
lab_cl <- colnames(label_score)[apply(label_score, 1, which.max)]
celltype <- unname(setNames(lab_cl, rownames(cl_means))[cl])
merged$celltype <- celltype
cat("\n细胞类型 x 分组:\n"); print(table(celltype, merged$group))
cat("\n聚类→标签:\n"); print(data.frame(cluster = rownames(cl_means), label = lab_cl))

# ---- PDK2/PDK4 每型统计 ----
grp <- merged$group; smp <- merged$sample
rows <- list()
for (g in pdk) {
  for (ct in sort(unique(celltype))) {
    idx <- which(celltype == ct)
    if (length(idx) < 20) next
    gi <- grp[idx]; si <- smp[idx]
    pb_dn <- sapply(unique(si[gi == "DN"]), function(s) mean(sub[g, idx][si == s]))
    pb_ct <- sapply(unique(si[gi == "Control"]), function(s) mean(sub[g, idx][si == s]))
    rows[[length(rows) + 1]] <- data.frame(
      gene = g, celltype = ct, n_nuclei = length(idx),
      pct_exp = round(mean(sub[g, idx] > 0) * 100, 1),
      mean_expr_all = round(mean(sub[g, idx]), 3),
      mean_expr_DN = round(mean(sub[g, idx][gi == "DN"]), 3),
      mean_expr_CTRL = round(mean(sub[g, idx][gi == "Control"]), 3),
      pseudobulk_log2FC = round(log2(mean(pb_dn) + 1e-3) - log2(mean(pb_ct) + 1e-3), 2),
      wilcox_p_3v3 = tryCatch(wilcox.test(pb_dn, pb_ct)$p.value, error = function(e) NA))
  }
}
pdk_stats <- setDT(rbindlist(rows)); setorder(pdk_stats, gene, -mean_expr_all)
fwrite(pdk_stats, "results/revision/M6_scrna_pdk_celltype.csv")
file.copy("results/revision/M6_scrna_pdk_celltype.csv",
          "manuscript/supplement/TableS12_scrna_pdk_celltype.csv", overwrite = TRUE)
print(pdk_stats)

# ---- Fig9 ----
Idents(merged) <- "celltype"
pA <- DimPlot(merged, reduction = "umap", label = TRUE, repel = TRUE, pt.size = .3) +
  ggtitle("Cell types (GSE131882 snRNA-seq)") + NoLegend()
pB <- FeaturePlot(merged, "PDK4", reduction = "umap", pt.size = .3, order = TRUE) +
  scale_colour_gradientn(colours = c("grey92", "#FDB863", "#D73027")) +
  ggtitle("PDK4") + NoLegend()
vtheme <- theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5))
pC <- VlnPlot(merged, "PDK4", group.by = "celltype", split.by = "group", split.plot = TRUE, pt.size = 0) & vtheme
pD <- VlnPlot(merged, "PDK2", group.by = "celltype", split.by = "group", split.plot = TRUE, pt.size = 0) & vtheme
png("results/figures/Fig9_scrna_pdk_localization.png", type = "windows",
    width = 3900, height = 3600, res = 300)
print(pA + pB + pC + pD + plot_layout(ncol = 2, heights = c(1, 1.25)))
dev.off()
saveRDS(merged, "results/revision/GSE131882_seurat_annotated.rds")
cat("\nSCRNA_DONE\n")

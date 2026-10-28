# 10_snrna_localization.R ----
# Single-nucleus RNA-seq localization of PDK2/PDK4 in early-DKD kidney
# (GSE131882: 3 control + 3 diabetic donors, Wilson PNAS 2019).
# Standard Seurat workflow, marker-gene cluster annotation, per-cell-type
# PDK2/PDK4 statistics with 3v3 pseudobulk comparison.
# Outputs:
#   results/snrna_pdk_celltype.csv (pipeline marker, also supplementary S12),
#   results/snrna/celltype_counts.csv, results/snrna/cluster_labels.csv,
#   results/snrna/seurat_annotated.rds (used by steps 11 and 13)

suppressWarnings(suppressMessages({
  library(Seurat); library(Matrix); library(org.Hs.eg.db); library(data.table)
}))
set.seed(42)
dir.create("results/snrna", showWarnings = FALSE, recursive = TRUE)
dir.create("results/supplement", showWarnings = FALSE, recursive = TRUE)

files <- c(
  CTRL_s1 = "data/scrna/GSM3823939_control.s1.dgecounts.rds.gz",
  CTRL_s2 = "data/scrna/GSM3823940_control.s2.dgecounts.rds.gz",
  CTRL_s3 = "data/scrna/GSM3823941_control.s3.dgecounts.rds.gz",
  DN_s1   = "data/scrna/GSM3823942_diabetes.s1.dgecounts.rds.gz",
  DN_s2   = "data/scrna/GSM3823943_diabetes.s2.dgecounts.rds.gz",
  DN_s3   = "data/scrna/GSM3823944_diabetes.s3.dgecounts.rds.gz")

sym_of <- function(genes) unname(mapIds(org.Hs.eg.db, keys = genes,
                                        column = "SYMBOL", keytype = "ENSEMBL"))

# Per-sample: ENSG -> symbol, collapse same-symbol UMIs (sparse group-sum),
# Seurat QC filters
objs <- list()
for (nm in names(files)) {
  con <- gzcon(gzfile(files[[nm]], "rb"))
  x <- readRDS(con); close(con)
  cm <- x$umicount$inex$all           # intron + exon counts
  rm(x); gc()
  cat(nm, "raw:", nrow(cm), "genes x", ncol(cm), "nuclei\n")
  sym <- sym_of(rownames(cm))
  keep <- !is.na(sym) & sym != ""
  cm <- cm[keep, ]; sym <- sym[keep]
  u <- unique(sym); gi <- match(sym, u)
  G <- Matrix::sparseMatrix(i = seq_along(sym), j = gi, x = 1,
                            dims = c(length(sym), length(u)))
  prod <- as(t(G) %*% cm, "dgCMatrix")
  dimnames(prod) <- list(make.unique(as.character(u)),
                         make.unique(as.character(colnames(cm))))
  storage.mode(prod@x) <- "numeric"
  s <- CreateSeuratObject(counts = prod, project = nm, min.cells = 3, min.features = 200)
  s$group <- if (grepl("^CTRL", nm)) "Control" else "DN"
  s$sample <- nm
  s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
  s <- subset(s, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
  objs[[nm]] <- s
  rm(cm, s, prod); gc()
  cat(nm, "QC-passed:", ncol(objs[[nm]]), "nuclei\n")
}
stopifnot(all(c("PDK2", "PDK4") %in% rownames(objs[[1]])))

# Merge + standard workflow (cached across re-runs)
cache <- "results/snrna/seurat_merged_preannot.rds"
if (file.exists(cache)) {
  merged <- readRDS(cache)
  cat("loaded merged object from cache\n")
} else {
  merged <- merge(objs[[1]], objs[-1], add.cell.ids = names(objs))
  merged <- JoinLayers(merged)
  merged <- NormalizeData(merged, normalization.method = "LogNormalize",
                          scale.factor = 1e4, verbose = FALSE)
  merged <- FindVariableFeatures(merged, nfeatures = 2000, verbose = FALSE)
  merged <- ScaleData(merged, vars.to.regress = "nCount_RNA", verbose = FALSE)
  merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
  merged <- FindNeighbors(merged, dims = 1:20, verbose = FALSE)
  merged <- FindClusters(merged, resolution = 0.6, verbose = FALSE)
  merged <- RunUMAP(merged, dims = 1:20, verbose = FALSE)
  cat("total nuclei:", ncol(merged), "| clusters:", nlevels(merged$seurat_clusters), "\n")
  saveRDS(merged, cache)
}

# Marker-gene annotation (cluster x label marker score, argmax)
markers <- c(Podo = "NPHS1,NPHS2,WT1,PLCE1", PT = "LRP2,CUBN,VCAM1,SLC34A1",
             `TAL/LOH` = "UMOD,SLC12A1", DCT = "SLC12A3,TRPV5,CALB1",
             `CD-PC` = "AQP2,AQP3", `CD-IC` = "ATP6V1B1,SLC26A4",
             Endo = "PECAM1,VWF,KDR", Mes = "PDGFRB,DES,NOTCH3",
             Fibro = "COL1A1,COL3A1,DCN", MyoFibro = "ACTA2,TAGLN,MYH11",
             Macro = "CD68,CSF1R,LST1", `T/NK` = "CD3D,NKG7,CD3E", B = "MS4A1,CD79A",
             PEC = "CLDN1,ITGA8,PAX8")
pdk <- c("PDK2", "PDK4")
keep_genes <- intersect(unique(c(unlist(strsplit(unlist(markers), ",")), pdk)),
                        rownames(merged))
sub <- as.matrix(GetAssayData(merged, layer = "data")[keep_genes, ])
cl <- as.character(Idents(merged))
cl_means <- t(sapply(split(seq_along(cl), cl), function(ii) rowMeans(sub[, ii, drop = FALSE])))
label_score <- sapply(names(markers), function(lb) {
  gs <- intersect(unlist(strsplit(markers[[lb]], ",")), colnames(cl_means))
  if (length(gs) == 0) return(rep(NA_real_, nrow(cl_means)))
  rowMeans(cl_means[, gs, drop = FALSE])
})
lab_cl <- colnames(label_score)[apply(label_score, 1, which.max)]
celltype <- unname(setNames(lab_cl, rownames(cl_means))[cl])
merged$celltype <- celltype
cat("cell types x group:\n"); print(table(celltype, merged$group))
fwrite(data.frame(cluster = rownames(cl_means), label = lab_cl),
       "results/snrna/cluster_labels.csv")
fwrite(as.data.frame(table(celltype = celltype, group = merged$group)),
       "results/snrna/celltype_counts.csv")

# Per-cell-type PDK2/PDK4 statistics (pseudobulk 3v3)
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
fwrite(pdk_stats, "results/snrna_pdk_celltype.csv")
file.copy("results/snrna_pdk_celltype.csv",
          "results/supplement/TableS12_scrna_pdk_celltype.csv", overwrite = TRUE)
print(pdk_stats)

saveRDS(merged, "results/snrna/seurat_annotated.rds")

cat("\nSNRNA_DONE\n")

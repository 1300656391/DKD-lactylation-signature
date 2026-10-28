# evidence_gain/2x_myeloid_ref/step1_annotate.R ----
# GSE209781 (human DKD scRNA atlas, 3 DKD + 3 normal kidneys, 6 10x samples)
# Goal: unpack -> standard Seurat QC/clustering (same parameters as the GSE131882 pipeline) -> cell-type annotation (13 kidney types + immune subtypes)
# Outputs are written only to evidence_gain/2x_myeloid_ref/
suppressWarnings(suppressMessages({ library(Seurat); library(Matrix) }))
set.seed(42)
OUT <- "evidence_gain/2x_myeloid_ref"
RAW <- file.path(OUT, "GSE209781_RAW.tar")
stopifnot(file.exists(RAW))
untar(RAW, exdir = file.path(OUT, "raw"))

members <- list.files(file.path(OUT, "raw"), pattern = "\\.tar\\.gz$", full.names = TRUE)
cat("samples:", length(members), "\n")
smeta <- sub("\\.tar\\.gz$", "", basename(members))
sgrp <- ifelse(grepl("DKD", smeta), "DN", "Control")

objs <- list()
for (k in seq_along(members)) {
  d <- file.path(OUT, "raw", paste0("x", k))
  untar(members[k], exdir = d)
  mfile <- list.files(d, pattern = "matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE)[1]
  stopifnot(!is.na(mfile))
  matdir <- dirname(mfile)
  cm <- Read10X(matdir)
  s <- CreateSeuratObject(counts = cm, project = smeta[k], min.cells = 3, min.features = 200)
  s$group <- sgrp[k]; s$sample <- smeta[k]
  s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
  s <- subset(s, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
  objs[[k]] <- s
  cat(smeta[k], "QC-passed:", ncol(s), "cells; genes:", nrow(s), "\n")
  unlink(d, recursive = TRUE)
}

merged <- merge(objs[[1]], objs[-1], add.cell.ids = smeta)
merged <- JoinLayers(merged)
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
merged <- FindVariableFeatures(merged, nfeatures = 2000, verbose = FALSE)
merged <- ScaleData(merged, vars.to.regress = "nCount_RNA", verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- FindNeighbors(merged, dims = 1:20, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.6, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:20, verbose = FALSE)
cat("total cells:", ncol(merged), "| clusters:", nlevels(merged$seurat_clusters), "\n")

# ---- Annotation: 13 kidney types (as in the main pipeline) + immune subtypes ----
markers <- c(Podo = "NPHS1,NPHS2,WT1,PLCE1", PT = "LRP2,CUBN,VCAM1,SLC34A1",
             `TAL/LOH` = "UMOD,SLC12A1", DCT = "SLC12A3,TRPV5,CALB1",
             `CD-PC` = "AQP2,AQP3", `CD-IC` = "ATP6V1B1,SLC26A4",
             Endo = "PECAM1,VWF,KDR", Mes = "PDGFRB,DES,NOTCH3",
             Fibro = "COL1A1,COL3A1,DCN", MyoFibro = "ACTA2,TAGLN,MYH11",
             Macro = "CD68,CSF1R,LST1,CD74", MonoNeu = "FCGR3A,S100A8,S100A9",
             DC = "FCER1A,CLEC10A,CD1C", `T/NK` = "CD3D,CD3E,NKG7,GNLY",
             B = "MS4A1,CD79A", Plasma = "MZB1,IGKC,JCHAIN", Mast = "TPSAB1,CPA3,KIT",
             PEC = "CLDN1,ITGA8,PAX8")
genes_all <- unique(unlist(strsplit(unlist(markers), ",")))
keep_genes <- intersect(genes_all, rownames(merged))
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
cat("\n== celltype x group ==\n"); print(table(celltype, merged$group))
cat("\n== cluster -> label ==\n"); print(data.frame(cluster = rownames(cl_means), label = lab_cl))

saveRDS(merged, file.path(OUT, "GSE209781_seurat_annotated.rds"))
write.csv(as.data.frame(table(celltype, group = merged$group)),
          file.path(OUT, "annot_type_counts.csv"), row.names = TRUE)
write.csv(data.frame(cluster = rownames(cl_means), label = lab_cl),
          file.path(OUT, "annot_cluster_map.csv"), row.names = FALSE)
cat("\nANNOTATE_DONE\n")

# evidence_gain/2_cross_deconv/step1_extract.R ----
# run with cwd = project root; relative paths only (ASCII-safe)
suppressWarnings(suppressMessages({ library(Seurat); library(Matrix) }))
set.seed(42)
OUT <- "evidence_gain/2_cross_deconv"

merged <- readRDS("results/revision/GSE131882_seurat_annotated.rds")
ct <- as.character(merged$celltype); gp <- as.character(merged$group)
cat("== snRNA celltype x group ==\n"); print(table(ct, gp))

dat <- GetAssayData(merged, layer = "data")
lin <- dat; lin@x <- expm1(lin@x)

tab  <- table(ct)
types <- names(tab)[tab >= 50]
cat("\nref types (>=50 nuclei):", paste(types, collapse = ", "), "\n")

S_lin <- sapply(types, function(t) Matrix::rowMeans(lin[, which(ct == t)]))
S_log <- sapply(types, function(t) Matrix::rowMeans(dat[, which(ct == t)]))
storage.mode(S_lin) <- "numeric"; storage.mode(S_log) <- "numeric"
colnames(S_lin) <- types; colnames(S_log) <- types

marker_list <- lapply(types, function(t) {
  idx  <- which(ct == t); oth <- which(ct != t)
  spec <- Matrix::rowMeans(dat[, idx]) - Matrix::rowMeans(dat[, oth])
  det  <- Matrix::rowMeans(dat[, idx] > 0)
  cand <- names(spec)[det > 0.10 & !grepl("^MT-", names(spec))]
  head(cand[order(-spec[cand])], 50)
})
names(marker_list) <- types
markers <- unique(unlist(marker_list))
cat("marker union:", length(markers), "\n")

write.csv(data.frame(celltype = types, n_nuclei = as.integer(tab[types])),
          file.path(OUT, "ref_type_counts.csv"), row.names = FALSE)
marker_map <- do.call(rbind, lapply(names(marker_list), function(t)
  data.frame(gene = marker_list[[t]], type = t, stringsAsFactors = FALSE)))
write.csv(marker_map, file.path(OUT, "ref_marker_list.csv"), row.names = FALSE)
write.csv(S_lin[markers, ], file.path(OUT, "ref_signature_marker_linear.csv"))
cat("ref matrix (marker):", nrow(S_lin[markers, ]), "genes x", ncol(S_lin[markers, ]), "types\n")

genes_need <- union(markers, c("PDK2", "PDK4", "PTGS2", "NFE2L2"))

cohorts <- list(
  list(id = "GSE96804",         expr = "data/GSE96804_expr_symbol.rds",           pheno = "data/GSE96804_pheno.csv"),
  list(id = "GSE30528",         expr = "data/GSE30528_expr_symbol.rds",           pheno = "data/GSE30528_pheno.csv"),
  list(id = "GSE30529",         expr = "data/GSE30529_expr_symbol.rds",           pheno = "data/GSE30529_pheno.csv"),
  list(id = "GSE142025",        expr = "data/GSE142025_expr.rds",                 pheno = "data/GSE142025_pheno.csv"),
  list(id = "GSE104948_A",      expr = "data/GSE104948-GPL22945_expr_symbol.rds", pheno = "data/GSE104948_pheno.csv"),
  list(id = "GSE104948_B",      expr = "data/GSE104948-GPL24120_expr_symbol.rds", pheno = "data/GSE104948_pheno.csv"),
  list(id = "GSE104948_pooled", expr = "data/GSE104948_expr_symbol.rds",          pheno = "data/GSE104948_pheno.csv")
)

for (co in cohorts) {
  e <- readRDS(co$expr)
  rn <- toupper(rownames(e))
  if (anyDuplicated(rn)) {
    mu <- rowMeans(as.matrix(e), na.rm = TRUE)
    keep <- tapply(seq_along(rn), rn, function(i) i[which.max(mu[i])])
    e <- e[as.integer(keep), ]; rn <- toupper(rownames(e))
  }
  rownames(e) <- rn
  ph <- read.csv(co$pheno, stringsAsFactors = FALSE)
  names(ph)[1] <- "sample"; names(ph)[2] <- "group"
  ph$group <- sub("^group:\\s*", "", ph$group)
  pref <- sub("_.*$", "", colnames(e))
  mi <- match(pref, ph$sample)
  if (all(is.na(mi))) {  # two-part column names (e.g., GSM4217781_A11A.txt): use the segment after the GSM id
    alt <- sub("\\.txt$", "", colnames(e))
    alt <- sub("^GSM[0-9]+_", "", alt)
    mi <- match(alt, ph$sample)
  }
  if (all(is.na(mi))) { cat("!! sample match failed:", co$id, "\n"); next }
  ok <- !is.na(mi)
  e <- e[, ok]; grp <- ph$group[mi[ok]]; smp <- ph$sample[mi[ok]]
  hit <- intersect(genes_need, rownames(e))
  cat(sprintf("[%s] samples=%d (DN=%d, Ctrl=%d) | marker hit=%d/%d | PDK2:%s PDK4:%s PTGS2:%s NFE2L2:%s | range[%.2f, %.2f]\n",
              co$id, ncol(e), sum(grp != "Control"), sum(grp == "Control"),
              length(intersect(markers, rownames(e))), length(markers),
              "PDK2" %in% rn, "PDK4" %in% rn, "PTGS2" %in% rn, "NFE2L2" %in% rn,
              min(as.matrix(e[hit, ])), max(as.matrix(e[hit, ]))))
  df <- as.data.frame(e[hit, ]); df$gene <- rownames(df)
  write.csv(df, file.path(OUT, paste0("bulk_", co$id, "_subset.csv")), row.names = FALSE)
  write.csv(data.frame(sample = smp, colname = colnames(e), group = grp),
            file.path(OUT, paste0("bulk_", co$id, "_pheno.csv")), row.names = FALSE)
}
cat("\nSTEP1_DONE\n")

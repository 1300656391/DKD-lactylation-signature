# evidence_gain/2x_myeloid_ref/step2_build_ref.R ----
# 合并参照：GSE131882 九个肾实质型（剔除其 B 型以避免共线）+ GSE209781 免疫型（>=50 细胞）
# 输出：合并 marker 参照矩阵、类型计数、（可行时）仅对照巨噬参照、以及按合并 marker 集重导出的 bulk 子集
suppressWarnings(suppressMessages({ library(Seurat); library(Matrix) }))
set.seed(42)
OUT <- "evidence_gain/2x_myeloid_ref"

kid <- readRDS("results/revision/GSE131882_seurat_annotated.rds")
imm <- readRDS(file.path(OUT, "GSE209781_seurat_annotated.rds"))

type_means <- function(obj, types_wanted) {
  ct <- as.character(obj$celltype)
  dat <- GetAssayData(obj, layer = "data")
  lin <- dat; lin@x <- expm1(lin@x)
  res_lin <- sapply(types_wanted, function(t) Matrix::rowMeans(lin[, which(ct == t)]))
  res_log <- sapply(types_wanted, function(t) Matrix::rowMeans(dat[, which(ct == t)]))
  storage.mode(res_lin) <- "numeric"; storage.mode(res_log) <- "numeric"
  list(lin = res_lin, log = res_log, n = table(ct)[types_wanted])
}

ktab <- table(as.character(kid$celltype))
kid_types <- names(ktab)[ktab >= 50]   # GSE131882 全部 10 型（含 B：GSE209781 无 B 簇，无共线）
IMM_WHITELIST <- c("Macro", "MonoNeu", "DC", "T/NK", "B", "Plasma", "Mast")  # 仅并入免疫型，肾型仍锚定 GSE131882
itab <- table(as.character(imm$celltype))
imm_types <- intersect(names(itab)[itab >= 50], IMM_WHITELIST)
cat("kidney types (GSE131882):", paste(kid_types, collapse = ", "), "\n")
cat("immune types (GSE209781):", paste(imm_types, collapse = ", "), "\n")

K <- type_means(kid, kid_types); I2 <- type_means(imm, imm_types)
stopifnot("Macro" %in% imm_types)

# 组合 log 均值（基因并集）
allg <- union(rownames(K$log), rownames(I2$log))
comb_log <- matrix(NA_real_, length(allg), length(kid_types) + length(imm_types),
                   dimnames = list(allg, c(kid_types, imm_types)))
comb_log[rownames(K$log), kid_types] <- K$log
comb_log[rownames(I2$log), imm_types] <- I2$log
comb_log[is.na(comb_log)] <- 0

# 每型 top50 特异标记（in - out，全类型联合计算），排除线粒体
marker_list <- lapply(colnames(comb_log), function(t) {
  spec <- comb_log[, t] - rowMeans(comb_log[, setdiff(colnames(comb_log), t), drop = FALSE])
  cand <- names(spec)[!grepl("^MT-", names(spec)) & !is.na(spec)]
  head(cand[order(-spec[cand])], 50)
})
names(marker_list) <- colnames(comb_log)
markers <- unique(unlist(marker_list))
cat("combined marker union:", length(markers), "\n")

# 线性参照矩阵
comb_lin <- matrix(0, length(markers), ncol(comb_log), dimnames = list(markers, colnames(comb_log)))
comb_lin[rownames(K$lin)[rownames(K$lin) %in% markers], kid_types] <-
  K$lin[rownames(K$lin)[rownames(K$lin) %in% markers], kid_types, drop = FALSE]
comb_lin[rownames(I2$lin)[rownames(I2$lin) %in% markers], imm_types] <-
  I2$lin[rownames(I2$lin)[rownames(I2$lin) %in% markers], imm_types, drop = FALSE]
write.csv(comb_lin, file.path(OUT, "ref_signature_combined_linear.csv"))

cnt <- rbind(data.frame(source = "GSE131882", type = kid_types, n = as.integer(K$n[kid_types])),
             data.frame(source = "GSE209781", type = imm_types, n = as.integer(I2$n[imm_types])))
write.csv(cnt, file.path(OUT, "ref_type_counts_combined.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(names(marker_list), function(t)
  data.frame(gene = marker_list[[t]], type = t, stringsAsFactors = FALSE))),
  file.path(OUT, "ref_marker_list_combined.csv"), row.names = FALSE)

# 宏噬的供体拆分（疾病/对照核数）+ 可行时仅对照巨噬 profile
igm <- imm$celltype == "Macro"
cat("\nMacro cells x group:\n"); print(table(imm$group[igm]))
n_macro_ctl <- sum(imm$group[igm] == "Control")
dat_i <- GetAssayData(imm, layer = "data"); lin_i <- dat_i; lin_i@x <- expm1(lin_i@x)
idx_ctl <- which(igm & imm$group == "Control")
if (n_macro_ctl >= 30) {
  prof <- Matrix::rowMeans(lin_i[, idx_ctl])
  write.csv(data.frame(gene = markers, Macro_control_only = as.numeric(prof[markers])),
            file.path(OUT, "ref_macro_control_profile.csv"), row.names = FALSE)
  cat("control-only macrophage profile exported (n=", n_macro_ctl, ")\n", sep = "")
} else {
  cat("control-only macrophage infeasible (n=", n_macro_ctl, " < 30) -> 仅用混合参照\n", sep = "")
}

# ---- 按合并 marker 集重导出 bulk 子集（只读 data/） ----
genes_need <- union(markers, c("PDK2", "PDK4", "PTGS2", "NFE2L2"))
cohorts <- list(
  list(id = "GSE96804",         expr = "data/GSE96804_expr_symbol.rds",           pheno = "data/GSE96804_pheno.csv"),
  list(id = "GSE30528",         expr = "data/GSE30528_expr_symbol.rds",           pheno = "data/GSE30528_pheno.csv"),
  list(id = "GSE30529",         expr = "data/GSE30529_expr_symbol.rds",           pheno = "data/GSE30529_pheno.csv"),
  list(id = "GSE142025",        expr = "data/GSE142025_expr.rds",                 pheno = "data/GSE142025_pheno.csv"),
  list(id = "GSE104948_A",      expr = "data/GSE104948-GPL22945_expr_symbol.rds", pheno = "data/GSE104948_pheno.csv"),
  list(id = "GSE104948_B",      expr = "data/GSE104948-GPL24120_expr_symbol.rds", pheno = "data/GSE104948_pheno.csv"),
  list(id = "GSE104948_pooled", expr = "data/GSE104948_expr_symbol.rds",          pheno = "data/GSE104948_pheno.csv"))
for (co in cohorts) {
  e <- readRDS(co$expr); rn <- toupper(rownames(e))
  if (anyDuplicated(rn)) {
    mu <- rowMeans(as.matrix(e), na.rm = TRUE)
    keep <- tapply(seq_along(rn), rn, function(i) i[which.max(mu[i])])
    e <- e[as.integer(keep), ]; rn <- toupper(rownames(e))
  }
  rownames(e) <- rn
  ph <- read.csv(co$pheno, stringsAsFactors = FALSE)
  names(ph)[1] <- "sample"; names(ph)[2] <- "group"; ph$group <- sub("^group:\\s*", "", ph$group)
  pref <- sub("_.*$", "", colnames(e)); mi <- match(pref, ph$sample)
  if (all(is.na(mi))) {
    alt <- sub("\\.txt$", "", colnames(e)); alt <- sub("^GSM[0-9]+_", "", alt); mi <- match(alt, ph$sample)
  }
  ok <- !is.na(mi); e <- e[, ok]
  hit <- intersect(genes_need, rownames(e))
  cat(sprintf("[%s] samples=%d | marker hit=%d/%d | Macro基因命中(CD68/CSF1R/LST1): %s\n",
              co$id, ncol(e), length(intersect(markers, rownames(e))), length(markers),
              paste(intersect(c("CD68","CSF1R","LST1"), rn), collapse = ",")))
  df <- as.data.frame(e[hit, ]); df$gene <- rownames(df)
  write.csv(df, file.path(OUT, paste0("bulk_", co$id, "_subset2.csv")), row.names = FALSE)
  write.csv(data.frame(sample = ph$sample[mi[ok]], colname = colnames(e), group = ph$group[mi[ok]]),
            file.path(OUT, paste0("bulk_", co$id, "_pheno2.csv")), row.names = FALSE)
}
cat("\nBUILD_REF_DONE\n")

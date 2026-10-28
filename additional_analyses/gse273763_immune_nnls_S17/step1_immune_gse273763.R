# evidence_gain/4_new_cohort/step1_immune_gse273763.R ----
# 执行①：免疫复现到 GSE273763（v6 方案）
#   ssGSEA+Bindea 18 细胞（主稿同法），双口径：提交者归一化矩阵（主）与 counts->log2CPM（敏感性）
#   预注册判定（锚点 B 节）：主终点=宏噬 DN↑ 且 BH<0.05（18 家族）；升格判据=主终点成立且全样本 PDK2x宏噬 rho<0、p<0.05
#   列名清洗：去除全部空白（含内部），assert 63 样本
#   附：导出 NNLS B 代理所需 marker 子集（与 2_cross_deconv 参照同基因集）
suppressWarnings(suppressMessages({ library(GSVA); library(data.table) }))
set.seed(42)
OUT <- "evidence_gain/4_new_cohort"

clean_cols <- function(m) { colnames(m) <- gsub("\\s+", "", colnames(m)); m }

# 两个矩阵
norm <- read.table(gzfile(file.path(OUT, "GSE273763_DKD_GRCH38_forGEO_July2026_v2.txt.gz")),
                   sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
cnt <- read.table(gzfile(file.path(OUT, "GSE273763_countsmatrix.txt.gz")),
                  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
norm <- clean_cols(norm); cnt <- clean_cols(cnt)
stopifnot(ncol(norm) == 63, ncol(cnt) == 63)
grp <- ifelse(grepl("^DKD", colnames(norm)) | grepl("_DKD$", colnames(norm)), "DKD", "Donor")
stopifnot(sum(grp == "DKD") == 50, sum(grp == "Donor") == 13)

# counts -> log2CPM（敏感性口径），ENSG -> symbol（用已映射文件）
sym_map <- cnt
# 已有符号化 counts（step: map_gse273763_symbols.R 产物）
cntsym <- read.csv(gzfile(file.path(OUT, "GSE273763_counts_symbol.csv.gz")), row.names = 1, check.names = FALSE)
cntsym <- clean_cols(cntsym)
cpm <- log2(t(t(cntsym) / colSums(cntsym)) * 1e6 + 1)
norm <- as.matrix(norm); storage.mode(norm) <- "numeric"
cpm <- as.matrix(cpm); storage.mode(cpm) <- "numeric"

# Bindea 18 签名（与 21_immune_replication.R 同法）
e <- new.env(); load("data/signature_collection.rda", envir = e)
sc <- get("signature_collection", envir = e)
cn <- grep("_Bindea_et_al$", names(sc), value = TRUE)
cn <- setdiff(cn, c("SW480_cancer_cells_Bindea_et_al"))
gsets <- sc[cn]; names(gsets) <- sub("_Bindea_et_al$", "", cn)
bindea_cells <- c("B_cells","T_cells","T_helper_cells","Tcm","Tem","Th1_cells","Th2_cells","TFH",
  "CD8_T_cells","Cytotoxic_cells","NK_cells","NK_CD56dim_cells","DC","iDC","Eosinophils",
  "Macrophages","Mast_cells","Neutrophils")
sig18 <- gsets[intersect(bindea_cells, names(gsets))]
sig18 <- lapply(sig18, function(g) unique(toupper(g)))
stopifnot(length(sig18) == 18)

run_ssgsea <- function(expr, tag) {
  rn <- toupper(rownames(expr)); rownames(expr) <- rn
  expr <- expr[!duplicated(rn), ]
  par <- ssgseaParam(exprData = expr, geneSets = sig18, normalize = TRUE)
  ss <- gsva(par, verbose = FALSE)
  ss <- as.data.frame(t(ss)); ss$group <- grp; ss$sample <- colnames(expr)
  fwrite(ss, file.path(OUT, paste0("gse273763_ssgsea_", tag, ".csv")))
  # 18 细胞组间检验（Wilcoxon + BH18），与主稿同法
  is_dn <- grp != "Donor"
  res <- do.call(rbind, lapply(bindea_cells, function(cc) {
    x <- ss[[cc]]
    data.frame(cell = cc, median_DN = median(x[is_dn]), median_Ctrl = median(x[!is_dn]),
               p = wilcox.test(x[is_dn], x[!is_dn])$p.value)
  }))
  res$p.adj_BH18 <- p.adjust(res$p, "BH")
  res$direction <- ifelse(res$median_DN > res$median_Ctrl, "up", "down")
  setorder(as.data.table(res), p.adj_BH18)
  fwrite(res, file.path(OUT, paste0("gse273763_immune_tests_", tag, ".csv")))
  cat("\n== ", tag, " 18-cell DN vs Donor (BH18) ==\n"); print(res)
  # PDK2/PDK4 x 全部细胞分数 Spearman（全样本 + DN 内）
  for (g in c("PDK2", "PDK4")) {
    ge <- as.numeric(expr[g, ])
    cr <- do.call(rbind, lapply(bindea_cells, function(cc)
      data.frame(gene = g, cell = cc,
                 rho_all = suppressWarnings(cor(ge, ss[[cc]], method = "spearman")),
                 rho_DN = suppressWarnings(cor(ge[is_dn], ss[[cc]][is_dn], method = "spearman")))))
    fwrite(cr, file.path(OUT, paste0("gse273763_pdk_corr_", tag, "_", g, ".csv")))
    cat("\n== ", tag, g, " x cells (rho_all / rho_DN) ==\n"); print(cr)
  }
}
run_ssgsea(norm, "submitter")
run_ssgsea(cpm, "log2cpm")

# 导出 NNLS B 代理子集（marker 基因 + PDK 基因，双口径）
# 注：step1b_immune_nnls_gate.py 在内存中完成等价计算（422-marker 参照 evidence_gain/2_cross_deconv/
#     ref_signature_marker_linear.csv）；本段导出的合并 marker 清单位于 2x_myeloid_ref。
markers <- read.csv("evidence_gain/2x_myeloid_ref/ref_marker_list_combined.csv", stringsAsFactors = FALSE)
mk <- unique(markers$gene)
for (mobj in list(list(nm = "submitter", m = norm), list(nm = "log2cpm", m = cpm))) {
  hit <- intersect(union(mk, c("PDK2", "PDK4")), rownames(mobj$m))
  df <- as.data.frame(mobj$m[hit, ]); df$gene <- rownames(df)
  fwrite(df, file.path(OUT, paste0("gse273763_bulk_subset_", mobj$nm, ".csv")))
  cat(sprintf("[%s NNLS subset] %d genes x %d samples\n", mobj$nm, length(hit), ncol(mobj$m)))
}
cat("\nSTEP1_IMMUNE_DONE\n")

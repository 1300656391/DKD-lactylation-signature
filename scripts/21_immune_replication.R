# 21_immune_replication.R ----
# R8 M-A1/M-B2：免疫分量补强
#   A. GSE142025（RNA-seq外部验证队列）免疫浸润复现：18细胞ssGSEA（同训练管线Bindea签名）DN vs 对照
#      + 候选基因/髓系评分Spearman相关复现
#   B. snRNA（GSE131882注释对象）免疫群比例 DN vs 对照
suppressWarnings(suppressMessages({
  library(data.table); library(GSVA); library(Seurat)
}))
set.seed(42)

# ---------- A. GSE142025 免疫复现 ----------
# 签名：复用训练管线的Bindea签名（suffix "_Bindea_et_al"），18免疫=剔除Normal_mucosa等非免疫
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
cat("18-cell set loaded:", length(sig18), "\n")
stopifnot(length(sig18) == 18)

e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])))
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
rownames(e14) <- toupper(rownames(e14))
keep <- p14$group %in% c("Advanced_DN", "Early_DN", "Control")
e14k <- as.matrix(e14[, keep])
grp <- p14$group[keep]
is_dn <- grp != "Control"
stage <- ifelse(grp == "Early_DN", "Early", ifelse(grp == "Advanced_DN", "Advanced", "Control"))

# ssGSEA（GSVA新API）
par <- ssgseaParam(exprData = e14k, geneSets = sig18, normalize = TRUE)
ss <- gsva(par, verbose = FALSE)
ss <- as.data.frame(t(ss)); ss$group <- grp; ss$stage <- stage

# DN vs Control Wilcoxon + BH18
res <- do.call(rbind, lapply(colnames(ss)[1:18], function(cc) {
  x <- ss[[cc]]; g <- is_dn
  data.frame(cell = cc, median_DN = median(x[g]), median_Ctrl = median(x[!g]),
             p = wilcox.test(x[g], x[!g])$p.value)
}))
res$p.adj_BH18 <- p.adjust(res$p, "BH")
res$direction <- ifelse(res$median_DN > res$median_Ctrl, "up", "down")
setorder(as.data.table(res), p.adj_BH18)
fwrite(res, "results/revision/M7_gse142025_immune.csv")
cat("\n== GSE142025 DN vs Control (18 cells) ==\n"); print(res)

# PDK2/PDK4 × 髓系/T细胞评分 Spearman（复现训练队列相关结构）
cand <- c("PDK2", "PDK4")
mye <- c("Macrophages", "Neutrophils", "DC", "iDC", "Eosinophils")
cor_res <- do.call(rbind, lapply(cand, function(g) {
  ge <- as.numeric(e14k[g, ])
  do.call(rbind, lapply(mye, function(cc) {
    s <- spearman_corr <- suppressWarnings(cor(ge, ss[[cc]], method = "spearman"))
    data.frame(gene = g, cell = cc, rho = round(s, 3))
  }))
}))
cat("\n== GSE142025 PDK x myeloid (Spearman) ==\n"); print(cor_res)
fwrite(cor_res, "results/revision/M7b_gse142025_pdk_myeloid.csv")

# 早期 vs 晚期亚组（巨噬细胞）
if (sum(stage != "Control") > 0) {
  ea <- ss[stage != "Control", ]
  w <- wilcox.test(ea$Macrophages[ea$stage == "Early"], ea$Macrophages[ea$stage == "Advanced"])
  cat("\nMacrophages Early vs Advanced: p =", signif(w$p.value, 3), "\n")
}

# ---------- B. snRNA免疫群比例 ----------
merged <- readRDS("results/revision/GSE131882_seurat_annotated.rds")
ct <- merged$celltype; gp <- merged$group
cat("\n== snRNA cell type x group (immune focus) ==\n")
tb <- table(ct, gp)
print(tb)
imm_types <- c("Macro", "T/NK", "B")
props <- do.call(rbind, lapply(imm_types, function(t) {
  dn <- sum(gp == "DN"); ctl <- sum(gp == "Control")
  data.frame(type = t, pct_DN = round(100 * sum(ct == t & gp == "DN") / dn, 2),
             pct_Ctrl = round(100 * sum(ct == t & gp == "Control") / ctl, 2))
}))
print(props)
fwrite(props, "results/revision/M7c_scrna_immune_proportion.csv")

cat("\nIMMUNE_REPLICATION_DONE\n")

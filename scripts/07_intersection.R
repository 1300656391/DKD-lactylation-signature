# 07_intersection.R ----
# S7: 候选基因确定 = LRG ∩ (DEGs ∪ WGCNA疾病模块基因)
# 输出: results/A_candidates.csv + 控制台Venn计数
# 运行: Rscript 07_intersection.R

suppressWarnings(suppressMessages(library(data.table)))

lrg  <- fread("data/LRG_geneset.csv")
degs <- fread("results/A_DEGs.csv")
wgc  <- fread("results/wgcna_sigmodule_genes.csv")

# 别名展开（新symbol→旧symbol，兼容芯片注释）
# R11-S1 fix: fread reads empty alias cells as "", not NA; na.omit kept "" as a phantom key (178 vs true 177)
lrg_alias <- lrg$alias[!is.na(lrg$alias) & nzchar(trimws(lrg$alias))]
lrg_keys <- unique(c(lrg$symbol, lrg_alias))
lrg_keys <- toupper(lrg_keys)
cat("LRG键数(含别名):", length(lrg_keys), "\n")

# LRG命中DEGs
deg_genes <- toupper(degs$gene)
is_lrg_deg <- deg_genes %in% lrg_keys
cat("DEGs总数:", length(deg_genes), "\n")
cat("LRG ∩ DEGs:", sum(is_lrg_deg), "\n")

# LRG命中WGCNA疾病模块
wgc_genes <- toupper(wgc$gene)
lrg_wgc <- intersect(lrg_keys, wgc_genes)
cat("LRG ∩ WGCNA疾病模块:", length(lrg_wgc), "\n")

# 候选 = LRG ∩ (DEGs ∪ 模块基因)
cand_keys <- union(deg_genes[is_lrg_deg], lrg_wgc)
cand_keys <- intersect(lrg_keys, cand_keys)
cat("候选基因数(LRG ∩ (DEGs ∪ 模块)):", length(cand_keys), "\n")

# 构建候选表：优先用DEG统计信息
deg_idx <- match(cand_keys, deg_genes)
wgc_flag <- cand_keys %in% wgc_genes
cand <- data.frame(
  gene = cand_keys,
  in_DEG = !is.na(deg_idx),
  in_WGCNA_module = wgc_flag,
  stringsAsFactors = FALSE
)
di <- deg_idx[!is.na(deg_idx)]
cand$logFC[!is.na(deg_idx)] <- degs$logFC[di]
cand$adj.P.Val[!is.na(deg_idx)] <- degs$adj.P.Val[di]
cand$aveExpr[!is.na(deg_idx)] <- degs$AveExpr[di]
# 标注LRG来源
src_map <- lrg$sources[match(cand_keys, toupper(lrg$symbol))]
cand$LRG_sources <- ifelse(is.na(src_map), "alias-hit", src_map)

cand <- cand[order(cand$in_DEG, cand$in_WGCNA_module, -abs(cand$logFC)), ]
fwrite(cand, "results/A_candidates.csv")
cat("\n候选基因列表:\n")
print(cand[, c("gene", "in_DEG", "in_WGCNA_module", "logFC", "adj.P.Val", "LRG_sources")])

# QC: 交集>=15
if (sum(is_lrg_deg) < 15) {
  cat("\n警告: LRG∩DEG < 15，按计划应扩充为'乳酸代谢+乳酸化'合并集（已含，无需动作）或评估备胎\n")
}
cat("\nINTERSECTION_DONE\n")

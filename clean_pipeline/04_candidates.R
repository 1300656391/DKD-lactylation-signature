# 04_candidates.R ----
# Candidate pool = LRG intersect (DEGs union WGCNA disease-module genes).
# Output: results/A_candidates.csv

suppressWarnings(suppressMessages(library(data.table)))
set.seed(42)

lrg  <- fread("data/LRG_geneset.csv")
degs <- fread("results/A_DEGs.csv")
wgc  <- fread("results/wgcna_sigmodule_genes.csv")

# Alias expansion (new symbol -> array-annotation legacy symbol).
# Empty alias cells read as "" (not NA) and must be dropped, otherwise ""
# survives as a phantom key.
lrg_alias <- lrg$alias[!is.na(lrg$alias) & nzchar(trimws(lrg$alias))]
lrg_keys <- toupper(unique(c(lrg$symbol, lrg_alias)))
cat("LRG keys (with aliases):", length(lrg_keys), "\n")

deg_genes <- toupper(degs$gene)
is_lrg_deg <- deg_genes %in% lrg_keys
cat("DEGs:", length(deg_genes), "| LRG hit DEGs:", sum(is_lrg_deg), "\n")

wgc_genes <- toupper(wgc$gene)
lrg_wgc <- intersect(lrg_keys, wgc_genes)
cat("LRG in WGCNA disease modules:", length(lrg_wgc), "\n")

cand_keys <- union(deg_genes[is_lrg_deg], lrg_wgc)
cand_keys <- intersect(lrg_keys, cand_keys)
cat("candidates (LRG intersect (DEG union module)):", length(cand_keys), "\n")

deg_idx <- match(cand_keys, deg_genes)
cand <- data.frame(
  gene = cand_keys,
  in_DEG = !is.na(deg_idx),
  in_WGCNA_module = cand_keys %in% wgc_genes,
  stringsAsFactors = FALSE
)
di <- deg_idx[!is.na(deg_idx)]
cand$logFC[!is.na(deg_idx)] <- degs$logFC[di]
cand$adj.P.Val[!is.na(deg_idx)] <- degs$adj.P.Val[di]
cand$aveExpr[!is.na(deg_idx)] <- degs$AveExpr[di]
src_map <- lrg$sources[match(cand_keys, toupper(lrg$symbol))]
cand$LRG_sources <- ifelse(is.na(src_map), "alias-hit", src_map)

cand <- cand[order(cand$in_DEG, cand$in_WGCNA_module, -abs(cand$logFC)), ]
fwrite(cand, "results/A_candidates.csv")
print(cand[, c("gene", "in_DEG", "in_WGCNA_module", "logFC", "adj.P.Val", "LRG_sources")])

if (sum(is_lrg_deg) < 15)
  cat("WARNING: LRG-hit DEGs < 15; candidate pool relies on module genes\n")

cat("\nCANDIDATES_DONE\n")

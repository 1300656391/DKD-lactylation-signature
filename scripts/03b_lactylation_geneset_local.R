# 03b_lactylation_geneset_local.R ----
# S6(本地版): 乳酸化相关基因集(LRG)构建 —— 三正交来源
#   源1: GO条目 (GO.db + org.Hs.eg.db, 关键词 lactate/lactyl/lactic) —— 离线
#   源2: Reactome通路 (reactome.db, 关键词 LACTATE/PYRUVATE/GLYCOLYSIS) —— 离线
#   源3: 文献策展 (修饰酶/代谢转运/已验证乳酸化底物, 逐基因出处)
# 输出: data/LRG_geneset.csv + data/LRG_geneset_detail.csv
# 运行: Rscript 03b_lactylation_geneset_local.R

suppressWarnings(suppressMessages({
  library(GO.db)
  library(org.Hs.eg.db)
  library(reactome.db)
  library(data.table)
}))

# ---------- 源1: GO ----------
go_terms_all <- as.list(GOTERM)
term_names <- vapply(go_terms_all, function(x) Term(x), character(1))
hit_go_ids <- names(term_names)[grepl("lactate|lactyl|lactic", term_names, ignore.case = TRUE)]
cat("GO命中条目数:", length(hit_go_ids), "\n")
cat(paste(sprintf("%s: %s", hit_go_ids, term_names[hit_go_ids]), collapse = "\n"), "\n\n")

# GO -> ENTREZ -> SYMBOL
go2eg <- as.list(org.Hs.egGO2ALLEGS)
egs <- unique(unlist(go2eg[intersect(hit_go_ids, names(go2eg))]))
eg2sym <- mapIds(org.Hs.eg.db, keys = egs, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
sym_go <- unique(na.omit(unname(eg2sym)))
src_go <- data.frame(symbol = sym_go, source = "GO",
                     detail = "GO terms: lactate/lactyl/lactic (GO.db)", stringsAsFactors = FALSE)
cat("GO来源基因数:", length(sym_go), "\n")

# ---------- 源2: Reactome ----------
p2name <- as.list(reactome.db::reactomePATHID2NAME)
hit_pw <- names(p2name)[grepl("lactate|pyruvate|glycolysis", unlist(p2name), ignore.case = TRUE)]
hit_pw <- hit_pw[grepl("^R-HSA-", hit_pw)]  # 仅人类通路
cat("\nReactome命中通路数:", length(hit_pw), "\n")
cat(paste(sprintf("%s: %s", hit_pw, unlist(p2name[hit_pw])), collapse = "\n"), "\n\n")

p2ext <- as.list(reactome.db::reactomePATHID2EXTID)
extids <- unique(unlist(p2ext[hit_pw]))
# EXTID识别：数字→Entrez；symbol样式→直接用；其余→尝试UniProt
num   <- grepl("^[0-9]+$", extids)
sym_re <- extids[!num & grepl("^[A-Z][A-Z0-9-]*$", extids)]
if (any(num)) {
  e2s <- mapIds(org.Hs.eg.db, keys = extids[num], column = "SYMBOL",
                keytype = "ENTREZID", multiVals = "first")
  sym_re <- c(sym_re, na.omit(unname(e2s)))
}
others <- setdiff(extids, c(sym_re, extids[num]))
if (length(others) > 0) {
  u2e <- tryCatch(mapIds(org.Hs.eg.db, keys = others, column = "ENTREZID",
                         keytype = "UNIPROT", multiVals = "first"),
                  error = function(e) NULL)
  if (!is.null(u2e)) {
    u2s <- mapIds(org.Hs.eg.db, keys = na.omit(unname(u2e)), column = "SYMBOL",
                  keytype = "ENTREZID", multiVals = "first")
    sym_re <- c(sym_re, na.omit(unname(u2s)))
  }
}
sym_re <- unique(sym_re)
src_re <- data.frame(symbol = sym_re, source = "Reactome",
                     detail = "Reactome: LACTATE/PYRUVATE/GLYCOLYSIS (reactome.db)", stringsAsFactors = FALSE)
cat("Reactome来源基因数:", length(sym_re), "\n")

# ---------- 源3: 文献策展 ----------
curated <- data.frame(
  symbol = c(
    "EP300", "AARS1", "LARS1", "KRT18",
    "HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "BRD4",
    "LDHA", "LDHB", "LDHC",
    "SLC16A1", "SLC16A3", "SLC16A7", "SLC5A8", "SLC5A12",
    "PDK1", "PDK2", "PDK3", "PDK4", "PDHA1", "PDHB", "PDP1", "PDP2",
    "MPC1", "MPC2", "ACSS1", "ACSS2",
    "KLF5", "STAT1", "RUNX1", "TFRC", "IGFBP5", "HIPK2", "FAS",
    "ACSL4", "ACSF2", "NAT10", "TRIM65", "GLIS1", "PFKFB3",
    "H3-3A", "H3C1", "H3C2", "H3C3", "H4C1", "H4C2"
  ),
  source = "Curated",
  detail = "文献策展：乳酸化修饰酶/乳酸代谢转运/DKD已验证乳酸化底物（Nature 2019; Redox Biol 2024; Cell Death Differ 2026; JBC 2026; Diabetologia 2024 等）",
  stringsAsFactors = FALSE
)
# 历史别名兼容（数组平台的旧symbol）
alias_map <- c("H3-3A" = "H3F3A", "H3C1" = "HIST1H3A", "H3C2" = "HIST1H3B", "H3C3" = "HIST1H3C",
               "H4C1" = "HIST1H4A", "H4C2" = "HIST1H4B")
curated$alias <- unname(alias_map[curated$symbol])

# ---------- 合并 ----------
all_src <- rbind(src_go, src_re, curated[, c("symbol", "source", "detail")])
all_src$symbol <- toupper(trimws(all_src$symbol))
all_src <- all_src[all_src$symbol != "" & !is.na(all_src$symbol), ]
all_src <- all_src[!duplicated(all_src[, c("symbol", "source")]), ]

lrg <- data.frame(symbol = sort(unique(all_src$symbol)), stringsAsFactors = FALSE)
lrg$n_sources <- vapply(lrg$symbol, function(s) length(unique(all_src$source[all_src$symbol == s])), integer(1))
lrg$sources <- vapply(lrg$symbol, function(s) paste(sort(unique(all_src$source[all_src$symbol == s])), collapse = "+"), character(1))

# 附上别名（用于旧注释平台映射）
lrg <- merge(lrg, curated[, c("symbol", "alias")], by = "symbol", all.x = TRUE)

fwrite(all_src, "data/LRG_geneset_detail.csv")
fwrite(lrg, "data/LRG_geneset.csv")

cat("\n==== LRG基因集统计 ====\n")
cat("GO来源基因:", length(sym_go), "\n")
cat("Reactome来源基因:", length(sym_re), "\n")
cat("策展来源基因:", nrow(curated), "\n")
cat("合并去重后总基因数:", nrow(lrg), "\n")
cat("多来源基因数:", sum(lrg$n_sources >= 2), "\n")
cat("\nLRG_GENASET_DONE\n")

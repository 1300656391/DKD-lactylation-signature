# 03_lactylation_geneset.R ----
# S6: 乳酸化相关基因集(LRG)构建 —— 三正交来源
#   源1: GO条目 (msigdbr, 关键词 lactate/lactyl/lactic)
#   源2: Reactome通路 (msigdbr CP集, 关键词 LACTATE/PYRUVATE/GLYCOLYSIS)
#   源3: 文献策展 (修饰酶/代谢转运 + 已验证乳酸化底物, 逐基因注明出处)
# 输出: data/LRG_geneset.csv (symbol, source, detail) + 控制台统计
# 运行: Rscript 03_lactylation_geneset.R

suppressWarnings(suppressMessages({
  library(msigdbr)
  library(data.table)
}))

get_msig <- function(...) {
  # 兼容msigdbr新旧API
  out <- tryCatch(msigdbr(species = "Homo sapiens", ...),
                  error = function(e) NULL)
  if (is.null(out)) out <- tryCatch(
    msigdbr(species = "Homo sapiens", category = "C5"),
    error = function(e) NULL)
  out
}

# ---------- 源1: GO ----------
go_db <- tryCatch(msigdbr(species = "Homo sapiens", db = "GO"),
                  error = function(e) msigdbr(species = "Homo sapiens", category = "C5"))
go_hits <- go_db[grepl("LACTATE|LACTYL|LACTIC", go_db$gs_name, ignore.case = TRUE), ]
go_terms <- unique(go_hits$gs_name)
cat("GO命中条目数:", length(go_terms), "\n")
print(go_terms)
src_go <- unique(data.frame(
  symbol = go_hits$gene_symbol,
  source = "GO",
  detail = paste0("GO: ", go_hits$gs_name),
  stringsAsFactors = FALSE
))

# ---------- 源2: Reactome ----------
cp_db <- tryCatch(msigdbr(species = "Homo sapiens", db = "CP"),
                  error = function(e) msigdbr(species = "Homo sapiens", category = "C2:CP:REACTOME"))
re_hits <- cp_db[grepl("REACTOME_", cp_db$gs_name) &
                 grepl("LACTATE|PYRUVATE|GLYCOLYSIS", cp_db$gs_name), ]
re_terms <- unique(re_hits$gs_name)
cat("\nReactome命中通路数:", length(re_terms), "\n")
print(re_terms)
src_re <- unique(data.frame(
  symbol = re_hits$gene_symbol,
  source = "Reactome",
  detail = paste0("Pathway: ", re_hits$gs_name),
  stringsAsFactors = FALSE
))

# ---------- 源3: 文献策展 ----------
curated <- data.frame(
  symbol = c(
    # 修饰写入酶 (writers)
    "EP300",                       # p300: 首个组蛋白乳酸化写入酶 (Nature 2019)
    "AARS1", "LARS1", "KRT18",     # 乳酸转移酶: AARS1(CDD 2026 DKD), LARS1(Cell Signal 2025 DKD), KRT18(Exp Mol Med 2026 DKD)
    # 去乳酸化酶 (erasers)
    "HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3",
    # 报道的读取酶 (readers)
    "BRD4",
    # 乳酸代谢与转运
    "LDHA", "LDHB", "LDHC",
    "SLC16A1", "SLC16A3", "SLC16A7",      # MCT1/4/2
    "SLC5A8", "SLC5A12",                  # 钠耦联 monocarboxylate 转运体
    "PDK1", "PDK2", "PDK3", "PDK4",       # 丙酮酸脱氢酶激酶
    "PDHA1", "PDHB", "PDP1", "PDP2",
    "MPC1", "MPC2",                       # 线粒体丙酮酸载体
    "ACSS1", "ACSS2",                     # 酰基CoA合成酶(短链酸代谢)
    # DKD/肾脏中已验证的乳酸化底物
    "KLF5",     # H3K14la/KLF5 轴 EMT (Redox Biol 2024)
    "STAT1",    # AARS1乳酸化STAT1促铁死亡 (Cell Death Differ 2026)
    "RUNX1",    # SIRT3-RUNX1反馈环 (Int J Biol Sci 2026)
    "TFRC",     # TFRC K374la 介导铁死亡 (JBC 2026)
    "IGFBP5",   # IGFBP5-H3K18la-EndoMT (Cell Mol Life Sci 2025)
    "HIPK2",    # H4K12la-HIPK2 (Int J Mol Med 2026)
    "FAS",      # K18乳酸转移酶调控Fas转录 (Exp Mol Med 2026)
    "ACSL4",    # H3K18la/LDHA-ACSL4铁死亡 (ABBS 2025)
    "ACSF2",    # ACSF2溶酶乳酸化肾小管损伤 (Diabetologia 2024)
    "NAT10",    # NAT10-TGFb/lactate (Int J Biol Sci 2026)
    "TRIM65",   # TRIM65铁死亡/糖酵解 (Cell Rep 2025)
    "GLIS1",    # Glis1抑制组蛋白乳酸化 (Life Sci 2025)
    "PFKFB3",   # PFKFB3-H4K12la-HIPK2 (Int J Mol Med 2026)
    # 代表性组蛋白基因(乳酸化主要发生载体)
    "H3-3A", "H3C1", "H3C2", "H3C3", "H4C1", "H4C2"
  ),
  stringsAsFactors = FALSE
)
curated$source <- "Curated"
curated$detail <- "文献策展：修饰酶/代谢转运/已验证乳酸化底物（详见论文补充材料引用）"

# ---------- 合并 ----------
all_src <- rbind(src_go, src_re, curated)
all_src$symbol <- toupper(trimws(all_src$symbol))
all_src <- all_src[all_src$symbol != "" & !is.na(all_src$symbol), ]
all_src <- all_src[!duplicated(all_src[, c("symbol", "source")]), ]

lrg <- data.frame(
  symbol = sort(unique(all_src$symbol)),
  stringsAsFactors = FALSE
)
lrg$n_sources <- vapply(lrg$symbol, function(s) length(unique(all_src$source[all_src$symbol == s])), integer(1))
lrg$sources <- vapply(lrg$symbol, function(s) paste(unique(all_src$source[all_src$symbol == s]), collapse = "+"), character(1))

fwrite(all_src, "data/LRG_geneset_detail.csv")
fwrite(lrg, "data/LRG_geneset.csv")

cat("\n==== LRG基因集统计 ====\n")
cat("GO来源基因:", length(unique(src_go$symbol)), "\n")
cat("Reactome来源基因:", length(unique(src_re$symbol)), "\n")
cat("策展来源基因:", length(unique(curated$symbol)), "\n")
cat("合并去重后总基因数:", nrow(lrg), "\n")
cat("多来源基因数:", sum(lrg$n_sources >= 2), "\n")
cat("\nLRG_GENASET_DONE\n")

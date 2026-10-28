# 02_lactylation_geneset.R ----
# Validate the provided lactylation-related gene (LRG) reference list and
# summarize its composition. The list itself (data/LRG_geneset.csv, 171
# unique genes; 177 symbol keys including aliases) was assembled from three
# sources - GO terms, Reactome pathways, and literature-curated lactylation
# machinery/substrates - and ships with the repository.
# Output: results/lrg_summary.txt

suppressWarnings(suppressMessages(library(data.table)))
set.seed(42)
dir.create("results", showWarnings = FALSE)

stopifnot(file.exists("data/LRG_geneset.csv"))
lrg <- fread("data/LRG_geneset.csv")
stopifnot(all(c("symbol", "alias", "sources") %in% colnames(lrg)))

lrg_alias <- lrg$alias[!is.na(lrg$alias) & nzchar(trimws(lrg$alias))]
keys <- toupper(unique(c(lrg$symbol, lrg_alias)))

src_tab <- table(lrg$sources)
lines <- c(
  "LRG reference list validation",
  sprintf("genes: %d unique symbols", length(unique(lrg$symbol))),
  sprintf("alias entries: %d", length(lrg_alias)),
  sprintf("lookup keys (symbols + aliases): %d", length(keys)),
  sprintf("multi-source genes: %d", sum(grepl("\\+|;|,", lrg$sources))),
  "",
  "composition by source:",
  capture.output(print(src_tab)),
  "",
  "expected values: 171 genes / 177 keys (used by steps 04 and 13)")
stopifnot(length(unique(lrg$symbol)) == 171, length(keys) == 177)
writeLines(lines, "results/lrg_summary.txt")
cat(paste(lines, collapse = "\n"), "\n")

cat("\nLRG_VALIDATE_DONE\n")

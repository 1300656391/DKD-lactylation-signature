# 01c_gse104948_diagnosis.R ----
# GSE104948 双平台diagnosis计数 + Bioc核心包加载自检

for (acc in c("GSE104948-GPL22945", "GSE104948-GPL24120")) {
  ml <- readLines(paste0("data/", acc, "_meta.txt"), warn = FALSE)
  ch <- ml[grepl("^!Sample_characteristics_ch1", ml)]
  vals <- lapply(ch, function(l) {
    s <- sub("^!Sample_characteristics_ch1\\d*\\s*=\\s*", "", l)
    trimws(gsub('^"|"$', "", strsplit(s, "\t", fixed = TRUE)[[1]]))
  })
  cat("\n=====", acc, "=====\n")
  for (i in seq_along(vals)) {
    v <- vals[[i]]
    if (any(grepl("diagnosis:", v))) {
      d <- grep("diagnosis:", v, value = TRUE)
      d <- sub("^diagnosis:\\s*", "", d)
      cat("diagnosis计数:\n"); print(table(d))
    }
    if (any(grepl("tissue:", v))) {
      t <- grep("tissue:", v, value = TRUE)
      cat("tissue计数:\n"); print(table(sub("^tissue:\\s*", "", t)))
    }
  }
}

cat("\n===== Bioc核心包加载自检 =====\n")
pkgs <- c("limma", "sva", "WGCNA", "clusterProfiler", "GSVA", "ComplexHeatmap",
          "GEOquery", "ReactomePA", "org.Hs.eg.db", "dcurves", "ggplot2")
for (p in pkgs) {
  ok <- suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
  cat(sprintf("%-16s %s\n", p, ifelse(ok, "OK", "MISSING")))
}
cat("\nDONE_01C\n")

# 02_probe_inspect.R ----
# S3前置：各数据集平台与探针格式探查

sets <- c("GSE96804", "GSE30528", "GSE30529",
          "GSE104948-GPL22945", "GSE104948-GPL24120")

for (acc in sets) {
  f <- paste0("data/", acc, "_expr.rds")
  if (!file.exists(f)) next
  e <- readRDS(f)
  ml <- readLines(paste0("data/", sub("-GPL.*$", "", acc), "_meta.txt"), warn = FALSE)
  plat <- ml[grepl("^!Sample_platform_id|^!Series_platform_id", ml)]
  cat("=====", acc, "=====\n")
  cat("维度:", nrow(e), "x", ncol(e), "\n")
  cat("平台行: ", paste(head(plat, 2), collapse = " | "), "\n")
  cat("行名样例: ", paste(head(rownames(e), 8), collapse = ", "), "\n")
  # 判断是否已经是基因symbol（启发式：大写字母开头、含字母且无下划线/点后缀数字）
  rn <- head(rownames(e), 500)
  sym_like <- mean(grepl("^[A-Z][A-Z0-9-]*$", rn))
  cat("疑似symbol比例(前500): ", round(sym_like, 2), "\n\n")
}
cat("INSPECT_DONE\n")

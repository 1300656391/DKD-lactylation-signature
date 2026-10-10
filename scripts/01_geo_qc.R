# 01_geo_qc.R ----
# S2: GEO series matrix 解析与QC
# 输入: data/GEO/GSE*_series_matrix.txt.gz
# 输出: data/GSE<acc>_expr.rds + data/GSE<acc>_meta.txt + results/qc_report.txt
# 运行: Rscript 01_geo_qc.R

suppressWarnings(suppressMessages({
  library(data.table)
}))

geo_dir  <- "data/GEO"
out_dir  <- "data"
res_dir  <- "results"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# 自动发现所有 series matrix 文件（GSE104948 为双平台，文件名带平台后缀）
sm_files <- list.files(geo_dir, pattern = "_series_matrix\\.txt\\.gz$", full.names = TRUE)
sm_files <- sm_files[file.size(sm_files) > 100000] # 剔除仅元数据的空壳文件
accs <- sub("_series_matrix\\.txt\\.gz$", "", basename(sm_files))

qc_lines <- character(0)
qadd <- function(...) qc_lines <<- c(qc_lines, paste0(...))

parse_series_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)

  keep_meta <- grepl(
    "^!Series_title|^!Series_summary|^!Series_overall_design|^!Sample_title|^!Sample_geo_accession|^!Sample_source_name|^!Sample_characteristics_ch1|^!Sample_platform_id|^!Platform_title|^!Sample_supplementary",
    lines
  )
  meta <- lines[keep_meta]

  beg <- grep("!series_matrix_table_begin", lines, fixed = TRUE)
  end <- grep("!series_matrix_table_end",   lines, fixed = TRUE)
  expr <- NULL
  if (length(beg) == 1 && length(end) >= 1) {
    expr <- read.delim(
      gzfile(path), skip = beg, nrows = end - beg - 1,
      row.names = 1, check.names = FALSE, na.strings = c("NA", "null", "")
    )
  }
  list(meta = meta, expr = expr)
}

for (acc in accs) {
  f <- file.path(geo_dir, paste0(acc, "_series_matrix.txt.gz"))
  qadd("================================================")
  qadd("数据集: ", acc)
  if (!file.exists(f)) { qadd("  !! 文件不存在，跳过"); next }

  parsed <- parse_series_matrix(f)
  expr <- parsed$expr
  meta <- parsed$meta

  writeLines(meta, file.path(out_dir, paste0(acc, "_meta.txt")))

  if (is.null(expr)) {
    qadd("  !! 无表达矩阵区块（RNA-seq常见），需改用原始文件/其它来源")
    next
  }

  qadd("  表达矩阵维度: ", nrow(expr), " 探针/特征 x ", ncol(expr), " 样本")
  qadd("  值域: [", round(min(expr, na.rm = TRUE), 2), ", ", round(max(expr, na.rm = TRUE), 2), "]")
  qadd("  中位数: ", round(median(as.matrix(expr), na.rm = TRUE), 2),
       " （若 >20 提示未log2，需转换）")
  na_frac <- mean(is.na(as.matrix(expr)))
  qadd("  NA比例: ", round(na_frac * 100, 2), "%")

  # 特征关键词汇总（每行characteristics去重）
  ch <- meta[grepl("^!Sample_characteristics_ch1", meta)]
  uniq_ch <- unique(unlist(strsplit(sub("^!Sample_characteristics_ch1\\d*\\s*=\\s*", "", ch), "\t")))
  qadd("  样本特征取值（前40个去重）:")
  for (v in utils::head(uniq_ch, 40)) qadd("    - ", v)

  saveRDS(expr, file.path(out_dir, paste0(acc, "_expr.rds")))
  qadd("  已保存: ", file.path(out_dir, paste0(acc, "_expr.rds")))
}

writeLines(qc_lines, file.path(res_dir, "qc_report.txt"))
cat(qc_lines, sep = "\n")
cat("\nQC_REPORT_DONE\n")

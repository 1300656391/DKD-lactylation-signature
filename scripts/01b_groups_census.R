# 01b_groups_census.R ----
# S2补充：各数据集分组标签普查 + GSE142025(RNA-seq)表达矩阵合并
# 输出: 控制台分组统计 + data/GSE142025_expr.rds + data/GSE142025_pheno.csv
# 运行: Rscript 01b_groups_census.R

suppressWarnings(suppressMessages(library(data.table)))

out_lines <- character(0)
qadd <- function(...) out_lines <<- c(out_lines, paste0(...))

# ---------- 1. 芯片数据集分组普查 ----------
meta_files <- list.files("data", pattern = "_meta\\.txt$", full.names = TRUE)

clean_v <- function(x) trimws(gsub('^"|"$', "", x))

for (mf in meta_files) {
  acc <- sub("_meta\\.txt$", "", basename(mf))
  if (acc == "GSE142025") next
  ml <- readLines(mf, warn = FALSE)
  qadd("================================================")
  qadd("数据集: ", acc)

  # characteristics 矩阵（每行一个特征通道）
  ch <- ml[grepl("^!Sample_characteristics_ch1", ml)]
  if (length(ch) > 0) {
    vals <- lapply(ch, function(l) {
      s <- sub("^!Sample_characteristics_ch1\\d*\\s*=\\s*", "", l)
      clean_v(strsplit(s, "\t", fixed = TRUE)[[1]])
    })
    cm <- do.call(rbind, vals)
    for (i in seq_len(nrow(cm))) {
      v <- cm[i, ]
      if (all(grepl(":", v)) && length(unique(v)) <= 30) {
        tab <- table(sub("^.*?:\\s*", "", v))
        qadd("  特征行", i, "(", unique(sub("\\s*:.*$", "", v))[1], "):")
        for (k in seq_along(tab)) qadd("    ", names(tab)[k], " = ", as.integer(tab[k]))
      }
    }
  }
  # 标题/来源前缀普查（用于GSE96804类标签藏在title的情况）
  for (tag in c("!Sample_title", "!Sample_source_name_ch1")) {
    l <- ml[grepl(paste0("^", tag), ml)]
    if (length(l) > 0) {
      s <- sub(paste0("^", tag, ".*?=\\s*"), "", l[1])
      v <- clean_v(strsplit(s, "\t", fixed = TRUE)[[1]])
      pref <- substr(v, 1, 4)
      qadd("  ", tag, " 前4字符频数: ",
           paste(names(table(pref)), table(pref), sep = "=", collapse = ", "))
    }
  }
}

# ---------- 2. GSE142025 RNA-seq 合并 ----------
qadd("================================================")
qadd("GSE142025 RNA-seq 合并中...")
raw_dir <- "data/GEO/GSE142025_raw"
sfiles <- list.files(raw_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
qadd("  样本文件数: ", length(sfiles))

mat_list <- lapply(sfiles, function(f) {
  d <- fread(f, header = TRUE)
  setnames(d, c("Symbol", "expr"))
  d
})
symbols <- mat_list[[1]]$Symbol
ok <- vapply(mat_list, function(d) identical(d$Symbol, symbols), logical(1))
qadd("  Symbol列完全一致: ", all(ok))
expr <- as.data.frame(do.call(cbind, lapply(mat_list, function(d) d$expr)))
rownames(expr) <- symbols
colnames(expr) <- tools::file_path_sans_ext(basename(sfiles))
qadd("  合并矩阵: ", nrow(expr), "基因 x ", ncol(expr), "样本")
qadd("  值域: [", round(min(expr), 2), ", ", round(max(expr), 2), "]")

# 分组：从series matrix元数据取
if (!file.exists("data/GSE142025_meta.txt")) {
  ml_all <- readLines(gzfile("data/GEO/GSE142025_series_matrix.txt.gz"), warn = FALSE)
  writeLines(ml_all[grepl("^!", ml_all)], "data/GSE142025_meta.txt")
}
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
title_l <- ml[grepl("^!Sample_title", ml)][1]
group_l <- ml[grepl("^!Sample_characteristics_ch1", ml)]
grp <- clean_v(strsplit(sub("^.*?=\\s*", "", group_l[1]), "\t", fixed = TRUE)[[1]])
ttl <- clean_v(strsplit(sub("^.*?=\\s*", "", title_l), "\t", fixed = TRUE)[[1]])
pheno <- data.frame(sample_id = ttl, group = grp, stringsAsFactors = FALSE)
qadd("  分组: ", paste(names(table(pheno$group)), table(pheno$group), sep = "=", collapse = ", "))

saveRDS(as.matrix(expr), "data/GSE142025_expr.rds")
fwrite(pheno, "data/GSE142025_pheno.csv")
qadd("  已保存: data/GSE142025_expr.rds + data/GSE142025_pheno.csv")

writeLines(out_lines, "results/groups_census.txt")
cat(out_lines, sep = "\n")
cat("\nCENSUS_DONE\n")

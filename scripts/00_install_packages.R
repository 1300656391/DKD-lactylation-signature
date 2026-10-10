# 00_install_packages.R ----
# SCI论文生产计划 · 环境搭建 S1
# 用途：安装全部分析所需R包（CRAN + Bioconductor）
# 运行：Rscript 00_install_packages.R
# 注：2026-10-08 首次执行；R 4.6.1

options(repos = c(CRAN = "https://cloud.r-project.org"))
N_CORES <- max(1, parallel::detectCores() - 1)

cran_pkgs <- c(
  # 数据处理
  "tidyverse", "data.table", "janitor",
  # 基因集与富集辅助
  "msigdbr", "VennDiagram", "ggpubr", "patchwork", "scales",
  # 机器学习与模型
  "glmnet", "caret", "e1071", "randomForest", "xgboost",
  "pROC", "rms", "ggDCA", "shapviz",
  # 网络与可视化
  "igraph", "ggraph",
  # 免疫浸润辅助
  "reshape2", "corrplot"
)

bioc_pkgs <- c(
  # 数据获取与预处理
  "GEOquery", "limma", "sva",
  # 共表达与富集
  "WGCNA", "clusterProfiler", "org.Hs.eg.db", "ReactomePA",
  # 免疫浸润/ssGSEA
  "GSVA", "immunedeconv",
  # 热图
  "ComplexHeatmap"
)

message("=== [1/3] 安装CRAN包（", length(cran_pkgs), "个） ===")
to_install <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) {
  install.packages(to_install, Ncpus = N_CORES)
}

message("=== [2/3] 安装Bioconductor包（", length(bioc_pkgs), "个） ===")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) {
  BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE, Ncpus = N_CORES)
}

message("=== [3/3] 加载自检 ===")
all_pkgs <- c(cran_pkgs, bioc_pkgs)
failed <- character(0)
for (p in all_pkgs) {
  ok <- suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
  if (!ok) failed <- c(failed, p)
}

if (length(failed) == 0) {
  cat("ALL_PACKAGES_OK\n")
} else {
  cat("FAILED_PACKAGES:", paste(failed, collapse = ", "), "\n")
  quit(status = 1)
}

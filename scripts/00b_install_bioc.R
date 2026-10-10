# 00b_install_bioc.R ----
# 修复脚本：Bioconductor包经清华TUNA镜像安装（bioconductor.org直连超时）
# 另：ggDCA在R4.6下架，改用dcurves
# 运行: Rscript 00b_install_bioc.R

options(timeout = 600)
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# ggDCA 替代
if (!requireNamespace("dcurves", quietly = TRUE)) install.packages("dcurves")

bioc_pkgs <- c(
  "GEOquery", "limma", "sva", "WGCNA",
  "clusterProfiler", "org.Hs.eg.db", "ReactomePA",
  "GSVA", "immunedeconv", "ComplexHeatmap"
)

to_install <- setdiff(bioc_pkgs, rownames(installed.packages()))
message("待安装Bioconductor包: ", paste(to_install, collapse = ", "))

if (length(to_install) > 0) {
  ok <- tryCatch({
    BiocManager::install(to_install, ask = FALSE, update = FALSE)
    TRUE
  }, error = function(e) { message("TUNA镜像失败: ", conditionMessage(e)); FALSE })

  if (!ok) {
    # 备用镜像：中科大
    options(BioC_mirror = "https://mirrors.ustc.edu.cn/bioconductor")
    BiocManager::install(to_install, ask = FALSE, update = FALSE)
  }
}

# 自检
all_pkgs <- c(bioc_pkgs, "dcurves")
failed <- character(0)
for (p in all_pkgs) {
  ok <- suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
  if (!ok) failed <- c(failed, p)
}
if (length(failed) == 0) cat("ALL_PACKAGES_OK\n") else cat("FAILED_PACKAGES:", paste(failed, collapse = ", "), "\n")

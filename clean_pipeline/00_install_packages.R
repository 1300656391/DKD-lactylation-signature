# 00_install_packages.R ----
# Install every R package used by steps 01-13. Run once before master.R:
#   Rscript 00_install_packages.R
# R >= 4.6 with Bioconductor >= 3.20; exact tested versions in sessionInfo_full.txt.

options(repos = c(CRAN = "https://cloud.r-project.org"))
N_CORES <- max(1, parallel::detectCores() - 1)

cran_pkgs <- c(
  "BiocManager", "remotes", "R.utils",
  "data.table", "ggplot2", "ggrepel", "patchwork", "viridis",
  "VennDiagram",
  "glmnet", "e1071", "xgboost", "pROC", "rms"
)
bioc_pkgs <- c(
  "GEOquery", "limma", "org.Hs.eg.db", "reactome.db",
  "WGCNA", "clusterProfiler", "ReactomePA", "GSVA",
  "ComplexHeatmap", "circlize", "Seurat"
)
github_pkgs <- c("MCPcounter" = "ebecht/MCPcounter")   # not on CRAN/Bioconductor

install_quietly <- function(pkgs, installer) {
  todo <- setdiff(pkgs, rownames(installed.packages()))
  if (length(todo) == 0) return(invisible(NULL))
  installer(todo)
}

install_quietly(cran_pkgs, function(p) install.packages(p, Ncpus = N_CORES))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install_quietly(bioc_pkgs, function(p) BiocManager::install(p, ask = FALSE, update = FALSE, Ncpus = N_CORES))
install_quietly(names(github_pkgs), function(p)
  remotes::install_github(github_pkgs[[p]], upgrade = "never"))

failed <- character(0)
for (p in c(cran_pkgs, bioc_pkgs, names(github_pkgs))) {
  if (!suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE))))
    failed <- c(failed, p)
}
if (length(failed) == 0) {
  cat("ALL_PACKAGES_OK\n")
} else {
  cat("FAILED_PACKAGES:", paste(failed, collapse = ", "), "\n")
  quit(status = 1)
}

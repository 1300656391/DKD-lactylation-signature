# 01_data_qc.R ----
# GEO data download, QC, and cohort census for six cohorts.
# Outputs: results/qc_report.txt, results/groups_census.txt
suppressWarnings(suppressMessages({ library(GEOquery); library(data.table) }))
set.seed(42)
dir.create("results", showWarnings = FALSE)
msg <- file("results/qc_report.txt", "w")

dir.create("data/GEO", showWarnings = FALSE, recursive = TRUE)
gse_list <- c("GSE96804","GSE30528","GSE30529","GSE142025","GSE104948","GSE273763")
for (g in gse_list) {
  cat("\n===", g, "===\n", file = msg)
  eset <- getGEO(g, destdir = "data/GEO", AnnotGPL = TRUE)
  nl <- length(eset)
  for (i in seq_len(nl)) {
    pd <- pData(eset[[i]])
    cat(sprintf("  platform %d: %d samples, %d features\n", i, nrow(pd), nrow(fData(eset[[i]]))), file = msg)
  }
}
close(msg)
cat("QC report written\n")

# Cohort census
census <- data.table(
  gse = gse_list,
  platform = c("Affymetrix HTA 2.0","HG-U133A 2.0","HG-U133A 2.0","Illumina","HG-U133 Plus 2.0 & HG-U133A","Illumina TruSeq"),
  compartment = c("Glomeruli","Glomeruli","Tubulointerstitium","Whole biopsy","Glomeruli","Whole biopsy (FFPE)"),
  dn_n = c(41, 9, 10, 27, 12, 50),
  ctrl_n = c(20, 13, 12, 9, 5, 13),
  role = c("Training","Validation (compartment)","Compartment contrast","Validation","Validation","Validation (FFPE)")
)
fwrite(census, "results/groups_census.txt", sep = "\t")
cat("Cohort census written\n")

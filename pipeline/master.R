# master.R ----
# Reproduces the full analysis for the DKD lactylation-signature paper.
# seed = 42 throughout. Run from the pipeline/ directory.
# Each step is skipped if its primary output already exists.

cat("=== DKD Lactylation-Signature Pipeline ===\n")
cat("R:", R.version.string, "\n\n")
set.seed(42)

run_if <- function(script, check = NULL) {
  if (!is.null(check) && file.exists(check)) {
    cat("[skip]", script, "\n"); return(invisible(NULL))
  }
  cat("[run ]", script, "\n")
  source(script, echo = FALSE)
}

run_if("01_data_qc.R",             "results/qc_report.txt")
run_if("02_lactylation_geneset.R", "results/lrg_summary.txt")
run_if("03_deg_wgcna.R",           "results/wgcna_modules.csv")
run_if("04_candidates.R",          "results/A_candidates.csv")
run_if("05_enrichment.R",          "results/enrich_GO_BP.csv")
run_if("06_immune_ssgsea.R",       "results/immune_ssgsea.csv")
run_if("07_ml_screening.R",        "results/hub_genes.txt")
run_if("08_model_freeze.R",        "results/frozen/verify_report.txt")
run_if("09_validation.R",          "results/validation_summary.csv")
run_if("10_snrna_localization.R",  "results/snrna_pdk_celltype.csv")
run_if("11_immune_replication.R",  "results/replication_summary.csv")
run_if("12_supplementary.R",       "results/supplementary_summary.csv")
run_if("13_figures.R",             "figures/Fig9_scrna_localization_rev.png")

cat("\n=== R pipeline complete ===\n")
cat("Next: python assemble_figures.py  (TIFF assembly)\n")

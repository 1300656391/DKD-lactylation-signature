# DKD Lactylation Signature — Analysis Code Repository

## Manuscript
Shi K, Li T. "Identification and validation of a lactate metabolism– and lactylation-related gene signature and its immune microenvironment correlates in diabetic kidney disease: a multi-cohort, machine-learning analysis." *Frontiers in Immunology*, submitted.

## Pipeline overview

| Stage | Script | Output |
|---|---|---|
| GEO QC & grouping | 01_geo_qc.R, 01b_groups_census.R, 01c_gse104948_diagnosis.R | groups_census.txt, qc_report.txt |
| Probe annotation & DEG | 02_probe_inspect.R, 04_preprocess_deg.R | A_DEGs.csv, A_all_genes_limma.csv |
| LRG construction | 03b_lactylation_geneset_local.R | LRG_geneset.csv |
| WGCNA | 05_wgcna.R | wgcna_*.csv |
| Intersection & enrichment | 07_intersection.R, 08_enrichment.R | intersection.txt, enrich_*.csv, gsea_Reactome.csv |
| Immune infiltration | 09_immune.R, 21_immune_replication.R | immune_*.csv, M7*.csv |
| ML model & validation | 10_ml_model.R | hub_genes.txt, AUC_summary.csv |
| Sensitivity & robustness analyses | 12*, 14*, 15*, 17*, 18*, 19* | results/revision/*.csv |
| Single-nucleus RNA-seq localization | 20_scrna_localization.R | M6_scrna_pdk_celltype.csv |
| Within-group immune correlations | 22_within_group_corr.R | M7d/M7e*.csv |
| Final figures (from frozen objects) | 13_figures_rev.R, 16_fig_round3.R | results/figures/*_rev.png |
| Submission figure fixes (panel tags, annotation placement, 300-dpi assembly) | 13b_fig2b_fix.R, 13c_fig4_regen.R, 16b_fig3_fix.R, 20b_fig9_tags.R, py_merge_figs.py, py_fig7_relayout.py, py_export_fig34.py | results/figures/tiff_300/* |
| Model freezing & verification | 15_freeze_ml.R, 15c_bootstrap_final.R, 17_verify_final.R | results/frozen/* |

## Frozen model objects (single source of truth)
`results/frozen/`: final_model_coefficients.csv (canonical deviance-criterion coefficients:
intercept 1.86, PDK2 -1.77, PDK4 -2.44), LASSO/SVM-RFE/XGBoost curves, verify_report.txt
(6/6 scenarios reproduce Table 2 AUC and Table 3 Youden statistics), sessionInfo.txt.

## Key findings reproduced by this code
- **Signature**: PDK2 + PDK4 (two-gene penalized logistic model, optimism-corrected AUC 0.998)
- **External validation**: 5 cohorts (GSE96804/GSE30528/GSE30529/GSE142025/GSE104948), AUCs 0.906-1.000
- **Immune replication** (GSE142025): 10/18 cell types significant; macrophage increase replicated (P=0.039, stage-dependent); PDK2-macrophage Spearman rho = -0.84 (all) / -0.78 (DN only)
- **snRNA localization** (GSE131882, 29,048 nuclei): PDK4 highest in podocytes (64.0% detection)
- **Control-source sensitivity** (GSE104948): living-donor AUC=1.000; hypertensive-nephropathy AUC=0.680

## Script execution order
01->...->11 (original pipeline) -> 12/12b/12c -> 13 (figures)
-> 14/14b (refinements) -> 15/15b/15c/15d (model freezing + bootstrap)
-> 16 (figures from frozen objects only) -> 17 (final verification + sessionInfo)
-> 18 (single-gene reference analyses) -> 19 (control-source sensitivity)
-> 20 (snRNA localization) -> 21 (immune replication) -> 22 (within-group correlations)
-> 13b/13c/16b/20b + py assembly scripts (submission figures: panel tags, titles, 300-dpi one-file-per-figure assembly)

## File roles in results/
- `results/*.csv` — primary analysis outputs (manuscript data sources)
- `results/revision/*.csv` — additional sensitivity and robustness analyses (M1-M7 series)
- `results/frozen/*.csv|txt` — frozen model objects; all reported statistics and figures derive from these
- Historical failed-run logs (run_freeze.txt, run_freeze_patch.txt, run_bootstrap.txt, run_bootstrap2.txt)
  document the debugging chain and are NOT part of the final pipeline; verify_report.txt is authoritative.

## NOTE on run_bootstrap2.txt
Its "Selection frequency top8" section uses an erroneous denominator (valid_B = 168, producing
PDK4 = 1.14). The correct per-gene frequencies (count and count/200) are in
results/revision/M4_selection_frequency.csv; see script 15d for the fix.

## Reproducibility
R 4.6.1 (see results/frozen/sessionInfo.txt), seed = 42, no internet required after data download.
GEO accessions: GSE96804, GSE30528, GSE30529, GSE142025, GSE104948, GSE131882.

## License
MIT License. See LICENSE file.

NOTE: GSE131882_seurat_annotated.rds (~650 MB) is NOT shipped; regenerate via
scripts/20_scrna_localization.R (~10 min).

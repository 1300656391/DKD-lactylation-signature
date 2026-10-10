# Code Repository README (Draft — 部署 GitHub/Zenodo 时使用)

## Manuscript
Shi K, Li T. "Identification and validation of a lactate metabolism– and lactylation-related gene signature for the diagnosis of diabetic kidney disease by integrating multiple cohorts and machine learning algorithms." Medicine (Baltimore), submitted.

## Pipeline overview
| Stage | Script | Output |
|---|---|---|
| GEO QC & grouping | 01_geo_qc.R, 01b_groups_census.R, 01c_gse104948_diagnosis.R | groups_census.txt, qc_report.txt |
| Probe annotation & DEG | 02_probe_inspect.R, 04_preprocess_deg.R | A_DEGs.csv, A_all_genes_limma.csv |
| LRG construction | 03_lactylation_geneset.R (online route, superseded), 03b_lactylation_geneset_local.R (used) | LRG_geneset.csv |
| WGCNA | 05_wgcna.R | wgcna_*.csv |
| Intersection & enrichment | 07_intersection.R, 08_enrichment.R | intersection.txt, enrich_*.csv, gsea_Reactome.csv |
| Immune infiltration | 09_immune.R | immune_*.csv |
| ML model & validation | 10_ml_model.R | hub_genes.txt, AUC_summary.csv |
| Figures (v1) | 11_figures.R | superseded by 13/16 |
| Revision analyses (R1-R5 review rounds) | 12*, 14*, 15*, 17*, 18* | results/revision/*.csv, results/frozen/* |
| Final figures (from frozen objects) | 13_figures_rev.R, 16_fig_round3.R | results/figures/*_rev.png |

## Frozen model objects (single source of truth)
`results/frozen/`: final_model_coefficients.csv (canonical deviance-criterion coefficients:
intercept 1.86, PDK2 −1.77, PDK4 −2.44), LASSO/SVM-RFE/XGBoost curves, verify_report.txt
(6/6 scenarios reproduce Table 2 AUC and Table 3 Youden statistics), sessionInfo.txt.

## Script execution order (revision chain)
01→…→11 (original pipeline) → 12/12b/12c (round-1 revision analyses) → 13 (round-1 figures)
→ 14/14b (round-2 fixes) → 15/15b/15c/15d (model freezing + bootstrap; see file headers for
known-issue notes on 15/15b, resolved in 15c) → 16 (round-3 figures, reads frozen objects
only) → 17 (final verification report + sessionInfo) → 18 (single-gene reference analyses).

## File roles in results/
- `results/*.csv` — primary analysis outputs (manuscript data sources)
- `results/revision/*.csv` — reviewer-requested supplementary analyses (M1–M9 series, incl. M5 control-source sensitivity)
- `results/frozen/*.csv|txt` — frozen model objects; all reported statistics and figures derive from these
- Historical failed-run logs (run_freeze.txt, run_freeze_patch.txt, run_bootstrap.txt, run_bootstrap2.txt) document
  the debugging chain and are NOT part of the final pipeline; verify_report.txt is authoritative.
- NOTE on run_bootstrap2.txt: its "Selection frequency top8" section uses an erroneous denominator
  (valid_B = 168, producing PDK4 = 1.14). The correct per-gene frequencies (count and count/200) are in
  results/revision/M4_selection_frequency.csv and supplementary Table S9b; see script 15d for the fix.

## Reproducibility
R 4.6.1 (see results/frozen/sessionInfo.txt), seed = 42, no internet required after data download.
GEO accessions: GSE96804, GSE30528, GSE30529, GSE142025, GSE104948.
NOTE: GSE131882_seurat_annotated.rds (~650 MB) is NOT shipped; regenerate via scripts/20_scrna_localization.R (~10 min).

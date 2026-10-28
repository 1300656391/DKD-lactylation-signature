# DKD Lactylation-Signature Analysis Pipeline

Code for: *A lactate–lactylation axis gene signature for diabetic kidney disease: diagnostic model development, six-cohort evaluation, and immune microenvironment correlates* (Frontiers in Immunology, submitted).

## Requirements

- R ≥ 4.6 with Bioconductor ≥ 3.20 — install all packages with `Rscript 00_install_packages.R` (exact tested versions in `sessionInfo_full.txt`)
- Python ≥ 3.10 with `Pillow` (`pip install -r python_requirements.txt`)

## Reproduction

```bash
Rscript 00_install_packages.R     # once
Rscript master.R                  # analysis + figures (~2-4 h)
python assemble_figures.py        # TIFF assembly (300 dpi, <= 179 mm)
```

Random seed = 42 throughout. Each step checks for its output before running; delete `results/` to force a full re-run.

## Pipeline

| Step | Script | Output |
|---|---|---|
| 0 | `00_install_packages.R` | R package installation |
| 1 | `01_data_qc.R` | Cohort QC and group census |
| 2 | `02_lactylation_geneset.R` | LRG list validation (n = 171 genes) |
| 3 | `03_deg_wgcna.R` | DEGs (limma) + WGCNA disease modules |
| 4 | `04_candidates.R` | Candidate LRG pool (n = 23) |
| 5 | `05_enrichment.R` | GO / Reactome / GSEA enrichment + per-module GO (S19) |
| 6 | `06_immune_ssgsea.R` | ssGSEA immune infiltration (training cohort) |
| 7 | `07_ml_screening.R` | LASSO + SVM-RFE + XGBoost consensus → 5 hub genes |
| 8 | `08_model_freeze.R` | Penalized model freeze + 200-sample bootstrap |
| 9 | `09_validation.R` | Six-cohort AUC / sens / spec / calibration / DCA |
| 10 | `10_snrna_localization.R` | Single-nucleus RNA-seq cell-type localization |
| 11 | `11_immune_replication.R` | GSE142025 replication + partial Spearman + MCP-counter (S21) |
| 12 | `12_supplementary.R` | Fig8B BH corrections (S20) + head-to-head (S18) |
| 13 | `13_figures.R` | All 9 main figures (PNG panels) |
| 14 | `assemble_figures.py` | TIFF assembly |

## Data

`data/LRG_geneset.csv`, `data/signature_collection.rda`, and `data/mcp_genes.tsv` ship with the repository. Everything else is downloaded from GEO and placed at:

| Path | GEO source |
|---|---|
| `data/GEO/<GSE>_series_matrix.txt.gz` | fetched automatically by step 1 (GSE104948 has two per-platform files) |
| `data/GEO/GPL571_family.soft.gz`, `GPL22945_family.soft.gz`, `GPL24120_family.soft.gz` | SOFT family files of the platform records |
| `data/GEO/GSE96804_family.soft.gz` | SOFT family file of the GSE96804 series record |
| `data/GEO/GSE142025_raw/` | extract `GSE142025_RAW.tar` here |
| `data/GSE273763/GSE273763_DKD_GRCH38_forGEO_July2026_v2.txt.gz`, `data/GSE273763/GSE273763_countsmatrix.txt.gz` | GSE273763 supplementary files |
| `data/GSE294519/GSE294519_Normalization_data_allsamples.txt.gz` | GSE294519 supplementary file |
| `data/scrna/GSM3823939–GSM3823944_*.dgecounts.rds.gz` | GSE131882 supplementary files |

## Key outputs

- `results/frozen/final_model_coefficients.csv` — frozen model (intercept, PDK2, PDK4)
- `results/frozen/verify_report.txt` — AUC / sensitivity / specificity verification
- `results/hub_genes.txt` — 5 hub genes
- `results/supplement/` — supplementary tables S5, S6, S12–S16, S18–S21 (S2–S4, S7–S9 are copies of primary outputs)
- `figures/*.tif` — 9 publication-ready figures

Exploratory sensitivity analyses (S10, S11, S15a/b, S17, DeLong comparisons, head-to-head paired bootstrap) ship with their frozen outputs in `additional_analyses/` (see its README); they are not part of the one-command pipeline.

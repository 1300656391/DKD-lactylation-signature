# code_repo_B — Orthogonal-layer analysis scripts and frozen outputs

Scripts supporting Sections 2.9/2.10 and 3.9/3.10/3.11 of the manuscript: two public
DKD kidney proteomes (PRIDE PXD041884; iProX IPX0003092000/PXD040617), proximal-tubule
DNA methylomes (GEO GSE233758 primary, GSE145745 replication), and perturbational
connectivity (iLINCS/L1000; queried 3 November 2026, saved API responses included for
the LIB_2 library). Python 3.12 with pandas/NumPy/SciPy; additionally openpyxl (B1),
requests (B3), rdata (B4).

Each module folder holds the analysis scripts together with the frozen output tables
they produced. `assemble_esm_tables.py` (package root) rebuilds the submission
supplementary tables S22–S27 (files ESM_28.csv–ESM_33.csv) byte-for-byte from those
outputs. The layers are analyzed as orthogonal evidence, not validation cohorts.

## B1_protein_PXD041884_iProX
| Script | Role | Output |
|---|---|---|
| analyze_pxd041884.py | Main PXD041884 FFPE proteome statistics (group medians, exact Wilcoxon, percentiles); also writes the PDH/pyruvate context gene table | TableB1_protein_pxd041884.csv (S22 core), TableB1_context_pdh_pyruvate.csv |
| map_uniprot2.py | UniProt accession -> gene symbol mapping via UniProt REST | uniprot_mapping.json |
| lactate_axis_panel.py | 36-node lactate-axis panel flags for the FFPE cohort | TableB1_lactate_axis_panel.csv (S22 panel flags) |
| acetylome_inventory.py | Acetyl-lysine peptide inventory, both deposited scales, exact Wilcoxon + BH | TableB1_acetylome_inventory.csv (= S23) |
| pepxml_axis_psms.py | pep.xml rank-1 PSM evidence sweep for candidate/axis nodes | pepxml_axis_psms.json (PSM tier evidence of Sec. 3.9) |
| iprox_analyze.py | Glomerular iProX matrix analysis: Fisher detection tests, detected-only Mann-Whitney, detection-dominated rule, median-centering sensitivity, global-context set | TableB1_iprox_panel.csv + TableB1_iprox_global_context.csv (= S24) |

## B4_methylation_GSE233758
| Script | Role | Output |
|---|---|---|
| s0_peek_header.py | Series-matrix header reconnaissance (sample annotation, batch tranches) | header_meta.txt |
| s1_parse_beta.py | Parse GSE233758 beta matrix (865,918 x 33) | beta_matrix.pkl |
| s2_extract_anno.py | EPIC ilm10b4 hg19 annotation extraction (RefGene groups, islands, SNPs) | anno_epic_full.pkl |
| s3_dnam_analysis.py | Per-probe OLS on M-values, region-level t tests, batch-adjusted sensitivity | probe_stats_raw.csv |
| s4_snp_qc.py | SNP-overlap QC (dbSNP147 common, MAF>5%) | snp_flags.csv |
| s5_make_table.py | Assemble TableB4_pdk_methylation (61 probes + SUMMARY rows); also writes region_level_stats.csv and probe_stats_snp_clean.csv | TableB4_pdk_methylation.csv (S25 core), region_level_stats.csv, probe_stats_snp_clean.csv |
| finalize_pdk_summary.py | Restructure the SUMMARY rows into the self-documenting best-probe scheme shipped with the submission (run after s5) | TableB4_pdk_methylation.csv (final) |
| s6_replication.py | GSE145745 untreated-arm replication pipeline (identical M-value/OLS, own 61-probe BH) | TableB4_replication.csv (S25 replication block) |
| s6b_region_repl.py | Region-level replication tests (console output; no file written) | — |
| lrg_methylation_sweep.py | Axis-wide 171-gene promoter/body sweep (per-gene averaged beta, Wilcoxon, BH within class) | TableB4_lrg_list_methylation.csv (= S26 primary block) |
| lrg_methylation_replication.py | LRG sweep replication in GSE145745 | TableB4_lrg_replication.csv (S26 replication block) |

## B3_ilincs_perturbation
| Script | Role | Output |
|---|---|---|
| s3_signature.py | Build the training-DEG query signature payload (top-150 up/down) | B3_signature_up.csv, B3_signature_down.csv, signature_payload.json |
| s22_final_queries.py | Definitive findConcordancesSC query vs LINCS L1000 LIB_5; builds the full-DEG signature (1,259 up / 1,386 down) inline | concordances_LIB_5_final_fullDEG.json etc. |
| s23_rank.py | Rank negative-similarity (reverser) rows | *_rankedNeg.json |
| s24_meta.py | Batched SignatureMeta metadata fetch (compound, cell line, dose, time, MoA) | signaturemeta_all.json |
| s28_build_table.py | Enrich ranked rows with metadata; per-compound aggregation | enriched_rows.json, aggregate_compounds.json |
| s34_build_table.py | Build the deliverable table from saved responses (no API calls) | TableB3_ilincs_reversers.csv (= S27) |

`downloads/concordances_LIB_2.json` and `downloads/concordances_LIB_2_rankedNeg.json`
(included, pre-placed) are the archived 3-November-2026 iLINCS API responses for the
LIB_2 library used by s23/s24/s28. clue.io retired January 2026; similarity = iLINCS weighted-KS
concordance, not tau.

## Rerunning from raw data

Place the inputs in each module folder (`downloads/` subdirectory unless noted), plus:

| Module | Inputs |
|---|---|
| B1 | PRIDE PXD041884: `downloads/Secocnd_Extraction_GPQ_Norm_TS_All_Peptides_FFPE_Protein_Rollup_Workbook.xlsx`, `downloads/Fritz_Acetyl_Quantitation_FFPE_Area2.xlsx`, 7 `downloads/*.pep.xml`; iProX IPX0003092000/PXD040617: `downloads_iprox/JTOM-2023-DN-Profiling.csv`, `downloads_iprox/iprox_IPX0003092000.xml` |
| B4 | GEO `downloads/GSE233758_series_matrix.txt.gz`, `downloads/GSE145745_series_matrix.txt.gz`; Illumina EPIC annotation (Bioconductor `IlluminaHumanMethylationEPICanno.ilm10b4.hg19`) unpacked under `downloads/IlluminaHumanMethylationEPICanno.ilm10b4.hg19/data/` (incl. `SNPs.147CommonSingle.rda`); the LRG gene list as `ESM_4.csv` (this submission's supplementary material) in the module folder |
| B3 | The full training DEG list as `ESM_3.csv` (this submission's supplementary material) in the module folder; live iLINCS API for s22/s24 (saved responses included for offline reruns) |

Run order: B1 map_uniprot2 -> analyze_pxd041884 -> lactate_axis_panel -> acetylome_inventory ->
pepxml_axis_psms -> iprox_analyze; B4 s0 -> s1 -> s2 -> s3 -> s4 -> s5 ->
finalize_pdk_summary -> s6 -> s6b -> lrg sweeps; B3 s3 -> s22 -> s23 -> s24 -> s28 ->
s34. Then `python assemble_esm_tables.py` rebuilds ESM_28-ESM_33.

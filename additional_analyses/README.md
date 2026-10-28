# Additional exploratory analyses

One-off sensitivity and exploratory analyses referenced in the manuscript
that are not part of the 14-step `clean_pipeline/`. Each sub-folder holds the
original analysis script (run from a project root that contains `data/`,
historical comments included) together with its frozen output table.

| Folder | Manuscript item | Script | Frozen output |
|---|---|---|---|
| `single_vs_signature_S10/` | Table S10 - single-gene versus signature AUCs with DeLong tests | `18_round5_fixes.R` (section M2) | `TableS10_single_vs_signature.csv` |
| `control_sensitivity_S11/` | Table S11 - GSE104948 control-source sensitivity (living-donor / hypertensive-nephropathy controls) | `19_ld_sensitivity.R` | `TableS11_control_sensitivity.csv` |
| `disease_arbitration_S15a/` | Table S15a - disease arbitration, proxy correlations, PDK2-PDK4 asymmetry (10,000-resample dependent-correlation bootstrap + Steiger z) | `24_disease_arbitration.R` | `TableS15a_disease_arbitration.csv` |
| `cytokine_podocyte_S15b/` | Table S15b - cytokine-axis and podocyte-module correlations | `26c_s15b_complete.R` | `TableS15b_cytokine_podocyte_module.csv` |
| `gse273763_immune_nnls_S17/` | Table S17 - GSE273763 ssGSEA (two normalizations) and NNLS B-cell fractions | `step1_immune_gse273763.R`, `step1b_immune_nnls_gate.py`, `step1c_s17_corr_full.py` | `TableS17_gse273763_immune.csv` + intermediate CSVs |
| `gse273763_delong/` | GSE273763 DeLong comparisons (PDK2/PDK4 versus signature) | `46_gse273763_score.py`, `47_gse273763_nnls.py`, `48_gse273763_s17.py` | `verify_report_gse273763.txt` |
| `head_to_head_bootstrap/` | Head-to-head paired bootstrap (signature versus PTGS2, B = 2000) | `head_to_head.py`, `49_head_to_head.py` | `signature_vs_PTGS2_bootstrap.csv` |

These analyses read the same frozen two-gene model as the main pipeline
(intercept 1.85981, PDK2 -1.76987, PDK4 -2.43690; see
`clean_pipeline/results/frozen/final_model_coefficients.csv` after running
step 8, or the frozen copy shipped with the repository).

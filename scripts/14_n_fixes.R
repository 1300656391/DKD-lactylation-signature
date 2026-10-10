# 14_n_fixes.R ----
# 按review_round2.md处理N1-N8：
#   N2: 生成 results/revision/M9_immune_tests_no_mucosa.csv（18行，BH在18检验族内重算）
#   N8预防: Fig2C 加严非编码/位点标签过滤后重绘（DANCR/-AS1/LOC/OTTHUM/NBPF/P\d伪基因）
#   m14: 补充表S1-S7实体组装到 manuscript/supplement/
#   N5: 摘要字数实测
suppressWarnings(suppressMessages({
  library(data.table); library(ComplexHeatmap); library(circlize)
}))
set.seed(42)

# ============ N2: M9 免疫检验表（剔除Normal_mucosa，18行，BH重算） ============
it <- read.csv("results/immune_tests.csv", stringsAsFactors = FALSE, check.names = FALSE)
it18 <- it[it$cell != "Normal_mucosa", ]
it18$p.adj_BH18 <- p.adjust(it18$p, method = "BH")
it18$sig_BH18 <- ifelse(it18$p.adj_BH18 < 0.05, "yes", "no")
dir.create("results/revision", showWarnings = FALSE, recursive = TRUE)
write.csv(it18, "results/revision/M9_immune_tests_no_mucosa.csv", row.names = FALSE)
cat("N2 done:", nrow(it18), "rows; significant (BH<0.05):", sum(it18$sig_BH18 == "yes"), "\n")
# 星号分档核对（供回应信）：图8A用19检验族BH，此表为18检验族BH
it18$star_19family <- ifelse(it18$p.adj < 0.001, "***",
                      ifelse(it18$p.adj < 0.01, "**",
                      ifelse(it18$p.adj < 0.05, "*", "ns")))
cat("Star tiers under 19-test BH (as plotted):\n"); print(it18[, c("cell","p.adj","star_19family")])

# ============ N8预防: Fig2C 加严过滤重绘 ============
train_expr <- readRDS("data/GSE96804_expr_symbol.rds")
train_ph   <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
deg <- fread("results/A_DEGs.csv")
nc_filter <- function(g) {
  grepl("^(MIR|SNOR|SNORA|SNORD|SCAR|RNU|LINC|RNA|DANCR|VTRNA|NPIPA|NBPF|LOC[0-9]|OTTHUM|ENSG)", g) |
    grepl("-AS[0-9]$", g) | grepl("P[0-9]$", g)
}
top_pool <- deg[order(deg$adj.P.Val)][1:120]
top50 <- top_pool[!nc_filter(top_pool$gene)][1:50]
excluded <- top_pool$gene[nc_filter(top_pool$gene)]
cat("Fig2C excluded artifact labels (top-120 pool):", paste(excluded, collapse = ", "), "\n")
mat <- t(scale(t(train_expr[top50$gene, ])))
ha <- HeatmapAnnotation(Group = train_ph$group,
                        col = list(Group = c(Control = "#2166AC", DN = "#B2182B")))
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
png("results/figures/Fig2C_heatmap_rev.png", type = "windows", width = 1800, height = 2400, res = 250)
draw(Heatmap(mat, name = "z-score", col = col_fun, top_annotation = ha,
             show_column_names = FALSE, row_names_gp = gpar(fontsize = 6),
             column_title = "Top 50 protein-coding DEGs (GSE96804)"))
dev.off()
cat("Fig2C redrawn:", nrow(top50), "rows\n")

# ============ m14: 补充表S1-S7实体组装 ============
supdir <- "manuscript/supplement"
dir.create(supdir, showWarnings = FALSE, recursive = TRUE)
# S1 队列汇总
s1 <- data.frame(
  cohort = c("GSE96804","GSE30528","GSE30529","GSE142025","GSE104948"),
  platform = c("Affy Human Transcriptome Array 2.0 (GPL17586)",
               "Affy HG-U133A 2.0 (GPL571)","Affy HG-U133A 2.0 (GPL571)",
               "RNA-seq","Affy U133 Plus2 & U133A, Brainarray ENTREZG CDF (GPL22945/GPL24120)"),
  compartment = c("Glomeruli","Glomeruli","Tubulointerstitium","Whole biopsy","Glomeruli"),
  group_n = c("DN 41 / Control 20","DKD 9 / Control 13","DKD 10 / Control 12",
              "Advanced DN 21 / Early DN 6 / Control 9","DN 12 / Control 5"),
  role = c("Training","Compartment-level validation (glomeruli)",
           "Compartment-level analysis (tubulointerstitium)",
           "External validation (+ staging)","External validation (platform-stratified)"),
  note = c("European biopsy biobank","Same LCM biopsy series as GSE30529",
           "Same LCM biopsy series as GSE30528","Includes early-stage DN",
           "Donor overlap with GSE96804 cannot be excluded")
)
write.csv(s1, file.path(supdir, "TableS1_cohorts.csv"), row.names = FALSE)
file.copy("results/A_DEGs.csv", file.path(supdir, "TableS2_DEGs_full.csv"), overwrite = TRUE)
file.copy("data/LRG_geneset.csv",  file.path(supdir, "TableS3_LRG_list.csv"),  overwrite = TRUE)
file.copy("results/A_candidates.csv", file.path(supdir, "TableS4_candidates.csv"), overwrite = TRUE)
# S5 各队列缺失基因
auc <- read.csv("results/AUC_summary.csv", stringsAsFactors = FALSE)
s5 <- auc[, c("cohort","n","missing_genes")]
write.csv(s5, file.path(supdir, "TableS5_missing_genes_by_cohort.csv"), row.names = FALSE)
file.copy("results/revision/M6_fixed_cutoff_performance.csv",
          file.path(supdir, "TableS6_training_cutoff_transfer.csv"), overwrite = TRUE)
file.copy("results/revision/M9_immune_tests_no_mucosa.csv",
          file.path(supdir, "TableS7_immune_tests_18cells.csv"), overwrite = TRUE)
cat("Supplement S1-S7 assembled in", supdir, "\n")

# ============ N5: 摘要字数实测 ============
ab <- c(
"Background: Diabetic kidney disease (DKD) is the leading cause of end-stage renal disease, and early recognition remains challenging. Lactate accumulation and protein lactylation link metabolic reprogramming to epigenetic regulation in renal fibrosis, yet lactylation-oriented biomarker discovery for DKD is lacking.",
"Methods: Five GEO cohorts of human kidney tissue (total n = 158) were analyzed. GSE96804 (diabetic nephropathy [DN], n = 41; controls, n = 20) served as training. A lactylation-related gene (LRG) list (171 genes; lactate-metabolic enzymes plus lactylation machinery/substrates) was compiled from Gene Ontology, Reactome, and literature. Differentially expressed genes (DEGs) were identified by limma and co-expression modules by WGCNA; candidate LRGs were screened by LASSO, SVM-RFE, and XGBoost; a penalized logistic signature was validated in four cohorts, and immune infiltration was quantified by single-sample GSEA.",
"Results: Of 2,645 DEGs, 23 candidate LRGs (11 differentially expressed) converged on five hub genes: ACSF2, HIST1H4B, ME1, PDK2, and PDK4. PDK2 was selected by all three algorithms; PDK4 by LASSO and XGBoost (highest gain). A two-gene (PDK2 + PDK4) signature showed complete class separation in training and AUCs of 0.906 in glomerular GSE30528, 0.992 in RNA-seq GSE142025 (early DN 0.889), and 1.000 in GSE104948, with attenuated performance in the tubulointerstitial compartment (GSE30529, 0.650). At cohort-specific Youden cutoffs, sensitivity/specificity were 0.78/1.00 (GSE30528) and 0.96/1.00 (GSE142025). Eleven of 18 immune cell types differed significantly; hub genes correlated with myeloid infiltration.",
"Conclusions: To our knowledge, this is the first multi-cohort, machine-learning-based lactylation-axis signature for DKD diagnosis. The PDK2/PDK4 pair warrants prospective evaluation as an early diagnostic biomarker."
)
cat("Abstract word count (incl. section labels):",
    sum(sapply(c("Background","Methods","Results","Conclusions"), function(x) length(strsplit(x," ")[[1]])) +
        sapply(ab, function(x) length(strsplit(x," ")[[1]]))), "\n")
cat("Abstract word count (body only):", sum(sapply(ab, function(x) length(strsplit(x," ")[[1]]))), "\n")

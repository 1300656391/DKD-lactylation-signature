# 22_within_group_corr.R ----
# R9 R3：ρ=−0.84的组结构混杂处理
#   DN组内 / 对照组内 / 按分期的 PDK2/PDK4×巨噬细胞Spearman（GSE142025）
#   + per-stage巨噬细胞中位数（补入S13）
#   + PDK2/PDK4不对称量化
suppressWarnings(suppressMessages({ library(data.table); library(GSVA) }))
set.seed(42)

# 重跑GSE142025 ssGSEA（复用21号脚本逻辑）
e <- new.env(); load("data/signature_collection.rda", envir = e)
sc <- get("signature_collection", envir = e)
cn <- grep("_Bindea_et_al$", names(sc), value = TRUE)
cn <- setdiff(cn, c("SW480_cancer_cells_Bindea_et_al"))
gsets <- sc[cn]; names(gsets) <- sub("_Bindea_et_al$", "", cn)
bindea18 <- c("B_cells","T_cells","T_helper_cells","Tcm","Tem","Th1_cells","Th2_cells","TFH",
  "CD8_T_cells","Cytotoxic_cells","NK_cells","NK_CD56dim_cells","DC","iDC","Eosinophils",
  "Macrophages","Mast_cells","Neutrophils")
sig18 <- lapply(gsets[intersect(bindea18, names(gsets))], function(g) unique(toupper(g)))

e14 <- readRDS("data/GSE142025_expr.rds")
ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
cl <- function(l) trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]]))
p14 <- data.frame(gsm = cl(ml[grepl("^!Sample_geo_accession", ml)][1]),
                  group = sub("^group:\\s*", "", cl(ml[grepl("^!Sample_characteristics_ch1", ml)][1])))
colnames(e14) <- p14$gsm[match(sub("_.*$", "", colnames(e14)), p14$gsm)]
rownames(e14) <- toupper(rownames(e14))
keep <- p14$group %in% c("Advanced_DN", "Early_DN", "Control")
e14k <- as.matrix(e14[, keep])
grp <- p14$group[keep]

par <- ssgseaParam(exprData = e14k, geneSets = sig18, normalize = TRUE)
ss <- gsva(par, verbose = FALSE)
ss <- as.data.frame(t(ss))
ss$group <- grp
ss$stage <- ifelse(grp == "Early_DN", "Early", ifelse(grp == "Advanced_DN", "Advanced", "Control"))

mac <- ss$Macrophages
pdk2 <- as.numeric(e14k["PDK2", ])
pdk4 <- as.numeric(e14k["PDK4", ])

# 组内相关
cors <- rbind(
  data.frame(subset = "All (DN+Ctrl, n=36)", PDK2_rho = suppressWarnings(cor(pdk2, mac, method="spearman")),
             PDK4_rho = suppressWarnings(cor(pdk4, mac, method="spearman"))),
  data.frame(subset = "DN only (n=27)", PDK2_rho = suppressWarnings(cor(pdk2[grp!="Control"], mac[grp!="Control"], method="spearman")),
             PDK4_rho = suppressWarnings(cor(pdk4[grp!="Control"], mac[grp!="Control"], method="spearman"))),
  data.frame(subset = "Control only (n=9)", PDK2_rho = suppressWarnings(cor(pdk2[grp=="Control"], mac[grp=="Control"], method="spearman")),
             PDK4_rho = suppressWarnings(cor(pdk4[grp=="Control"], mac[grp=="Control"], method="spearman"))),
  data.frame(subset = "Early DN (n=6)", PDK2_rho = suppressWarnings(cor(pdk2[grp=="Early_DN"], mac[grp=="Early_DN"], method="spearman")),
             PDK4_rho = suppressWarnings(cor(pdk4[grp=="Early_DN"], mac[grp=="Early_DN"], method="spearman"))),
  data.frame(subset = "Advanced DN (n=21)", PDK2_rho = suppressWarnings(cor(pdk2[grp=="Advanced_DN"], mac[grp=="Advanced_DN"], method="spearman")),
             PDK4_rho = suppressWarnings(cor(pdk4[grp=="Advanced_DN"], mac[grp=="Advanced_DN"], method="spearman")))
)
cors$PDK2_rho <- round(cors$PDK2_rho, 3)
cors$PDK4_rho <- round(cors$PDK4_rho, 3)
fwrite(cors, "results/revision/M7d_within_group_macro_corr.csv")
cat("== Within-group Spearman (GSE142025) ==\n"); print(cors)

# Per-stage巨噬细胞中位数
stg <- data.frame(
  group = c("Control", "Early DN", "Advanced DN"),
  n = c(sum(grp=="Control"), sum(grp=="Early_DN"), sum(grp=="Advanced_DN")),
  median_macrophage = c(round(median(mac[grp=="Control"]), 3),
                        round(median(mac[grp=="Early_DN"]), 3),
                        round(median(mac[grp=="Advanced_DN"]), 3)))
fwrite(stg, "results/revision/M7e_perstage_macrophage.csv")
cat("\n== Per-stage macrophage medians ==\n"); print(stg)

# PDK2 vs PDK4不对称：Fisher z检验
z1 <- atanh(cors$PDK2_rho[1]); z2 <- atanh(cors$PDK4_rho[1])
zdiff <- (z1 - z2) / sqrt(1/(36-3) + 1/(36-3))
p_asym <- 2 * pnorm(-abs(zdiff))
cat(sprintf("\nPDK2 vs PDK4 rho asymmetry (all-sample): z=%.2f, p=%.3f\n", zdiff, p_asym))

cat("\nWITHIN_GROUP_DONE\n")

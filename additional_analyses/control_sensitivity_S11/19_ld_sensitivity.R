# 19_ld_sensitivity.R ----
# Round-6 A5/B7 修复：GSE104948 对照来源敏感性分析
#   DN vs LD（活体供肾，21例）——直接检验"肿瘤切除对照"混杂：活体供肾无肿瘤/热缺血混杂
#   DN vs HT（高血压肾病，15例）——疾病特异性对照：签名是否区分 DKD 与另一肾脏病
#   全部使用冻结系数（canonical deviance），z标准化在各自分析子集内（与正文口径一致）
suppressWarnings(suppressMessages({ library(data.table); library(pROC) }))
hub <- readLines("results/hub_genes.txt")
fm <- fread("results/frozen/final_model_coefficients.csv")
bh <- setNames(fm$canonical_deviance, fm$term)
zrows <- function(m) t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
sig_of <- function(e) {  # e: genes x samples (DN+对照子集), 内部z标准化; 缺失基因系数为0自动跳过
  hg <- intersect(hub, rownames(e))
  Xv <- t(zrows(e[hg, , drop = FALSE]))
  as.vector(Xv %*% bh[hg] + bh[["(Intercept)"]])
}
p4_of <- function(e) as.numeric(zrows(as.matrix(e["PDK4", , drop = FALSE]))[1, ])

grp_from_title <- function(p) {
  ml <- readLines(paste0("data/GSE104948-", p, "_meta.txt"), warn = FALSE)
  gsm <- trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", ml[grepl("^!Sample_geo_accession", ml)][1]), "\t", fixed = TRUE)[[1]]))
  ttl <- trimws(gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", ml[grepl("^!Sample_title", ml)][1]), "\t", fixed = TRUE)[[1]]))
  data.frame(gsm = gsm, grp = sub("^[A-Z][0-9]+-Glom-", "", sub("[0-9]+$", "", ttl)), stringsAsFactors = FALSE)
}
load_expr <- function(p) readRDS(paste0("data/GSE104948-", p, "_expr_symbol.rds"))

res <- list(); i <- 0
add <- function(label, e, case_idx, ctrl_idx, note) {
  es <- as.matrix(e[, c(case_idx, ctrl_idx), drop = FALSE])
  y <- c(rep(1, length(case_idx)), rep(0, length(ctrl_idx)))
  r <- roc(y, sig_of(es), quiet = TRUE)
  r4 <- roc(y, p4_of(es) * ifelse(median(es["PDK4", y == 1]) < median(es["PDK4", y == 0]), -1, 1), quiet = TRUE)
  ci <- coords(r, "best", ret = c("sensitivity", "specificity"), best.method = "youden")
  i <<- i + 1
  res[[i]] <<- data.frame(
    comparison = label, n_case = length(case_idx), n_ctrl = length(ctrl_idx),
    AUC_signature = round(as.numeric(auc(r)), 3),
    CI_lower = round(ci.auc(r)[1], 3), CI_upper = round(ci.auc(r)[3], 3),
    AUC_PDK4 = round(as.numeric(auc(r4)), 3),
    youden_sens = round(ci$sensitivity, 3), youden_spec = round(ci$specificity, 3),
    note = note)
}

for (p in c("GPL22945", "GPL24120")) {
  g <- grp_from_title(p); e <- load_expr(p)
  stopifnot(all(colnames(e) == g$gsm))
  dn <- which(g$grp == "DN"); ld <- which(g$grp == "LD_"); ht <- which(g$grp == "HT"); tn <- which(g$grp == "TN")
  cat(p, ": DN", length(dn), "LD", length(ld), "HT", length(ht), "TN", length(tn), "\n")
  if (length(ld) > 0) add(paste0("DN_vs_LD_", p), e, dn, ld, "living-donor control (per-platform)")
  if (length(ht) > 0) add(paste0("DN_vs_HT_", p), e, dn, ht, "hypertensive-nephropathy disease control (per-platform)")
}
# 池化（两平台各自z标准化后合并分数）
g1 <- grp_from_title("GPL22945"); e1 <- load_expr("GPL22945")
g2 <- grp_from_title("GPL24120"); e2 <- load_expr("GPL24120")
pool <- function(grp1, grp2, label, note) {
  s_list <- list(); y_list <- list()
  for (cfg in list(list(g1, e1), list(g2, e2))) {
    g <- cfg[[1]]; e <- cfg[[2]]
    idx_c <- which(g$grp == grp1); idx_n <- which(g$grp == grp2)
    if (length(idx_c) == 0 || length(idx_n) == 0) next
    es <- e[, c(idx_c, idx_n), drop = FALSE]
    s_list[[length(s_list) + 1]] <- sig_of(es)
    y_list[[length(y_list) + 1]] <- c(rep(1, length(idx_c)), rep(0, length(idx_n)))
  }
  if (length(s_list) < 2) { cat(label, ": 单平台，跳过池化\n"); return(invisible()) }
  s <- unlist(s_list); y <- unlist(y_list)
  r <- roc(y, s, quiet = TRUE)
  ci <- coords(r, "best", ret = c("sensitivity", "specificity"), best.method = "youden")
  i <<- i + 1
  res[[i]] <<- data.frame(comparison = label, n_case = sum(y == 1), n_ctrl = sum(y == 0),
    AUC_signature = round(as.numeric(auc(r)), 3), CI_lower = round(ci.auc(r)[1], 3),
    CI_upper = round(ci.auc(r)[3], 3), AUC_PDK4 = NA,
    youden_sens = round(ci$sensitivity, 3), youden_spec = round(ci$specificity, 3), note = note)
}
pool("DN", "LD_", "DN_vs_LD_pooled", "living-donor control (pooled, per-platform z)")
pool("DN", "HT", "DN_vs_HT_pooled", "hypertensive disease control (pooled, per-platform z)")
out <- rbindlist(res)
fwrite(out, "results/revision/M5_gse104948_control_sensitivity.csv")
file.copy("results/revision/M5_gse104948_control_sensitivity.csv",
          "manuscript/supplement/TableS11_control_sensitivity.csv", overwrite = TRUE)
print(out)
cat("\nLD_SENSITIVITY_DONE\n")

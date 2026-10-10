# 04_preprocess_deg.R ----
# S3: 探针→基因映射 + 按基因折叠 + 各队列pheno
# S4: GSE96804训练集 limma DEGs
# 输出: data/<acc>_expr_symbol.rds + data/<acc>_pheno.csv + results/A_DEGs.csv
#       + results/figures/Fig2_qc_deg_draft.png
# 运行: Rscript 04_preprocess_deg.R

suppressWarnings(suppressMessages({
  library(data.table)
  library(limma)
  library(org.Hs.eg.db)
  library(ggplot2)
}))

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. GPL探针注释表解析
# ============================================================
find_table_line <- function(f) {
  con <- file(f, "r"); on.exit(close(con))
  n <- 0L
  repeat {
    l <- readLines(con, n = 1L, warn = FALSE)
    if (length(l) == 0L) stop("未找到!platform_table_begin: ", f)
    n <- n + 1L
    if (grepl("!platform_table_begin", l, fixed = TRUE)) return(n)
  }
}

read_gpl_table <- function(gpl) {
  f <- paste0("data/GEO/", gpl, "_full.txt")
  skip_n <- find_table_line(f)
  dt <- fread(f, skip = skip_n, sep = "\t", header = TRUE, quote = "",
              fill = TRUE, showProgress = FALSE)
  dt <- dt[!grepl("!platform_table_end", as.character(dt[[1]]), fixed = TRUE)]
  if (ncol(dt) < 2) stop(gpl, ": 平台表解析异常，仅", ncol(dt), "列")
  dt
}

first_token <- function(x) trimws(sub("\\/\\/\\/.*$", "", x))
first_field <- function(x, k) {          # gene_assignment第k个"//"字段
  a <- sub("\\/\\/\\/.*$", "", x)        # 只取第一个assignment
  f <- strsplit(a, "//", fixed = TRUE)
  vapply(f, function(p) {
    p <- trimws(p)
    if (length(p) >= k && p[k] != "") p[k] else NA_character_
  }, character(1))
}

build_map <- function(gpl) {
  dt <- read_gpl_table(gpl)
  id <- as.character(dt[[1]])
  sym <- rep(NA_character_, length(id))
  cn <- colnames(dt)
  if ("Gene Symbol" %in% cn)       sym <- first_token(dt[["Gene Symbol"]])
  else if ("Symbol" %in% cn)       sym <- first_token(dt[["Symbol"]])
  else if ("gene_assignment" %in% cn) sym <- first_field(dt[["gene_assignment"]], 2)  # 第2字段为symbol
  else if ("ENTREZ_GENE_ID" %in% cn) {
    ez <- trimws(as.character(dt[["ENTREZ_GENE_ID"]]))
    ez[!grepl("^[0-9]+$", ez)] <- NA
    sym <- mapIds(org.Hs.eg.db, keys = ez, column = "SYMBOL",
                  keytype = "ENTREZID", multiVals = "first")
  }
  m <- data.frame(probe = id, symbol = unname(sym), stringsAsFactors = FALSE)
  m <- m[!is.na(m$symbol) & m$symbol != "" & m$symbol != "---", ]
  m[!duplicated(m$probe), ]
}

# GSE96804专用：从family.soft内嵌平台表解析完整注释
# （acc.cgi平台表被截断7416行；hta20sttranscriptcluster.db已从Bioc除名）
build_map_familysoft <- function() {
  f <- "data/GEO/GSE96804_family.soft.gz"
  con <- gzfile(f, "rt"); on.exit(close(con))
  n <- 0L; begin_n <- NA_integer_; end_n <- NA_integer_
  repeat {
    l <- readLines(con, n = 1L, warn = FALSE)
    if (length(l) == 0L) break
    n <- n + 1L
    if (is.na(begin_n) && grepl("platform_table_begin", l, ignore.case = TRUE)) {
      begin_n <- n
    } else if (!is.na(begin_n) && is.na(end_n) &&
               grepl("platform_table_end", l, ignore.case = TRUE)) {
      end_n <- n; break
    }
  }
  if (is.na(begin_n) || is.na(end_n)) stop("family.soft中未找到平台表边界")
  cat("  平台表区间: 行", begin_n, "-", end_n, "\n")
  dt <- fread(f, skip = begin_n, nrows = end_n - begin_n - 1L, sep = "\t",
              header = TRUE, quote = "", fill = TRUE, showProgress = FALSE)
  cat("  family.soft平台表行数:", nrow(dt), "\n")
  sym <- first_field(as.character(dt$gene_assignment), 2)
  m <- data.frame(probe = sub("_st$", "", as.character(dt[[1]])),
                  symbol = sym, stringsAsFactors = FALSE)
  m <- m[!is.na(m$symbol) & m$symbol != "" & m$symbol != "---", ]
  m[!duplicated(m$probe), ]
}

cat("解析GPL注释...\n")
map_list <- list(
  GSE96804          = build_map_familysoft(),
  GSE30528          = build_map("GPL571"),
  GSE30529          = build_map("GPL571"),
  `GSE104948-GPL22945` = build_map("GPL22945"),
  `GSE104948-GPL24120` = build_map("GPL24120")
)
for (nm in names(map_list)) cat(" ", nm, "注释映射:", nrow(map_list[[nm]]), "探针\n")

# ============================================================
# 2. 表达矩阵映射+按基因折叠（多探针取均值最高者）
# ============================================================
collapse_to_symbol <- function(expr, map) {
  rn <- rownames(expr)
  rn_key <- sub("_st$", "", rn)                    # HTA数组探针后缀归一化
  keep <- rn_key %in% map$probe
  e <- expr[keep, , drop = FALSE]
  sym <- map$symbol[match(rn_key[keep], map$probe)]
  mu <- rowMeans(e, na.rm = TRUE)
  ord <- order(sym, -mu)
  e <- e[ord, , drop = FALSE]; sym <- sym[ord]
  first_idx <- !duplicated(sym)
  out <- e[first_idx, , drop = FALSE]
  rownames(out) <- sym[first_idx]
  out
}

map_and_save <- function(acc) {
  expr <- readRDS(paste0("data/", acc, "_expr.rds"))
  m <- collapse_to_symbol(expr, map_list[[acc]])
  saveRDS(m, paste0("data/", acc, "_expr_symbol.rds"))
  cat(" ", acc, ":", nrow(m), "基因 x", ncol(m), "样本\n")
  m
}

for (acc in names(map_list)) map_and_save(acc)

# ============================================================
# 3. 各队列pheno
# ============================================================
clean_v <- function(x) trimws(gsub('^"|"$', "", x))
split_line <- function(l) {
  if (grepl("=", l, fixed = TRUE)) {
    s <- sub("^[^=]*=", "", l)
  } else {
    s <- sub("^[^\t]*\t", "", l)
  }
  clean_v(strsplit(s, "\t", fixed = TRUE)[[1]])
}

# GSE96804: source_name前缀 diab/unaf
ml <- readLines("data/GSE96804_meta.txt", warn = FALSE)
src <- split_line(ml[grepl("^!Sample_source_name_ch1", ml)][1])
pheno96804 <- data.frame(
  sample = colnames(readRDS("data/GSE96804_expr_symbol.rds")),
  stringsAsFactors = FALSE)
pheno96804$group <- ifelse(grepl("^diab", src, ignore.case = TRUE), "DN", "Control")
write.csv(pheno96804, "data/GSE96804_pheno.csv", row.names = FALSE)
cat("GSE96804分组:", paste(names(table(pheno96804$group)), table(pheno96804$group), sep = "=", collapse = ", "), "\n")

# GSE30528/30529: title前缀 Diab/Cont
for (acc in c("GSE30528", "GSE30529")) {
  ml <- readLines(paste0("data/", acc, "_meta.txt"), warn = FALSE)
  ttl <- split_line(ml[grepl("^!Sample_title", ml)][1])
  ph <- data.frame(sample = colnames(readRDS(paste0("data/", acc, "_expr_symbol.rds"))),
                   stringsAsFactors = FALSE)
  ph$group <- ifelse(grepl("^Diab", ttl, ignore.case = TRUE), "DKD", "Control")
  write.csv(ph, paste0("data/", acc, "_pheno.csv"), row.names = FALSE)
  cat(acc, "分组:", paste(names(table(ph$group)), table(ph$group), sep = "=", collapse = ", "), "\n")
}

# GSE104948: diagnosis特征行 → 取DN与Tumor Nephrectomy
build_104948 <- function(acc) {
  ml <- readLines(paste0("data/", acc, "_meta.txt"), warn = FALSE)
  ch <- ml[grepl("^!Sample_characteristics_ch1", ml)]
  vals <- lapply(ch, function(l) split_line(l))
  expr <- readRDS(paste0("data/", acc, "_expr_symbol.rds"))
  n <- ncol(expr)
  maxlen <- max(vapply(vals, length, integer(1)))
  mat <- do.call(rbind, lapply(vals, function(v) { length(v) <- maxlen; v }))
  diagv <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    hits <- mat[, i][grepl("diagnosis:", mat[, i])]
    if (length(hits) > 0) diagv[i] <- sub("^diagnosis:\\s*", "", hits[1])
  }
  stopifnot(length(diagv) == n)
  cat(" ", acc, "diagnosis分布:", paste(names(table(diagv)), table(diagv), sep = "=", collapse = ", "), "\n")
  keep <- diagv %in% c("Diabetic Nephropathy", "Tumor Nephrectomy")
  list(expr = expr[, keep, drop = FALSE],
       grp = ifelse(diagv[keep] == "Diabetic Nephropathy", "DN", "Control"))
}
a <- build_104948("GSE104948-GPL22945")
b <- build_104948("GSE104948-GPL24120")
common <- intersect(rownames(a$expr), rownames(b$expr))
e104948 <- cbind(a$expr[common, , drop = FALSE], b$expr[common, , drop = FALSE])
colnames(e104948) <- make.unique(colnames(e104948))
ph104948 <- data.frame(sample = colnames(e104948),
                       group = c(a$grp, b$grp), stringsAsFactors = FALSE)
saveRDS(e104948, "data/GSE104948_expr_symbol.rds")
write.csv(ph104948, "data/GSE104948_pheno.csv", row.names = FALSE)
cat("GSE104948合并分组:", paste(names(table(ph104948$group)), table(ph104948$group), sep = "=", collapse = ", "), "\n")

# ============================================================
# 4. S4: 训练集DEGs（GSE96804, DN vs Control）
# ============================================================
expr <- readRDS("data/GSE96804_expr_symbol.rds")
ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
stopifnot(identical(ph$sample, colnames(expr)))
group <- factor(ph$group, levels = c("Control", "DN"))
design <- model.matrix(~ group)
fit <- lmFit(expr, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef = "groupDN", number = Inf, adjust.method = "BH")
tt$gene <- rownames(tt)
deg <- tt[abs(tt$logFC) >= 0.5 & tt$adj.P.Val < 0.05, ]
deg <- deg[order(deg$adj.P.Val), ]
setcolorder <- NULL
fwrite(as.data.frame(tt), "results/A_all_genes_limma.csv")
fwrite(as.data.frame(deg), "results/A_DEGs.csv")

cat("\n==== DEG统计 ====\n")
cat("检测基因数:", nrow(tt), "\n")
cat("显著DEG数(|logFC|>=0.5 & adj.P<0.05):", nrow(deg), "\n")
cat("上调:", sum(deg$logFC > 0), " 下调:", sum(deg$logFC < 0), "\n")

# ============================================================
# 5. QC图（草稿版）
# ============================================================
pca_df <- data.frame(
  PC1 = prcomp(t(expr))$x[, 1],
  PC2 = prcomp(t(expr))$x[, 2],
  group = ph$group)
p1 <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 2.5, alpha = .8) + theme_bw(base_size = 12) +
  stat_ellipse(level = .9, linetype = 2) +
  ggtitle("PCA GSE96804 (training)")
vol <- data.frame(logFC = tt$logFC, neglog10 = -log10(tt$adj.P.Val))
vol$col <- ifelse(vol$logFC >= .5 & vol$neglog10 > -log10(.05), "Up",
           ifelse(vol$logFC <= -.5 & vol$neglog10 > -log10(.05), "Down", "NS"))
p2 <- ggplot(vol, aes(logFC, neglog10, color = col)) +
  geom_point(size = .6, alpha = .5) + theme_bw(base_size = 12) +
  scale_color_manual(values = c(Up = "#D6604D", Down = "#2166AC", NS = "grey70")) +
  geom_hline(yintercept = -log10(.05), linetype = 3) +
  geom_vline(xintercept = c(-.5, .5), linetype = 3) +
  ggtitle(sprintf("Volcano (DEG n=%d)", nrow(deg)))
ggsave("results/figures/Fig2_qc_deg_draft.png",
       patchwork::wrap_plots(p1, p2), width = 10, height = 4.5, dpi = 200)

cat("\nPREPROCESS_DEG_DONE\n")

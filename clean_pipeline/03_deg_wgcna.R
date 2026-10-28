# 03_deg_wgcna.R ----
# Parse GEO series matrices, merge the GSE142025 RNA-seq cohort, map probes to
# gene symbols, then run limma DEGs and WGCNA on the training cohort (GSE96804).
# Outputs:
#   data/<acc>_expr.rds, data/<acc>_meta.txt          (series-matrix parse)
#   data/GSE142025_expr.rds, data/GSE142025_pheno.csv (RNA-seq merge)
#   data/<acc>_expr_symbol.rds, data/<acc>_pheno.csv  (probe-mapped matrices)
#   results/A_all_genes_limma.csv, results/A_DEGs.csv
#   results/wgcna_modules.csv, results/wgcna_module_trait.csv,
#   results/wgcna_sigmodule_genes.csv

suppressWarnings(suppressMessages({
  library(data.table); library(limma); library(org.Hs.eg.db); library(WGCNA)
}))
set.seed(42)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Series-matrix parse -> data/<acc>_expr.rds + <acc>_meta.txt
# ============================================================
parse_series_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  keep_meta <- grepl(
    "^!Series_title|^!Series_summary|^!Series_overall_design|^!Sample_title|^!Sample_geo_accession|^!Sample_source_name|^!Sample_characteristics_ch1|^!Sample_platform_id|^!Platform_title|^!Sample_supplementary",
    lines)
  beg <- grep("!series_matrix_table_begin", lines, fixed = TRUE)
  end <- grep("!series_matrix_table_end",   lines, fixed = TRUE)
  expr <- NULL
  if (length(beg) == 1 && length(end) >= 1) {
    # beg = begin marker, beg+1 = header, data rows end one line above the
    # end marker; a stub with no data rows yields a zero-row table
    expr <- read.delim(gzfile(path), skip = beg, nrows = end - beg - 2,
                       row.names = 1, check.names = FALSE,
                       na.strings = c("NA", "null", ""))
    if (nrow(expr) == 0 ||
        any(grepl("!series_matrix_table_end", rownames(expr), fixed = TRUE))) {
      expr <- NULL
    }
  }
  list(meta = lines[keep_meta], expr = expr)
}

# All series matrices are parsed (metadata-only stubs such as GSE142025 and
# GSE273763 still carry the phenotype lines needed downstream); the
# expression block is only saved when present.
sm_files <- list.files("data/GEO", pattern = "_series_matrix\\.txt\\.gz$",
                       full.names = TRUE)
for (f in sm_files) {
  acc <- sub("_series_matrix\\.txt\\.gz$", "", basename(f))
  parsed <- parse_series_matrix(f)
  writeLines(parsed$meta, paste0("data/", acc, "_meta.txt"))
  if (is.null(parsed$expr)) {
    cat(acc, ": no expression block in series matrix\n"); next
  }
  saveRDS(parsed$expr, paste0("data/", acc, "_expr.rds"))
  cat(acc, ":", nrow(parsed$expr), "features x", ncol(parsed$expr), "samples\n")
}

# ============================================================
# 2. GSE142025 RNA-seq: merge per-sample Illumina files
# ============================================================
if (!file.exists("data/GSE142025_expr.rds")) {
  raw_dir <- "data/GEO/GSE142025_raw"
  sfiles <- list.files(raw_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
  mat_list <- lapply(sfiles, function(f) {
    d <- fread(f, header = TRUE); setnames(d, c("Symbol", "expr")); d
  })
  symbols <- mat_list[[1]]$Symbol
  stopifnot(all(vapply(mat_list, function(d) identical(d$Symbol, symbols), logical(1))))
  expr14 <- as.matrix(do.call(cbind, lapply(mat_list, function(d) d$expr)))
  rownames(expr14) <- symbols
  colnames(expr14) <- tools::file_path_sans_ext(basename(sfiles))
  saveRDS(expr14, "data/GSE142025_expr.rds")
  clean_v <- function(x) trimws(gsub('^"|"$', "", x))
  strip_fields <- function(l) clean_v(strsplit(sub("^[^\t]*\t", "", l), "\t", fixed = TRUE)[[1]])
  ml <- readLines("data/GSE142025_meta.txt", warn = FALSE)
  ttl <- strip_fields(ml[grepl("^!Sample_title", ml)][1])
  grp <- sub("^group:\\s*", "",
             strip_fields(ml[grepl("^!Sample_characteristics_ch1", ml)][1]))
  ph14 <- data.frame(sample_id = ttl, group = grp, stringsAsFactors = FALSE)
  fwrite(ph14, "data/GSE142025_pheno.csv")
  cat("GSE142025 merged:", nrow(expr14), "genes x", ncol(expr14), "samples;",
      paste(names(table(grp)), table(grp), sep = "=", collapse = ", "), "\n")
}

# ============================================================
# 3. GPL probe annotation
# ============================================================
# Platform annotation: the SOFT family file from the GEO platform record
# (GPLxxx_family.soft.gz, used as downloaded). A plain-text copy named
# GPLxxx_full.txt is accepted as a fallback.
gpl_file <- function(gpl) {
  for (f in c(file.path("data/GEO", paste0(gpl, "_family.soft.gz")),
              file.path("data/GEO", paste0(gpl, "_full.txt")))) {
    if (file.exists(f)) return(f)
  }
  stop(gpl, ": platform annotation file not found (download ", gpl,
       "_family.soft.gz from the GEO platform record into data/GEO/)")
}

find_table_line <- function(f) {
  con <- if (grepl("\\.gz$", f)) gzfile(f, "rt") else file(f, "r")
  on.exit(close(con))
  n <- 0L
  repeat {
    l <- readLines(con, n = 1L, warn = FALSE)
    if (length(l) == 0L) stop("platform table not found: ", f)
    n <- n + 1L
    if (grepl("!platform_table_begin", l, fixed = TRUE)) return(n)
  }
}

read_gpl_table <- function(gpl) {
  f <- gpl_file(gpl)
  dt <- fread(f, skip = find_table_line(f), sep = "\t", header = TRUE, quote = "",
              fill = TRUE, showProgress = FALSE)
  dt <- dt[!grepl("!platform_table_end", as.character(dt[[1]]), fixed = TRUE)]
  if (ncol(dt) < 2) stop(gpl, ": platform table parse failed, only ", ncol(dt), " columns")
  dt
}

first_token <- function(x) trimws(sub("\\/\\/\\/.*$", "", x))
first_field <- function(x, k) {          # k-th "//" field of gene_assignment
  a <- sub("\\/\\/\\/.*$", "", x)        # first assignment only
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
  if ("Gene Symbol" %in% cn)            sym <- first_token(dt[["Gene Symbol"]])
  else if ("Symbol" %in% cn)            sym <- first_token(dt[["Symbol"]])
  else if ("gene_assignment" %in% cn)   sym <- first_field(dt[["gene_assignment"]], 2)
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

# GSE96804: parse the full platform table embedded in the family.soft file
# (the GPL17586 acc.cgi table is truncated; no Bioconductor db exists)
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
  if (is.na(begin_n) || is.na(end_n)) stop("platform table boundaries not found in family.soft")
  dt <- fread(f, skip = begin_n, nrows = end_n - begin_n - 1L, sep = "\t",
              header = TRUE, quote = "", fill = TRUE, showProgress = FALSE)
  sym <- first_field(as.character(dt$gene_assignment), 2)
  m <- data.frame(probe = sub("_st$", "", as.character(dt[[1]])),
                  symbol = sym, stringsAsFactors = FALSE)
  m <- m[!is.na(m$symbol) & m$symbol != "" & m$symbol != "---", ]
  m[!duplicated(m$probe), ]
}

map_list <- list(
  GSE96804             = build_map_familysoft(),
  GSE30528             = build_map("GPL571"),
  GSE30529             = build_map("GPL571"),
  `GSE104948-GPL22945` = build_map("GPL22945"),
  `GSE104948-GPL24120` = build_map("GPL24120")
)
for (nm in names(map_list)) cat(" ", nm, "probe map:", nrow(map_list[[nm]]), "probes\n")

# ============================================================
# 4. Probe -> symbol collapse (multi-probe genes keep max-mean probe)
# ============================================================
collapse_to_symbol <- function(expr, map) {
  rn <- rownames(expr)
  rn_key <- sub("_st$", "", rn)          # normalize HTA array probe suffix
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

for (acc in names(map_list)) {
  expr <- readRDS(paste0("data/", acc, "_expr.rds"))
  m <- collapse_to_symbol(expr, map_list[[acc]])
  saveRDS(m, paste0("data/", acc, "_expr_symbol.rds"))
  cat(" ", acc, ":", nrow(m), "genes x", ncol(m), "samples\n")
}

# ============================================================
# 5. Cohort phenotypes
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

# GSE96804: source_name prefix diab/unaf
ml <- readLines("data/GSE96804_meta.txt", warn = FALSE)
src <- split_line(ml[grepl("^!Sample_source_name_ch1", ml)][1])
pheno96804 <- data.frame(sample = colnames(readRDS("data/GSE96804_expr_symbol.rds")),
                         stringsAsFactors = FALSE)
pheno96804$group <- ifelse(grepl("^diab", src, ignore.case = TRUE), "DN", "Control")
write.csv(pheno96804, "data/GSE96804_pheno.csv", row.names = FALSE)
cat("GSE96804 groups:", paste(names(table(pheno96804$group)), table(pheno96804$group),
                              sep = "=", collapse = ", "), "\n")

# GSE30528 / GSE30529: title prefix Diab/Cont
for (acc in c("GSE30528", "GSE30529")) {
  ml <- readLines(paste0("data/", acc, "_meta.txt"), warn = FALSE)
  ttl <- split_line(ml[grepl("^!Sample_title", ml)][1])
  ph <- data.frame(sample = colnames(readRDS(paste0("data/", acc, "_expr_symbol.rds"))),
                   stringsAsFactors = FALSE)
  ph$group <- ifelse(grepl("^Diab", ttl, ignore.case = TRUE), "DKD", "Control")
  write.csv(ph, paste0("data/", acc, "_pheno.csv"), row.names = FALSE)
  cat(acc, "groups:", paste(names(table(ph$group)), table(ph$group),
                            sep = "=", collapse = ", "), "\n")
}

# GSE104948: diagnosis characteristic -> keep DN and tumor-nephrectomy controls,
# then pool both platforms on shared genes
build_104948 <- function(acc) {
  ml <- readLines(paste0("data/", acc, "_meta.txt"), warn = FALSE)
  ch <- ml[grepl("^!Sample_characteristics_ch1", ml)]
  vals <- lapply(ch, split_line)
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
  grp_all <- ifelse(diagv == "Diabetic Nephropathy", "DN",
             ifelse(diagv == "Tumor Nephrectomy", "Control", NA))
  # per-platform pheno covers every sample of the platform matrix; samples
  # outside the DN / tumor-nephrectomy contrast carry an empty group
  write.csv(data.frame(sample = colnames(expr), group = grp_all),
            paste0("data/", acc, "_pheno.csv"), row.names = FALSE)
  keep <- !is.na(grp_all)
  list(expr = expr[, keep, drop = FALSE], grp = grp_all[keep])
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
cat("GSE104948 pooled groups:", paste(names(table(ph104948$group)),
                                      table(ph104948$group), sep = "=", collapse = ", "), "\n")

# ============================================================
# 6. Training-cohort DEGs (GSE96804, DN vs Control, limma)
# ============================================================
expr <- readRDS("data/GSE96804_expr_symbol.rds")
ph <- read.csv("data/GSE96804_pheno.csv", stringsAsFactors = FALSE)
stopifnot(identical(ph$sample, colnames(expr)))
group <- factor(ph$group, levels = c("Control", "DN"))
design <- model.matrix(~ group)
fit <- eBayes(lmFit(expr, design))
tt <- topTable(fit, coef = "groupDN", number = Inf, adjust.method = "BH")
tt$gene <- rownames(tt)
deg <- tt[abs(tt$logFC) >= 0.5 & tt$adj.P.Val < 0.05, ]
deg <- deg[order(deg$adj.P.Val), ]
fwrite(as.data.frame(tt), "results/A_all_genes_limma.csv")
fwrite(as.data.frame(deg), "results/A_DEGs.csv")
cat("DEG stats: tested =", nrow(tt),
    "| significant (|logFC|>=0.5 & adj.P<0.05) =", nrow(deg),
    "| up =", sum(deg$logFC > 0), "| down =", sum(deg$logFC < 0), "\n")

# ============================================================
# 7. WGCNA (training cohort)
# ============================================================
mads <- apply(expr, 1, mad)
dat <- t(expr[order(-mads)[1:min(5000, nrow(expr))], ])
cat("WGCNA input:", nrow(dat), "samples x", ncol(dat), "genes\n")

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(dat, powerVector = powers, networkType = "signed",
                         verbose = 0, corFnc = "cor")
sp <- sft$powerEstimate
if (is.na(sp)) sp <- 6
cat("soft-threshold power =", sp, "\n")

net <- blockwiseModules(
  dat, power = sp, networkType = "signed", TOMType = "signed",
  minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  maxBlockSize = 6000, verbose = 0)
moduleColors <- labels2colors(net$colors)
cat("modules:", length(unique(moduleColors)), "\n")

trait <- as.numeric(ph$group == "DN")
MEs <- orderMEs(moduleEigengenes(dat, moduleColors)$eigengenes)
mod_trait <- as.data.frame(cor(MEs, trait, use = "p"))
names(mod_trait)[1] <- "corr_r"
mod_trait$p <- corPvalueStudent(mod_trait$corr_r, nrow(dat))
mod_trait$p_BH <- signif(p.adjust(mod_trait$p, method = "BH"), 3)
mod_trait$module <- sub("^ME", "", rownames(mod_trait))
fwrite(mod_trait, "results/wgcna_module_trait.csv")

kme <- signedKME(dat, MEs)
wgc_out <- data.frame(gene = colnames(dat), module = moduleColors)
for (m in colnames(kme)) wgc_out[[m]] <- kme[, m]
fwrite(wgc_out, "results/wgcna_modules.csv")

sig_mods <- mod_trait$module[abs(mod_trait[, 1]) > 0.3 & mod_trait$p < 0.05]
sig_mods <- setdiff(sig_mods, "grey")
sig_genes <- wgc_out$gene[wgc_out$module %in% sig_mods]
fwrite(data.frame(gene = sig_genes), "results/wgcna_sigmodule_genes.csv")
cat("disease-associated modules (|r|>0.3, p<0.05):", paste(sig_mods, collapse = ", "),
    "| genes:", length(sig_genes), "\n")

cat("\nDEG_WGCNA_DONE\n")

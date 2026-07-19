#!/usr/bin/env Rscript
# fetch_preprocess_breast.R
#
# Downloads and preprocesses three breast cancer GEO datasets:
#   GSE2034  Wang et al. 2005, Lancet       GPL96  Affymetrix HGU133A
#   GSE2990  Sotiriou et al. 2006, JNCI     GPL96  Affymetrix HGU133A
#   GSE9893  Chanrion et al. 2008, CCR      custom two-colour cDNA array
#
# Output: survival_NPDR/data/{accession}/dat.rds for each dataset
#   dat.rds is a data.frame with columns: time, status, <gene/probe features>
#   Features are filtered to the top 5,000 by variance across samples.
#   Samples with missing time or status are excluded.
#
# Dependencies: GEOquery, Biobase, dplyr
# Install if needed:
#   if (!requireNamespace("BiocManager")) install.packages("BiocManager")
#   BiocManager::install(c("GEOquery", "Biobase"))
#   install.packages("dplyr")

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
root      <- tryCatch(rprojroot::find_root(rprojroot::has_file("README.md")),
                      error = function(e) getwd())
data_root <- file.path(root, "data")
N_TOP     <- 5000L   # variance filter threshold (shared across all three datasets)

# ---------------------------------------------------------------------------
# Helper: variance-filter + save
# ---------------------------------------------------------------------------
save_dat <- function(dat, out_path) {
  outcome_cols <- c("time", "status")
  attr_mat     <- dat[, setdiff(names(dat), outcome_cols)]
  attr_mat     <- attr_mat[, sapply(attr_mat, is.numeric), drop = FALSE]

  var_vals  <- apply(attr_mat, 2, var, na.rm = TRUE)
  top_genes <- names(sort(var_vals, decreasing = TRUE)[seq_len(min(N_TOP, ncol(attr_mat)))])
  attr_mat  <- attr_mat[, top_genes, drop = FALSE]

  out <- cbind(time = dat$time, status = dat$status, attr_mat)
  saveRDS(out, out_path)
  invisible(out)
}

# ---------------------------------------------------------------------------
# Helper: deduplicate probes → keep probe with highest mean expression per gene
# ---------------------------------------------------------------------------
dedup_by_mean <- function(expr_df, gene_col = "gene") {
  sample_cols      <- setdiff(names(expr_df), c("ID", gene_col))
  expr_df$row_mean <- rowMeans(expr_df[, sample_cols, drop = FALSE], na.rm = TRUE)
  expr_df          <- expr_df[order(expr_df[[gene_col]], -expr_df$row_mean), ]
  expr_df          <- expr_df[!duplicated(expr_df[[gene_col]]), ]
  expr_df$row_mean <- NULL
  rownames(expr_df) <- expr_df[[gene_col]]
  expr_df[, sample_cols, drop = FALSE]
}

# ===========================================================================
# GSE2034 — node-negative breast cancer, time to relapse (months)
# ===========================================================================
message("\n=== GSE2034 ===")
out2034 <- file.path(data_root, "GSE2034", "dat.rds")
dir.create(dirname(out2034), recursive = TRUE, showWarnings = FALSE)

gset2034 <- getGEO("GSE2034", GSEMatrix = TRUE, getGPL = TRUE)[[1]]

# Expression matrix
expr2034  <- exprs(gset2034)
fdata2034 <- fData(gset2034)

# Map probes → gene symbols (GPL96 'Gene Symbol' column)
sym_col <- grep("gene.symbol|symbol", names(fdata2034), ignore.case = TRUE, value = TRUE)[1]
if (is.na(sym_col)) stop("GSE2034: cannot find Gene Symbol column in fData")

probe_gene <- fdata2034[, sym_col, drop = FALSE]
colnames(probe_gene) <- "gene"
probe_gene$gene <- trimws(as.character(probe_gene$gene))

# Remove probes without a gene symbol
probe_gene <- probe_gene[nchar(probe_gene$gene) > 0 & !is.na(probe_gene$gene), , drop = FALSE]
expr2034   <- expr2034[rownames(probe_gene), , drop = FALSE]

expr_df2034       <- as.data.frame(t(expr2034))   # samples × probes
expr_df2034_t     <- as.data.frame(t(expr2034))   # keep samples as rows for dedup
expr_probe_df     <- as.data.frame(expr2034)       # probes × samples
expr_probe_df$ID  <- rownames(expr_probe_df)
expr_probe_df$gene <- probe_gene[rownames(expr_probe_df), "gene"]
expr_by_gene2034  <- dedup_by_mean(expr_probe_df)  # genes × samples

# Phenotype: time (months), relapse status
# Prefer local phenotype_data.txt if present; fall back to GEO pData
pheno_file2034 <- file.path(data_root, "GSE2034", "phenotype_data.txt")
if (file.exists(pheno_file2034)) {
  pheno2034 <- read.delim(pheno_file2034, comment.char = "#",
                          header = TRUE, check.names = FALSE,
                          stringsAsFactors = FALSE)
  # Column names from phenotype_data.txt:
  #   PID, GEO accession number, lymph node status,
  #   time to relapse or last follow-up (months), relapse (1=True), ...
  time_col   <- grep("time", names(pheno2034), ignore.case = TRUE, value = TRUE)[1]
  status_col <- grep("relapse.*1|^relapse", names(pheno2034), ignore.case = TRUE, value = TRUE)[1]
  geo_col    <- grep("geo|accession|asscession", names(pheno2034), ignore.case = TRUE, value = TRUE)[1]
  pheno2034$time   <- suppressWarnings(as.numeric(pheno2034[[time_col]]))
  pheno2034$status <- suppressWarnings(as.numeric(pheno2034[[status_col]]))
  rownames(pheno2034) <- trimws(as.character(pheno2034[[geo_col]]))
} else {
  pd2034     <- pData(gset2034)
  time_col   <- grep("time", names(pd2034), ignore.case = TRUE, value = TRUE)[1]
  status_col <- grep("relapse|event|status", names(pd2034), ignore.case = TRUE, value = TRUE)[1]
  pheno2034  <- data.frame(
    time   = suppressWarnings(as.numeric(pd2034[[time_col]])),
    status = suppressWarnings(as.numeric(pd2034[[status_col]])),
    row.names = rownames(pd2034)
  )
}

# Align samples
common2034 <- intersect(colnames(expr_by_gene2034), rownames(pheno2034))
if (length(common2034) == 0) stop("GSE2034: no overlapping sample IDs between expression and phenotype")

expr_aligned2034  <- as.data.frame(t(expr_by_gene2034[, common2034, drop = FALSE]))
pheno_aligned2034 <- pheno2034[common2034, c("time", "status"), drop = FALSE]
dat2034 <- cbind(pheno_aligned2034, expr_aligned2034)
dat2034 <- dat2034[complete.cases(dat2034[, c("time", "status")]), ]

message(sprintf("  Samples after QC: %d  |  Events: %d  |  Genes before var-filter: %d",
                nrow(dat2034), sum(dat2034$status, na.rm = TRUE),
                ncol(dat2034) - 2))

saved2034 <- save_dat(dat2034, out2034)
message(sprintf("  Saved: %s  [%d samples × %d features]", out2034,
                nrow(saved2034), ncol(saved2034) - 2))

# ===========================================================================
# GSE2990 — tamoxifen-treated breast cancer, relapse-free survival (years)
# ===========================================================================
message("\n=== GSE2990 ===")
out2990 <- file.path(data_root, "GSE2990", "dat.rds")
dir.create(dirname(out2990), recursive = TRUE, showWarnings = FALSE)

gset2990 <- getGEO("GSE2990", GSEMatrix = TRUE, getGPL = TRUE)[[1]]

# Expression — same GPL96 platform, same probe→gene mapping approach
expr2990   <- exprs(gset2990)
fdata2990  <- fData(gset2990)
sym_col2   <- grep("gene.symbol|symbol", names(fdata2990), ignore.case = TRUE, value = TRUE)[1]
if (is.na(sym_col2)) stop("GSE2990: cannot find Gene Symbol column in fData")

probe_gene2 <- fdata2990[, sym_col2, drop = FALSE]
colnames(probe_gene2) <- "gene"
probe_gene2$gene <- trimws(as.character(probe_gene2$gene))
probe_gene2 <- probe_gene2[nchar(probe_gene2$gene) > 0 & !is.na(probe_gene2$gene), , drop = FALSE]
expr2990    <- expr2990[rownames(probe_gene2), , drop = FALSE]

expr_probe_df2 <- as.data.frame(expr2990)
expr_probe_df2$ID   <- rownames(expr_probe_df2)
expr_probe_df2$gene <- probe_gene2[rownames(expr_probe_df2), "gene"]
expr_by_gene2990 <- dedup_by_mean(expr_probe_df2)

# Phenotype: prefer local GSE2990_suppl_info.txt
# Columns: geo_accn, sample_name, treatment, dataset, grade, node, size, age, time, status
suppl2990 <- file.path(data_root, "GSE2990", "GSE2990_suppl_info.txt")
if (file.exists(suppl2990)) {
  pheno2990_raw <- tryCatch(
    read.delim(suppl2990, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (!is.null(pheno2990_raw)) {
    time_col2   <- grep("^time|follow", names(pheno2990_raw), ignore.case = TRUE, value = TRUE)[1]
    status_col2 <- grep("status|relapse|event", names(pheno2990_raw), ignore.case = TRUE, value = TRUE)[1]
    geo_col2    <- grep("geo_accn|geo.acc|gsm", names(pheno2990_raw), ignore.case = TRUE, value = TRUE)[1]
    pheno2990 <- data.frame(
      time   = suppressWarnings(as.numeric(pheno2990_raw[[time_col2]])),
      status = suppressWarnings(as.numeric(pheno2990_raw[[status_col2]])),
      row.names = trimws(as.character(pheno2990_raw[[geo_col2]])),
      stringsAsFactors = FALSE
    )
  } else {
    pheno2990 <- NULL
  }
}

if (is.null(pheno2990) || nrow(pheno2990) == 0) {
  # Fall back to GEO pData
  pd2990     <- pData(gset2990)
  time_col2  <- grep("time|follow|survival", names(pd2990), ignore.case = TRUE, value = TRUE)[1]
  stat_col2  <- grep("relapse|event|status", names(pd2990), ignore.case = TRUE, value = TRUE)[1]
  pheno2990  <- data.frame(
    time   = suppressWarnings(as.numeric(pd2990[[time_col2]])),
    status = suppressWarnings(as.numeric(pd2990[[stat_col2]])),
    row.names = rownames(pd2990),
    stringsAsFactors = FALSE
  )
}

common2990 <- intersect(colnames(expr_by_gene2990), rownames(pheno2990))
if (length(common2990) == 0) stop("GSE2990: no overlapping sample IDs between expression and phenotype")

expr_aligned2990  <- as.data.frame(t(expr_by_gene2990[, common2990, drop = FALSE]))
pheno_aligned2990 <- pheno2990[common2990, c("time", "status"), drop = FALSE]
dat2990 <- cbind(pheno_aligned2990, expr_aligned2990)
dat2990 <- dat2990[complete.cases(dat2990[, c("time", "status")]), ]

message(sprintf("  Samples after QC: %d  |  Events: %d  |  Genes before var-filter: %d",
                nrow(dat2990), sum(dat2990$status, na.rm = TRUE),
                ncol(dat2990) - 2))

saved2990 <- save_dat(dat2990, out2990)
message(sprintf("  Saved: %s  [%d samples × %d features]", out2990,
                nrow(saved2990), ncol(saved2990) - 2))

# ===========================================================================
# GSE9893 — ER+ tamoxifen-treated breast cancer, custom two-colour cDNA array
#           Chanrion et al. 2008 (Clin Cancer Res)
#
# Expression: series matrix from GEO (pre-normalised log-ratios).
# Clinical:   GSE9893_clinicalData.txt (local); event = "deceased".
# Features:   cDNA clone IDs (no HGNC gene-symbol mapping available on this platform).
# ===========================================================================
message("\n=== GSE9893 ===")
out9893 <- file.path(data_root, "GSE9893", "dat.rds")
dir.create(dirname(out9893), recursive = TRUE, showWarnings = FALSE)

gset9893 <- getGEO("GSE9893", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
expr9893 <- exprs(gset9893)   # features × samples (clone IDs as rownames)

# Clinical data from local file
clin_file9893 <- file.path(data_root, "GSE9893_clinicalData.txt")
if (!file.exists(clin_file9893)) {
  clin_file9893 <- file.path(data_root, "GSE9893", "GSE9893_clinicalData.txt")
}
if (!file.exists(clin_file9893)) stop("GSE9893: cannot find GSE9893_clinicalData.txt")

clin9893 <- read.delim(clin_file9893, sep = "\t", header = TRUE,
                        check.names = FALSE, stringsAsFactors = FALSE)

# Expected columns: "Follow-up period (months)", "State of health"
# (from gse9893_survival_curve.R)
time_col9   <- "Follow-up period (months)"
status_col9 <- "State of health"

if (!all(c(time_col9, status_col9) %in% names(clin9893))) {
  # Try to auto-detect
  time_col9   <- grep("follow|time", names(clin9893), ignore.case = TRUE, value = TRUE)[1]
  status_col9 <- grep("health|status|vital|state", names(clin9893), ignore.case = TRUE, value = TRUE)[1]
}

time9893   <- suppressWarnings(as.numeric(trimws(clin9893[[time_col9]])))
status_raw <- tolower(trimws(clin9893[[status_col9]]))
status9893 <- ifelse(grepl("^deceased", status_raw), 1,
               ifelse(status_raw == "" | is.na(status_raw), NA, 0))

# GEO GSM IDs: typically the first column or rownames of the clinical file
# GSE9893 clinical file uses tumor sample IDs (EB5012, etc.) — match via pData
pd9893     <- pData(gset9893)
gsm_ids    <- rownames(pd9893)

# The clinical file row order matches the GEO sample order for GSE9893
# (confirmed by n = 154 in both)
if (nrow(clin9893) != ncol(expr9893)) {
  warning(sprintf("GSE9893: clinical rows (%d) != expression columns (%d); attempting order match",
                  nrow(clin9893), ncol(expr9893)))
}
n9893 <- min(nrow(clin9893), ncol(expr9893))

pheno9893 <- data.frame(
  time   = time9893[seq_len(n9893)],
  status = status9893[seq_len(n9893)],
  row.names = gsm_ids[seq_len(n9893)],
  stringsAsFactors = FALSE
)

expr9893_sub <- expr9893[, seq_len(n9893), drop = FALSE]
colnames(expr9893_sub) <- gsm_ids[seq_len(n9893)]

expr_aligned9893  <- as.data.frame(t(expr9893_sub))
dat9893 <- cbind(pheno9893, expr_aligned9893)
dat9893 <- dat9893[complete.cases(dat9893[, c("time", "status")]), ]

message(sprintf("  Samples after QC: %d  |  Events: %d  |  Features before var-filter: %d",
                nrow(dat9893), sum(dat9893$status, na.rm = TRUE),
                ncol(dat9893) - 2))
message("  NOTE: features are cDNA clone IDs (e.g. 101F6, 13CDNA73), not HGNC gene symbols")

saved9893 <- save_dat(dat9893, out9893)
message(sprintf("  Saved: %s  [%d samples × %d features]", out9893,
                nrow(saved9893), ncol(saved9893) - 2))

# ===========================================================================
# Summary
# ===========================================================================
message("\n=== Preprocessing complete ===")
message(sprintf("  GSE2034: %d samples, %d features, time unit = months (relapse)",
                nrow(saved2034), ncol(saved2034) - 2))
message(sprintf("  GSE2990: %d samples, %d features, time unit = years (relapse-free survival)",
                nrow(saved2990), ncol(saved2990) - 2))
message(sprintf("  GSE9893: %d samples, %d features, time unit = months (overall survival)",
                nrow(saved9893), ncol(saved9893) - 2))

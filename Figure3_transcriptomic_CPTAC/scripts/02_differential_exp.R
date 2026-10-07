#!/usr/bin/env Rscript
# ======================================================================
# 02_differential_exp.R
# Differential expression (S1 vs S2) using official edgeR + limma
# ----------------------------------------------------------------------
# Pipeline (Law et al., 2014; Ritchie et al., 2015; Robinson & Oshlack, 2010):
#   1. edgeR::filterByExpr  - automatic low-count filtering by group
#   2. edgeR::calcNormFactors (TMM)
#   3. limma::voom          - log2-CPM + LOWESS precision weights
#   4. limma::lmFit         - weighted least squares
#   5. limma::eBayes        - empirical Bayes moderation (d0 ESTIMATED, not fixed)
#   6. limma::topTable      - S1 vs S2 statistics
#
# Output: data/processed/deg_limma_voom.csv
#   columns: gene, logFC, AveExpr, t, p_value, FDR, category
# ======================================================================

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- sub("^--file=", "", file_arg)
if (length(script_path) == 0 || script_path == "") {
  script_path <- "scripts/upstream/R/02_differential_exp.R"
}
REPO <- dirname(dirname(normalizePath(script_path)))
DATA_RAW  <- file.path(REPO, "data", "raw")
DATA_PROC <- file.path(REPO, "data", "processed")
dir.create(DATA_PROC, showWarnings = FALSE, recursive = TRUE)

counts_path <- file.path(DATA_RAW, "KIRC_counts_for_limma.csv")
sample_path <- file.path(DATA_RAW, "sample_info_for_limma.csv")
stopifnot(file.exists(counts_path), file.exists(sample_path))

library(edgeR)
library(limma)

cat("==== edgeR + limma-voom: S1 vs S2 DEG ====\n")

counts <- read.csv(counts_path, row.names = 1)
sample_info <- read.csv(sample_path)
subtype <- as.character(sample_info$Subtype)
stopifnot(length(subtype) == ncol(counts))
cat(sprintf("  %d genes x %d samples (S1=%d, S2=%d)\n",
            nrow(counts), ncol(counts),
            sum(subtype == "Subtype_1"), sum(subtype == "Subtype_2")))

group <- factor(subtype, levels = c("Subtype_2", "Subtype_1"))  # S2 = reference
y <- DGEList(counts = counts, group = group)
keep <- filterByExpr(y)
cat(sprintf("  filterByExpr: %d genes retained (of %d)\n", sum(keep), nrow(y)))
y <- y[keep, , keep.lib.sizes = FALSE]

y <- calcNormFactors(y, method = "TMM")
cat(sprintf("  TMM factors range: [%.3f, %.3f]\n",
            min(y$samples$norm.factor), max(y$samples$norm.factor)))

design <- model.matrix(~ group)
v <- voom(y, design = design, plot = FALSE)

fit <- lmFit(v, design)
fit <- eBayes(fit, trend = FALSE, robust = FALSE)

tt <- topTable(fit, coef = "groupSubtype_1", n = Inf, sort.by = "none")

res <- data.frame(
  gene    = rownames(tt),
  logFC   = tt$logFC,
  AveExpr = tt$AveExpr,
  t       = tt$t,
  p_value = tt$P.Value,
  FDR     = tt$adj.P.Val,
  stringsAsFactors = FALSE
)
res$category <- "NS"
res$category[res$FDR < 0.05 & res$logFC >  0.5] <- "UP_in_S1"
res$category[res$FDR < 0.05 & res$logFC < -0.5] <- "DOWN_in_S1"

n_up   <- sum(res$category == "UP_in_S1")
n_down <- sum(res$category == "DOWN_in_S1")
cat(sprintf("  UP_in_S1: %d | DOWN_in_S1: %d | NS: %d\n",
            n_up, n_down, nrow(res) - n_up - n_down))

out_csv <- file.path(DATA_PROC, "deg_limma_voom.csv")
write.csv(res, out_csv, row.names = FALSE)
cat("  Written:", out_csv, "\n")
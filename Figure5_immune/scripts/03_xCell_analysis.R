#!/usr/bin/env Rscript
# 03_xCell_analysis.R â€?True xCell analysis with spillover correction
# Input: data/tcga_kirc/KIRC_expr_tpm.csv, data/processed/subtype_assignment_balanced.csv
# Output: data/xcell/xcell_scores_raw.csv, xcell_scores_zscore.csv, xcell_subtype_diff.csv
#
# Repo-relative paths: detects this script's location and reads/writes under
# <repo>/data/... (same locations the Python figure scripts use).
# Honours MEL_DATA_ROOT env var when set (matches the Python scripts).
# Provenance: originally run from C:/Users/Administrator/WorkBuddy/2026-06-04-09-29-28
# (subtype file was read from a `consensus_cluster/` subdir; now uses data/processed/).

library(xCell)

# ---- Paths (auto-detect; honours MEL_DATA_ROOT like the Python scripts) ----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  SCRIPT_DIR <- getwd()
}
PROJ_ROOT  <- dirname(SCRIPT_DIR)
DATA_ROOT  <- Sys.getenv("MEL_DATA_ROOT", PROJ_ROOT)
XCELL_DIR  <- file.path(DATA_ROOT, "data", "xcell")
OUT_DIR    <- XCELL_DIR
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load linear TPM
cat("Loading TPM data...\n")
expr <- read.csv(file.path(DATA_ROOT, "data", "tcga_kirc", "KIRC_expr_tpm.csv"),
                 row.names = 1, check.names = FALSE)
expr_mat <- as.matrix(expr)
cat(sprintf("  %d genes x %d samples\n", nrow(expr_mat), ncol(expr_mat)))

# Run xCell (with spillover correction if available)
cat("Running xCell analysis...\n")
spill_rds <- file.path(OUT_DIR, "xcell_spill.rds")
if (file.exists(spill_rds)) {
  my_spill <- readRDS(spill_rds)
  cat(sprintf("  Spillover loaded: K %dx%d, fv %dx%d\n",
              nrow(my_spill$K), ncol(my_spill$K),
              nrow(my_spill$fv), ncol(my_spill$fv)))
  xcell_raw <- xCellAnalysis(expr_mat, spill = my_spill, rnaseq = TRUE)
} else {
  cat("  No spillover data, running without correction\n")
  xcell_raw <- xCellAnalysis(expr_mat, spill = FALSE, rnaseq = TRUE)
}
cat(sprintf("  xCell complete: %d cell types x %d samples\n",
            nrow(xcell_raw), ncol(xcell_raw)))

# Save raw scores
write.csv(xcell_raw, file.path(OUT_DIR, "xcell_scores_raw.csv"))
cat("  Saved: xcell_scores_raw.csv\n")

# Z-score normalize
cat("Z-score normalizing...\n")
xcell_z <- t(scale(t(xcell_raw)))
write.csv(xcell_z, file.path(OUT_DIR, "xcell_scores_zscore.csv"))
cat("  Saved: xcell_scores_zscore.csv\n")

# S1 vs S2 stats
cat("Computing S1 vs S2 stats...\n")
sub <- read.csv(file.path(DATA_ROOT, "data", "processed",
                          "subtype_assignment_balanced.csv"))
sub_map <- setNames(sub$Subtype, sub$Patient)
s1 <- names(sub_map)[sub_map == "Subtype_1"]
s2 <- names(sub_map)[sub_map == "Subtype_2"]
common <- intersect(colnames(xcell_z), names(sub_map))
s1_idx <- intersect(common, s1)
s2_idx <- intersect(common, s2)

results <- data.frame()
for (ct in rownames(xcell_z)) {
  s1_vals <- as.numeric(xcell_z[ct, s1_idx])
  s2_vals <- as.numeric(xcell_z[ct, s2_idx])
  if (length(s1_vals) == 0 || length(s2_vals) == 0) next
  test <- wilcox.test(s1_vals, s2_vals)
  results <- rbind(results, data.frame(
    cell_type = ct,
    S1_mean   = mean(s1_vals),
    S2_mean   = mean(s2_vals),
    mean_zscore_diff = mean(s2_vals) - mean(s1_vals),
    pvalue    = test$p.value,
    stringsAsFactors = FALSE
  ))
}

results$fdr <- p.adjust(results$pvalue, method = "BH")
results <- results[order(results$fdr), ]
write.csv(results, file.path(OUT_DIR, "xcell_subtype_diff.csv"), row.names = FALSE)

cat(sprintf("  Significant (FDR<0.05): %d/%d\n",
            sum(results$fdr < 0.05), nrow(results)))
cat("\nDone.\n")

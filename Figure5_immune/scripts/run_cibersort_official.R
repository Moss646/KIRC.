#!/usr/bin/env Rscript
# run_cibersort_official.R â€?Official CIBERSORT R package analysis
# Input: KIRC_expr_tpm.csv (raw TPM, non-log, no negative, no missing)
# Output: cibersort_results_official_full.csv (22 cell types + P-value + Correlation + RMSE)
# Method: CIBERSORT with LM22 signature matrix (547 genes, nu-SVR)
#
# Repo-relative paths: detects this script's location and writes the CSV into
# <repo>/data/processed/ (same place the Python figure scripts read from).
# Honours MEL_DATA_ROOT env var when set (matches the Python scripts).
# Provenance: originally run from C:/Users/Administrator/WorkBuddy/2026-06-04-09-29-28
# (outputs went to a `precluster/` subdir; that step is now consolidated here).

library(CIBERSORT)
cat("CIBERSORT R package loaded (v0.1.0)\n")

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
OUT_DIR    <- file.path(DATA_ROOT, "data", "processed")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Load raw TPM expression matrix ----
cat("Loading raw TPM expression matrix...\n")
expr <- read.csv(file.path(DATA_ROOT, "data", "tcga_kirc", "KIRC_expr_tpm.csv"),
                 row.names = 1, check.names = FALSE)
expr_mat <- as.matrix(expr)
cat(sprintf("  %d genes x %d samples\n", nrow(expr_mat), ncol(expr_mat)))
cat(sprintf("  Min=%.2f, Max=%.2f\n", min(expr_mat), max(expr_mat)))
cat(sprintf("  Any negative: %s, Any NA: %s\n",
            any(expr_mat < 0), any(is.na(expr_mat))))

# ---- Load LM22 signature matrix ----
cat("Loading LM22 signature matrix...\n")
data(LM22)
cat(sprintf("  LM22: %d genes x %d cell types\n", nrow(LM22), ncol(LM22)))

# ---- Run CIBERSORT ----
cat("Running CIBERSORT (perm=1000, QN=TRUE)...\n")
cat("  This may take 10-30 minutes for 533 samples...\n")

results <- cibersort(sig_matrix = LM22,
                     mixture_file = expr_mat,
                     perm = 1000,
                     QN = TRUE)

cat(sprintf("  CIBERSORT complete: %d samples x %d columns\n",
            nrow(results), ncol(results)))

# ---- Check P-value filtering ----
p_col <- which(colnames(results) == "P-value")
if (length(p_col) > 0) {
  n_sig <- sum(results[, p_col] < 0.05)
  n_total <- nrow(results)
  cat(sprintf("  P < 0.05: %d/%d samples (%.1f%%)\n",
              n_sig, n_total, n_sig/n_total*100))

  # Filter by P < 0.05 (standard CIBERSORT practice)
  results_filtered <- results[results[, p_col] < 0.05, ]
  cat(sprintf("  After P-value filtering: %d samples\n", nrow(results_filtered)))
} else {
  cat("  No P-value column found\n")
  results_filtered <- results
}

# ---- Save results ----
# Full results (all samples, with quality metrics)
out_full <- file.path(OUT_DIR, "cibersort_results_official_full.csv")
write.csv(results, out_full)
cat(sprintf("  Saved full results: %s\n", out_full))

# Filtered results (P < 0.05 only, with quality metrics)
out_filtered <- file.path(OUT_DIR, "cibersort_results_official_filtered.csv")
write.csv(results_filtered, out_filtered)
cat(sprintf("  Saved filtered results: %s\n", out_filtered))

# Cell proportions only (for compatibility with fig4_immune_merged script)
# Extract just the 22 cell type columns, no P-value/Correlation/RMSE
cell_cols <- colnames(results)[1:22]
results_cells <- results[, cell_cols]
out_cells <- file.path(OUT_DIR, "cibersort_results.csv")
write.csv(results_cells, out_cells)
cat(sprintf("  Saved cell proportions (all %d samples): %s\n",
            nrow(results_cells), out_cells))

cat("\nDone!\n")

#!/usr/bin/env Rscript
# ============================================================
# run_ssgsea_immune.R
# ssGSEA for Hallmark immune pathways using the R GSVA package.
#
# Paths are auto-detected from the script location. Override the data
# root by setting the MEL_DATA_ROOT environment variable (same
# convention used by the Python scripts).
#
# Input:
#   - Expression: <DATA_ROOT>/data/tcga_kirc/KIRC_expr_log2_tpm.csv
#       (log2(TPM+1) transformed, genes x samples)
#   - Gene sets:  <DATA_ROOT>/data/gene_sets/h.all.v2024.1.Hs.symbols.gmt
#       (MSigDB Hallmark collection, 50 gene sets)
#   - Subtypes:   <DATA_ROOT>/data/processed/subtype_assignment_balanced.csv
#
# Method:
#   - ssGSEA implemented in GSVA R package v2.4.9 (R v4.5.0)
#   - Parameters: ssgseaParam with default settings
#
# Output (written to <DATA_ROOT>/data/processed/):
#   - ssgsea_immune_scores_R.csv (per-sample NES, pathways x samples)
#   - ssgsea_immune_stats_R.csv  (S1 vs S2 stats)
# ============================================================

suppressPackageStartupMessages({
    library(GSVA)
    library(limma)
})

cat("R version:", R.version.string, "\n")
cat("GSVA package version:", as.character(packageVersion("GSVA")), "\n")

# ---- Paths (auto-detect from this script; MEL_DATA_ROOT override) ----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
    SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
    SCRIPT_DIR <- getwd()
}
PROJ_ROOT <- dirname(SCRIPT_DIR)
DATA_ROOT <- Sys.getenv("MEL_DATA_ROOT", unset = PROJ_ROOT)

EXPR_PATH <- file.path(DATA_ROOT, "data", "tcga_kirc", "KIRC_expr_log2_tpm.csv")
GMT_PATH  <- file.path(DATA_ROOT, "data", "gene_sets", "h.all.v2024.1.Hs.symbols.gmt")
SUB_PATH  <- file.path(DATA_ROOT, "data", "processed", "subtype_assignment_balanced.csv")
OUT_DIR   <- file.path(DATA_ROOT, "data", "processed")

# ---- Load expression matrix ----
cat("\n[1/5] Loading expression matrix...\n")
expr <- read.csv(EXPR_PATH, row.names = 1, check.names = FALSE)
expr_mat <- as.matrix(expr)
cat(sprintf("  Genes: %d, Samples: %d\n", nrow(expr_mat), ncol(expr_mat)))
cat(sprintf("  Expression range: [%f, %f]\n", min(expr_mat), max(expr_mat)))

# ---- Load subtype assignments ----
cat("\n[2/5] Loading subtype assignments...\n")
sub <- read.csv(SUB_PATH)
sub_map <- setNames(sub$Subtype, sub$Patient)
s1_samples <- names(sub_map[sub_map == "Subtype_1"])
s2_samples <- names(sub_map[sub_map == "Subtype_2"])
cat(sprintf("  S1: %d samples, S2: %d samples\n",
            length(s1_samples), length(s2_samples)))

# ---- Load Hallmark gene sets from GMT ----
cat("\n[3/5] Loading Hallmark gene sets from GMT...\n")
gene_sets <- list()
con <- file(GMT_PATH, "r")
n_total <- 0
while (length(line <- readLines(con, n = 1)) > 0) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
        gs_name <- parts[1]
        gs_genes <- parts[3:length(parts)]
        gs_genes <- gs_genes[gs_genes != ""]
        gene_sets[[gs_name]] <- gs_genes
        n_total <- n_total + 1
    }
}
close(con)
cat(sprintf("  Total gene sets in GMT: %d\n", n_total))

genes_available <- rownames(expr_mat)
gene_sets_filtered <- lapply(gene_sets, function(gs) {
    intersect(gs, genes_available)
})
gene_sets_filtered <- gene_sets_filtered[sapply(gene_sets_filtered, length) >= 1]
cat(sprintf("  Gene sets after filtering (>=1 gene in expr): %d\n",
            length(gene_sets_filtered)))

# ---- Run ssGSEA ----
cat("\n[4/5] Running ssGSEA (GSVA R package v2.4.9)...\n")
cat("  This may take a few minutes...\n")

param <- ssgseaParam(
    exprData = expr_mat,
    geneSets = gene_sets_filtered,
    normalize = TRUE
)

ssgsea_scores <- gsva(param)

cat(sprintf("  Done! Output: %d pathways x %d samples\n",
            nrow(ssgsea_scores), ncol(ssgsea_scores)))
cat(sprintf("  Score range: [%f, %f]\n",
            min(ssgsea_scores), max(ssgsea_scores)))

# ---- S1 vs S2 comparison (Wilcoxon + BH) ----
cat("\n[5/5] Computing S1 vs S2 differences (Wilcoxon + BH)...\n")

s1_cols <- intersect(s1_samples, colnames(ssgsea_scores))
s2_cols <- intersect(s2_samples, colnames(ssgsea_scores))
cat(sprintf("  S1 samples in output: %d, S2 samples in output: %d\n",
            length(s1_cols), length(s2_cols)))

stats_rows <- list()
for (pwy in rownames(ssgsea_scores)) {
    label <- gsub("HALLMARK_", "", pwy)
    label <- gsub("_", " ", label)
    label <- tools::toTitleCase(label)

    v1 <- ssgsea_scores[pwy, s1_cols]
    v2 <- ssgsea_scores[pwy, s2_cols]

    s1_mean <- mean(v1)
    s2_mean <- mean(v2)
    diff_val <- s1_mean - s2_mean

    wt <- wilcox.test(v1, v2)
    p_raw <- wt$p.value

    stats_rows[[length(stats_rows) + 1]] <- data.frame(
        pathway = pwy,
        label = label,
        S1_mean = s1_mean,
        S2_mean = s2_mean,
        diff = diff_val,
        p_raw = p_raw
    )
}
stats_df <- do.call(rbind, stats_rows)
stats_df$p_bh <- p.adjust(stats_df$p_raw, method = "BH")

cat("\nS1 vs S2 pathway differences:\n")
for (i in 1:nrow(stats_df)) {
    r <- stats_df[i, ]
    sig <- if (r$p_bh < 0.001) "***" else if (r$p_bh < 0.01) "**" else
           if (r$p_bh < 0.05) "*" else ""
    cat(sprintf("  %s  diff=%+.4f  p_raw=%.2e  p_bh=%.2e  %s\n",
                r$label, r$diff, r$p_raw, r$p_bh, sig))
}

# ---- Save results ----
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

scores_df <- as.data.frame(ssgsea_scores)
scores_df$pathway <- rownames(scores_df)
scores_df <- scores_df[, c("pathway", setdiff(colnames(scores_df), "pathway"))]
write.csv(scores_df, file.path(OUT_DIR, "ssgsea_immune_scores_R.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved scores: %s\n",
            file.path(OUT_DIR, "ssgsea_immune_scores_R.csv")))

write.csv(stats_df, file.path(OUT_DIR, "ssgsea_immune_stats_R.csv"),
          row.names = FALSE)
cat(sprintf("Saved stats: %s\n",
            file.path(OUT_DIR, "ssgsea_immune_stats_R.csv")))

# ---- Summary ----
cat("\n=== Summary ===\n")
cat(sprintf("Method: ssGSEA via GSVA R package v%s\n",
            as.character(packageVersion("GSVA"))))
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Expression: log2(TPM+1), %d genes x %d samples\n",
            nrow(expr_mat), ncol(expr_mat)))
cat(sprintf("Gene sets: %d (from h.all.v2024.1.Hs.symbols.gmt)\n",
            length(gene_sets_filtered)))
cat(sprintf("Output: %d pathways x %d samples\n",
            nrow(ssgsea_scores), ncol(ssgsea_scores)))
cat(sprintf("Score range: [%f, %f]\n",
            min(ssgsea_scores), max(ssgsea_scores)))
cat(sprintf("BH correction: across all %d pathways\n", nrow(stats_df)))
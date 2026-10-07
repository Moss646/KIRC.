# ============================================================
# gsva_R.R
# GSVA pathway activity scores using the R GSVA package.
# Input:
#   - data/raw/KIRC_expr_log2_tpm.csv (log2(TPM+1), genes x samples)
#   - data/gene_sets/msigdb_lipid_dedup_236.gmt
# Output:
#   - data/processed/gsva_scores.csv (pathways x samples)
# ============================================================

suppressPackageStartupMessages({
    library(GSVA)
})

cat("R version:", R.version.string, "\n")
cat("GSVA package version:", as.character(packageVersion("GSVA")), "\n")

# ---- Resolve project root from script location ----
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0) script_path <- getwd()
R_DIR <- dirname(normalizePath(script_path, winslash = "/"))
PROJECT_ROOT <- normalizePath(file.path(R_DIR, ".."),
                              winslash = "/")
DATA_RAW  <- file.path(PROJECT_ROOT, "data", "raw")
DATA_PROC <- file.path(PROJECT_ROOT, "data", "processed")
GENE_SETS <- file.path(PROJECT_ROOT, "data", "gene_sets")
dir.create(DATA_PROC, showWarnings = FALSE, recursive = TRUE)

EXPR_PATH <- file.path(DATA_RAW, "KIRC_expr_log2_tpm.csv")
GMT_PATH  <- file.path(GENE_SETS, "msigdb_lipid_dedup_236.gmt")
if (!file.exists(GMT_PATH)) {
    GMT_PATH <- file.path(DATA_PROC, "msigdb_lipid_dedup_236.gmt")
}
OUT_PATH  <- file.path(DATA_PROC, "gsva_scores.csv")

# ---- Load expression matrix ----
cat("\n[1/4] Loading expression matrix...\n")
expr <- read.csv(EXPR_PATH, row.names = 1, check.names = FALSE)
expr_mat <- as.matrix(expr)
cat(sprintf("  Genes: %d, Samples: %d\n", nrow(expr_mat), ncol(expr_mat)))

# ---- Load gene sets from GMT ----
cat("\n[2/4] Loading gene sets from GMT...\n")
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

gs_sizes <- sapply(gene_sets_filtered, length)
cat(sprintf("  Gene set size: min=%d, median=%d, max=%d\n",
            min(gs_sizes), median(gs_sizes), max(gs_sizes)))

# ---- Run GSVA ----
cat("\n[3/4] Running GSVA (kcdf=Gaussian)...\n")
cat("  This may take a few minutes for 14k genes x 533 samples...\n")

param <- gsvaParam(
    exprData = expr_mat,
    geneSets = gene_sets_filtered,
    kcdf = "Gaussian"
)

gsva_scores <- gsva(param)

cat(sprintf("  Done! Output: %d pathways x %d samples\n",
            nrow(gsva_scores), ncol(gsva_scores)))

# ---- Save ----
cat("\n[4/4] Saving results...\n")
gsva_df <- as.data.frame(gsva_scores)
write.csv(gsva_df, OUT_PATH)
cat(sprintf("  Saved: %s\n", OUT_PATH))

cat("\n=== Summary ===\n")
cat(sprintf("Method: R GSVA package v%s\n",
            as.character(packageVersion("GSVA"))))
cat(sprintf("KCDF: Gaussian\n"))
cat(sprintf("Expression: log2(TPM+1), %d genes x %d samples\n",
            nrow(expr_mat), ncol(expr_mat)))
cat(sprintf("Gene sets: %d (from msigdb_lipid_dedup_236.gmt)\n",
            length(gene_sets_filtered)))
cat(sprintf("Output: %d pathways x %d samples\n",
            nrow(gsva_scores), ncol(gsva_scores)))
cat(sprintf("Score range: [%.4f, %.4f]\n",
            min(gsva_scores), max(gsva_scores)))
cat(sprintf("Score mean: %.4f, SD: %.4f\n",
            mean(gsva_scores), sd(gsva_scores)))
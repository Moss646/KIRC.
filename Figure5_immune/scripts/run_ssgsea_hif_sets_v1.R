#!/usr/bin/env Rscript
# ============================================================
# ssGSEA for recognized HIF / hypoxia target gene sets
# (sensitivity analysis for the HALLMARK_HYPOXIA null result)
# ============================================================
# Run from repo root:  Rscript scripts/R/run_ssgsea_hif_sets_v1.R
#
# Purpose
#   The manuscript reports no S1/S2 difference for the HALLMARK_HYPOXIA
#   ssGSEA score (S1 0.392 vs S2 0.395, BH P = 0.096) while individual
#   HIF pathway genes (EPAS1, CA9, VEGFA, EGLN3, HIF1A, SLC2A1) are all
#   higher in S2 (BH P < 0.01). This script re-evaluates hypoxia/HIF
#   activity with an a priori panel of recognized HIF/hypoxia target
#   gene sets from MSigDB v2024.1 to test whether the null depends on
#   the choice of gene set.
#
# Pre-committed set selection rule (fixed before running; results are
# reported for ALL selected sets, none dropped post hoc):
#   (a) Hallmark reference      : HALLMARK_HYPOXIA (same gmt as the
#                                 original pipeline)
#   (b) GO:BP annotation        : GOBP_CELLULAR_RESPONSE_TO_DECREASED_
#                                 OXYGEN_LEVELS (GO:0036294). The GO term
#                                 "response to hypoxia" (GO:0001666) was
#                                 obsoleted in the 2023 GO revision; this
#                                 is its surviving successor term, taken
#                                 from MSigDB v2023.2 (the release still
#                                 carrying the pre-revision annotations).
#   (c) Canonical pathways      : REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA,
#                                 PID_HIF1_TFPATHWAY, PID_HIF2PATHWAY,
#                                 PID_HIF1A_PATHWAY, BIOCARTA_HIF_PATHWAY
#   (d) Literature signatures   : SEMENZA_HIF1_TARGETS,
#                                 ELVIDGE_HYPOXIA_UP, MANALO_HYPOXIA_UP,
#                                 WINTER_HYPOXIA_UP,
#                                 BUFFA_HYPOXIA_METAGENE
#   (e) VHL axis (ccRCC-relevant): JIANG_HYPOXIA_VIA_VHL
#   Excluded by the pre-committed minimum set size of 15 genes
#   (standard ssGSEA practice): MAINA_HYPOXIA_VHL_TARGETS_UP (6 genes),
#   REACTOME_REGULATION_OF_GENE_EXPRESSION_BY_HYPOXIA_INDUCIBLE_FACTOR
#   (11 genes).
#
# Method (identical to the published immune pipeline):
#   - ssGSEA via GSVA R package, ssgseaParam(normalize = TRUE)
#   - Expression: log2(TPM+1), genes x samples (same matrix as Fig 5)
#   - S1 vs S2: two-sided Wilcoxon rank-sum (= Mann-Whitney U)
#   - BH correction across the selected HIF/hypoxia sets
#   - Effect size: rank-biserial correlation
#     r_rb = 2*U/(n1*n2) - 1  (positive = S1 higher)
#
# Outputs
#   data/processed/ssgsea_hif_sets_scores_v1.csv   (sets x samples)
#   results/tables/hif_gene_sets_stats_v1.csv      (stats incl. effect size)
#
# Gene set files (MSigDB v2024.1, same release as h.all used by the
# immune pipeline):
#   data/gene_sets/h.all.v2024.1.Hs.symbols.gmt
#   data/gene_sets/c2.all.v2024.1.Hs.symbols.gmt
#   data/gene_sets/c5.all.v2023.2.Hs.symbols.gmt
#   (GO:BP from release 2023.2: the "response to hypoxia" GO terms were
#   retired from the ontology in late 2023 and are absent from MSigDB
#   v2024.1; 2023.2 is the most recent release still containing
#   GOBP_RESPONSE_TO_HYPOXIA / GOBP_CELLULAR_RESPONSE_TO_HYPOXIA)
# ============================================================

suppressPackageStartupMessages({
    library(GSVA)
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
SUB_PATH  <- file.path(DATA_ROOT, "data", "processed",
                       "subtype_assignment_balanced.csv")
GMT_H     <- file.path(DATA_ROOT, "data", "gene_sets",
                       "h.all.v2024.1.Hs.symbols.gmt")
GMT_C2    <- file.path(DATA_ROOT, "data", "gene_sets",
                       "c2.all.v2024.1.Hs.symbols.gmt")
GMT_C5    <- file.path(DATA_ROOT, "data", "gene_sets",
                       "c5.all.v2023.2.Hs.symbols.gmt")
OUT_SCORES <- file.path(DATA_ROOT, "data", "processed",
                        "ssgsea_hif_sets_scores_v1.csv")
OUT_STATS  <- file.path(DATA_ROOT, "results", "tables",
                        "hif_gene_sets_stats_v1.csv")

MIN_SIZE <- 15

# ---- Pre-committed panel ---------------------------------------------------
SET_PANEL <- data.frame(
    set_name = c(
        "HALLMARK_HYPOXIA",
        "GOBP_CELLULAR_RESPONSE_TO_DECREASED_OXYGEN_LEVELS",
        "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",
        "PID_HIF1_TFPATHWAY",
        "PID_HIF2PATHWAY",
        "PID_HIF1A_PATHWAY",
        "BIOCARTA_HIF_PATHWAY",
        "SEMENZA_HIF1_TARGETS",
        "ELVIDGE_HYPOXIA_UP",
        "MANALO_HYPOXIA_UP",
        "WINTER_HYPOXIA_UP",
        "BUFFA_HYPOXIA_METAGENE",
        "JIANG_HYPOXIA_VIA_VHL"
    ),
    collection = c(
        "Hallmark", "GO:BP", "Reactome", "PID", "PID", "PID",
        "BioCarta", "C2:CGP", "C2:CGP", "C2:CGP", "C2:CGP", "C2:CGP",
        "C2:CGP"
    ),
    stringsAsFactors = FALSE
)

# ---- Load expression matrix ------------------------------------------------
cat("\n[1/6] Loading expression matrix...\n")
expr <- read.csv(EXPR_PATH, row.names = 1, check.names = FALSE)
expr_mat <- as.matrix(expr)
cat(sprintf("  Genes: %d, Samples: %d\n", nrow(expr_mat), ncol(expr_mat)))

# ---- Load subtype assignments ----------------------------------------------
cat("\n[2/6] Loading subtype assignments...\n")
sub <- read.csv(SUB_PATH)
sub_map <- setNames(sub$Subtype, sub$Patient)
s1_samples <- names(sub_map[sub_map == "Subtype_1"])
s2_samples <- names(sub_map[sub_map == "Subtype_2"])
cat(sprintf("  S1: %d samples, S2: %d samples\n",
            length(s1_samples), length(s2_samples)))

# ---- Load gene sets from GMT files -----------------------------------------
cat("\n[3/6] Loading gene sets from GMT files...\n")
read_gmt_lines <- function(path) {
    gene_sets <- list()
    con <- file(path, "r")
    while (length(line <- readLines(con, n = 1)) > 0) {
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) >= 3) {
            gs_name <- parts[1]
            gs_genes <- parts[3:length(parts)]
            gs_genes <- gs_genes[gs_genes != ""]
            gene_sets[[gs_name]] <- gs_genes
        }
    }
    close(con)
    gene_sets
}

all_sets <- list()
all_sets <- c(all_sets, read_gmt_lines(GMT_H))
all_sets <- c(all_sets, read_gmt_lines(GMT_C2))
if (file.exists(GMT_C5)) {
    all_sets <- c(all_sets, read_gmt_lines(GMT_C5))
} else {
    cat("  WARNING: C5 gmt not found; GO:BP sets will be skipped\n")
}

# Intersect panel sets with genes present in the expression matrix
genes_available <- rownames(expr_mat)
panel <- lapply(SET_PANEL$set_name, function(s) {
    if (!s %in% names(all_sets)) return(NULL)
    intersect(all_sets[[s]], genes_available)
})
names(panel) <- SET_PANEL$set_name

n_missing <- sum(sapply(panel, is.null))
if (n_missing > 0) {
    cat(sprintf("  WARNING: %d panel sets absent from gmt files: %s\n",
                n_missing,
                paste(SET_PANEL$set_name[sapply(panel, is.null)],
                      collapse = ", ")))
}
too_small <- names(panel)[!sapply(panel, is.null) &
                          sapply(panel, length) < MIN_SIZE]
if (length(too_small) > 0) {
    stop(sprintf("Panel sets below MIN_SIZE=%d after intersection: %s\n",
                 MIN_SIZE, paste(too_small, collapse = ", ")))
}
panel <- panel[!sapply(panel, is.null)]
cat(sprintf("  Panel sets available (>= %d genes in matrix): %d/%d\n",
            MIN_SIZE, length(panel), nrow(SET_PANEL)))
for (s in names(panel)) {
    cat(sprintf("    %-64s %4d genes\n", s, length(panel[[s]])))
}

# ---- Run ssGSEA --------------------------------------------------------------
cat("\n[4/6] Running ssGSEA (GSVA)...\n")
param <- ssgseaParam(
    exprData = expr_mat,
    geneSets = panel,
    normalize = TRUE
)
ss <- gsva(param)
cat(sprintf("  Output: %d sets x %d samples\n", nrow(ss), ncol(ss)))

# ---- S1 vs S2 statistics -----------------------------------------------------
cat("\n[5/6] Wilcoxon rank-sum S1 vs S2 + BH...\n")
s1_cols <- intersect(s1_samples, colnames(ss))
s2_cols <- intersect(s2_samples, colnames(ss))
cat(sprintf("  S1 in scores: %d, S2 in scores: %d\n",
            length(s1_cols), length(s2_cols)))

rows <- list()
for (s in rownames(ss)) {
    v1 <- ss[s, s1_cols]
    v2 <- ss[s, s2_cols]
    wt <- wilcox.test(v1, v2)
    n1 <- length(v1); n2 <- length(v2)
    r_rb <- 2 * as.numeric(wt$statistic) / (n1 * n2) - 1
    rows[[length(rows) + 1]] <- data.frame(
        set_name = s,
        collection = SET_PANEL$collection[match(s, SET_PANEL$set_name)],
        n_genes_in_matrix = length(panel[[s]]),
        S1_mean = mean(v1),
        S2_mean = mean(v2),
        diff_S1_minus_S2 = mean(v1) - mean(v2),
        rank_biserial = r_rb,
        higher = if (mean(v1) > mean(v2)) "S1" else "S2",
        p_raw = wt$p.value,
        stringsAsFactors = FALSE
    )
}
stats_df <- do.call(rbind, rows)
stats_df$p_bh <- p.adjust(stats_df$p_raw, method = "BH")
stats_df <- stats_df[match(SET_PANEL$set_name[
    SET_PANEL$set_name %in% stats_df$set_name], stats_df$set_name), ]

cat("\nResults (S1 vs S2):\n")
for (i in seq_len(nrow(stats_df))) {
    r <- stats_df[i, ]
    sig <- if (r$p_bh < 0.001) "***" else if (r$p_bh < 0.01) "**" else
           if (r$p_bh < 0.05) "*" else "ns"
    cat(sprintf("  %-64s S1=%.4f S2=%.4f  r_rb=%+.3f  p_bh=%.3g  %s\n",
                r$set_name, r$S1_mean, r$S2_mean, r$rank_biserial,
                r$p_bh, sig))
}

# ---- Save --------------------------------------------------------------------
cat("\n[6/6] Saving outputs...\n")
dir.create(dirname(OUT_SCORES), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(OUT_STATS), showWarnings = FALSE, recursive = TRUE)

scores_df <- as.data.frame(ss)
scores_df$set_name <- rownames(scores_df)
scores_df <- scores_df[, c("set_name",
                           setdiff(colnames(scores_df), "set_name"))]
write.csv(scores_df, OUT_SCORES, row.names = FALSE)
write.csv(stats_df, OUT_STATS, row.names = FALSE)
cat(sprintf("Saved scores: %s\n", OUT_SCORES))
cat(sprintf("Saved stats:  %s\n", OUT_STATS))

cat("\n=== Summary ===\n")
cat(sprintf("Method: ssGSEA via GSVA v%s, ssgseaParam(normalize=TRUE)\n",
            as.character(packageVersion("GSVA"))))
cat(sprintf("Expression: log2(TPM+1), %d genes x %d samples\n",
            nrow(expr_mat), ncol(expr_mat)))
cat(sprintf("Sets: %d (pre-committed panel, MIN_SIZE=%d)\n",
            nrow(stats_df), MIN_SIZE))
cat("BH correction: across the selected HIF/hypoxia sets\n")

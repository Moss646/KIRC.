#!/usr/bin/env Rscript
# ============================================================
# GSVA (default method) for the HIF/hypoxia target gene set panel
# ============================================================
# Run from repo root:  Rscript scripts/R/run_gsva_hif_sets_v1.R
#
# Why a second scoring method
#   ssGSEA scores of the same 13 sets split in direction (some S1-higher,
#   some S2-higher) although 82-100% of significant genes within EVERY
#   set are higher in S2 at the gene level. ssGSEA is rank-based and its
#   normalize=TRUE step rescales per sample ACROSS the sets in the call
#   (the HALLMARK_HYPOXIA score alone shifted from ~0.39 to ~0.71 when
#   the panel changed from 50 sets to 13), so the set-level split may be
#   a scoring artifact. GSVA scores are computed per gene set
#   independently (no cross-set normalization) and are the method already
#   used in the manuscript's Figure 3 pipeline (gsvaParam kcdf="Gaussian").
#
# Method
#   - GSVA via GSVA R package, gsvaParam(exprData, geneSets,
#     kcdf = "Gaussian")  -- identical to the Figure 3 pipeline
#   - Expression: log2(TPM+1), genes x samples (same matrix as Fig 5)
#   - S1 vs S2: two-sided Wilcoxon rank-sum (= Mann-Whitney U)
#   - BH correction across the 13 panel sets
#   - Effect size: rank-biserial correlation (positive = S1 higher)
#
# Output
#   data/processed/gsva_hif_sets_scores_v1.csv    (sets x samples)
#   results/tables/hif_gene_sets_gsva_stats_v1.csv
# ============================================================

suppressPackageStartupMessages({
    library(GSVA)
})

cat("R version:", R.version.string, "\n")
cat("GSVA package version:", as.character(packageVersion("GSVA")), "\n")

# ---- Paths ----
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
GMT_H  <- file.path(DATA_ROOT, "data", "gene_sets",
                    "h.all.v2024.1.Hs.symbols.gmt")
GMT_C2 <- file.path(DATA_ROOT, "data", "gene_sets",
                    "c2.all.v2024.1.Hs.symbols.gmt")
GMT_C5 <- file.path(DATA_ROOT, "data", "gene_sets",
                    "c5.all.v2023.2.Hs.symbols.gmt")
OUT_SCORES <- file.path(DATA_ROOT, "data", "processed",
                        "gsva_hif_sets_scores_v1.csv")
OUT_STATS  <- file.path(DATA_ROOT, "results", "tables",
                        "hif_gene_sets_gsva_stats_v1.csv")

PANEL_DEF <- data.frame(
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

# ---- Load data ----
cat("\n[1/4] Loading expression and subtypes...\n")
expr <- read.csv(EXPR_PATH, row.names = 1, check.names = FALSE)
expr_mat <- as.matrix(expr)
sub <- read.csv(SUB_PATH)
sub_map <- setNames(sub$Subtype, sub$Patient)
s1_samples <- names(sub_map[sub_map == "Subtype_1"])
s2_samples <- names(sub_map[sub_map == "Subtype_2"])
cat(sprintf("  Genes: %d, Samples: %d (S1 %d / S2 %d)\n",
            nrow(expr_mat), ncol(expr_mat),
            length(s1_samples), length(s2_samples)))

# ---- Load gene sets ----
cat("\n[2/4] Loading gene sets...\n")
read_gmt_lines <- function(path) {
    gene_sets <- list()
    con <- file(path, "r")
    while (length(line <- readLines(con, n = 1)) > 0) {
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) >= 3) {
            gene_sets[[parts[1]]] <- parts[3:length(parts)]
        }
    }
    close(con)
    gene_sets
}
all_sets <- c(read_gmt_lines(GMT_H), read_gmt_lines(GMT_C2),
              read_gmt_lines(GMT_C5))
panel <- lapply(PANEL_DEF$set_name, function(s) {
    intersect(all_sets[[s]], rownames(expr_mat))
})
names(panel) <- PANEL_DEF$set_name
stopifnot(!any(sapply(panel, function(g) length(g) < 15)))
cat(sprintf("  %d panel sets ready\n", length(panel)))

# ---- GSVA ----
cat("\n[3/4] Running GSVA (gsvaParam, kcdf='Gaussian')...\n")
param <- gsvaParam(
    exprData = expr_mat,
    geneSets = panel,
    kcdf = "Gaussian"
)
gs <- gsva(param)
cat(sprintf("  Output: %d sets x %d samples\n", nrow(gs), ncol(gs)))

# ---- S1 vs S2 statistics ----
cat("\nWilcoxon rank-sum S1 vs S2 + BH...\n")
s1_cols <- intersect(s1_samples, colnames(gs))
s2_cols <- intersect(s2_samples, colnames(gs))
rows <- list()
for (s in rownames(gs)) {
    v1 <- gs[s, s1_cols]
    v2 <- gs[s, s2_cols]
    wt <- wilcox.test(v1, v2)
    n1 <- length(v1); n2 <- length(v2)
    r_rb <- 2 * as.numeric(wt$statistic) / (n1 * n2) - 1
    rows[[length(rows) + 1]] <- data.frame(
        set_name = s,
        collection = PANEL_DEF$collection[match(s, PANEL_DEF$set_name)],
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

cat("\nResults (S1 vs S2; positive rank_biserial = S1 higher):\n")
for (i in seq_len(nrow(stats_df))) {
    r <- stats_df[i, ]
    sig <- if (r$p_bh < 0.001) "***" else if (r$p_bh < 0.01) "**" else
           if (r$p_bh < 0.05) "*" else "ns"
    cat(sprintf("  %-64s S1=%+.3f S2=%+.3f  r_rb=%+.3f  p_bh=%.3g  %s\n",
                r$set_name, r$S1_mean, r$S2_mean, r$rank_biserial,
                r$p_bh, sig))
}

# ---- Save ----
dir.create(dirname(OUT_SCORES), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(OUT_STATS), showWarnings = FALSE, recursive = TRUE)
scores_df <- as.data.frame(gs)
scores_df$set_name <- rownames(scores_df)
scores_df <- scores_df[, c("set_name",
                           setdiff(colnames(scores_df), "set_name"))]
write.csv(scores_df, OUT_SCORES, row.names = FALSE)
write.csv(stats_df, OUT_STATS, row.names = FALSE)
cat(sprintf("\nSaved scores: %s\n", OUT_SCORES))
cat(sprintf("Saved stats:  %s\n", OUT_STATS))

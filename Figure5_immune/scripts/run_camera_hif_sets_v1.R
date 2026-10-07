#!/usr/bin/env Rscript
# ============================================================
# camera (limma) competitive gene-set test for the HIF panel
# ============================================================
# Run from repo root:  Rscript scripts/R/run_camera_hif_sets_v1.R
#
# Why a second method
#   The ssGSEA set-level scores of the same 13 HIF/hypoxia sets split in
#   direction (some S1-higher, some S2-higher), although at the gene level
#   82-100% of significant genes within EVERY set are higher in S2
#   (hif_gene_sets_per_gene_check_v1.csv). ssGSEA scores are rank-based
#   and per-sample normalized across the panel (changing the panel alone
#   shifted the HALLMARK_HYPOXIA score from ~0.39 to ~0.71), so the
#   set-level split may be a scoring artifact rather than biology.
#   camera() (limma) is a recognized competitive gene-set test on the
#   absolute log2 expression scale that avoids per-sample normalization.
#
# Method
#   - Expression: log2(TPM+1), genes x samples (same matrix as Fig 5)
#   - Design: ~ Subtype with Subtype_1 as reference (coefficient = S2 vs S1)
#   - camera(index = panel set, design) per set. Newer limma returns only
#     a two-sided PValue plus a Direction column; the one-sided PUp and
#     PDown are reconstructed from the standard relation
#     p_two_sided = 2 * min(p_up, p_down), so the output keeps the same
#     columns as older limma versions.
#
# Output
#   results/tables/hif_gene_sets_camera_stats_v1.csv
# ============================================================

suppressPackageStartupMessages({
    library(limma)
})

cat("R version:", R.version.string, "\n")
cat("limma version:", as.character(packageVersion("limma")), "\n")

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
OUT_STATS <- file.path(DATA_ROOT, "results", "tables",
                       "hif_gene_sets_camera_stats_v1.csv")

# same pre-committed panel as run_ssgsea_hif_sets_v1.R
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
common <- intersect(colnames(expr_mat), names(sub_map))
expr_mat <- expr_mat[, common]
grp <- factor(sub_map[common], levels = c("Subtype_1", "Subtype_2"))
cat(sprintf("  Genes: %d, Samples: %d (S1 %d / S2 %d)\n",
            nrow(expr_mat), ncol(expr_mat),
            sum(grp == "Subtype_1"), sum(grp == "Subtype_2")))

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

# ---- camera per set ----
cat("\n[3/4] camera per set (design ~ Subtype; coefficient = S2 vs S1)...\n")
design <- model.matrix(~ grp)
colnames(design)[2] <- "S2_vs_S1"

rows <- list()
for (i in seq_len(nrow(PANEL_DEF))) {
    s <- PANEL_DEF$set_name[i]
    cm <- camera(expr_mat, index = panel[[s]], design = design,
                 contrast = 2, inter.gene.cor = NA)

    p_two <- cm$PValue
    dirn  <- cm$Direction
    if (identical(dirn, "Up")) {
        p_up   <- p_two / 2
        p_down <- 1 - p_two / 2
    } else {
        p_down <- p_two / 2
        p_up   <- 1 - p_two / 2
    }

    rows[[length(rows) + 1]] <- data.frame(
        set_name = s,
        collection = PANEL_DEF$collection[i],
        n_genes = length(panel[[s]]),
        mean_logFC_S2_vs_S1 = mean(
            rowMeans(expr_mat[panel[[s]], grp == "Subtype_2", drop = FALSE]) -
            rowMeans(expr_mat[panel[[s]], grp == "Subtype_1", drop = FALSE])),
        direction = dirn,
        p_up_S2 = p_up,
        p_down_S1 = p_down,
        p_two_sided = p_two,
        stringsAsFactors = FALSE
    )
}
stats_df <- do.call(rbind, rows)
stats_df$p_bh <- p.adjust(stats_df$p_two_sided, method = "BH")

cat("\nResults (camera; positive direction = set genes higher in S2):\n")
for (i in seq_len(nrow(stats_df))) {
    r <- stats_df[i, ]
    sig <- if (r$p_bh < 0.001) "***" else if (r$p_bh < 0.01) "**" else
           if (r$p_bh < 0.05) "*" else "ns"
    cat(sprintf("  %-64s logFC=%+.3f  dir=%-4s  p_bh=%.3g  %s\n",
                r$set_name, r$mean_logFC_S2_vs_S1, r$direction,
                r$p_bh, sig))
}

# ---- Save ----
cat("\n[4/4] Saving...\n")
dir.create(dirname(OUT_STATS), showWarnings = FALSE, recursive = TRUE)
write.csv(stats_df, OUT_STATS, row.names = FALSE)
cat(sprintf("Saved: %s\n", OUT_STATS))
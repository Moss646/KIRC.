#!/usr/bin/env Rscript
# 03_fgsea.R - GSEA via fgsea R package
# Input:  data/processed/deg_limma_voom_with_symbols.csv
#         data/gene_sets/msigdb_lipid_dedup_236.gmt
# Output: data/processed/gsea_full_results.csv
#
# Paths are resolved relative to this script's location
# (scripts/upstream/R -> project root).

suppressPackageStartupMessages(library(fgsea))

# ---- Resolve paths (scripts/upstream/R -> project root) ----
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0) script_path <- getwd()
R_DIR <- dirname(normalizePath(script_path, winslash = "/"))
PROJECT_ROOT <- normalizePath(file.path(R_DIR, ".."), winslash = "/")
PROC  <- file.path(PROJECT_ROOT, "data", "processed")
GENE_SETS <- file.path(PROJECT_ROOT, "data", "gene_sets")

# ---- Load DEG results (limma-voom) ----
cat("[1/4] Loading DEG results (limma-voom)...\n")
deg <- read.csv(file.path(PROC, "deg_limma_voom_with_symbols.csv"))
cat(sprintf("  DEG rows: %d\n", nrow(deg)))

# Build ranked list from logFC (remove NAs, empty symbols, duplicates)
deg <- deg[!is.na(deg$gene_symbol) & deg$gene_symbol != "", ]
deg <- deg[order(deg$p_value), ]
deg <- deg[!duplicated(deg$gene_symbol), ]
ranks <- deg$logFC
names(ranks) <- deg$gene_symbol
ranks <- sort(ranks, decreasing = TRUE)
cat(sprintf("  Ranked: %d genes, range = %.3f to %.3f\n",
            length(ranks), min(ranks), max(ranks)))

# ---- Load GMT (236 gene sets, GO:BP-deduplicated) ----
cat("[2/4] Loading GMT (236 gene sets)...\n")
gmt_file <- file.path(GENE_SETS, "msigdb_lipid_dedup_236.gmt")
if (!file.exists(gmt_file)) {
    gmt_file <- file.path(PROC, "msigdb_lipid_dedup_236.gmt")
}
read_gmt <- function(file) {
    lines <- readLines(file)
    gsets <- list()
    for (l in lines) {
        parts <- strsplit(l, "\t")[[1]]
        gsets[[parts[1]]] <- parts[-(1:2)]
    }
    gsets
}
gmt <- read_gmt(gmt_file)
cat(sprintf("  GMT gene sets: %d\n", length(gmt)))

# ---- Run fgsea ----
cat("[3/4] Running fgsea (10,000 permutations)...\n")
set.seed(42)
gsea_res <- fgsea(pathways = gmt, stats = ranks,
                  minSize = 5, maxSize = 2000, nPermSimple = 10000)
gsea_res <- gsea_res[order(gsea_res$pval)]
cat(sprintf("  GSEA complete: %d pathways (p<0.05: %d, FDR<0.05: %d)\n",
            nrow(gsea_res), sum(gsea_res$pval < 0.05),
            sum(gsea_res$padj < 0.05)))

# ---- Format output ----
cat("[4/4] Saving results...\n")
out_df <- data.frame(
    pathway = gsea_res$pathway,
    ES      = gsea_res$ES,
    NES     = gsea_res$NES,
    p_value = gsea_res$pval,
    FDR     = gsea_res$padj,
    n_genes = sapply(gsea_res$leadingEdge, length),
    stringsAsFactors = FALSE
)
write.csv(out_df, file.path(PROC, "gsea_full_results.csv"), row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(PROC, "gsea_full_results.csv")))
cat(sprintf("Method: fgsea R package v%s, nPermSimple = 10000, ranked by limma-voom logFC\n",
            as.character(packageVersion("fgsea"))))
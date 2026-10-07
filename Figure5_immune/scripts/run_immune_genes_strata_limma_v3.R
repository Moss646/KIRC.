#!/usr/bin/env Rscript
# ======================================================================
# run_immune_genes_strata_limma_v3.R
# Stratum-specific limma-voom (S1 vs S2) for the 23 immune marker genes.
#
# Inputs:
#   <repo>/data/raw/KIRC_counts_for_limma.csv
#   <repo>/data/processed/immune_purity_strata_v3.csv
#   <repo>/data/processed/subtype_assignment_balanced.csv
#   <repo>/data/processed/deg_limma_voom_with_symbols.csv
#
# Outputs:
#   results/tables/immune_genes_stratified_limma_v3.csv
#   data/processed/immune_genes_voom_log2cpm_all_v3.csv
# ======================================================================

suppressMessages({library(edgeR); library(limma)})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
    SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
    SCRIPT_DIR <- getwd()
}
PROJ_ROOT <- dirname(SCRIPT_DIR)
DATA_ROOT <- Sys.getenv("MEL_DATA_ROOT", unset = PROJ_ROOT)

counts_path <- file.path(DATA_ROOT, "data", "raw", "KIRC_counts_for_limma.csv")
sym_path    <- file.path(DATA_ROOT, "data", "processed",
                         "deg_limma_voom_with_symbols.csv")
strata_path <- file.path(DATA_ROOT, "data", "processed",
                         "immune_purity_strata_v3.csv")
sub_path    <- file.path(DATA_ROOT, "data", "processed",
                         "subtype_assignment_balanced.csv")
out_stats   <- file.path(DATA_ROOT, "results", "tables",
                         "immune_genes_stratified_limma_v3.csv")
out_voom    <- file.path(DATA_ROOT, "data", "processed",
                         "immune_genes_voom_log2cpm_all_v3.csv")

targets <- c("B2M","HLA-A","HLA-B","HLA-C","HLA-E","TAP1",
             "HLA-DRA","HLA-DRB1","HLA-DQA1","HLA-DQB1","CD74",
             "PRF1","GZMB","GZMA","GNLY","NKG7","IFNG",
             "PDCD1","LAG3","CTLA4","HAVCR2","CD274","TIGIT")

counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
sub    <- read.csv(sub_path)
strata_df <- read.csv(strata_path)
sym    <- read.csv(sym_path, colClasses = c("character","character"))

sym_map <- setNames(sym$gene_symbol, sym$gene)
cat(sprintf("[imm-limma] counts %d x %d | strata rows %d\n",
            nrow(counts), ncol(counts), nrow(strata_df)))

out_list <- list()
voom_all <- NULL
for (st in c("All", "High", "Low")) {
    pts <- strata_df$patient[strata_df$stratum == st]
    keep_cols <- intersect(colnames(counts), pts)
    stopifnot(length(keep_cols) == length(pts))
    info <- sub[match(keep_cols, sub$Patient), , drop = FALSE]
    grp  <- factor(info$Subtype, levels = c("Subtype_2", "Subtype_1"))
    y    <- DGEList(counts = counts[, keep_cols], group = grp)
    keep <- filterByExpr(y)
    y    <- y[keep, , keep.lib.sizes = FALSE]
    y    <- calcNormFactors(y, method = "TMM")
    design <- model.matrix(~ grp)
    v    <- voom(y, design = design, plot = FALSE)
    fit  <- lmFit(v, design)
    fit  <- eBayes(fit, trend = FALSE, robust = FALSE)
    tt   <- topTable(fit, coef = "grpSubtype_1", n = Inf, sort.by = "none")

    tt$gene_symbol <- sym_map[rownames(tt)]
    sub_tt <- tt[tt$gene_symbol %in% targets, ]
    sub_tt <- sub_tt[order(-abs(sub_tt$t)), ]
    sub_tt <- sub_tt[!duplicated(sub_tt$gene_symbol), ]
    out_list[[st]] <- data.frame(
        stratum = st, gene_symbol = sub_tt$gene_symbol,
        logFC = sub_tt$logFC, AveExpr = sub_tt$AveExpr,
        t = sub_tt$t, p_value = sub_tt$P.Value, FDR = sub_tt$adj.P.Val,
        row.names = NULL, stringsAsFactors = FALSE)

    if (st == "All") {
        E <- as.data.frame(v$E)
        E$gene_symbol <- sym_map[rownames(E)]
        E <- E[E$gene_symbol %in% targets, ]
        E <- E[order(-apply(E[, keep_cols, drop = FALSE], 1,
                            function(x) abs(mean(x)))), ]
        E <- E[!duplicated(E$gene_symbol), ]
        voom_all <- data.frame(gene_symbol = E$gene_symbol,
                               E[, keep_cols, drop = FALSE],
                               check.names = FALSE, stringsAsFactors = FALSE)
        rownames(voom_all) <- NULL
    }
    cat(sprintf("[imm-limma] %-4s n=%d (S1=%d, S2=%d): %d genes kept\n",
                st, length(keep_cols), sum(grp == "Subtype_1"),
                sum(grp == "Subtype_2"), sum(keep)))
}

res <- do.call(rbind, out_list)
res <- res[order(match(res$gene_symbol, targets)), ]
res <- res[order(match(res$stratum, c("All","High","Low"))), ]

dir.create(dirname(out_stats), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_voom), showWarnings = FALSE, recursive = TRUE)

write.csv(res, out_stats, row.names = FALSE)
write.csv(voom_all, out_voom, row.names = FALSE)
cat("[imm-limma] written:", out_stats, "\n")
cat("[imm-limma] written:", out_voom, "\n")
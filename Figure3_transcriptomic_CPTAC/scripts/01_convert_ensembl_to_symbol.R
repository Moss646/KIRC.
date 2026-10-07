# Script to convert Ensembl IDs to gene symbols in limma-voom results.
#
# Gene symbols are obtained from the local Bioconductor annotation package
# org.Hs.eg.db. If a previous symbol table is present (either the v9 backup
# file or an earlier deg_limma_voom_with_symbols.csv), it is used as a
# secondary lookup so that overlapping genes keep their original labels.
# The fallback is optional: on a fresh checkout neither file exists and the
# script proceeds with org.Hs.eg.db only.
#
# Usage: Rscript 01_convert_ensembl_to_symbol.R
#   (run from scripts/upstream/R/ or project root)

suppressWarnings({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ---- Resolve paths (scripts/upstream/R -> project root) ----
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0) script_path <- getwd()
R_DIR <- dirname(normalizePath(script_path, winslash = "/"))
PROJECT_ROOT <- normalizePath(file.path(R_DIR, ".."), winslash = "/")
OUT <- file.path(PROJECT_ROOT, "data", "processed")

cat("Input dir:", OUT, "\n")

# Read the (corrected) limma-voom results
deg <- read.csv(file.path(OUT, "deg_limma_voom.csv"), stringsAsFactors = FALSE)

# ---- Optional fallback table ----
fallback <- character(0)

ws_bak <- file.path(OUT, "deg_limma_voom_with_symbols.csv.bak_v9")
ws_cur <- file.path(OUT, "deg_limma_voom_with_symbols.csv")

ws_old <- NULL
if (file.exists(ws_bak)) {
  ws_old <- read.csv(ws_bak, stringsAsFactors = FALSE)
  cat("Fallback source: ", basename(ws_bak), "\n", sep = "")
} else if (file.exists(ws_cur)) {
  ws_old <- read.csv(ws_cur, stringsAsFactors = FALSE)
  cat("Fallback source: ", basename(ws_cur), "\n", sep = "")
} else {
  cat("No previous symbol table found; proceeding with org.Hs.eg.db only\n")
}

if (!is.null(ws_old) && all(c("gene", "gene_symbol") %in% colnames(ws_old))) {
  ws_old_sym <- ws_old[, c("gene", "gene_symbol")]
  ws_old_sym$gene <- sub("\\.[0-9]+$", "", ws_old_sym$gene)
  ws_old_sym <- ws_old_sym[!is.na(ws_old_sym$gene_symbol) &
                           ws_old_sym$gene_symbol != "", ]
  fallback <- setNames(ws_old_sym$gene_symbol, ws_old_sym$gene)
}

# Strip Ensembl version numbers (ENSG00000000003.15 -> ENSG00000000003)
ensembl_ids <- sub("\\.[0-9]+$", "", deg$gene)

# Map Ensembl -> SYMBOL via local org.Hs.eg.db
cat("Mapping Ensembl IDs to symbols via org.Hs.eg.db...\n")
sym_db <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                 column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
gene_symbol <- as.character(sym_db)
names(gene_symbol) <- names(sym_db)

# Fill from the fallback table when org.Hs.eg.db returns NA/empty
n_fallback <- 0
if (length(fallback) > 0) {
  for (i in seq_along(ensembl_ids)) {
    if (is.na(gene_symbol[i]) || gene_symbol[i] == "") {
      f <- fallback[ensembl_ids[i]]
      if (length(f) > 0 && !is.na(f) && f != "") {
        gene_symbol[i] <- f
        n_fallback <- n_fallback + 1
      }
    }
  }
}
cat("Filled via previous symbol table:", n_fallback, "\n")

# Assemble output, preserving column order expected by Fig3_main_*.py
out <- data.frame(
  gene        = deg$gene,
  gene_symbol = gene_symbol,
  logFC       = deg$logFC,
  AveExpr     = deg$AveExpr,
  t           = deg$t,
  p_value     = deg$p_value,
  FDR         = deg$FDR,
  category    = deg$category,
  stringsAsFactors = FALSE
)

output_file <- file.path(OUT, "deg_limma_voom_with_symbols.csv")
write.csv(out, output_file, row.names = FALSE)
cat("Results saved to:", output_file, "\n")

n_total    <- nrow(out)
n_with_sym <- sum(!is.na(out$gene_symbol) & out$gene_symbol != "")
n_without  <- n_total - n_with_sym
cat("Total genes:", n_total, "\n")
cat("Genes with symbol:", n_with_sym, "\n")
cat("Genes without symbol:", n_without, "\n")
cat("\nFirst 10 rows:\n")
print(head(out, 10))
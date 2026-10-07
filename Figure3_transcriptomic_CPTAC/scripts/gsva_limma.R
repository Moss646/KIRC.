# ============================================================
# gsva_limma.R
# GSVA pathway-level differential analysis using R limma eBayes.
# Reads pre-computed GSVA scores (from the GSVA R package, Gaussian
# KCDF) and subtype assignments, fits linear models per pathway, and
# applies empirical Bayes moderation via limma::eBayes().
# Output: data/processed/gsva_subtype_statistics.csv
# ============================================================

suppressPackageStartupMessages(library(limma))

# ---- Resolve project root from script location ----
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0) script_path <- getwd()
R_DIR <- dirname(normalizePath(script_path, winslash = "/"))
PROJECT_ROOT <- normalizePath(file.path(R_DIR, ".."),
                              winslash = "/")
DATA_PROC <- file.path(PROJECT_ROOT, "data", "processed")
dir.create(DATA_PROC, showWarnings = FALSE, recursive = TRUE)

# ---- Input paths ----
gsva_path <- file.path(DATA_PROC, "gsva_scores.csv")
sub_path  <- file.path(DATA_PROC, "subtype_assignment_balanced.csv")

cat("Loading GSVA scores and subtype assignments...\n")
gsva <- read.csv(gsva_path, row.names = 1, check.names = FALSE)
gsva_mat <- as.matrix(gsva)

sub <- read.csv(sub_path)
sub_map <- setNames(sub$Subtype, sub$Patient)

common_samples <- intersect(colnames(gsva_mat), names(sub_map))
cat(sprintf("  Samples: %d (S1=%d, S2=%d)\n",
            length(common_samples),
            sum(sub_map[common_samples] == "Subtype_1"),
            sum(sub_map[common_samples] == "Subtype_2")))
cat(sprintf("  Pathways: %d\n", nrow(gsva_mat)))

subtype <- factor(sub_map[common_samples],
                  levels = c("Subtype_1", "Subtype_2"))
design <- model.matrix(~ subtype)
colnames(design) <- c("Intercept", "Subtype_2")
cat(sprintf("  Design matrix: %d x %d\n", nrow(design), ncol(design)))

cat("\nFitting linear models with limma...\n")
fit <- lmFit(gsva_mat[, common_samples], design)
fit <- eBayes(fit)

cat(sprintf("  Fit pathways: %d\n", nrow(fit)))
cat(sprintf("  Prior df: %.4f\n", fit$df.prior))
cat(sprintf("  Prior variance (s2.prior): %.6f\n", fit$s2.prior))

results_table <- topTable(fit, coef = "Subtype_2", number = Inf,
                          sort.by = "P")
cat(sprintf("  topTable rows: %d\n", nrow(results_table)))

pathway_names <- rownames(results_table)
gene_idx <- match(pathway_names, rownames(gsva_mat))
n <- nrow(results_table)

df_prior_val <- if (length(fit$df.prior) > 1) fit$df.prior[1] else fit$df.prior
df_resid_val <- if (length(fit$df.residual) > 1) fit$df.residual[1] else fit$df.residual

results_out <- data.frame(
    pathway     = pathway_names,
    mean_diff   = -results_table$logFC,
    p_value     = results_table$P.Value,
    FDR         = results_table$adj.P.Val,
    logFC       = -results_table$logFC,
    t           = results_table$t,
    df_prior    = rep(df_prior_val, n),
    df_residual = rep(df_resid_val, n),
    df_total    = rep(df_prior_val + df_resid_val, n),
    s2_ols      = (fit$sigma[gene_idx])^2,
    s2_post     = fit$s2.post[gene_idx],
    stringsAsFactors = FALSE
)
rownames(results_out) <- NULL

n_sig <- sum(results_out$FDR < 0.05)
n_up  <- sum(results_out$mean_diff > 0 & results_out$FDR < 0.05)
n_dn  <- sum(results_out$mean_diff < 0 & results_out$FDR < 0.05)
cat(sprintf("\n  Results: %d pathways, FDR<0.05: %d (%d up in S1, %d down in S1)\n",
            nrow(results_out), n_sig, n_up, n_dn))

cat("\n  Top 10 pathways:\n")
cat(sprintf("  %-55s %12s %10s %10s\n", "Pathway", "FDR", "mean_diff", "t"))
for (i in 1:min(10, nrow(results_out))) {
    pw <- results_out$pathway[i]
    pw_short <- substr(gsub("GOBP_|REACTOME_|HALLMARK_|KEGG_", "", pw), 1, 55)
    cat(sprintf("  %-55s %12.2e %10.4f %10.2f\n",
                pw_short, results_out$FDR[i],
                results_out$mean_diff[i], results_out$t[i]))
}

out_path <- file.path(DATA_PROC, "gsva_subtype_statistics.csv")
write.csv(results_out, out_path, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_path))
cat("Method: R limma eBayes (lmFit + eBayes, Smyth 2004)\n")
cat(sprintf("limma version: %s\n", as.character(packageVersion("limma"))))
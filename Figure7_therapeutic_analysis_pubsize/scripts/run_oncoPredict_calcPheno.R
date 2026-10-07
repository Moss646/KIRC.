# =============================================================================
# Figure 7 upstream: OFFICIAL oncoPredict::calcPhenotype
# GDSC2 -> TCGA-KIRC drug IC50 prediction (9 target drugs)
# =============================================================================
# Output:  data/processed/oncoPredict_calcPheno_results.csv
#          data/processed/oncoPredict_calcPheno_ic50.csv
#          data/processed/oncoPredict_calcPheno_<Drug>_ic50.csv
#
# Input (data/raw/):
#   GDSC2_Expr.rds          GDSC2 cell-line expression (oncoPredict training)
#   GDSC2_Res.rds           GDSC2 drug sensitivity (IC50)
#   kirc_expr_matched.csv   TCGA-KIRC log-expression matrix (genes x samples)
#   kirc_subtypes.csv       Subtype assignment (columns: X=sample, Subtype)
#
# R packages: oncoPredict (with its Bioconductor deps), glmnet, ridge
#
# Path strategy: repo root is auto-detected from the script location, and
#   the MEL_DATA_ROOT environment variable can override it (same convention
#   as the other Figure5 / Figure7 R scripts).
# =============================================================================

suppressMessages(library(oncoPredict))
suppressMessages(library(glmnet))
suppressMessages(library(ridge))

set.seed(42)

# ----- Paths -----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    SCRIPT_DIR <- dirname(normalizePath(script_path, winslash = "/"))
} else {
    SCRIPT_DIR <- getwd()
}
PROJ_ROOT <- Sys.getenv("MEL_DATA_ROOT",
                        unset = dirname(dirname(SCRIPT_DIR)))
DATA_ROOT <- file.path(PROJ_ROOT, "data")

RAW  <- file.path(DATA_ROOT, "raw")
PROC <- file.path(DATA_ROOT, "processed")
dir.create(PROC, showWarnings = FALSE, recursive = TRUE)

# calcPhenotype writes a temporary calcPhenotype_Output/ in the working
# directory; lock that to PROC so nothing appears elsewhere.
setwd(PROC)

# ----- Load data -----
cat("\n=== Loading data ===\n")
gdsc_expr <- readRDS(file.path(RAW, "GDSC2_Expr.rds"))
gdsc_res  <- readRDS(file.path(RAW, "GDSC2_Res.rds"))
kirc_all  <- read.csv(file.path(RAW, "kirc_expr_matched.csv"), row.names = 1)
kirc_sub  <- read.csv(file.path(RAW, "kirc_subtypes.csv"))

cat("GDSC2_Expr:", dim(gdsc_expr), "\n")
cat("GDSC2_Res:",  dim(gdsc_res),  "\n")
cat("KIRC expr:",  dim(kirc_all),  "\n\n")

# ----- Match genes -----
common_genes <- intersect(rownames(gdsc_expr), rownames(kirc_all))
cat("Common genes:", length(common_genes), "\n")
gdsc_expr_sub <- gdsc_expr[common_genes, , drop = FALSE]
kirc_expr_sub <- as.matrix(kirc_all[common_genes, , drop = FALSE])

# ----- Subtype mapping (dots -> dashes) -----
s1_samples <- kirc_sub$X[kirc_sub$Subtype == "Subtype_1"]
s2_samples <- kirc_sub$X[kirc_sub$Subtype == "Subtype_2"]
kirc_ids <- gsub("\\.", "-", colnames(kirc_expr_sub))
s1_idx <- which(kirc_ids %in% s1_samples)
s2_idx <- which(kirc_ids %in% s2_samples)
cat("S1:", length(s1_idx), "S2:", length(s2_idx), "\n\n")

# ----- Target drugs (GDSC2 internal column name -> display name) -----
target_drugs <- c(
  "Axitinib_1021", "Sorafenib_1085", "Rapamycin_1084",
  "Cediranib_1922", "Crizotinib_1083", "Dasatinib_1079",
  "Trametinib_1372", "AZD2014_1441", "AZD8055_1059"
)
target_names <- c(
  "Axitinib", "Sorafenib", "Rapamycin",
  "Cediranib", "Crizotinib", "Dasatinib",
  "Trametinib", "AZD2014", "AZD8055"
)

# ----- Results storage -----
all_ic50 <- matrix(NA, nrow = ncol(kirc_expr_sub), ncol = length(target_drugs))
colnames(all_ic50) <- target_names
rownames(all_ic50) <- colnames(kirc_expr_sub)

results <- data.frame(
  Drug = character(),
  S1_mean = numeric(), S2_mean = numeric(),
  S1_sd = numeric(), S2_sd = numeric(),
  log2FC = numeric(), p = numeric(), FDR = numeric(),
  n_cell_lines = integer(),
  stringsAsFactors = FALSE
)

for (i in 1:length(target_drugs)) {
  drug_col <- target_drugs[i]
  drug_name <- target_names[i]

  cat(sprintf("=== [%d/9] %s ===\n", i, drug_name))

  if (!drug_col %in% colnames(gdsc_res)) {
    cat("  SKIP: drug not found\n"); next
  }

  train_y <- gdsc_res[, drug_col]
  valid_idx <- which(!is.na(train_y))
  if (length(valid_idx) < 10) {
    cat(sprintf("  SKIP: only %d valid cell lines\n", length(valid_idx))); next
  }

  # Two-column trick: keep the training phenotype matrix 2-column so that
  # calcPhenotype's internal trainingPtype[, a] indexing does not collapse
  # to a vector.
  train_y_clean <- cbind(train_y[valid_idx], train_y[valid_idx])
  colnames(train_y_clean) <- c(drug_name, paste0(drug_name, "_dup"))
  rownames(train_y_clean) <- colnames(gdsc_expr_sub[, valid_idx, drop = FALSE])
  train_x_clean <- gdsc_expr_sub[, valid_idx, drop = FALSE]

  cat(sprintf("  N cells=%d, N genes=%d\n",
              ncol(train_x_clean), nrow(train_x_clean)))

  tryCatch({
    unlink("./calcPhenotype_Output", recursive = TRUE)

    preds <- calcPhenotype(
      trainingExprData = train_x_clean,
      trainingPtype   = train_y_clean,
      testExprData    = kirc_expr_sub,
      batchCorrect    = "standardize",
      powerTransformPhenotype = TRUE,
      removeLowVaryingGenes = 0.2,
      minNumSamples = 10,
      selection = 1,
      printOutput = FALSE,
      removeLowVaringGenesFrom = "homogenizeData"
    )

    if (is.null(preds) || !(drug_name %in% colnames(preds))) {
      cat("  WARNING: predictions not returned for this drug\n"); next
    }
    ic50 <- preds[, drug_name]
    names(ic50) <- rownames(preds)

    all_ic50[, i] <- ic50

    write.csv(
      data.frame(sample = names(ic50), predicted_IC50 = ic50),
      file.path(PROC, paste0("oncoPredict_calcPheno_", drug_name, "_ic50.csv")),
      row.names = FALSE
    )

    s1_vals <- ic50[s1_idx]; s2_vals <- ic50[s2_idx]
    s1_m  <- mean(s1_vals, na.rm = TRUE); s2_m <- mean(s2_vals, na.rm = TRUE)
    s1_sd <- sd(s1_vals, na.rm = TRUE);   s2_sd <- sd(s2_vals, na.rm = TRUE)

    # IC50 is already on a log scale; log2FC = S2 - S1
    log2fc <- s2_m - s1_m
    p_val <- t.test(s1_vals, s2_vals)$p.value

    cat(sprintf("  S1=%.3f  S2=%.3f  log2FC(S2-S1)=%+.3f  p=%.2e\n",
                s1_m, s2_m, log2fc, p_val))

    results <- rbind(results, data.frame(
      Drug = drug_name,
      S1_mean = s1_m, S2_mean = s2_m,
      S1_sd = s1_sd, S2_sd = s2_sd,
      log2FC = log2fc, p = p_val, FDR = NA,
      n_cell_lines = ncol(train_x_clean),
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
  })
}

unlink("./calcPhenotype_Output", recursive = TRUE)

results$FDR <- p.adjust(results$p, method = "BH")

write.csv(results, file.path(PROC, "oncoPredict_calcPheno_results.csv"),
          row.names = FALSE)
write.csv(all_ic50, file.path(PROC, "oncoPredict_calcPheno_ic50.csv"),
          row.names = TRUE)

cat("\n============================================\n")
cat("=== Official oncoPredict::calcPhenotype results ===\n")
cat("============================================\n")
cat(sprintf("%-12s %10s %10s %10s %10s\n",
            "Drug", "S1_mean", "S2_mean", "log2FC", "FDR"))
cat(strrep("-", 60), "\n")
for (i in 1:nrow(results)) {
  r <- results[i, ]
  sens <- if (r$S1_mean < r$S2_mean) "S1_sens" else "S2_sens"
  cat(sprintf("%-12s %10.3f %10.3f %+10.3f %10.2e  %s\n",
              r$Drug, r$S1_mean, r$S2_mean, r$log2FC, r$FDR, sens))
}
cat("\nS1_mean < S2_mean: S1 predicted IC50 lower = S1 more sensitive\n")
cat("\n=== DONE ===\n")
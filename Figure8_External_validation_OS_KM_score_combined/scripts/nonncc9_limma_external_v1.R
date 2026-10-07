## limma for the 9 non-NCC FAO genes in E-MTAB / ICGC / CPTAC.
## Called by nonncc9_fao_cross_cohort_v1.py.
## First command-line argument is the relay directory.
## Groups = NCC S1-like / S2-like. BH within cohort.
## logFC > 0 = higher in S1-like.

args <- commandArgs(trailingOnly = TRUE)
relay <- if (length(args) >= 1) args[1] else tempdir()
out_csv <- file.path(relay, "nonncc9_limma_external_v1.csv")

specs <- list(
  list(cohort = "E-MTAB", expr = file.path(relay, "emtab_expr.csv"),
       grp = file.path(relay, "emtab_groups.csv")),
  list(cohort = "ICGC", expr = file.path(relay, "icgc_expr.csv"),
       grp = file.path(relay, "icgc_groups.csv")),
  list(cohort = "CPTAC", expr = file.path(relay, "cptac_expr.csv"),
       grp = file.path(relay, "cptac_groups.csv"))
)

library(limma)

rows <- list()
for (sp in specs) {
  expr <- read.csv(sp$expr, row.names = 1, check.names = FALSE)
  grp <- read.csv(sp$grp, stringsAsFactors = FALSE)
  stopifnot(identical(colnames(expr), grp$sample))

  g <- factor(grp$group, levels = c("S2like", "S1like"))
  design <- model.matrix(~ g)
  fit <- lmFit(as.matrix(expr), design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")

  n1 <- sum(grp$group == "S1like")
  n2 <- sum(grp$group == "S2like")

  for (gn in rownames(tt)) {
    rows[[length(rows) + 1]] <- data.frame(
      cohort = sp$cohort,
      gene = gn,
      logFC = tt[gn, "logFC"],
      t = tt[gn, "t"],
      P.Value = tt[gn, "P.Value"],
      adj.P.Val = tt[gn, "adj.P.Val"],
      n_s1 = n1,
      n_s2 = n2,
      stringsAsFactors = FALSE
    )
  }

  cat(sprintf("%s: %d genes tested, %d BH<0.05\n",
              sp$cohort, nrow(tt), sum(tt$adj.P.Val < 0.05)))
}

out <- do.call(rbind, rows)
write.csv(out, out_csv, row.names = FALSE)
cat(sprintf("written: %s (%d rows)\n", out_csv, nrow(out)))
# R package installation for the KIRC lipid-metabolism subtyping study.
# Run once before executing any R script in the repository.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

bioc_packages <- c(
    "limma",
    "edgeR",
    "GSVA",
    "fgsea",
    "org.Hs.eg.db",
    "AnnotationDbi"
)

cran_packages <- c(
    "glmnet",
    "ridge",
    "oncoPredict",
    "CIBERSORT",
    "xCell"
)

for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
}

for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
    }
}

cat("\nAll R packages installed.\n")
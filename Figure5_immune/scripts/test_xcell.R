# Diagnostic probe for the xCell package's internal spillover data.
# Not part of the analysis pipeline.
#
# Usage:
#   Rscript test_xcell.R                    # read from installed xCell namespace
#   Rscript test_xcell.R <path_to_xCell.data.rda>

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
    e <- new.env()
    load(args[1], envir = e)
    xd <- e$xCell.data
    cat("Loaded from file:", args[1], "\n")
} else {
    suppressPackageStartupMessages(library(xCell))
    xd <- xCell.data
    cat("Loaded from installed xCell package\n")
}

cat("class:", class(xd), "\n")
cat("names:", paste(names(xd), collapse = ", "), "\n")
if (is.list(xd$spill)) {
    cat("spill sub-names:", paste(names(xd$spill), collapse = ", "), "\n")
    cat("K  dim:", paste(dim(xd$spill$K),  collapse = "x"), "\n")
    cat("fv dim:", paste(dim(xd$spill$fv), collapse = "x"), "\n")
}
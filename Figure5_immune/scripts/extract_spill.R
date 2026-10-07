# Extract xCell spill data to RDS (avoids .rda gzip issue).
# Run this ONCE before 03_xCell_analysis.R so xCell uses spillover correction.
# Output: <repo>/data/xcell/xcell_spill.rds

library(xCell)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  SCRIPT_DIR <- getwd()
}
PROJ_ROOT <- dirname(SCRIPT_DIR)
DATA_ROOT <- Sys.getenv("MEL_DATA_ROOT", PROJ_ROOT)

OUT_DIR <- file.path(DATA_ROOT, "data", "xcell")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

spill <- xCell.data$spill
saveRDS(spill, file.path(OUT_DIR, "xcell_spill.rds"))

cat("Spill data saved to", file.path(OUT_DIR, "xcell_spill.rds"), "\n")
cat("names:", paste(names(spill), collapse = ", "), "\n")
cat("K dim:", paste(dim(spill$K), collapse = "x"), "\n")
cat("fv dim:", paste(dim(spill$fv), collapse = "x"), "\n")
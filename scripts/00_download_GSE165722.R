# 00_download_GSE165722.R
# Downloads supplementary files for GSE165722 from GEO and extracts to 'data/' folder

if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
library(GEOquery)

# Create output directory
dir.create("data", recursive = FALSE, showWarnings = FALSE)

# Download supplementary files
getGEOSuppFiles("GSE165722", baseDir = "data")

# Extract all tar files
tar_files <- list.files("data/GSE165722", pattern = "\\.tar$", full.names = TRUE)
for (tarfile in tar_files) {
  message("Extracting ", tarfile)
  untar(tarfile, exdir = "data/GSE165722/extracted")
}

message("Done. Extracted matrices are in 'data/GSE165722/extracted/'")
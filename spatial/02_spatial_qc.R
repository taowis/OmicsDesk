# spatial/spatial_qc.R
# Quality control for spatial transcriptomics datasets (e.g., Visium, Slide-seq).
# Dependencies: Seurat (>=5), SeuratData, dplyr, ggplot2
# Usage:
#   obj <- spatial_qc(load_path = "data/visium/section1/")
#   saveRDS(obj, file = "results/rds/spatial_qc.rds")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

#' spatial_qc
#' @param load_path Path to Visium/Spatial data (folder with 'filtered_feature_bc_matrix' and 'spatial/')
#' @param min_counts Minimum counts per spot
#' @param max_mito  Maximum percent mitochondrial RNA allowed
#' @return Seurat object after QC
spatial_qc <- function(load_path, min_counts = 500, max_mito = 20) {
  obj <- Load10X_Spatial(data.dir = load_path)
  # Calculate mito percent if MT- genes present
  if (any(grepl("^MT-", rownames(obj)))) {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  } else {
    obj[["percent.mt"]] <- 0
  }
  # Basic filters
  obj <- subset(obj, subset = nCount_Spatial >= min_counts & percent.mt <= max_mito)
  obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)
  return(obj)
}

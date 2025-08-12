# spatial/spatial_analysis.R
# Spatially variable gene detection and spatial clustering.
# Dependencies: Seurat, sctransform, spatstat (optional), SeuratDisk
# Usage:
#   st <- readRDS("results/rds/spatial_qc.rds")
#   st2 <- spatial_analysis(st)
#   saveRDS(st2, "results/rds/spatial_analysis.rds")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

#' spatial_analysis
#' @param spatial_obj Seurat object with Spatial assay
#' @param n_pcs Number of principal components for clustering/UMAP
#' @return Spatial object with clustering and SVGs
spatial_analysis <- function(spatial_obj, n_pcs = 30) {
  spatial_obj <- RunPCA(spatial_obj)
  spatial_obj <- FindNeighbors(spatial_obj, dims = 1:n_pcs)
  spatial_obj <- FindClusters(spatial_obj, resolution = 0.5)
  spatial_obj <- RunUMAP(spatial_obj, dims = 1:n_pcs)
  # Identify spatially variable features (SVGs) via Moran's I proxy using SpatiallyVariableFeatures
  if ("Spatial" %in% names(spatial_obj@assays)) {
    svg <- SpatiallyVariableFeatures(spatial_obj, assay = "SCT", selection.method = "moransi", nfeatures = 200)
    VariableFeatures(spatial_obj) <- svg
  }
  return(spatial_obj)
}

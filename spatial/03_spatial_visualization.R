# spatial/spatial_visualization.R
# Plotting helpers for spatial data.
# Dependencies: Seurat, ggplot2, patchwork
# Usage:
#   st <- readRDS("results/rds/spatial_mapped.rds")
#   p <- plot_spatial_labels(st, label_col = "predicted.id")
#   ggsave("results/plots/spatial_labels.png", p, width=7, height=5)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

#' plot_spatial_labels
#' @param spatial_obj Spatial object with predicted labels (e.g., predicted.id)
#' @param label_col Column in meta.data with labels to display
#' @return ggplot object
plot_spatial_labels <- function(spatial_obj, label_col = "predicted.id") {
  if (!label_col %in% colnames(spatial_obj@meta.data)) {
    stop(sprintf("Column '%s' not found in spatial object meta.data", label_col))
  }
  p1 <- SpatialDimPlot(spatial_obj, group.by = label_col, pt.size.factor = 1.6, label = TRUE)
  return(p1)
}

#' plot_svg_feature
#' @param spatial_obj Spatial object
#' @param features Vector of gene names to plot
#' @return patchwork plot
plot_svg_feature <- function(spatial_obj, features) {
  p <- SpatialFeaturePlot(spatial_obj, features = features, ncol = 2)
  return(p)
}

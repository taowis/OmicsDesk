# integration/integration_visualization.R
# Compare embeddings and clustering before/after integration.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

#' plot_integration_comparison
#' @param obj_raw Seurat object before integration
#' @param obj_int Seurat object after integration
#' @param group Column to color by (e.g., 'orig.ident' or 'celltype')
#' @return patchwork plot with side-by-side UMAPs
plot_integration_comparison <- function(obj_raw, obj_int, group = "orig.ident") {
  p1 <- DimPlot(obj_raw, reduction = "umap", group.by = group) + ggtitle("Before integration")
  p2 <- DimPlot(obj_int, reduction = "umap", group.by = group) + ggtitle("After integration")
  return(p1 + p2)
}

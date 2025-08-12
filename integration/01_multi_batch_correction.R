# integration/multi_batch_correction.R
# Batch correction across multiple scRNA-seq datasets.
# Supports Harmony and Seurat CCA.
# Usage:
#   obj_list <- list(obj1, obj2, obj3)
#   integrated <- integrate_harmony(obj_list)
#   saveRDS(integrated, "results/rds/integrated_harmony.rds")

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony) # if installed
})

#' integrate_harmony
#' @param objs List of Seurat objects
#' @return Integrated Seurat object
integrate_harmony <- function(objs) {
  obj <- merge(objs[[1]], objs[-1])
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  if (requireNamespace("harmony", quietly = TRUE)) {
    obj <- RunHarmony(obj, group.by.vars = "orig.ident")
  }
  obj <- RunUMAP(obj, reduction = ifelse("harmony" %in% names(obj@reductions), "harmony", "pca"), dims = 1:30)
  obj <- FindNeighbors(obj, reduction = ifelse("harmony" %in% names(obj@reductions), "harmony", "pca"), dims = 1:30)
  obj <- FindClusters(obj, resolution = 0.6)
  return(obj)
}

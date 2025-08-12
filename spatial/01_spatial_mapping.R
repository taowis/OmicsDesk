# spatial/spatial_mapping.R
# Map scRNA-seq references to spatial spots via label transfer (Seurat anchors).
# Dependencies: Seurat (>=5)
# Usage:
#   st <- readRDS("results/rds/spatial_qc.rds")
#   sc <- readRDS("results/rds/scrna_annotated.rds")
#   mapped <- map_scrna_to_spatial(sc, st)
#   saveRDS(mapped, "results/rds/spatial_mapped.rds")

suppressPackageStartupMessages({
  library(Seurat)
})

#' map_scrna_to_spatial
#' @param scrna_obj Annotated scRNA-seq reference (with celltype labels in Idents or meta.data$celltype)
#' @param spatial_obj Spatial Seurat object (after preprocessing/QC)
#' @param label_col Column name with cell type labels in scrna_obj@meta.data
#' @return Spatial object with transferred labels and prediction scores
map_scrna_to_spatial <- function(scrna_obj, spatial_obj, label_col = "celltype") {
  if (!label_col %in% colnames(scrna_obj@meta.data)) {
    stop(sprintf("Label column '%s' not found in scRNA object meta.data", label_col))
  }
  scrna_obj <- UpdateSeuratObject(scrna_obj)
  spatial_obj <- UpdateSeuratObject(spatial_obj)

  scrna_obj <- NormalizeData(scrna_obj)
  scrna_obj <- FindVariableFeatures(scrna_obj)
  spatial_obj <- NormalizeData(spatial_obj)
  spatial_obj <- FindVariableFeatures(spatial_obj)

  anchors <- FindTransferAnchors(reference = scrna_obj, query = spatial_obj,
                                 dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
  preds <- TransferData(anchorset = anchors, refdata = scrna_obj[[label_col, drop=TRUE]], dims = 1:30)
  spatial_obj <- AddMetaData(spatial_obj, metadata = preds)
  return(spatial_obj)
}

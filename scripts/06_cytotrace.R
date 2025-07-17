# 05_cytotrace.R - CytoTRACE analysis

library(CytoTRACE)
expr_matrix <- as.matrix(seurat_subset@assays$RNA@data)
cyto_out <- CytoTRACE(expr_matrix)
seurat_subset$CytoTRACE <- cyto_out$CytoTRACE[Cells(seurat_subset)]
FeaturePlot(seurat_subset, features = "CytoTRACE", cols = c("blue", "yellow", "red"))
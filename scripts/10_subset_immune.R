# 10_subset_immune.R - Immune subcluster analysis

library(Seurat)

immune_clusters <- c("Macrophage", "T cell", "G-MDSC")  # Update cluster names as needed
seurat_immune <- subset(seurat_obj, idents = immune_clusters)

seurat_immune <- NormalizeData(seurat_immune)
seurat_immune <- FindVariableFeatures(seurat_immune)
seurat_immune <- ScaleData(seurat_immune)
seurat_immune <- RunPCA(seurat_immune)
seurat_immune <- RunUMAP(seurat_immune, dims = 1:15)
seurat_immune <- FindNeighbors(seurat_immune, dims = 1:15)
seurat_immune <- FindClusters(seurat_immune, resolution = 0.5)

DimPlot(seurat_immune, reduction = "umap", label = TRUE)
saveRDS(seurat_immune, file = "seurat_immune.rds")
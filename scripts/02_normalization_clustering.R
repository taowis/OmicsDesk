# 02_normalization_clustering.R
# Normalization, dimensionality reduction, clustering, UMAP and tSNE visualization

library(Seurat)
library(dplyr)
library(SingleR)
library(celldex)
library(ggplot2)

# Load the filtered and merged Seurat object
seurat_obj <- readRDS("results/seurat_qc.rds")

# Normalize the data
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)

# Run PCA
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

# Find neighbors and clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Set identity as clusters (for downstream reference)
Idents(seurat_obj) <- "seurat_clusters"

# Extract reference
ref <- celldex::HumanPrimaryCellAtlasData()

# Extract test data matrix
data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

# Run SingleR (this may take a minute)
singleR.results <- SingleR(test = data.input, ref = ref, labels = ref$label.main)

# Attach results back to Seurat
seurat_obj$celltype <- singleR.results$labels

# Run UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
umap_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "celltype", label = TRUE)
ggsave("results/plots/umap_clusters_default.png", plot = umap_plot, width = 7, height = 5)

# Run tSNE (For GSE165722 dataset)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:10)
tsne_plot <- DimPlot(seurat_obj, reduction = "tsne", group.by = "celltype", label = TRUE)
ggsave("results/plots/tsne_clusters_celltype.png", plot = tsne_plot, width = 7, height = 5)

# Save clustered object
saveRDS(seurat_obj, file = "results/seurat_clustered.rds")
message("✅ Normalization, clustering, and plotting complete.")
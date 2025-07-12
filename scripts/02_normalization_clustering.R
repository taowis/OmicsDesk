# 02_normalization_clustering.R
# Normalization, dimensionality reduction, clustering, UMAP and tSNE visualization

library(Seurat)
library(dplyr)
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

# --- Manual Annotation (example) ---
# Modify the mapping below as needed for your dataset
# Define cluster-to-type mapping
seurat_obj$celltype <- plyr::mapvalues(
  seurat_obj$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23"),
  to = c("Effector", "HT-CLNP", "HT-CLNP", "Effector", "Fibro", "Homeostatic", "HT-CLNP", 
         "HT-CLNP", "Effector", "Fibro", "Adhesion", "Adhesion", "Regulatory", "Fibro", 
         "Effector", "Adhesion", "Regulatory", "Homeostatic", "Regulatory", "Regulatory", 
         "Regulatory", "Regulatory", "Regulatory", "Regulatory")
)

# Run UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
ggsave("results/plots/umap_clusters_default.png", plot = umap_plot, width = 7, height = 5)

# Run tSNE (For GSE165722 dataset)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:10)
tsne_plot <- DimPlot(seurat_obj, reduction = "tsne", group.by = "celltype", label = TRUE)
ggsave("results/plots/tsne_clusters_celltype.png", plot = tsne_plot, width = 7, height = 5)

# Save clustered object
saveRDS(seurat_obj, file = "results/seurat_clustered.rds")
message("✅ Normalization, clustering, and plotting complete.")
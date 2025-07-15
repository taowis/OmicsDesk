# 02_normalization_clustering.R
# Normalization, dimensionality reduction, clustering, UMAP and tSNE visualization

library(Seurat)
library(dplyr)
library(SingleR)
library(celldex)
library(ggplot2)

# Load the filtered and merged Seurat object
seurat_obj <- readRDS("results/rds/seurat_qc.rds")

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
data.input <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")

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

# Run AUCell or GSVA for Signature Score Enrichment

## ✅ Step 1: Load AUCell and create NPC gene sets
library(AUCell)

npc_gene_sets <- list(
  "Effector_NPCs"    = c("CHI3L1", "TNFAIP6", "MMP13"),
  "HT_CLNP"          = c("CD24", "CD44", "POU5F1"),
  "Adhesion_NPCs"    = c("ITGA6", "ITGB1", "COL6A1"),
  "Homeostatic_NPCs" = c("SOX9", "FOXF1"),
  "Fibro_NPCs"       = c("COL1A1", "FAP", "ACTA2"),
  "Regulatory_NPCs"  = c("IL10", "FOXP3", "TGFBR2")
)

## ✅ Step 2: Extract the expression matrix (make sure to use the correct assay/layer)
### Because seurat_obj have been normalized by SCTransform, the assay shoul use "SCT", not "RNA"
exprMatrix <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")

## ✅ Step 3: Build cell rankings and compute AUCell enrichment scores
cells_rankings <- AUCell::AUCell_buildRankings(exprMatrix, plotStats = TRUE, nCores = 1)

### Calculate enrichment for each gene set (by default using the top 5% of expressed genes as the enrichment threshold)
cells_AUC <- AUCell_calcAUC(geneSets = npc_gene_sets, rankings = cells_rankings)

## ✅ Step 4: Integrate the enrichment scores back into the Seurat object as metadata
auc_matrix <- as.data.frame(t(getAUC(cells_AUC)))
colnames(auc_matrix) <- paste0("AUC_", colnames(auc_matrix))

## Add to metadata
seurat_obj <- AddMetaData(seurat_obj, metadata = auc_matrix)

## ✅ Step 5: Visualization (FeaturePlot, RidgePlot)
### 5.1 Visualization of a single subtype
png("results/plots/feature_plot.png", width=800, height=800)
FeaturePlot(seurat_obj, features = "AUC_Effector_NPCs")
dev.off()

### 5.2 Visualization of multiple subtypes together
library(patchwork)
p1 <- FeaturePlot(seurat_obj, features = "AUC_Effector_NPCs") + ggtitle("Effector")
p2 <- FeaturePlot(seurat_obj, features = "AUC_HT_CLNP") + ggtitle("HT-CLNP")
p3 <- FeaturePlot(seurat_obj, features = "AUC_Adhesion_NPCs") + ggtitle("Adhesion")
p4 <- FeaturePlot(seurat_obj, features = "AUC_Homeostatic_NPCs") + ggtitle("Homeostatic")
p5 <- FeaturePlot(seurat_obj, features = "AUC_Fibro_NPCs") + ggtitle("Fibro")
p6 <- FeaturePlot(seurat_obj, features = "AUC_Regulatory_NPCs") + ggtitle("Regulatory")
png("results/plots/ multiple_subtypes_feature_plot.png", width=800, height=800)
(p1 | p2 | p3) / (p4 | p5 | p6)
dev.off()

## ✅ Step 6: (Optional) Aggregate scores at the cluster level and perform automated annotation
Idents(seurat_obj) <- "seurat_clusters"

npc_aggregated <- FetchData(seurat_obj, vars = colnames(auc_matrix)) %>%
  mutate(cluster = Idents(seurat_obj)) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean))

print(npc_aggregated)

npc_cluster_annotation <- npc_aggregated %>%
  rowwise() %>%
  mutate(predicted_npc = names(.)[which.max(c_across(starts_with("AUC_")))]) %>%
  select(cluster, predicted_npc)

print(npc_cluster_annotation)

### Map clusters to NPC subtype names
seurat_obj$npc_subtype_auc <- plyr::mapvalues(
  x = as.character(seurat_obj$seurat_clusters),
  from = npc_cluster_annotation$cluster,
  to   = npc_cluster_annotation$predicted_npc
)

# UMAP plot，group.by = npc_subtype_auc
umap_npc <- DimPlot(seurat_obj, reduction = "umap", group.by = "npc_subtype_auc", label = TRUE) +
  ggtitle("UMAP: NPC Subtypes (AUCell-based)")

# Save to image file
ggsave("results/plots/umap_npc_subtypes_auc.png", plot = umap_npc, width = 7, height = 5)

# tSNE plot
tsne_npc <- DimPlot(seurat_obj, reduction = "tsne", group.by = "npc_subtype_auc", label = TRUE) +
  ggtitle("tSNE: NPC Subtypes (AUCell-based)")

# Save images to result folder
ggsave("results/plots/tsne_npc_subtypes_auc.png", plot = tsne_npc, width = 7, height = 5)

# Save clustered object
saveRDS(seurat_obj, file = "results/rds/seurat_clustered.rds")
message("✅ Normalization, clustering, and plotting complete.")
# 04_trajectory_monocle.R
# Pseudotime inference using Monocle3 with Seurat integration

library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)

# Load annotated Seurat object
seurat_obj <- readRDS("results/seurat_annotated.rds")

# Convert to CellDataSet
cds <- as.cell_data_set(seurat_obj)

# Cluster and learn trajectory graph
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# Optionally set root cells based on a known starting cluster (e.g., cluster 0)
# Replace "0" with another starting cluster as appropriate
root_cells <- colnames(cds)[colData(cds)$seurat_clusters == "0"]
cds <- order_cells(cds, root_cells = root_cells)

# Plot and save pseudotime trajectory
p <- plot_cells(cds,
                color_cells_by = "pseudotime",
                label_groups_by_cluster = FALSE,
                label_leaves = TRUE,
                label_branch_points = TRUE)

ggsave("results/pseudotime_plot.png", plot = p, width = 6, height = 5)

# Optionally embed pseudotime back into Seurat object and save
seurat_obj$pseudotime <- pseudotime(cds)
saveRDS(seurat_obj, "results/seurat_with_pseudotime.rds")

# Save Monocle object
saveRDS(cds, "results/cds_pseudotime.rds")

message("✅ Pseudotime analysis complete.")
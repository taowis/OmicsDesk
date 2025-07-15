# 04_trajectory_monocle.R (Robust version)
# Pseudotime inference using Monocle3 with fallback for automatic root assignment

library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)

message("🔁 Loading Seurat object...")
seurat_obj <- readRDS("results/rds/seurat_annotated.rds")

# Convert to CellDataSet
cds <- as.cell_data_set(seurat_obj)

message("🧠 Running clustering and trajectory learning...")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# Attempt to choose root cells from early cluster (e.g., 0)
if ("seurat_clusters" %in% colnames(colData(cds))) {
  root_cluster <- "0"
  root_cells <- colnames(cds)[colData(cds)$seurat_clusters == root_cluster]

  if (length(root_cells) > 0) {
    cds <- order_cells(cds, root_cells = root_cells)
    message(paste0("✅ Root cells set from seurat_cluster ", root_cluster))
  } else {
    cds <- order_cells(cds)
    message("⚠️ No root cells found in cluster '0'. Using Monocle3 automatic root.")
  }
} else {
  cds <- order_cells(cds)
  message("⚠️ seurat_clusters not found. Using Monocle3 automatic root.")
}

# Plot pseudotime trajectory if values are valid
p <- NULL
if (!all(is.na(pseudotime(cds)))) {
  p <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  label_groups_by_cluster = FALSE,
                  label_leaves = TRUE,
                  label_branch_points = TRUE)
  ggsave("results/plots/pseudotime_plot.png", plot = p, width = 6, height = 5)
  message("✅ Pseudotime plot saved to results/pseudotime_plot.png")
} else {
  message("❌ Pseudotime inference failed: pseudotime values are all NA.")
}

# Save results
seurat_obj$pseudotime <- pseudotime(cds)
saveRDS(seurat_obj, "results/rds/seurat_with_pseudotime.rds")
saveRDS(cds, "results/rds/cds_pseudotime.rds")

message("✅ Trajectory analysis complete.")
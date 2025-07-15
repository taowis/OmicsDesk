# 03_markers_annotation_figures.R
# Marker detection, cluster annotation, and figure generation

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggplot2)

# Load clustered Seurat object
seurat_obj <- readRDS("results/rds/seurat_clustered.rds")
seurat_obj <- PrepSCTFindMarkers(seurat_obj)

# --- Marker Detection ---
markers <- FindAllMarkers(seurat_obj,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)

write.csv(markers, "results/tables/marker_genes_all_clusters.csv", row.names = FALSE)

# --- DotPlot of top markers ---
top10 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

dotplot <- DotPlot(seurat_obj, features = unique(top10$gene), group.by = "celltype") +
  RotatedAxis() +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
  )
ggsave("results/plots/dotplot_top_markers.png", plot = dotplot, width = 14, height = 6, bg = "white")

# --- Cell type proportions per sample ---
prop_df <- as.data.frame(table(seurat_obj$orig.ident, seurat_obj$celltype))
colnames(prop_df) <- c("Sample", "CellType", "Freq")

barplot <- ggplot(prop_df, aes(x = Sample, y = Freq, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Proportion", title = "Cell type distribution per sample") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )
ggsave("results/plots/barplot_celltype_distribution.png", plot = barplot, width = 8, height = 5)

# Save final annotated Seurat object
saveRDS(seurat_obj, file = "results/rds/seurat_annotated.rds")
message("✅ Marker detection, annotation, and figure export complete.")
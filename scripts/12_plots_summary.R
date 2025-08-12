# 06_plots_summary.R
# Generate summary plots

library(Seurat)
library(ggplot2)
library(patchwork)

seurat_obj <- readRDS("results/rds/seurat_annotated.rds")

p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "celltype", label = TRUE) + ggtitle("Cell Type")
p2 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "percent.mt"), ncol = 2)

ggsave("results/plots/umap_celltypes.png", plot = p1, width = 8, height = 6)
ggsave("results/plots/violin_qc.png", plot = p2, width = 8, height = 4)
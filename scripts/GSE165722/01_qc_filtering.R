# 01_qc_filtering_merged.R
# QC script that merges all samples first, detects mitochondrial genes, and handles missing cases

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

# Load previously saved list of Seurat objects
seurat_list <- readRDS("results/seurat_list_all_samples.rds")

# Merge all samples into one Seurat object
seurat_obj <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))

# Try to detect mitochondrial genes
mito.genes <- grep("^MT-|^mt-", rownames(seurat_obj), ignore.case = TRUE, value = TRUE)
if (length(mito.genes) == 0) {
  message("⚠️ No mitochondrial genes found. Assigning percent.mt = 0")
  seurat_obj[["percent.mt"]] <- 0
} else {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features = mito.genes)
}

# Plot QC violin plots (using raw counts layer explicitly)
vln_plot <- VlnPlot(seurat_obj,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3,
                    slot = "counts")
ggsave("results/plots/qc_vlnplot.png", plot = vln_plot, width = 10, height = 4)

# Apply filtering thresholds
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

# Save result
saveRDS(seurat_obj, file = "results/seurat_qc.rds")
message("✅ QC complete. Filtered Seurat object saved to results/seurat_qc_merged.rds")
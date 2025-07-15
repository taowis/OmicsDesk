# 01_qc_filtering.R
# Load and QC scRNA-seq count matrix

library(Seurat)
library(Matrix)
library(dplyr)

# For 10x Cell Ranger format
# Load matrix (modify path to your input file)
counts <- Read10X(data.dir = "data/")

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, project = "SingleCellDesk", min.cells = 3, min.features = 200)

# Add mitochondrial gene percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Visual QC
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

# Save object
saveRDS(seurat_obj, file = "results/rds/seurat_qc.rds")
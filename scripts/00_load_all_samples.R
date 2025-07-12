# 00_load_all_samples.R
# Loads all extracted GSE165722 count matrices into Seurat objects

library(Seurat)
library(dplyr)

# Directory where extracted sample folders are located
sample_base <- "data/GSE165722/extracted"

# Find subfolders that contain 'filtered_feature_bc_matrix'
sample_dirs <- list.dirs(sample_base, recursive = TRUE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("filtered_feature_bc_matrix", sample_dirs)]

# Load each sample into a named list
seurat_list <- list()
for (dir in sample_dirs) {
  sample_name <- basename(dirname(dir))
  message("Loading sample: ", sample_name)
  counts <- Read10X(data.dir = dir)
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
  seurat_obj$orig.ident <- sample_name
  seurat_list[[sample_name]] <- seurat_obj
}

# Optionally save list of Seurat objects
saveRDS(seurat_list, file = "results/seurat_list_all_samples.rds")

message("Loaded ", length(seurat_list), " samples into Seurat objects.")
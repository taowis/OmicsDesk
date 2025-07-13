# 00_load_all_samples_for_GSE165722_matrix_format.R
# Loads count matrices from Ji Tu's GSE165722 format using CellIndex and CellName mapping

library(Seurat)
library(Matrix)
library(data.table)

data_dir <- "data/GSE165722/extracted"
sample_files <- list.files(data_dir, pattern = "counts.tsv.gz$", full.names = TRUE)
seurat_list <- list()

for (count_file in sample_files) {
  sample_id <- gsub("_Sample.*", "", basename(count_file))
  message("Loading sample: ", basename(count_file))

  # Load count matrix and extract gene names
  mat <- fread(count_file)
  gene_names <- mat[[1]]
  mat[[1]] <- NULL
  rownames(mat) <- make.unique(gene_names)  # preserve original gene names

  # Load barcode mapping
  map_file <- gsub("counts.tsv.gz", "cellname.txt.gz", count_file)
  if (file.exists(map_file)) {
    cell_map <- fread(map_file)
    colnames(mat) <- cell_map$CellName
  } else {
    warning("Missing cellname.txt.gz for ", sample_id)
  }

  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = mat, project = sample_id)
  seurat_list[[basename(count_file)]] <- seurat_obj
}

saveRDS(seurat_list, file = "results/seurat_list_all_samples.rds")
message("✅ Loaded ", length(seurat_list), " samples.")
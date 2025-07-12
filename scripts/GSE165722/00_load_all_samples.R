# 00_load_all_samples_for_GSE165722_matrix_format.R
# Loads count matrices from Ji Tu's GSE165722 format using CellIndex and CellName mapping

library(Seurat)
library(Matrix)
library(data.table)

input_dir <- "data/GSE165722/extracted"
files <- list.files(input_dir, pattern = "counts.tsv.gz$", full.names = TRUE)

seurat_list <- list()

for (count_file in files) {
  sample_id <- sub("\\.counts\\.tsv\\.gz$", "", basename(count_file))
  map_file <- file.path(input_dir, paste0(sample_id, ".cellname.txt.gz"))

  message("Loading sample: ", sample_id)

  # Load matrix and cell name mapping
  mat <- fread(count_file)
  rownames(mat) <- mat[[1]]
  mat[[1]] <- NULL
  col_ids <- colnames(mat)

  cell_map <- fread(map_file)
  colnames(cell_map) <- c("CellName", "CellIndex")
  colname_map <- setNames(cell_map$CellName, cell_map$CellIndex)

  # Replace column names (e.g., C1 -> real barcode)
  real_colnames <- colname_map[col_ids]
  if (any(is.na(real_colnames))) {
    stop("Missing cell name mapping for some columns in ", sample_id)
  }
  colnames(mat) <- real_colnames

  # Convert to sparse and create Seurat object
  sparse_mat <- Matrix(as.matrix(mat), sparse = TRUE)
  seurat_obj <- CreateSeuratObject(counts = sparse_mat, project = sample_id)
  seurat_obj$orig.ident <- sample_id
  seurat_list[[sample_id]] <- seurat_obj
}

saveRDS(seurat_list, file = "results/seurat_list_all_samples.rds")
message("✅ Loaded ", length(seurat_list), " samples.")
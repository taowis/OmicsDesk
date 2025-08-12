# integration/scRNA_ST_integration.R
# End-to-end scRNA-seq ↔ Spatial integration using Seurat anchors.
# Usage:
#   sc <- readRDS("results/rds/scrna_annotated.rds")
#   st <- readRDS("results/rds/spatial_qc.rds")
#   st_mapped <- map_scrna_to_spatial(sc, st, label_col="celltype")
#   saveRDS(st_mapped, "results/rds/spatial_mapped.rds")

source(file.path("spatial", "spatial_mapping.R"))

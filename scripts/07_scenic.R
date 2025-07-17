# 06_scenic.R - SCENIC transcription factor analysis

library(SCENIC)
exprMat <- as.matrix(seurat_subset@assays$RNA@counts)
cellInfo <- data.frame(cell_id = colnames(exprMat), cluster = Idents(seurat_subset))
dir.create("SCENIC_output", showWarnings = FALSE)
saveRDS(exprMat, "SCENIC_output/exprMat.Rds")
saveRDS(cellInfo, "SCENIC_output/cellInfo.Rds")

scenicOptions <- initializeScenic(org = "hgnc", dbDir = "cisTarget_databases", datasetTitle = "NPC_SCENIC")
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat)

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
AUCell::plotTSNE(AUC = getAUC(regulonAUC), cellInfo = cellInfo, aucType = "AUC")
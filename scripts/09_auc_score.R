# 08_auc_score.R - AUCell/GSVA-based signature scoring

library(AUCell)
geneSets <- list(SASP = c("IL6", "CXCL8", "MMP3", "SERPINE1", "CCL2"))
exprMat <- as.matrix(seurat_subset@assays$RNA@data)
cells_rankings <- AUCell_buildRankings(exprMat, plotStats = TRUE, nCores = 4)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
seurat_subset$AUC_SASP <- as.vector(getAUC(cells_AUC)["SASP", ])
FeaturePlot(seurat_subset, features = "AUC_SASP", cols = c("lightblue", "red"))
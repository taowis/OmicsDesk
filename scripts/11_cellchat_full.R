# 09_cellchat_full.R - CellChat on all clusters

library(CellChat)
library(patchwork)
library(Seurat)

data.input <- GetAssayData(seurat_obj, slot = "data", assay = "RNA")
meta <- data.frame(group = Idents(seurat_obj))
rownames(meta) <- colnames(seurat_obj)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

netVisual_circle(cellchat@net$weight, vertex.weight = as.numeric(table(cellchat@idents)))
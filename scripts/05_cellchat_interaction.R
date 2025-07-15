# 05_cellchat_interaction.R
# Cell-cell communication analysis using CellChat

library(CellChat)
library(patchwork)
library(Seurat)

seurat_obj <- readRDS("results/rds/seurat_annotated.rds")
data.input <- GetAssayData(seurat_obj, assay = "SCT", layer = "data")
meta <- seurat_obj@meta.data

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "celltype")
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)), weight.scale = T)

saveRDS(cellchat, file = "results/rds/cellchat.rds")
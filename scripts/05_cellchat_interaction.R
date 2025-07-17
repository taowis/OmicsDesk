# 05_cellchat_interaction.R
# Cell-cell communication analysis using CellChat

library(CellChat)
library(patchwork)
library(Seurat)
library(ggplot2)

#seurat_obj <- readRDS("results/rds/seurat_annotated.rds")
seurat_obj <- readRDS("results/rds/seurat_clustered.rds")
data.input <- GetAssayData(seurat_obj, assay = "SCT", layer = "data")
meta <- seurat_obj@meta.data

# 1. Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")

# 2. Set ligand-receptor database (choose human or mouse)
CellChatDB <- CellChatDB.human  # or CellChatDB.mouse
cellchat@DB <- CellChatDB

# 3. Preprocessing steps (these are **critical**)
cellchat <- subsetData(cellchat)  # Error if not ready
# 🔴 So we need to run these first:
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 4. (Optional) PPI projection
# cellchat <- projectData(cellchat, PPI.human)

# 5. Now subsetData will work
cellchat <- subsetData(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

png("results/plots/cellchat_circle_plot.png", width = 2000, height = 2000, res = 300)
netVisual_circle(
  cellchat@net$count,
  vertex.weight = as.numeric(table(cellchat@idents)),
  weight.scale = TRUE
)
dev.off()

saveRDS(cellchat, file = "results/rds/cellchat.rds")
# 07_gsea.R - GSEA pathway analysis using fgsea

library(fgsea)
library(msigdbr)
library(dplyr)

deg <- FindMarkers(seurat_subset, ident.1 = "Grade_IV", ident.2 = "Grade_II")
ranks <- deg %>% arrange(desc(avg_log2FC)) %>%
  mutate(gene = rownames(.)) %>% select(gene, avg_log2FC) %>% deframe()

msigdb <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(msigdb$gene_symbol, msigdb$gs_name)
gsea_res <- fgsea(pathways = pathways, stats = ranks, nperm = 1000)
plotEnrichment(pathways[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]], stats = ranks)
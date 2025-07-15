# SingleCellDesk

**SingleCellDesk** is a modular and user-friendly R toolkit for downstream analysis of single-cell RNA-seq (scRNA-seq) datasets. It enables fast and customizable exploration of count matrix–level data using standard workflows, including quality control, clustering, trajectory inference, and cell-cell interaction analysis.

## 📄 Project Links
- 📂 [Source Code](https://github.com/biosciences/SingleCellDesk): Explore the full repository
- 🔗 [Live Demo Report](https://biosciences.github.io/SingleCellDesk/SingleCellDesk_Analysis.html): View the interactive HTML output


## 🔍 Features

- Quality control with Seurat
- Clustering and UMAP embedding
- Cell type annotation via marker genes or SingleR
- Pseudotime trajectory analysis using Monocle3
- Ligand-receptor signaling inference with CellChat
- Modular script-based structure suitable for any 10x-style input

# The Cell Type Annotation

This project compares two approaches for annotating cell types in nucleus pulposus (NPC) single-cell RNA-seq data.

---

## 🌐 Global Reference Annotation (SingleR)

We used **SingleR** to assign global cell type identities based on reference datasets such as the Human Primary Cell Atlas. This method helps identify broad cell categories (e.g., MSCs, T cells) based on transcriptional similarity across tissues.

The result is visualized in `tsne_clusters_celltype.png`.

---

## 📌 NPC Subtype Annotation via Signature Score Enrichment (AUCell)

To refine cell type annotation beyond global references, we applied **AUCell** to perform signature score enrichment using NPC subtype marker genes reported by Ji Tu *et al.*. This approach quantifies the activation level of each subtype signature per cell. Unlike tsne_clusters_celltype.png, which shows general cell types inferred from reference datasets (e.g., SingleR), tsne_npc_subtypes_auc.png reflects functional NPC subtypes based on marker gene enrichment, offering tissue-specific insights.

The result is visualized in `tsne_npc_subtypes_auc.png`.

---

## 🔍 Key Differences Between Annotation Strategies

| Aspect | 🌐 **Global Reference Annotation (SingleR)** | 📌 **NPC Subtype Annotation via Signature Score Enrichment (AUCell)** |
|--------|---------------------------------------------|------------------------------------------------------------------------|
| **Goal** | Assign broad cell type labels using a global reference dataset | Quantify enrichment of known NPC subtype signatures |
| **Reference** | Global references (e.g., Human Primary Cell Atlas, Blueprint) | NPC subtype marker genes from Ji Tu *et al.* |
| **Resolution** | General (e.g., MSCs, T cells, monocytes) | Fine-grained (e.g., HT-CLNP, FibroNPCs, Effector NPCs) |
| **Method** | Correlation with reference expression profiles | AUCell AUC scoring over ranked gene expression |
| **Output** | Discrete cell type labels per cell | Continuous enrichment scores per cell per signature |
| **Best for** | Cross-tissue or pan-cell type comparisons | Tissue-specific functional annotation of IVD cells |


## 📁 Folder Structure

```
SingleCellDesk/
├── data/                     # Input data (10x count matrix)
├── results/                  # Output Seurat objects, plots, and tables
├── scripts/                  # Modular R scripts for each analysis step
├── paper/                    # JOSS paper and references
├── tests/                    # Unit tests using testthat
├── install.R                 # R dependency installer
├── DESCRIPTION               # Package metadata
├── .github/                  # GitHub Actions CI workflow
└── SingleCellDesk.Rproj      # RStudio project file
```

## 🚀 Getting Started

1. Clone the repository:
```bash
git clone https://github.com/your_username/SingleCellDesk.git
cd SingleCellDesk
```

2. Install R packages:
```r
source("install.R")
```

3. Run the full analysis workflow:
```r
source("scripts/01_qc_filtering.R")
source("scripts/02_normalization_clustering.R")
source("scripts/03_annotation_markers.R")
source("scripts/04_trajectory_monocle.R")
source("scripts/05_cellchat_interaction.R")
source("scripts/06_plots_summary.R")
```

## 🧪 Unit Testing

Run tests to ensure each script performs as expected:
```r
devtools::test()
```

## 📖 Example Dataset

To demonstrate the pipeline, we provide compatibility with GSE165722 (human NP scRNA-seq) and include a mock 10x output folder in `/data`.

## 📥 Loading Your Data

We provide two dedicated loader scripts depending on your data format:

### A. For 10x Genomics Data (e.g., Cell Ranger output folders)

Use:
```r
source("scripts/00_load_all_samples_for_10x_Genomics_data.R")
```

Expected folder layout (one per sample):
```
data/
└── sample1/
    └── filtered_feature_bc_matrix/
        ├── matrix.mtx.gz
        ├── barcodes.tsv.gz
        └── features.tsv.gz
```

### B. For Matrix Format like GSE165722 (TSV + cell name mapping)

Use:
```r
source("scripts/00_load_all_samples_for_GSE165722_matrix_format.R")
```

Expected files per sample:
- `SampleX.counts.tsv.gz`: Count matrix with short column names (e.g. C1, C2, ...)
- `SampleX.cellname.txt.gz`: Two-column mapping file with real cell barcodes

Data source citation:

> J Tu, W Li, S Yang, P Yang, Q Yan, S Wang, K Lai, X Bai, C Wu, W Ding, J Cooper‐White, A Diwan, C Yang, H Yang, J Zou (2021). Single‐Cell Transcriptome Profiling Reveals Multicellular Ecosystem of Nucleus Pulposus during Degeneration Progression. Advanced Science 9(3). Impact Factor: 15.1 (2022)

Output: a list of Seurat objects saved to `results/seurat_list_all_samples.rds`.

## 📄 License

MIT License © Kaitao Lai

## 📬 Citation

If you use SingleCellDesk in your research, please cite the associated JOSS paper (under review). You can also cite the software repository:

> K Lai (2024). *SingleCellDesk: A Modular Toolkit for Downstream Single-Cell RNA-seq Analysis in R*. Journal of Open Source Software (under review). https://github.com/biosciences/SingleCellDesk
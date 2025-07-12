# Title
SingleCellDesk: A Modular Toolkit for Downstream Single-Cell RNA-seq Analysis in R

# Authors
Kaitao Lai  
University of Sydney  
Corresponding author: kaitao.lai@sydney.edu.au

# Summary
We present **SingleCellDesk**, a modular R-based toolkit designed to facilitate downstream analysis of single-cell RNA-seq (scRNA-seq) datasets. This streamlined and structured framework enables analysis of count matrix-level data for quality control, clustering, annotation, trajectory inference, and cell-cell interaction exploration. We illustrate the utility of SingleCellDesk using GSE165722, a human nucleus pulposus (NP) scRNA-seq dataset across various stages of intervertebral disc degeneration (IVDD). SingleCellDesk simplifies scRNA-seq workflows and enables flexible interpretation of complex datasets.

# Statement of Need
Although many scRNA-seq datasets are publicly available, downstream analysis workflows are often inconsistent and fragmented. **SingleCellDesk** addresses this need by offering:
- A structured and user-friendly approach for processing count matrix data
- Integration of standard analytical steps including clustering, annotation, and pseudotime analysis
- Support for visualization and interpretation of cell-cell communication patterns
- Easy adaptation to various biological systems and study designs

This toolkit is valuable for both computational and experimental researchers seeking efficient and transparent ways to explore single-cell transcriptomics.

# Implementation
The toolkit is implemented in R and available as a modular project with well-organized scripts:

- **01_qc_filtering.R**: Performs Seurat-based QC and filtering
- **02_normalization_clustering.R**: Normalizes, clusters, and visualizes cells via PCA and UMAP
- **03_annotation_markers.R**: Identifies cluster markers and supports manual/automated annotation
- **04_trajectory_monocle.R**: Infers pseudotime using Monocle3
- **05_cellchat_interaction.R**: Computes ligand-receptor signaling networks via CellChat
- **06_plots_summary.R**: Assembles summary figures and publication-ready plots

Seurat objects and figures are saved in `results/`. A script `install.R` installs all required packages.

# Example Usage
```r
source("scripts/01_qc_filtering.R")
source("scripts/02_normalization_clustering.R")
source("scripts/03_annotation_markers.R")
source("scripts/04_trajectory_monocle.R")
source("scripts/05_cellchat_interaction.R")
```

# Workflow Capabilities
The pipeline includes support for fixed seed settings and R environment management (e.g., via `renv`) to aid in consistent analysis output. Optional Docker and GitHub Actions configurations are available for collaborative project environments.

# Demonstration Highlights
- CellChat uncovers signaling changes across IVDD stages
- Monocle3 trajectories show alternative bifurcation patterns of progenitor states
- Cell-type annotation using SingleR complements manual marker-based labeling

# Software Repository
https://github.com/biosciences/SingleCellDesk

# Acknowledgements
This work demonstrates analysis using data from Tu et al. (2022), *Advanced Science*. GEO accession: GSE165722.

# References
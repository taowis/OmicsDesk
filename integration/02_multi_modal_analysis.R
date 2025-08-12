# integration/multi_modal_analysis.R
# Multi-modal integration stubs for MOFA+ or scVI/totalVI (Python).
# This file provides R-side wrappers and guidance.
# For scVI/totalVI, consider using reticulate to call Python modules.

#' notes_mofa_scvi
#' This function prints guidance on running MOFA+ or scVI/totalVI.
notes_mofa_scvi <- function() {
  cat("MOFA+/scVI guidance:\n")
  cat("- MOFA+: export matrices and run MOFA+ in Python/R, then import factors for visualization.\n")
  cat("- scVI/totalVI: use reticulate, build AnnData, run model in Python, and import embeddings to Seurat.\n")
}

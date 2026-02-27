# 03_features.R — Feature selection, scaling, and PCA
#
# Identifies highly variable genes (HVGs), scales the data, runs PCA,
# and produces an elbow plot to help choose the number of principal components.
#
# If SCTransform was used in module 02, FindVariableFeatures and ScaleData
# are already done — this module detects that and skips to PCA.

source("R/utils.R")

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

#' Run feature selection, scaling, and PCA.
#'
#' @param seu A Seurat object (output of `run_normalize()`).
#' @param config Named list from `load_config()`.
#' @return A list with components:
#'   - `seu`: Seurat object with PCA computed.
#'   - `plots`: Named list of ggplot2 objects (elbow plot, HVG plot).
#'   - `summary`: Character string summarising results.
run_features <- function(seu, config) {

  set.seed(config$seed)
  n_hvg  <- config$features$n_variable_features %||% 2000
  norm_method <- config$normalization$method
  out_dir <- file.path(config$output$dir, "clustering")  # shared with module 04
  plots <- list()

  # --- HVG + Scaling (skip if SCTransform already handled it) -----------------
  active_assay <- Seurat::DefaultAssay(seu)

  if (tolower(norm_method) == "sctransform" && active_assay == "SCT") {
    log_message("SCTransform detected — HVG selection and scaling already done. Skipping to PCA.")
  } else {
    log_message(sprintf("Finding %d highly variable features (vst method)...", n_hvg))
    seu <- Seurat::FindVariableFeatures(
      seu,
      selection.method = "vst",
      nfeatures = n_hvg,
      verbose = FALSE
    )

    # Variable feature plot
    top10 <- utils::head(Seurat::VariableFeatures(seu), 10)
    plots$hvg <- Seurat::VariableFeaturePlot(seu)
    plots$hvg <- Seurat::LabelPoints(
      plot = plots$hvg,
      points = top10,
      repel = TRUE
    )
    save_plot(plots$hvg, file.path(out_dir, "variable_features"), width = 10, height = 6)

    log_message("Scaling data (variable features only)...")
    seu <- Seurat::ScaleData(seu, verbose = FALSE)
  }

  # --- PCA --------------------------------------------------------------------
  log_message("Running PCA (50 components)...")
  seu <- Seurat::RunPCA(seu, npcs = 50, verbose = FALSE)

  # Elbow plot
  plots$elbow <- Seurat::ElbowPlot(seu, ndims = 50) +
    ggplot2::theme_classic() +
    ggplot2::geom_vline(
      xintercept = config$clustering$n_pcs,
      linetype = "dashed",
      colour = "firebrick"
    ) +
    ggplot2::labs(
      title = "PCA Elbow Plot",
      subtitle = sprintf("Dashed line = %d PCs (config setting)", config$clustering$n_pcs)
    )
  save_plot(plots$elbow, file.path(out_dir, "elbow_plot"))

  # PC loadings for first 2 PCs
  plots$pca_loadings <- Seurat::VizDimLoadings(seu, dims = 1:2, reduction = "pca") +
    ggplot2::theme_classic()
  save_plot(plots$pca_loadings, file.path(out_dir, "pca_loadings"), width = 10, height = 6)

  # --- Save intermediate object -----------------------------------------------
  out_path <- file.path(config$output$dir, "03_seurat_pca.rds")
  tryCatch(
    saveRDS(seu, out_path),
    error = function(e) log_message(
      paste0("Could not save PCA object: ", e$message), "ERROR"
    )
  )
  log_message(sprintf("PCA Seurat object saved to %s", out_path))

  # --- Summary ----------------------------------------------------------------
  n_var <- length(Seurat::VariableFeatures(seu))
  summary_text <- paste0(
    "=== Feature Selection & PCA Summary ===\n",
    sprintf("Variable features: %d\n", n_var),
    "PCA: 50 components computed\n",
    sprintf("Recommended PCs for downstream: %d (set in config.yaml > clustering > n_pcs)\n",
            config$clustering$n_pcs),
    "Check the elbow plot to verify this is a good choice."
  )

  list(seu = seu, plots = plots, summary = summary_text)
}

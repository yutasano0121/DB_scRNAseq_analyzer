# 06_integration.R — Multi-sample / batch correction integration
#
# Supported methods (selected via config.yaml):
#   - harmony:  Fast, scales well to >10 samples (default)
#   - rpca:     Seurat RPCA-based integration
#   - cca:      Seurat CCA-based integration (slower, better for divergent batches)
#   - sketch:   Seurat v5 Sketch integration (for very large datasets)
#
# This module should be run AFTER normalisation (module 02) and PCA (module 03)
# but BEFORE clustering (module 04). If integration is enabled, clustering
# should use the corrected reduction (e.g. "harmony") instead of "pca".

source("R/utils.R")

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

#' Integrate multiple samples / correct batch effects.
#'
#' @param seu A Seurat object with PCA computed (output of `run_features()`).
#'   Must contain multiple samples in `orig.ident`.
#' @param config Named list from `load_config()`.
#' @return A list with components:
#'   - `seu`: Integrated Seurat object with corrected reduction.
#'   - `plots`: Named list of ggplot2 objects (before/after UMAP).
#'   - `summary`: Character string summarising integration.
run_integration <- function(seu, config) {

  set.seed(config$seed)
  method     <- tolower(config$integration$method %||% "harmony")
  n_pcs      <- config$clustering$n_pcs
  norm_method <- config$normalization$method
  out_dir    <- file.path(config$output$dir, "integration")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plots      <- list()

  # --- Check if integration is enabled ----------------------------------------
  if (!isTRUE(config$integration$enabled)) {
    log_message("Integration is disabled in config (integration > enabled: false). Skipping.")
    return(list(
      seu = seu,
      plots = list(),
      summary = "Integration skipped (disabled in config)."
    ))
  }

  # Check that multiple samples exist
  n_samples <- length(unique(seu$orig.ident))
  if (n_samples < 2) {
    log_message(
      "Only one sample found — integration requires multiple samples. Skipping.",
      level = "WARN"
    )
    return(list(
      seu = seu,
      plots = list(),
      summary = "Integration skipped (only 1 sample present)."
    ))
  }

  log_message(sprintf(
    "Integrating %d samples using method: %s", n_samples, method
  ))

  # --- Pre-integration UMAP (for comparison) ----------------------------------
  if (!"umap" %in% Seurat::Reductions(seu)) {
    seu <- Seurat::RunUMAP(seu, dims = 1:n_pcs, verbose = FALSE)
  }

  plots$umap_pre <- Seurat::DimPlot(
    seu, reduction = "umap", group.by = "orig.ident"
  ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Before Integration")
  save_plot(plots$umap_pre, file.path(out_dir, "umap_pre_integration"))

  # --- Run integration --------------------------------------------------------
  reduction_to_use <- switch(method,
    harmony = integrate_harmony(seu, n_pcs),
    rpca    = integrate_seurat(seu, n_pcs, norm_method, reduction = "rpca"),
    cca     = integrate_seurat(seu, n_pcs, norm_method, reduction = "cca"),
    sketch  = integrate_sketch(seu, config),
    stop(sprintf(
      "Unknown integration method '%s'. Use 'harmony', 'rpca', 'cca', or 'sketch'.",
      method
    ))
  )

  seu <- reduction_to_use$seu
  red_name <- reduction_to_use$reduction_name

  # --- Post-integration UMAP -------------------------------------------------
  log_message(sprintf("Recomputing UMAP from '%s' reduction...", red_name))
  seu <- Seurat::RunUMAP(seu, reduction = red_name, dims = 1:n_pcs, verbose = FALSE)
  seu <- Seurat::FindNeighbors(seu, reduction = red_name, dims = 1:n_pcs, verbose = FALSE)

  plots$umap_post <- Seurat::DimPlot(
    seu, reduction = "umap", group.by = "orig.ident"
  ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("After Integration (%s)", method))
  save_plot(plots$umap_post, file.path(out_dir, "umap_post_integration"))

  # Side-by-side comparison
  plots$comparison <- patchwork::wrap_plots(plots$umap_pre, plots$umap_post, ncol = 2)
  save_plot(plots$comparison, file.path(out_dir, "integration_comparison"), width = 14, height = 6)

  # --- Save intermediate object -----------------------------------------------
  out_path <- file.path(config$output$dir, "06_seurat_integrated.rds")
  tryCatch(
    saveRDS(seu, out_path),
    error = function(e) log_message(
      paste0("Could not save integrated object: ", e$message), "ERROR"
    )
  )
  log_message(sprintf("Integrated Seurat object saved to %s", out_path))

  # --- Summary ----------------------------------------------------------------
  summary_text <- paste0(
    "=== Integration Summary ===\n",
    sprintf("Method: %s\n", method),
    sprintf("Samples integrated: %d\n", n_samples),
    sprintf("PCs used: %d\n", n_pcs),
    sprintf("Corrected reduction: '%s'\n", red_name),
    "UMAP and neighbour graph recomputed from corrected reduction.\n",
    "Proceed to clustering (module 04) — it will use the corrected embeddings."
  )

  list(seu = seu, plots = plots, summary = summary_text)
}

# ---------------------------------------------------------------------------
# Integration method implementations
# ---------------------------------------------------------------------------

#' Harmony integration.
#'
#' @param seu Seurat object with PCA.
#' @param n_pcs Number of PCs.
#' @return List with `seu` and `reduction_name`.
integrate_harmony <- function(seu, n_pcs) {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("harmony package is required for Harmony integration. Install with: install.packages('harmony')")
  }
  log_message("Running Harmony...")
  seu <- harmony::RunHarmony(seu, group.by.vars = "orig.ident")
  list(seu = seu, reduction_name = "harmony")
}

#' Seurat CCA or RPCA integration.
#'
#' @param seu Seurat object.
#' @param n_pcs Number of PCs.
#' @param norm_method Normalisation method ("LogNormalize" or "SCTransform").
#' @param reduction "rpca" or "cca".
#' @return List with `seu` and `reduction_name`.
integrate_seurat <- function(seu, n_pcs, norm_method, reduction = "rpca") {
  log_message(sprintf("Running Seurat %s integration...", toupper(reduction)))

  seu_list <- Seurat::SplitObject(seu, split.by = "orig.ident")

  if (tolower(norm_method) == "sctransform") {
    seu_list <- lapply(seu_list, function(x) {
      Seurat::SCTransform(x, verbose = FALSE)
    })
    features <- Seurat::SelectIntegrationFeatures(seu_list, nfeatures = 3000)
    seu_list <- Seurat::PrepSCTIntegration(seu_list, anchor.features = features)

    anchors <- Seurat::FindIntegrationAnchors(
      seu_list,
      normalization.method = "SCT",
      anchor.features = features,
      dims = 1:n_pcs,
      reduction = reduction
    )
    seu_int <- Seurat::IntegrateData(
      anchors,
      normalization.method = "SCT",
      dims = 1:n_pcs
    )
  } else {
    seu_list <- lapply(seu_list, function(x) {
      x <- Seurat::NormalizeData(x, verbose = FALSE)
      Seurat::FindVariableFeatures(x, verbose = FALSE)
    })
    anchors <- Seurat::FindIntegrationAnchors(
      seu_list,
      dims = 1:n_pcs,
      reduction = reduction
    )
    seu_int <- Seurat::IntegrateData(anchors, dims = 1:n_pcs)
  }

  Seurat::DefaultAssay(seu_int) <- "integrated"
  seu_int <- Seurat::ScaleData(seu_int, verbose = FALSE)
  seu_int <- Seurat::RunPCA(seu_int, npcs = n_pcs, verbose = FALSE)

  list(seu = seu_int, reduction_name = "pca")
}

#' Sketch integration (for very large datasets).
#'
#' @param seu Seurat object.
#' @param config Full config list.
#' @return List with `seu` and `reduction_name`.
integrate_sketch <- function(seu, config) {
  n_sketch <- config$integration$sketch_ncells %||% 5000
  n_pcs    <- config$clustering$n_pcs

  log_message(sprintf("Running Sketch integration (ncells = %d)...", n_sketch))

  seu <- Seurat::SketchData(
    seu,
    ncells = n_sketch,
    method = "LeverageScore",
    sketched.assay = "sketch"
  )
  Seurat::DefaultAssay(seu) <- "sketch"
  seu <- Seurat::FindVariableFeatures(seu, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  seu <- Seurat::RunPCA(seu, npcs = n_pcs, verbose = FALSE)

  list(seu = seu, reduction_name = "pca")
}

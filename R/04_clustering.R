# 04_clustering.R — Nearest-neighbour graph, clustering, and UMAP / tSNE
#
# Builds a k-NN graph from PCA, clusters cells using Leiden (preferred) or
# Louvain, and computes UMAP (and optionally tSNE) for 2D visualisation.

source("R/utils.R")

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

#' Cluster cells and compute dimensionality reduction embeddings.
#'
#' @param seu A Seurat object with PCA computed (output of `run_features()`).
#' @param config Named list from `load_config()`.
#' @return A list with components:
#'   - `seu`: Seurat object with clusters and UMAP.
#'   - `plots`: Named list of ggplot2 objects.
#'   - `summary`: Character string summarising clustering results.
run_clustering <- function(seu, config) {

  set.seed(config$seed)
  n_pcs      <- config$clustering$n_pcs
  resolution <- config$clustering$resolution
  algorithm  <- config$clustering$algorithm %||% 4
  run_tsne   <- isTRUE(config$clustering$run_tsne)
  out_dir    <- file.path(config$output$dir, "clustering")
  plots      <- list()

  # --- Nearest-neighbour graph ------------------------------------------------
  log_message(sprintf("Building k-NN graph using %d PCs...", n_pcs))
  seu <- Seurat::FindNeighbors(seu, dims = 1:n_pcs, verbose = FALSE)

  # --- Clustering -------------------------------------------------------------
  algo_name <- if (algorithm == 4) "Leiden" else "Louvain"

  # Leiden requires leidenAlg or igraph with Leiden support; fall back to Louvain
  if (algorithm == 4 && !requireNamespace("leidenAlg", quietly = TRUE)) {
    log_message(
      "leidenAlg package not found — falling back to Louvain (algorithm = 1). Install leidenAlg for Leiden clustering.",
      level = "WARN"
    )
    algorithm <- 1
    algo_name <- "Louvain (fallback)"
  }

  log_message(sprintf(
    "Clustering with %s, resolution = %.2f...", algo_name, resolution
  ))
  seu <- Seurat::FindClusters(
    seu,
    resolution = resolution,
    algorithm  = algorithm,
    verbose    = FALSE
  )

  n_clusters <- length(levels(Seurat::Idents(seu)))
  log_message(sprintf("Found %d clusters.", n_clusters))

  # --- UMAP ------------------------------------------------------------------
  log_message(sprintf("Computing UMAP (%d PCs)...", n_pcs))
  seu <- Seurat::RunUMAP(seu, dims = 1:n_pcs, verbose = FALSE)

  plots$umap <- Seurat::DimPlot(
    seu,
    reduction = "umap",
    label = TRUE,
    label.size = 4
  ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf(
      "UMAP — %d clusters (%s, res=%.2f)", n_clusters, algo_name, resolution
    ))
  save_plot(plots$umap, file.path(out_dir, "umap_clusters"))

  # UMAP coloured by sample (useful for multi-sample)
  if ("orig.ident" %in% colnames(seu@meta.data)) {
    plots$umap_sample <- Seurat::DimPlot(
      seu,
      reduction = "umap",
      group.by = "orig.ident"
    ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "UMAP — coloured by sample")
    save_plot(plots$umap_sample, file.path(out_dir, "umap_by_sample"))
  }

  # --- tSNE (optional) -------------------------------------------------------
  if (run_tsne) {
    log_message("Computing tSNE...")
    seu <- Seurat::RunTSNE(seu, dims = 1:n_pcs, verbose = FALSE)

    plots$tsne <- Seurat::DimPlot(
      seu,
      reduction = "tsne",
      label = TRUE,
      label.size = 4
    ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "tSNE — clusters")
    save_plot(plots$tsne, file.path(out_dir, "tsne_clusters"))
  }

  # --- Save intermediate object -----------------------------------------------
  out_path <- file.path(config$output$dir, "04_seurat_clustered.rds")
  tryCatch(
    saveRDS(seu, out_path),
    error = function(e) log_message(
      paste0("Could not save clustered object: ", e$message), "ERROR"
    )
  )
  log_message(sprintf("Clustered Seurat object saved to %s", out_path))

  # --- Summary ----------------------------------------------------------------
  cluster_sizes <- table(Seurat::Idents(seu))
  size_str <- paste(
    sprintf("  Cluster %s: %d cells", names(cluster_sizes), as.integer(cluster_sizes)),
    collapse = "\n"
  )
  summary_text <- paste0(
    "=== Clustering Summary ===\n",
    sprintf("Algorithm: %s\n", algo_name),
    sprintf("Resolution: %.2f\n", resolution),
    sprintf("PCs used: %d\n", n_pcs),
    sprintf("Clusters found: %d\n", n_clusters),
    "\nCluster sizes:\n", size_str
  )

  list(seu = seu, plots = plots, summary = summary_text)
}

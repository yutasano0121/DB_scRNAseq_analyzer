# 08_trajectory.R — Trajectory inference and pseudotime
#
# Method:
#   1. Monocle3 — Pure R trajectory inference with pseudotime ordering

source("R/utils.R")

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

#' Run trajectory inference.
#'
#' @param seu A Seurat object with clusters and UMAP (output of `run_clustering()`).
#' @param config Named list from `load_config()`.
#' @return A list with components:
#'   - `seu`: Seurat object (may have pseudotime added to metadata).
#'   - `plots`: Named list of ggplot2 objects.
#'   - `summary`: Character string summarising results.
run_trajectory <- function(seu, config) {

  set.seed(config$seed)
  out_dir <- file.path(config$output$dir, "trajectory")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plots   <- list()
  summary_parts <- c("=== Trajectory Summary ===")

  # Monocle3 is the only trajectory method supported in this project.
  result <- run_monocle3(seu, config, out_dir)
  seu <- result$seu
  plots <- c(plots, result$plots)
  summary_parts <- c(summary_parts, result$summary)

  # --- Save intermediate object -----------------------------------------------
  out_path <- file.path(config$output$dir, "08_seurat_trajectory.rds")
  tryCatch(
    saveRDS(seu, out_path),
    error = function(e) log_message(
      paste0("Could not save trajectory object: ", e$message), "ERROR"
    )
  )

  list(
    seu = seu,
    plots = plots,
    summary = paste(summary_parts, collapse = "\n\n")
  )
}

# ---------------------------------------------------------------------------
# Monocle3
# ---------------------------------------------------------------------------

#' Run Monocle3 trajectory analysis.
#'
#' Converts the Seurat object to a Monocle3 CellDataSet, learns the trajectory
#' graph, orders cells in pseudotime, and produces plots.
#'
#' @param seu Seurat object.
#' @param config Config list.
#' @param out_dir Output directory.
#' @return List with `seu`, `plots`, and `summary`.
run_monocle3 <- function(seu, config, out_dir) {

  if (!requireNamespace("monocle3", quietly = TRUE) ||
      !requireNamespace("SeuratWrappers", quietly = TRUE)) {
    log_message(
      "monocle3 or SeuratWrappers not installed — skipping Monocle3 trajectory.",
      "WARN"
    )
    return(list(
      seu = seu,
      plots = list(),
      summary = "Monocle3: skipped (packages not installed)."
    ))
  }

  log_message("Running Monocle3 trajectory inference...")
  plots <- list()

  # Convert Seurat → CellDataSet (transfers UMAP coordinates)
  cds <- tryCatch(
    SeuratWrappers::as.cell_data_set(seu),
    error = function(e) {
      log_message(paste0("Seurat → Monocle3 conversion failed: ", e$message), "ERROR")
      NULL
    }
  )

  if (is.null(cds)) {
    return(list(
      seu = seu, plots = list(),
      summary = "Monocle3: conversion from Seurat failed."
    ))
  }

  # Cluster cells (Monocle3's own clustering for trajectory)
  cds <- monocle3::cluster_cells(cds)

  # Learn trajectory graph
  cds <- monocle3::learn_graph(cds, use_partition = TRUE)

  # Order cells — use the cluster with most cells as root (heuristic)
  # A proper root selection would be interactive or marker-based
  root_cells <- select_root_cells(cds, seu)
  cds <- monocle3::order_cells(cds, root_cells = root_cells)

  # Transfer pseudotime back to Seurat
  seu$monocle3_pseudotime <- monocle3::pseudotime(cds)

  # Plot: pseudotime on trajectory
  plots$monocle3_pseudotime <- monocle3::plot_cells(
    cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE
  )
  save_plot(plots$monocle3_pseudotime, file.path(out_dir, "monocle3_pseudotime"))

  # Plot: clusters on trajectory
  plots$monocle3_clusters <- monocle3::plot_cells(
    cds,
    color_cells_by = "cluster",
    show_trajectory_graph = TRUE,
    label_cell_groups = TRUE
  )
  save_plot(plots$monocle3_clusters, file.path(out_dir, "monocle3_clusters"))

  # UMAP coloured by pseudotime via ggplot2
  plots$umap_pseudotime <- ggplot2::ggplot(
    seu@meta.data,
    ggplot2::aes(
      x = seu@reductions$umap@cell.embeddings[, 1],
      y = seu@reductions$umap@cell.embeddings[, 2],
      colour = monocle3_pseudotime
    )
  ) +
    ggplot2::geom_point(size = 0.3, alpha = 0.6) +
    viridis::scale_color_viridis(option = "magma", na.value = "grey80") +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = "Pseudotime (Monocle3)",
      x = "UMAP 1", y = "UMAP 2",
      colour = "Pseudotime"
    )
  save_plot(plots$umap_pseudotime, file.path(out_dir, "umap_pseudotime"))

  summary_text <- sprintf(
    "Monocle3: trajectory learned, pseudotime computed for %d cells.\nPseudotime range: [%.2f, %.2f]",
    sum(!is.na(seu$monocle3_pseudotime)),
    min(seu$monocle3_pseudotime, na.rm = TRUE),
    max(seu$monocle3_pseudotime, na.rm = TRUE)
  )

  log_message("Monocle3 trajectory complete.")
  list(seu = seu, plots = plots, summary = summary_text)
}

#' Heuristic root cell selection for Monocle3.
#'
#' Selects cells in the largest cluster's earliest partition as root.
#' In a real analysis, the user should select root cells based on biological
#' knowledge (e.g., stem cell markers).
#'
#' @param cds A Monocle3 CellDataSet.
#' @param seu The original Seurat object (for metadata).
#' @return Character vector of cell IDs to use as root.
select_root_cells <- function(cds, seu) {
  # Use cells from the largest cluster as a default root heuristic
  cluster_counts <- table(Seurat::Idents(seu))
  largest_cluster <- names(which.max(cluster_counts))
  root_cells <- colnames(seu)[Seurat::Idents(seu) == largest_cluster]

  log_message(sprintf(
    "Root cells: using cluster %s (%d cells) as heuristic root. Override with biological markers for best results.",
    largest_cluster, length(root_cells)
  ))

  root_cells
}

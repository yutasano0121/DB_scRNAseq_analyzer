# 01_qc.R — Quality control metrics, visualisation, and cell filtering
#
# Computes per-cell QC metrics (nFeature_RNA, nCount_RNA, percent.mt,
# percent.ribo), generates violin/scatter plots, and filters cells based on
# thresholds from config.yaml.

source("R/utils.R")

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

#' Run QC metrics computation, visualisation, and filtering.
#'
#' @param seu A Seurat object (output of `run_load_data()`).
#' @param config Named list from `load_config()`.
#' @return A list with components:
#'   - `seu`: Filtered Seurat object.
#'   - `plots`: Named list of ggplot2 objects.
#'   - `summary`: Character string summarising cells before/after filtering.
run_qc <- function(seu, config) {

  set.seed(config$seed)
  species   <- config$species
  qc_config <- config$qc
  out_dir   <- file.path(config$output$dir, "qc")

  cells_before <- ncol(seu)
  log_message(sprintf("Starting QC: %d cells, %d genes.", cells_before, nrow(seu)))

  # --- Compute QC metrics -----------------------------------------------------
  mt_pattern   <- mito_prefix(species)
  ribo_pattern <- "^RP[SL]"

  seu[["percent.mt"]]   <- Seurat::PercentageFeatureSet(seu, pattern = mt_pattern)
  seu[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seu, pattern = ribo_pattern)

  log_message("QC metrics computed: nFeature_RNA, nCount_RNA, percent.mt, percent.ribo.")

  # --- Visualise before filtering ---------------------------------------------
  plots <- list()

  # Violin plots
  plots$violin <- Seurat::VlnPlot(
    seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
    ncol = 4,
    pt.size = 0
  ) + ggplot2::theme_classic()
  save_plot(plots$violin, file.path(out_dir, "violin_qc_pre_filter"),
            width = 14, height = 5)

  # Scatter: nCount vs nFeature, coloured by percent.mt
  plots$scatter_mt <- ggplot2::ggplot(
    seu@meta.data,
    ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.mt)
  ) +
    ggplot2::geom_point(size = 0.3, alpha = 0.5) +
    viridis::scale_color_viridis(option = "magma") +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = "UMI Counts vs Genes Detected",
      x = "UMI Counts (nCount_RNA)",
      y = "Genes Detected (nFeature_RNA)",
      colour = "% Mito"
    )
  save_plot(plots$scatter_mt, file.path(out_dir, "scatter_count_vs_feature"))

  # --- Filter cells -----------------------------------------------------------
  log_message(sprintf(
    "Filtering thresholds: nFeature [%d, %d], nCount [%d, %d], max_mito=%.1f%%",
    qc_config$min_features, qc_config$max_features,
    qc_config$min_counts,   qc_config$max_counts,
    qc_config$max_mito
  ))

  seu <- subset(
    seu,
    subset = nFeature_RNA > qc_config$min_features &
             nFeature_RNA < qc_config$max_features &
             nCount_RNA   > qc_config$min_counts   &
             nCount_RNA   < qc_config$max_counts   &
             percent.mt   < qc_config$max_mito
  )

  cells_after <- ncol(seu)
  cells_removed <- cells_before - cells_after
  log_message(sprintf(
    "Filtering complete: %d → %d cells (%d removed, %.1f%%).",
    cells_before, cells_after, cells_removed,
    100 * cells_removed / cells_before
  ))

  # --- Visualise after filtering -----------------------------------------------
  plots$violin_post <- Seurat::VlnPlot(
    seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
    ncol = 4,
    pt.size = 0
  ) + ggplot2::theme_classic()
  save_plot(plots$violin_post, file.path(out_dir, "violin_qc_post_filter"),
            width = 14, height = 5)

  # --- Save intermediate object -----------------------------------------------
  out_path <- file.path(config$output$dir, "01_seurat_qc.rds")
  tryCatch(
    saveRDS(seu, out_path),
    error = function(e) log_message(
      paste0("Could not save QC object: ", e$message), "ERROR"
    )
  )
  log_message(sprintf("QC Seurat object saved to %s", out_path))

  # --- Summary ----------------------------------------------------------------
  summary_text <- paste0(
    "=== QC Summary ===\n",
    sprintf("Cells before filtering: %d\n", cells_before),
    sprintf("Cells after filtering:  %d\n", cells_after),
    sprintf("Cells removed:          %d (%.1f%%)\n",
            cells_removed, 100 * cells_removed / cells_before),
    sprintf("Genes:                  %d\n", nrow(seu)),
    "\nThresholds applied:\n",
    sprintf("  nFeature_RNA: [%d, %d]\n", qc_config$min_features, qc_config$max_features),
    sprintf("  nCount_RNA:   [%d, %d]\n", qc_config$min_counts, qc_config$max_counts),
    sprintf("  percent.mt:   < %.1f%%\n", qc_config$max_mito)
  )

  list(seu = seu, plots = plots, summary = summary_text)
}

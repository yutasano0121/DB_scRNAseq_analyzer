# 05_annotation.R — Cell type annotation
#
# Two approaches:
#   1. Marker-based (always run): FindAllMarkers per cluster, heatmap of top markers
#   2. SingleR (optional): Automated reference-based annotation using celldex

source("R/utils.R")

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

#' Annotate clusters with cell type labels.
#'
#' @param seu A Seurat object with clusters (output of `run_clustering()`).
#' @param config Named list from `load_config()`.
#' @return A list with components:
#'   - `seu`: Seurat object (with SingleR labels added if run).
#'   - `markers`: Data frame of marker genes per cluster.
#'   - `plots`: Named list of ggplot2 objects.
#'   - `summary`: Character string summarising results.
run_annotation <- function(seu, config) {

  set.seed(config$seed)
  species <- config$species
  out_dir <- file.path(config$output$dir, "annotation")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plots <- list()

  # --- Find marker genes per cluster -----------------------------------------
  log_message("Finding marker genes for each cluster (Wilcoxon, positive only)...")
  markers <- Seurat::FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = FALSE
  )

  n_markers <- nrow(markers)
  log_message(sprintf("Found %d marker genes across %d clusters.",
                       n_markers, length(unique(markers$cluster))))

  # Save markers table
  markers_path <- file.path(out_dir, "cluster_markers.csv")
  tryCatch(
    utils::write.csv(markers, markers_path, row.names = FALSE),
    error = function(e) log_message(
      paste0("Could not save markers CSV: ", e$message), "ERROR"
    )
  )

  # Top markers heatmap
  top_markers <- markers |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(avg_log2FC, n = 5)

  plots$heatmap <- tryCatch(
    {
      p <- Seurat::DoHeatmap(seu, features = top_markers$gene) +
        ggplot2::theme(text = ggplot2::element_text(size = 6))
      save_plot(p, file.path(out_dir, "marker_heatmap"), width = 14, height = 10)
      p
    },
    error = function(e) {
      log_message(paste0("Heatmap generation failed: ", e$message), "WARN")
      NULL
    }
  )

  # Dot plot of top 3 markers per cluster
  top3 <- markers |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(avg_log2FC, n = 3)

  plots$dotplot <- tryCatch(
    {
      p <- Seurat::DotPlot(seu, features = unique(top3$gene)) +
        ggplot2::coord_flip() +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      save_plot(p, file.path(out_dir, "marker_dotplot"), width = 12, height = 8)
      p
    },
    error = function(e) {
      log_message(paste0("Dot plot generation failed: ", e$message), "WARN")
      NULL
    }
  )

  # --- SingleR (optional) ----------------------------------------------------
  singler_ran <- FALSE

  if (requireNamespace("SingleR", quietly = TRUE) &&
      requireNamespace("celldex", quietly = TRUE)) {
    log_message("Running SingleR automated annotation...")

    ref <- tryCatch(
      {
        if (tolower(species) == "human") {
          celldex::HumanPrimaryCellAtlasData()
        } else {
          celldex::MouseRNAseqData()
        }
      },
      error = function(e) {
        log_message(
          paste0("Could not load SingleR reference: ", e$message,
                 ". Skipping automated annotation."),
          "WARN"
        )
        NULL
      }
    )

    if (!is.null(ref)) {
      sce <- Seurat::as.SingleCellExperiment(seu)
      pred <- SingleR::SingleR(
        test = sce,
        ref = ref,
        labels = ref$label.main
      )
      seu$singler_label <- pred$labels
      singler_ran <- TRUE

      log_message(sprintf(
        "SingleR assigned %d unique cell types.",
        length(unique(pred$labels))
      ))

      # UMAP coloured by SingleR labels
      plots$umap_singler <- Seurat::DimPlot(
        seu,
        reduction = "umap",
        group.by = "singler_label",
        label = TRUE,
        label.size = 3,
        repel = TRUE
      ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "UMAP — SingleR Annotations")
      save_plot(plots$umap_singler, file.path(out_dir, "umap_singler"), width = 10, height = 8)
    }
  } else {
    log_message(
      "SingleR or celldex not installed — skipping automated annotation. Install both for reference-based labelling.",
      level = "WARN"
    )
  }

  # --- Save intermediate object -----------------------------------------------
  out_path <- file.path(config$output$dir, "05_seurat_annotated.rds")
  tryCatch(
    saveRDS(seu, out_path),
    error = function(e) log_message(
      paste0("Could not save annotated object: ", e$message), "ERROR"
    )
  )
  log_message(sprintf("Annotated Seurat object saved to %s", out_path))

  # --- Summary ----------------------------------------------------------------
  summary_text <- paste0(
    "=== Annotation Summary ===\n",
    sprintf("Marker genes found: %d\n", n_markers),
    sprintf("Clusters with markers: %d\n", length(unique(markers$cluster))),
    sprintf("Markers saved to: %s\n", markers_path),
    sprintf("SingleR automated annotation: %s\n",
            if (singler_ran) "yes" else "skipped (packages not available)")
  )

  list(seu = seu, markers = markers, plots = plots, summary = summary_text)
}

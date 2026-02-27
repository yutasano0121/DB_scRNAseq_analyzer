# 09_interactome.R — Cell-cell communication analysis
#
# Uses CellChat to infer ligand-receptor interactions between cell types.
# Requires normalised data and cell type annotations (from module 05).

source("R/utils.R")

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

#' Run cell-cell communication analysis with CellChat.
#'
#' @param seu A Seurat object with cell type annotations.
#' @param config Named list from `load_config()`.
#' @param group_by Character. Metadata column with cell type labels.
#'   Defaults to "singler_label" if available, otherwise "seurat_clusters".
#' @return A list with components:
#'   - `cellchat`: The CellChat object.
#'   - `plots`: Named list of ggplot2 / CellChat plot objects.
#'   - `summary`: Character string summarising results.
run_interactome <- function(seu, config, group_by = NULL) {

  set.seed(config$seed)
  species <- config$species
  out_dir <- file.path(config$output$dir, "interactome")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plots <- list()

  # --- Check CellChat availability -------------------------------------------
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    log_message(
      "CellChat package not installed — skipping interactome analysis. Install with: remotes::install_github('sqjin/CellChat')",
      "WARN"
    )
    return(list(
      cellchat = NULL,
      plots = list(),
      summary = "Interactome: skipped (CellChat not installed)."
    ))
  }

  # --- Determine cell type column ---------------------------------------------
  if (is.null(group_by)) {
    if ("singler_label" %in% colnames(seu@meta.data)) {
      group_by <- "singler_label"
    } else {
      group_by <- "seurat_clusters"
    }
  }

  log_message(sprintf("Running CellChat (grouping by '%s')...", group_by))

  # --- Select database -------------------------------------------------------
  cellchat_db <- if (tolower(species) == "human") {
    CellChat::CellChatDB.human
  } else {
    CellChat::CellChatDB.mouse
  }

  # --- Create CellChat object ------------------------------------------------
  cellchat <- tryCatch(
    {
      cc <- CellChat::createCellChat(
        object = seu,
        group.by = group_by,
        assay = "RNA"
      )
      cc@DB <- cellchat_db
      cc
    },
    error = function(e) {
      log_message(paste0("CellChat object creation failed: ", e$message), "ERROR")
      NULL
    }
  )

  if (is.null(cellchat)) {
    return(list(
      cellchat = NULL,
      plots = list(),
      summary = "Interactome: CellChat object creation failed."
    ))
  }

  # --- Run CellChat pipeline -------------------------------------------------
  log_message("Preprocessing CellChat data...")
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)

  log_message("Computing communication probabilities (this may take several minutes)...")
  cellchat <- CellChat::computeCommunProb(cellchat)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)

  log_message("Computing pathway-level communication...")
  cellchat <- CellChat::computeCommunProbPathway(cellchat)

  log_message("Aggregating cell-cell communication network...")
  cellchat <- CellChat::aggregateNet(cellchat)

  # --- Visualisations ---------------------------------------------------------
  n_celltypes <- length(unique(seu@meta.data[[group_by]]))

  # Network heatmap
  tryCatch(
    {
      pdf(file.path(out_dir, "interaction_heatmap.pdf"), width = 10, height = 8)
      CellChat::netVisual_heatmap(cellchat)
      dev.off()

      png(file.path(out_dir, "interaction_heatmap.png"), width = 10, height = 8,
          units = "in", res = 300)
      CellChat::netVisual_heatmap(cellchat)
      dev.off()

      log_message("Interaction heatmap saved.")
    },
    error = function(e) log_message(
      paste0("Heatmap visualisation failed: ", e$message), "WARN"
    )
  )

  # Bubble plot (top pathways)
  tryCatch(
    {
      sources <- seq_len(min(3, n_celltypes))
      targets <- seq_len(min(n_celltypes, 6))

      pdf(file.path(out_dir, "bubble_plot.pdf"), width = 12, height = 8)
      CellChat::netVisual_bubble(
        cellchat,
        sources.use = sources,
        targets.use = targets,
        remove.isolate = TRUE
      )
      dev.off()

      png(file.path(out_dir, "bubble_plot.png"), width = 12, height = 8,
          units = "in", res = 300)
      CellChat::netVisual_bubble(
        cellchat,
        sources.use = sources,
        targets.use = targets,
        remove.isolate = TRUE
      )
      dev.off()

      log_message("Bubble plot saved.")
    },
    error = function(e) log_message(
      paste0("Bubble plot failed: ", e$message), "WARN"
    )
  )

  # Aggregate network circle plot
  tryCatch(
    {
      pdf(file.path(out_dir, "circle_plot.pdf"), width = 8, height = 8)
      CellChat::netVisual_circle(
        cellchat@net$count,
        vertex.weight = table(seu@meta.data[[group_by]]),
        weight.scale = TRUE,
        title.name = "Number of Interactions"
      )
      dev.off()

      png(file.path(out_dir, "circle_plot.png"), width = 8, height = 8,
          units = "in", res = 300)
      CellChat::netVisual_circle(
        cellchat@net$count,
        vertex.weight = table(seu@meta.data[[group_by]]),
        weight.scale = TRUE,
        title.name = "Number of Interactions"
      )
      dev.off()

      log_message("Circle plot saved.")
    },
    error = function(e) log_message(
      paste0("Circle plot failed: ", e$message), "WARN"
    )
  )

  # --- Save CellChat object --------------------------------------------------
  cc_path <- file.path(out_dir, "cellchat_object.rds")
  tryCatch(
    saveRDS(cellchat, cc_path),
    error = function(e) log_message(
      paste0("Could not save CellChat object: ", e$message), "ERROR"
    )
  )
  log_message(sprintf("CellChat object saved to %s", cc_path))

  # --- Summary ----------------------------------------------------------------
  n_interactions <- tryCatch(
    nrow(CellChat::subsetCommunication(cellchat)),
    error = function(e) NA
  )

  n_pathways <- tryCatch(
    length(cellchat@netP$pathways),
    error = function(e) NA
  )

  summary_text <- paste0(
    "=== Interactome Summary (CellChat) ===\n",
    sprintf("Species database: %s\n", species),
    sprintf("Cell types (grouped by '%s'): %d\n", group_by, n_celltypes),
    sprintf("Significant interactions: %s\n",
            if (is.na(n_interactions)) "unknown" else as.character(n_interactions)),
    sprintf("Active signalling pathways: %s\n",
            if (is.na(n_pathways)) "unknown" else as.character(n_pathways)),
    sprintf("CellChat object saved to: %s\n", cc_path)
  )

  list(cellchat = cellchat, plots = plots, summary = summary_text)
}

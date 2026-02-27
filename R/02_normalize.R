# 02_normalize.R — Normalisation of gene expression
#
# Supports two methods (selected via config.yaml):
#   1. LogNormalize — NormalizeData (fast, standard)
#   2. SCTransform  — Variance-stabilising transformation (preferred for
#                     datasets with variable sequencing depth)

source("R/utils.R")

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

#' Normalise a Seurat object.
#'
#' @param seu A Seurat object (output of `run_qc()`).
#' @param config Named list from `load_config()`.
#' @return A list with components:
#'   - `seu`: Normalised Seurat object.
#'   - `summary`: Character string describing what was done.
run_normalize <- function(seu, config) {

  set.seed(config$seed)
  method       <- config$normalization$method
  scale_factor <- config$normalization$scale_factor %||% 10000
  vars_to_regress <- config$normalization$vars_to_regress

  log_message(sprintf("Normalisation method: %s", method))

  if (tolower(method) == "sctransform") {
    # SCTransform replaces NormalizeData + FindVariableFeatures + ScaleData
    log_message("Running SCTransform (this may take a few minutes)...")

    regress <- if (length(vars_to_regress) > 0) vars_to_regress else NULL
    seu <- Seurat::SCTransform(
      seu,
      vars.to.regress = regress,
      verbose = FALSE
    )

    summary_text <- paste0(
      "=== Normalisation Summary ===\n",
      "Method: SCTransform\n",
      sprintf("Variables regressed: %s\n",
              if (is.null(regress)) "none" else paste(regress, collapse = ", ")),
      "Active assay set to: SCT\n",
      "Note: SCTransform replaces NormalizeData + FindVariableFeatures + ScaleData."
    )

  } else if (tolower(method) == "lognormalize") {
    log_message(sprintf("Running LogNormalize (scale.factor = %d)...", scale_factor))

    seu <- Seurat::NormalizeData(
      seu,
      normalization.method = "LogNormalize",
      scale.factor = scale_factor,
      verbose = FALSE
    )

    summary_text <- paste0(
      "=== Normalisation Summary ===\n",
      "Method: LogNormalize\n",
      sprintf("Scale factor: %d\n", scale_factor),
      "Active assay: RNA\n",
      "Next step: Feature selection & PCA (module 03)."
    )

  } else {
    stop(sprintf(
      "Unknown normalisation method '%s'.\n  Use 'LogNormalize' or 'SCTransform' in config.yaml > normalization > method.",
      method
    ))
  }

  # --- Save intermediate object -----------------------------------------------
  out_path <- file.path(config$output$dir, "02_seurat_normalized.rds")
  tryCatch(
    saveRDS(seu, out_path),
    error = function(e) log_message(
      paste0("Could not save normalised object: ", e$message), "ERROR"
    )
  )
  log_message(sprintf("Normalised Seurat object saved to %s", out_path))

  list(seu = seu, summary = summary_text)
}

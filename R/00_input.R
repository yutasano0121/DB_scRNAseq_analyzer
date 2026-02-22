# 00_input.R — Load 10X Genomics data into a BPCells-backed Seurat object
#
# Supports two input formats:
#   1. 10X MEX directory (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
#   2. 10X HDF5 (.h5) matrix file
#
# All count matrices are written to BPCells on-disk format before creating the
# Seurat object, ensuring memory-efficient downstream analysis.

source("R/utils.R")

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

#' Load 10X data and create a BPCells-backed Seurat object.
#'
#' Reads one or more samples from 10X matrix output (MEX or HDF5), writes
#' each to a BPCells on-disk directory, creates individual Seurat objects,
#' and merges them if multiple samples are provided.
#'
#' @param config Named list from `load_config()`.
#' @return A Seurat object with BPCells-backed count matrix.
run_load_data <- function(config) {

  input_type  <- config$input$type
  sample_paths <- config$input$paths
  bpcells_dir  <- file.path(config$output$dir, "..", "data", "bpcells")

  # --- Validate inputs -------------------------------------------------------
  if (is.null(sample_paths) || length(sample_paths) == 0) {
    stop(
      "No input paths found in config.yaml.\n",
      "  Add at least one sample under 'input > paths', for example:\n",
      "  input:\n",
      "    paths:\n",
      "      sample1: \"data/raw/sample1/mex/\""
    )
  }

  if (!input_type %in% c("mex", "h5")) {
    stop(sprintf(
      "Unknown input type '%s'. Use 'mex' (MEX folder) or 'h5' (HDF5 file).",
      input_type
    ))
  }

  log_message(sprintf(
    "Loading %d sample(s) as '%s' format.",
    length(sample_paths), input_type
  ))

  # --- Load each sample -------------------------------------------------------
  seu_list <- list()

  for (sample_name in names(sample_paths)) {
    path <- sample_paths[[sample_name]]
    bp_dir <- normalizePath(file.path(bpcells_dir, sample_name), mustWork = FALSE)

    log_message(sprintf("Processing sample '%s': %s", sample_name, path))

    # Validate path exists
    if (!file.exists(path)) {
      stop(sprintf(
        "Input path not found for sample '%s': %s\n  Check config.yaml > input > paths.",
        sample_name, path
      ))
    }

    # Read raw matrix
    raw_mat <- read_10x_matrix(path, input_type, sample_name)

    # Write to BPCells on-disk format
    mat <- write_bpcells_matrix(raw_mat, bp_dir, sample_name)

    # Create Seurat object
    seu <- Seurat::CreateSeuratObject(
      counts   = mat,
      project  = sample_name,
      min.cells    = 3,
      min.features = 200
    )

    log_message(sprintf(
      "  Sample '%s': %d cells, %d genes after initial filter (min.cells=3, min.features=200).",
      sample_name, ncol(seu), nrow(seu)
    ))

    seu_list[[sample_name]] <- seu
  }

  # --- Merge if multiple samples ----------------------------------------------
  if (length(seu_list) == 1) {
    seu_merged <- seu_list[[1]]
  } else {
    log_message("Merging samples...")
    seu_merged <- merge(
      seu_list[[1]],
      y = seu_list[-1],
      add.cell.ids = names(seu_list)
    )
    log_message(sprintf(
      "Merged object: %d cells, %d genes.",
      ncol(seu_merged), nrow(seu_merged)
    ))
  }

  # --- Save intermediate object -----------------------------------------------
  out_path <- file.path(config$output$dir, "00_seurat_raw.rds")
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    saveRDS(seu_merged, out_path),
    error = function(e) log_message(
      paste0("Could not save intermediate object: ", e$message), "ERROR"
    )
  )
  log_message(sprintf("Raw Seurat object saved to %s", out_path))

  seu_merged
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Read a 10X count matrix from MEX directory or HDF5 file.
#'
#' @param path Character. Path to MEX directory or .h5 file.
#' @param input_type Character. "mex" or "h5".
#' @param sample_name Character. Sample label for error messages.
#' @return A sparse matrix (dgCMatrix or similar).
read_10x_matrix <- function(path, input_type, sample_name) {
  tryCatch(
    {
      if (input_type == "mex") {
        validate_mex_directory(path, sample_name)
        Seurat::Read10X(data.dir = path)
      } else {
        mat <- Seurat::Read10X_h5(filename = path)
        # Read10X_h5 may return a list for multi-modal data; take Gene Expression
        if (is.list(mat)) {
          if ("Gene Expression" %in% names(mat)) {
            mat[["Gene Expression"]]
          } else {
            mat[[1]]
          }
        } else {
          mat
        }
      }
    },
    error = function(e) {
      stop(sprintf(
        "Failed to read sample '%s' from '%s'.\n  Error: %s\n  %s",
        sample_name, path, e$message,
        if (input_type == "mex") {
          "Check that the folder contains barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz."
        } else {
          "Check that the .h5 file is a valid 10X HDF5 matrix."
        }
      ))
    }
  )
}

#' Validate that a 10X MEX directory contains the expected files.
#'
#' @param path Character. Path to the MEX directory.
#' @param sample_name Character. Sample label for error messages.
validate_mex_directory <- function(path, sample_name) {
  expected <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  # Also accept older format (genes.tsv.gz instead of features.tsv.gz)
  found <- list.files(path)
  missing <- expected[!expected %in% found]

  # Allow genes.tsv.gz as a substitute for features.tsv.gz (older 10X format)
  if ("features.tsv.gz" %in% missing && "genes.tsv.gz" %in% found) {
    missing <- setdiff(missing, "features.tsv.gz")
  }

  if (length(missing) > 0) {
    stop(sprintf(
      "Sample '%s': MEX directory is missing files: %s\n  Path: %s",
      sample_name, paste(missing, collapse = ", "), path
    ))
  }
}

#' Write a sparse matrix to BPCells on-disk format.
#'
#' If the BPCells directory already exists, it is reopened (not overwritten)
#' to avoid redundant writes on re-runs.
#'
#' @param raw_mat A sparse matrix.
#' @param bp_dir Character. Absolute path for the BPCells directory.
#' @param sample_name Character. Sample label for log messages.
#' @return A BPCells IterableMatrix.
write_bpcells_matrix <- function(raw_mat, bp_dir, sample_name) {
  if (dir.exists(bp_dir)) {
    log_message(sprintf(
      "  BPCells directory already exists for '%s', reopening: %s",
      sample_name, bp_dir
    ))
    return(BPCells::open_matrix_dir(dir = bp_dir))
  }

  log_message(sprintf("  Writing BPCells matrix for '%s' to %s", sample_name, bp_dir))
  dir.create(bp_dir, recursive = TRUE, showWarnings = FALSE)

  tryCatch(
    {
      BPCells::write_matrix_dir(mat = raw_mat, dir = bp_dir)
      BPCells::open_matrix_dir(dir = bp_dir)
    },
    error = function(e) {
      stop(sprintf(
        "Failed to write BPCells matrix for sample '%s'.\n  Directory: %s\n  Error: %s",
        sample_name, bp_dir, e$message
      ))
    }
  )
}

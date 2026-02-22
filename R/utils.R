# utils.R — Shared helpers for the scRNA-seq analysis app
#
# Sources all utilities needed by module scripts and the Shiny app.
# Load at the top of each module with: source("R/utils.R")
#
# Dependencies used directly here (must be installed):
#   yaml, ggplot2

# ---------------------------------------------------------------------------
# 1. Package source registry
# ---------------------------------------------------------------------------

#' Registry of all required packages and their installation sources.
#'
#' Used by `check_and_install_dependencies()`. Edit here when adding packages.
PACKAGE_SOURCES <- list(
  cran = c(
    "Seurat", "SeuratObject",
    "ggplot2", "dplyr",
    "yaml",
    "harmony",
    "bslib", "shiny",
    "promises", "future",
    "ggsci", "viridis",
    "patchwork", "ggrepel",
    "rmarkdown",
    "remotes"
  ),
  bioc = c(
    "SingleR", "celldex",
    "scuttle",
    "DESeq2",
    "BiocParallel"
  ),
  github = list(
    "bnprks/BPCells"               = "BPCells",        # GitHub-only (not on Bioc)
    "satijalab/seurat-wrappers"    = "SeuratWrappers",
    "cole-trapnell-lab/monocle3"   = "monocle3",
    "sqjin/CellChat"               = "CellChat"
  )
)

# ---------------------------------------------------------------------------
# 2. Dependency checking & auto-installation
# ---------------------------------------------------------------------------

#' Check for missing packages and optionally install them.
#'
#' Called at app startup. When `auto_install = FALSE` (default), returns a
#' character vector of missing package names rather than stopping, so the
#' Shiny UI can display a friendly message. When `auto_install = TRUE`,
#' installs all missing packages from the appropriate source.
#'
#' @param auto_install Logical. If TRUE, install missing packages automatically.
#'   Default FALSE (report-only mode for Shiny startup).
#' @return Invisibly returns a character vector of packages that were missing
#'   (empty if all are present).
check_and_install_dependencies <- function(auto_install = FALSE) {
  all_pkgs <- c(
    PACKAGE_SOURCES$cran,
    PACKAGE_SOURCES$bioc,
    unlist(PACKAGE_SOURCES$github)
  )
  all_pkgs <- unique(all_pkgs)

  missing <- all_pkgs[!sapply(all_pkgs, requireNamespace, quietly = TRUE)]

  if (length(missing) == 0) {
    log_message("All required packages are installed.")
    return(invisible(character(0)))
  }

  msg <- paste0("Missing packages: ", paste(missing, collapse = ", "))

  if (!auto_install) {
    log_message(msg, level = "WARN")
    return(invisible(missing))
  }

  log_message(paste0(msg, " — installing now..."), level = "WARN")

  # --- CRAN ---
  cran_missing <- intersect(missing, PACKAGE_SOURCES$cran)
  if (length(cran_missing) > 0) {
    log_message(paste0("Installing from CRAN: ", paste(cran_missing, collapse = ", ")))
    tryCatch(
      utils::install.packages(cran_missing),
      error = function(e) log_message(paste0("CRAN install error: ", e$message), "ERROR")
    )
  }

  # --- Bioconductor ---
  bioc_missing <- intersect(missing, PACKAGE_SOURCES$bioc)
  if (length(bioc_missing) > 0) {
    log_message(paste0("Installing from Bioconductor: ", paste(bioc_missing, collapse = ", ")))
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      utils::install.packages("BiocManager")
    }
    tryCatch(
      BiocManager::install(bioc_missing, ask = FALSE),
      error = function(e) log_message(paste0("Bioc install error: ", e$message), "ERROR")
    )
  }

  # --- GitHub ---
  for (repo in names(PACKAGE_SOURCES$github)) {
    pkg <- PACKAGE_SOURCES$github[[repo]]
    if (pkg %in% missing) {
      log_message(paste0("Installing from GitHub: ", repo))
      if (!requireNamespace("remotes", quietly = TRUE)) {
        utils::install.packages("remotes")
      }
      tryCatch(
        remotes::install_github(repo, upgrade = "never"),
        error = function(e) log_message(
          paste0("GitHub install error (", repo, "): ", e$message), "ERROR"
        )
      )
    }
  }

  log_message("Installation complete. Please restart R and re-run the app.")
  invisible(missing)
}

# ---------------------------------------------------------------------------
# 3. Config loader
# ---------------------------------------------------------------------------

#' Load the YAML configuration file.
#'
#' Returns a named list. Validates that the file exists and is parseable.
#'
#' @param path Character. Path to the config YAML file.
#' @return Named list of configuration values.
load_config <- function(path = "config.yaml") {
  if (!file.exists(path)) {
    stop(sprintf(
      "Configuration file not found: '%s'.\n  Make sure you are running the app from the project root directory.",
      path
    ))
  }
  tryCatch(
    yaml::read_yaml(path),
    error = function(e) stop(sprintf(
      "Could not parse config file '%s'.\n  Check for YAML syntax errors.\n  Details: %s",
      path, e$message
    ))
  )
}

# ---------------------------------------------------------------------------
# 4. Logging
# ---------------------------------------------------------------------------

#' Print a timestamped log message to the console.
#'
#' @param msg Character. The message text.
#' @param level Character. One of "INFO", "WARN", "ERROR". Default "INFO".
log_message <- function(msg, level = "INFO") {
  message(sprintf("[%s] %s %s", level, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

# ---------------------------------------------------------------------------
# 5. Plot saving
# ---------------------------------------------------------------------------

#' Save a ggplot2 object as both PDF and PNG.
#'
#' Creates the output directory if it does not exist.
#'
#' @param plot A ggplot2 object.
#' @param path Character. Full output path (extension stripped and replaced).
#'   E.g. "output/qc/violin.pdf" — will save violin.pdf AND violin.png.
#' @param width Numeric. Plot width in inches. Default 8.
#' @param height Numeric. Plot height in inches. Default 6.
#' @param dpi Integer. PNG resolution. Default 300.
#' @return Invisibly returns the plot object.
save_plot <- function(plot, path, width = 8, height = 6, dpi = 300) {
  base <- tools::file_path_sans_ext(path)
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    ggplot2::ggsave(
      filename = paste0(base, ".pdf"),
      plot = plot, width = width, height = height
    )
    ggplot2::ggsave(
      filename = paste0(base, ".png"),
      plot = plot, width = width, height = height, dpi = dpi
    )
    log_message(paste0("Plot saved: ", base, ".{pdf,png}"))
  }, error = function(e) {
    log_message(paste0("Failed to save plot '", path, "': ", e$message), "ERROR")
  })

  invisible(plot)
}

# ---------------------------------------------------------------------------
# 6. Mito gene prefix helper
# ---------------------------------------------------------------------------

#' Return the mitochondrial gene prefix for a given species.
#'
#' @param species Character. "human" or "mouse".
#' @return Character scalar: "^MT-" (human) or "^mt-" (mouse).
mito_prefix <- function(species) {
  switch(tolower(species),
    human = "^MT-",
    mouse = "^mt-",
    stop(sprintf("Unknown species '%s'. Use 'human' or 'mouse'.", species))
  )
}

# ---------------------------------------------------------------------------
# 7. Human-readable error wrapper
# ---------------------------------------------------------------------------

#' Run an expression and return a user-friendly error message on failure.
#'
#' Intended for use in Shiny `observeEvent` blocks so raw stack traces are
#' never shown to end users.
#'
#' @param expr An expression to evaluate.
#' @param context Character. Short description of the operation for the error message.
#' @return The result of `expr`, or a character string describing the error.
safe_run <- function(expr, context = "operation") {
  tryCatch(
    expr,
    error = function(e) {
      log_message(paste0(context, " failed: ", e$message), "ERROR")
      paste0(
        "Something went wrong during: ", context, "\n\n",
        "Error: ", e$message, "\n\n",
        "What to check:\n",
        "  - Is the input file/folder correct?\n",
        "  - Are all required packages installed? (see Settings tab)\n",
        "  - Does config.yaml have the right settings for this step?"
      )
    }
  )
}

# 10_report.R — Auto-generate summary report
#
# Collects outputs from all completed modules and renders a self-contained
# HTML report using R Markdown (or Quarto if available). The report includes:
#   - config.yaml parameters used
#   - QC metrics and filtering summary
#   - Clustering and annotation results
#   - All saved plots (embedded as PNGs)
#   - Session info for reproducibility

source("R/utils.R")

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

#' Generate a summary report of all analysis steps.
#'
#' @param seu A Seurat object (the final state after all modules).
#' @param config Named list from `load_config()`.
#' @param module_results Named list of results from each module (optional).
#'   Keys: "qc", "features", "clustering", "annotation", "integration",
#'         "deg", "trajectory", "interactome"
#' @return A list with components:
#'   - `report_path`: Path to the rendered report.
#'   - `summary`: Character string confirming the report was generated.
run_report <- function(seu, config, module_results = list()) {

  out_dir <- file.path(config$output$dir, "report")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Collect all PNG files from output subdirectories -----------------------
  plot_files <- list.files(
    config$output$dir,
    pattern = "\\.png$",
    recursive = TRUE,
    full.names = TRUE
  )

  log_message(sprintf("Found %d plot files to include in the report.", length(plot_files)))

  # --- Write the Rmd template ------------------------------------------------
  rmd_path <- file.path(out_dir, "analysis_report.Rmd")
  write_report_template(rmd_path, config, plot_files, seu, module_results)

  # --- Render -----------------------------------------------------------------
  report_path <- file.path(out_dir, "analysis_report.html")

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    log_message(
      "rmarkdown package not installed — cannot render the report. Install with: install.packages('rmarkdown')",
      "WARN"
    )
    return(list(
      report_path = rmd_path,
      summary = paste0(
        "Report template written to: ", rmd_path, "\n",
        "Install rmarkdown to render: install.packages('rmarkdown')"
      )
    ))
  }

  log_message("Rendering report (this may take a moment)...")

  tryCatch(
    {
      rmarkdown::render(
        input = rmd_path,
        output_file = basename(report_path),
        output_dir = out_dir,
        quiet = TRUE
      )
      log_message(sprintf("Report rendered: %s", report_path))
    },
    error = function(e) {
      log_message(paste0("Report rendering failed: ", e$message), "ERROR")
      report_path <- rmd_path
    }
  )

  list(
    report_path = report_path,
    summary = paste0("=== Report ===\nGenerated at: ", report_path)
  )
}

# ---------------------------------------------------------------------------
# Report template writer
# ---------------------------------------------------------------------------

#' Write the R Markdown template for the analysis report.
#'
#' @param rmd_path Path to write the .Rmd file.
#' @param config Config list.
#' @param plot_files Character vector of PNG paths.
#' @param seu Seurat object (for metadata summaries).
#' @param module_results Named list of module outputs.
write_report_template <- function(rmd_path, config, plot_files, seu, module_results) {

  # Build plot sections
  plot_sections <- ""
  if (length(plot_files) > 0) {
    # Group plots by subdirectory
    plot_dirs <- unique(dirname(plot_files))
    for (pd in plot_dirs) {
      section_name <- basename(pd)
      section_plots <- plot_files[dirname(plot_files) == pd]
      plot_sections <- paste0(
        plot_sections,
        sprintf("\n### %s\n\n", tools::toTitleCase(section_name))
      )
      for (pf in section_plots) {
        plot_label <- tools::file_path_sans_ext(basename(pf))
        plot_sections <- paste0(
          plot_sections,
          sprintf("**%s**\n\n", gsub("_", " ", plot_label)),
          sprintf("![%s](%s){width=90%%}\n\n", plot_label, normalizePath(pf, mustWork = FALSE))
        )
      }
    }
  }

  # Build module summary sections
  summary_sections <- ""
  for (mod_name in names(module_results)) {
    res <- module_results[[mod_name]]
    if (!is.null(res$summary)) {
      summary_sections <- paste0(
        summary_sections,
        sprintf("\n### %s\n\n```\n%s\n```\n\n", tools::toTitleCase(mod_name), res$summary)
      )
    }
  }

  # Compose the full Rmd
  rmd_content <- paste0(
    '---
title: "', config$project_name %||% "scRNA-seq Analysis", ' — Analysis Report"
date: "', format(Sys.time(), "%Y-%m-%d %H:%M"), '"
output:
  html_document:
    toc: true
    toc_float: true
    theme: flatly
    self_contained: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
```

## Project Overview

| Parameter | Value |
|-----------|-------|
| Project name | ', config$project_name %||% "—", ' |
| Species | ', config$species %||% "—", ' |
| Normalisation | ', config$normalization$method %||% "—", ' |
| Clustering PCs | ', config$clustering$n_pcs %||% "—", ' |
| Clustering resolution | ', config$clustering$resolution %||% "—", ' |
| Integration | ', if (isTRUE(config$integration$enabled)) config$integration$method else "disabled", ' |
| DEG test | ', config$deg$test %||% "—", ' |
| Random seed | ', config$seed %||% "—", ' |

## Dataset Summary

| Metric | Value |
|--------|-------|
| Total cells | ', ncol(seu), ' |
| Total genes | ', nrow(seu), ' |
| Samples | ', length(unique(seu$orig.ident)), ' |
| Clusters | ', length(levels(Seurat::Idents(seu))), ' |

## Module Summaries

', summary_sections, '

## Figures

', plot_sections, '

## Session Info

```{r session-info}
sessionInfo()
```
')

  writeLines(rmd_content, rmd_path)
  log_message(sprintf("Report template written to %s", rmd_path))
}

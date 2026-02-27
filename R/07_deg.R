# 07_deg.R — Differential expression analysis
#
# Supports:
#   - Wilcoxon rank-sum test (default, fast)
#   - MAST (single-cell aware, optional)
#   - DESeq2 pseudobulk (for multi-sample comparisons with biological replicates)
#
# Produces volcano plots and exports results as CSV.

source("R/utils.R")

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

#' Run differential expression analysis.
#'
#' If `ident_1` and `ident_2` are provided in the config or as parameters,
#' runs a pairwise comparison. Otherwise, runs FindAllMarkers across all
#' clusters/cell types.
#'
#' @param seu A Seurat object with clusters or cell type annotations.
#' @param config Named list from `load_config()`.
#' @param ident_1 Character. Identity class 1 for pairwise comparison (optional).
#' @param ident_2 Character. Identity class 2 for pairwise comparison (optional).
#' @param group_by Character. Metadata column to group by (default: active idents).
#' @return A list with components:
#'   - `deg_results`: Data frame of DEG results.
#'   - `plots`: Named list of ggplot2 objects.
#'   - `summary`: Character string summarising results.
run_deg <- function(seu, config, ident_1 = NULL, ident_2 = NULL, group_by = NULL) {

  set.seed(config$seed)
  deg_config <- config$deg
  test_method <- deg_config$test %||% "wilcox"
  logfc_thresh <- deg_config$logfc_threshold %||% 0.25
  padj_thresh  <- deg_config$padj_threshold %||% 0.05
  min_pct      <- deg_config$min_pct %||% 0.1
  out_dir      <- file.path(config$output$dir, "deg")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plots <- list()

  # Set identity class if group_by is specified
  if (!is.null(group_by)) {
    Seurat::Idents(seu) <- group_by
  }

  # Make sure we're using the RNA assay for DEG
  Seurat::DefaultAssay(seu) <- "RNA"

  # --- Pseudobulk DESeq2 path ------------------------------------------------
  if (tolower(test_method) == "deseq2") {
    return(run_pseudobulk_deg(seu, config, group_by, out_dir))
  }

  # --- Standard Seurat DEG path -----------------------------------------------
  log_message(sprintf(
    "DEG test: %s, logFC threshold: %.2f, min.pct: %.2f",
    test_method, logfc_thresh, min_pct
  ))

  if (!is.null(ident_1)) {
    # Pairwise comparison
    log_message(sprintf("Pairwise comparison: %s vs %s",
                         ident_1, ident_2 %||% "all others"))

    deg_results <- Seurat::FindMarkers(
      seu,
      ident.1 = ident_1,
      ident.2 = ident_2,
      test.use = test_method,
      min.pct = min_pct,
      logfc.threshold = logfc_thresh,
      verbose = FALSE
    )
    deg_results$gene <- rownames(deg_results)
    deg_results$cluster <- ident_1

    label <- sprintf("%s_vs_%s", ident_1, ident_2 %||% "rest")
  } else {
    # All markers
    log_message("Running FindAllMarkers across all clusters...")

    deg_results <- Seurat::FindAllMarkers(
      seu,
      test.use = test_method,
      min.pct = min_pct,
      logfc.threshold = logfc_thresh,
      verbose = FALSE
    )
    label <- "all_clusters"
  }

  n_sig <- sum(deg_results$p_val_adj < padj_thresh, na.rm = TRUE)
  log_message(sprintf("Found %d genes with adj. p < %.2g", n_sig, padj_thresh))

  # Save results
  csv_path <- file.path(out_dir, paste0("deg_", label, ".csv"))
  utils::write.csv(deg_results, csv_path, row.names = FALSE)
  log_message(sprintf("DEG results saved to %s", csv_path))

  # --- Volcano plot -----------------------------------------------------------
  if ("avg_log2FC" %in% colnames(deg_results)) {
    plots$volcano <- make_volcano_plot(deg_results, logfc_thresh, padj_thresh, label)
    save_plot(plots$volcano, file.path(out_dir, paste0("volcano_", label)))
  }

  # --- Summary ----------------------------------------------------------------
  summary_text <- paste0(
    "=== DEG Summary ===\n",
    sprintf("Test: %s\n", test_method),
    sprintf("Comparison: %s\n", label),
    sprintf("Total genes tested: %d\n", nrow(deg_results)),
    sprintf("Significant (adj. p < %.2g): %d\n", padj_thresh, n_sig),
    sprintf("Results saved to: %s\n", csv_path)
  )

  list(deg_results = deg_results, plots = plots, summary = summary_text)
}

# ---------------------------------------------------------------------------
# Volcano plot
# ---------------------------------------------------------------------------

#' Create a volcano plot from DEG results.
#'
#' @param deg Data frame with `avg_log2FC` and `p_val_adj` columns.
#' @param logfc_thresh Numeric. Log2FC significance threshold.
#' @param padj_thresh Numeric. Adjusted p-value threshold.
#' @param title Character. Plot title suffix.
#' @return A ggplot2 object.
make_volcano_plot <- function(deg, logfc_thresh, padj_thresh, title = "") {
  deg$significance <- ifelse(
    deg$p_val_adj < padj_thresh & abs(deg$avg_log2FC) > logfc_thresh,
    "Significant", "Not significant"
  )

  # Label top genes
  top_genes <- deg |>
    dplyr::filter(significance == "Significant") |>
    dplyr::arrange(p_val_adj) |>
    utils::head(10)

  p <- ggplot2::ggplot(
    deg,
    ggplot2::aes(x = avg_log2FC, y = -log10(p_val_adj), colour = significance)
  ) +
    ggplot2::geom_point(alpha = 0.5, size = 0.8) +
    ggplot2::scale_colour_manual(
      values = c("Significant" = "firebrick", "Not significant" = "grey70")
    ) +
    ggplot2::geom_vline(xintercept = c(-logfc_thresh, logfc_thresh),
                         linetype = "dashed", colour = "grey50") +
    ggplot2::geom_hline(yintercept = -log10(padj_thresh),
                         linetype = "dashed", colour = "grey50") +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste("Volcano Plot", title),
      x = "log2 Fold Change",
      y = "-log10(adj. p-value)",
      colour = NULL
    )

  # Add labels for top genes if ggrepel is available
  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(top_genes) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = top_genes,
      ggplot2::aes(label = gene),
      size = 3,
      max.overlaps = 15,
      colour = "black"
    )
  }

  p
}

# ---------------------------------------------------------------------------
# Pseudobulk DESeq2
# ---------------------------------------------------------------------------

#' Run pseudobulk differential expression using DESeq2.
#'
#' Aggregates counts per sample per cell type, then runs DESeq2 for proper
#' statistical treatment of biological replicates.
#'
#' @param seu Seurat object.
#' @param config Config list.
#' @param group_by Metadata column for cell types.
#' @param out_dir Output directory.
#' @return A list with `deg_results`, `plots`, and `summary`.
run_pseudobulk_deg <- function(seu, config, group_by, out_dir) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 package is required for pseudobulk analysis. Install with: BiocManager::install('DESeq2')")
  }

  log_message("Running pseudobulk DEG analysis with DESeq2...")

  group_col <- group_by %||% "seurat_clusters"
  padj_thresh <- config$deg$padj_threshold %||% 0.05
  plots <- list()

  # Aggregate expression per sample per cell type
  pseudo <- Seurat::AggregateExpression(
    seu,
    assays = "RNA",
    group.by = c(group_col, "orig.ident"),
    return.seurat = FALSE
  )

  pseudo_mat <- pseudo$RNA

  # Build sample metadata from column names
  col_info <- strsplit(colnames(pseudo_mat), "_")
  sample_meta <- data.frame(
    celltype = sapply(col_info, `[`, 1),
    sample   = sapply(col_info, `[`, 2),
    row.names = colnames(pseudo_mat)
  )

  # Run DESeq2 per cell type
  all_results <- list()
  celltypes <- unique(sample_meta$celltype)

  for (ct in celltypes) {
    ct_cols <- rownames(sample_meta)[sample_meta$celltype == ct]
    if (length(ct_cols) < 2) next

    tryCatch({
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = round(pseudo_mat[, ct_cols, drop = FALSE]),
        colData   = sample_meta[ct_cols, , drop = FALSE],
        design    = ~1  # simple comparison; extend for condition contrasts
      )
      dds <- DESeq2::DESeq(dds, quiet = TRUE)
      res <- DESeq2::results(dds)
      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)
      res_df$celltype <- ct
      all_results[[ct]] <- res_df
    }, error = function(e) {
      log_message(sprintf("DESeq2 failed for celltype '%s': %s", ct, e$message), "WARN")
    })
  }

  if (length(all_results) == 0) {
    log_message("No pseudobulk results produced. Check that you have multiple samples.", "WARN")
    return(list(
      deg_results = data.frame(),
      plots = list(),
      summary = "Pseudobulk DESeq2: no results (need multiple samples per cell type)."
    ))
  }

  deg_results <- do.call(rbind, all_results)

  csv_path <- file.path(out_dir, "deg_pseudobulk.csv")
  utils::write.csv(deg_results, csv_path, row.names = FALSE)
  log_message(sprintf("Pseudobulk DEG results saved to %s", csv_path))

  summary_text <- paste0(
    "=== Pseudobulk DEG Summary (DESeq2) ===\n",
    sprintf("Cell types analysed: %d\n", length(all_results)),
    sprintf("Total genes: %d\n", nrow(deg_results)),
    sprintf("Results saved to: %s\n", csv_path)
  )

  list(deg_results = deg_results, plots = plots, summary = summary_text)
}

# app.R — scRNA-seq Analysis App (Shiny + bslib Bootstrap 5)
#
# Launch with: shiny::runApp("app.R")
# or from RStudio: click "Run App" button.
#
# All analysis logic lives in R/00_input.R ... R/10_report.R.
# This file only wires the UI and dispatches to those modules.

library(shiny)
library(bslib)
library(future)

source("R/utils.R")
source("R/00_input.R")
source("R/01_qc.R")
source("R/02_normalize.R")
source("R/03_features.R")
source("R/04_clustering.R")
source("R/05_annotation.R")
source("R/06_integration.R")
source("R/07_deg.R")
source("R/08_trajectory.R")
source("R/09_interactome.R")
source("R/10_report.R")

# Allow package installation to run in a background process so the UI stays
# responsive. Must be set before any future::future() call.
future::plan(future::multisession)

# ---------------------------------------------------------------------------
# Startup checks
# ---------------------------------------------------------------------------

# Load config (fail fast with a clear message if config.yaml is missing)
APP_CONFIG <- tryCatch(
  load_config("config.yaml"),
  error = function(e) {
    message("STARTUP ERROR: ", e$message)
    NULL
  }
)

# Check for missing packages (report-only; do not install automatically on startup)
MISSING_PKGS <- check_and_install_dependencies(auto_install = FALSE)

# ---------------------------------------------------------------------------
# Theme
# ---------------------------------------------------------------------------

app_theme <- bslib::bs_theme(
  version    = 5,
  bootswatch = "flatly",
  base_font  = bslib::font_google("Inter"),
  heading_font = bslib::font_google("Inter")
)

# ---------------------------------------------------------------------------
# UI helpers
# ---------------------------------------------------------------------------

#' Build a standard analysis tab panel.
#'
#' @param id        Character. Short identifier used to name input/output IDs.
#' @param title     Character. Tab label shown in the navbar.
#' @param step_num  Character. Module number, e.g. "00".
#' @param subtitle  Character. One-sentence scientific description of this step.
#' @param has_plot  Logical. Whether to show a plot output area.
module_tab <- function(id, title, step_num, subtitle, has_plot = TRUE) {
  bslib::nav_panel(
    title = title,
    bslib::layout_columns(
      col_widths = c(3, 9),

      # ---- Left sidebar: controls ----------------------------------------
      bslib::card(
        bslib::card_header(paste0("Step ", step_num, ": ", title)),
        bslib::card_body(
          tags$p(class = "text-muted small", subtitle),
          tags$hr(),
          shiny::actionButton(
            inputId = paste0("run_", id),
            label   = paste0("Run Step ", step_num),
            class   = "btn-primary w-100"
          ),
          tags$br(), tags$br(),
          shiny::actionButton(
            inputId = paste0("reset_", id),
            label   = "Reset",
            class   = "btn-outline-secondary btn-sm w-100"
          )
        )
      ),

      # ---- Right panel: output -------------------------------------------
      bslib::layout_column_wrap(
        width  = 1,
        heights_equal = "row",

        if (has_plot) {
          bslib::card(
            bslib::card_header("Output"),
            bslib::card_body(
              shiny::plotOutput(paste0("plot_", id), height = "450px")
            )
          )
        },

        bslib::card(
          bslib::card_header("Log"),
          bslib::card_body(
            shiny::verbatimTextOutput(paste0("log_", id))
          )
        )
      )
    )
  )
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ui <- bslib::page_navbar(
  title  = "scRNA-seq Analysis",
  theme  = app_theme,
  id     = "main_nav",

  # -- Module tabs ----------------------------------------------------------

  module_tab(
    id       = "input",
    title    = "Load Data",
    step_num = "00",
    subtitle = paste(
      "Import 10X count matrices (MEX folder or HDF5 file) and create an",
      "on-disk Seurat object backed by BPCells to minimise RAM usage."
    ),
    has_plot = FALSE
  ),

  module_tab(
    id       = "qc",
    title    = "QC & Filtering",
    step_num = "01",
    subtitle = paste(
      "Compute per-cell quality metrics (gene count, UMI count, % mitochondrial reads),",
      "visualise distributions, and remove low-quality cells and potential doublets."
    )
  ),

  module_tab(
    id       = "normalize",
    title    = "Normalization",
    step_num = "02",
    subtitle = paste(
      "Normalise gene expression to correct for sequencing depth differences.",
      "Choose LogNormalize (fast) or SCTransform (recommended for heterogeneous depth)."
    ),
    has_plot = FALSE
  ),

  module_tab(
    id       = "features",
    title    = "Features & PCA",
    step_num = "03",
    subtitle = paste(
      "Identify highly variable genes, scale data, and perform PCA.",
      "An elbow plot helps you choose the optimal number of principal components."
    )
  ),

  module_tab(
    id       = "clustering",
    title    = "Clustering & UMAP",
    step_num = "04",
    subtitle = paste(
      "Build a k-nearest-neighbour graph, cluster cells (Leiden algorithm),",
      "and project them into 2D with UMAP for visual exploration."
    )
  ),

  module_tab(
    id       = "annotation",
    title    = "Annotation",
    step_num = "05",
    subtitle = paste(
      "Identify cluster marker genes and assign cell type labels.",
      "Optionally run SingleR for automated reference-based annotation."
    )
  ),

  module_tab(
    id       = "integration",
    title    = "Integration",
    step_num = "06",
    subtitle = paste(
      "Correct for batch effects when analysing multiple samples.",
      "Supported methods: Harmony (fast), RPCA, CCA, and Sketch integration."
    )
  ),

  module_tab(
    id       = "deg",
    title    = "DEG Analysis",
    step_num = "07",
    subtitle = paste(
      "Find differentially expressed genes between cell types or conditions.",
      "Wilcoxon test (default) or pseudobulk DESeq2 for multi-sample comparisons."
    )
  ),

  module_tab(
    id       = "trajectory",
    title    = "Trajectory",
    step_num = "08",
    subtitle = paste(
      "Infer developmental trajectories and pseudotime ordering with Monocle3."
    )
  ),

  module_tab(
    id       = "interactome",
    title    = "Cell-Cell Communication",
    step_num = "09",
    subtitle = paste(
      "Predict ligand-receptor interactions between cell types using CellChat.",
      "Visualise communication networks as chord diagrams and bubble plots."
    )
  ),

  module_tab(
    id       = "report",
    title    = "Report",
    step_num = "10",
    subtitle = paste(
      "Generate a self-contained HTML or PDF summary report of all analysis steps,",
      "including figures, parameters used, and session info."
    ),
    has_plot = FALSE
  ),

  # -- Settings tab ---------------------------------------------------------

  bslib::nav_spacer(),

  bslib::nav_panel(
    title = shiny::icon("gear"),   # gear icon
    value = "settings",

    bslib::layout_columns(
      col_widths = c(6, 6),

      # Dependency status card
      bslib::card(
        bslib::card_header("Package Status"),
        bslib::card_body(
          shiny::uiOutput("pkg_status"),
          tags$hr(),
          shiny::actionButton(
            "install_pkgs",
            "Install Missing Packages",
            class = "btn-warning"
          ),
          tags$p(
            class = "text-muted small mt-2",
            "Installs missing packages from CRAN, Bioconductor, and GitHub.",
            "Restart R after installation."
          )
        )
      ),

      # Config card
      bslib::card(
        bslib::card_header("Configuration"),
        bslib::card_body(
          shiny::verbatimTextOutput("config_display"),
          tags$hr(),
          shiny::actionButton("reload_config", "Reload config.yaml",
                              class = "btn-outline-secondary btn-sm")
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  # ---- Reactive config (reloaded on demand) --------------------------------
  config <- shiny::reactiveVal(APP_CONFIG)

  shiny::observeEvent(input$reload_config, {
    new_cfg <- tryCatch(
      load_config("config.yaml"),
      error = function(e) {
        shiny::showNotification(
          paste0("Could not reload config: ", e$message),
          type = "error", duration = 8
        )
        NULL
      }
    )
    if (!is.null(new_cfg)) {
      config(new_cfg)
      shiny::showNotification("config.yaml reloaded.", type = "message", duration = 3)
    }
  })

  # ---- Display config -------------------------------------------------------
  output$config_display <- shiny::renderText({
    cfg <- config()
    if (is.null(cfg)) return("config.yaml could not be loaded.")
    paste(capture.output(str(cfg, max.level = 2)), collapse = "\n")
  })

  # ---- Package status -------------------------------------------------------
  output$pkg_status <- shiny::renderUI({
    if (length(MISSING_PKGS) == 0) {
      tags$p(class = "text-success", shiny::icon("check-circle"),
             " All required packages are installed.")
    } else {
      shiny::tagList(
        tags$p(class = "text-danger fw-bold",
               shiny::icon("exclamation-triangle"),
               " Missing packages detected:"),
        tags$ul(lapply(MISSING_PKGS, function(p) tags$li(p))),
        tags$p(class = "text-muted small",
               "Click 'Install Missing Packages' below, then restart R.")
      )
    }
  })

  shiny::observeEvent(input$install_pkgs, {
    shiny::showNotification(
      "Starting installation — this may take several minutes. Check the R console for progress.",
      type = "message", duration = 10
    )
    future::future({
      source("R/utils.R")
      check_and_install_dependencies(auto_install = TRUE)
    }, seed = TRUE)
  })

  # ---- Shared pipeline state -------------------------------------------------
  # The Seurat object flows through the pipeline. Each module updates it.
  seu <- shiny::reactiveVal(NULL)

  # Per-step log text and plots
  step_log <- shiny::reactiveValues(
    input = "", qc = "", normalize = "", features = "",
    clustering = "", annotation = "", integration = "",
    deg = "", trajectory = "", interactome = "", report = ""
  )
  step_plots <- shiny::reactiveValues(
    qc = NULL, features = NULL, clustering = NULL, annotation = NULL,
    integration = NULL, deg = NULL, trajectory = NULL, interactome = NULL
  )

  # Module results (passed to report generator)
  module_results <- shiny::reactiveValues()

  # ---- Step 00: Load Data ---------------------------------------------------
  shiny::observeEvent(input$run_input, {
    step_log$input <- safe_run({
      result <- run_load_data(config())
      seu(result)
      sprintf(
        "[%s] Data loaded successfully.\n  Cells: %d\n  Genes: %d\n  Samples: %d",
        format(Sys.time(), "%H:%M:%S"),
        ncol(result), nrow(result), length(unique(result$orig.ident))
      )
    }, context = "Load Data")
  })

  # ---- Step 01: QC & Filtering -----------------------------------------------
  shiny::observeEvent(input$run_qc, {
    shiny::req(seu())
    step_log$qc <- safe_run({
      result <- run_qc(seu(), config())
      seu(result$seu)
      step_plots$qc <- result$plots$violin_post %||% result$plots$violin
      module_results$qc <- result
      result$summary
    }, context = "QC & Filtering")
  })

  # ---- Step 02: Normalization ------------------------------------------------
  shiny::observeEvent(input$run_normalize, {
    shiny::req(seu())
    step_log$normalize <- safe_run({
      result <- run_normalize(seu(), config())
      seu(result$seu)
      module_results$normalize <- result
      result$summary
    }, context = "Normalization")
  })

  # ---- Step 03: Features & PCA ----------------------------------------------
  shiny::observeEvent(input$run_features, {
    shiny::req(seu())
    step_log$features <- safe_run({
      result <- run_features(seu(), config())
      seu(result$seu)
      step_plots$features <- result$plots$elbow
      module_results$features <- result
      result$summary
    }, context = "Features & PCA")
  })

  # ---- Step 04: Clustering & UMAP -------------------------------------------
  shiny::observeEvent(input$run_clustering, {
    shiny::req(seu())
    step_log$clustering <- safe_run({
      result <- run_clustering(seu(), config())
      seu(result$seu)
      step_plots$clustering <- result$plots$umap
      module_results$clustering <- result
      result$summary
    }, context = "Clustering & UMAP")
  })

  # ---- Step 05: Annotation ---------------------------------------------------
  shiny::observeEvent(input$run_annotation, {
    shiny::req(seu())
    step_log$annotation <- safe_run({
      result <- run_annotation(seu(), config())
      seu(result$seu)
      step_plots$annotation <- result$plots$dotplot %||% result$plots$heatmap
      module_results$annotation <- result
      result$summary
    }, context = "Annotation")
  })

  # ---- Step 06: Integration --------------------------------------------------
  shiny::observeEvent(input$run_integration, {
    shiny::req(seu())
    step_log$integration <- safe_run({
      result <- run_integration(seu(), config())
      seu(result$seu)
      step_plots$integration <- result$plots$umap_post %||% result$plots$comparison
      module_results$integration <- result
      result$summary
    }, context = "Integration")
  })

  # ---- Step 07: DEG Analysis -------------------------------------------------
  shiny::observeEvent(input$run_deg, {
    shiny::req(seu())
    step_log$deg <- safe_run({
      result <- run_deg(seu(), config())
      step_plots$deg <- result$plots$volcano
      module_results$deg <- result
      result$summary
    }, context = "DEG Analysis")
  })

  # ---- Step 08: Trajectory ---------------------------------------------------
  shiny::observeEvent(input$run_trajectory, {
    shiny::req(seu())
    step_log$trajectory <- safe_run({
      result <- run_trajectory(seu(), config())
      seu(result$seu)
      step_plots$trajectory <- result$plots$umap_pseudotime %||%
                               result$plots$monocle3_pseudotime
      module_results$trajectory <- result
      result$summary
    }, context = "Trajectory")
  })

  # ---- Step 09: Cell-Cell Communication -------------------------------------
  shiny::observeEvent(input$run_interactome, {
    shiny::req(seu())
    step_log$interactome <- safe_run({
      result <- run_interactome(seu(), config())
      module_results$interactome <- result
      result$summary
    }, context = "Cell-Cell Communication")
  })

  # ---- Step 10: Report -------------------------------------------------------
  shiny::observeEvent(input$run_report, {
    shiny::req(seu())
    step_log$report <- safe_run({
      mod_res <- shiny::reactiveValuesToList(module_results)
      result <- run_report(seu(), config(), module_results = mod_res)
      result$summary
    }, context = "Report")
  })

  # ---- Reset handlers -------------------------------------------------------
  step_ids <- c("input", "qc", "normalize", "features", "clustering",
                "annotation", "integration", "deg", "trajectory",
                "interactome", "report")

  for (sid in step_ids) {
    local({
      id <- sid
      shiny::observeEvent(input[[paste0("reset_", id)]], {
        step_log[[id]] <- ""
        if (id %in% names(step_plots)) step_plots[[id]] <- NULL
      })
    })
  }

  # ---- Log outputs -----------------------------------------------------------
  for (sid in step_ids) {
    local({
      id <- sid
      output[[paste0("log_", id)]] <- shiny::renderText({
        step_log[[id]]
      })
    })
  }

  # ---- Plot outputs ----------------------------------------------------------
  plot_ids <- c("qc", "features", "clustering", "annotation",
                "integration", "deg", "trajectory", "interactome")

  for (sid in plot_ids) {
    local({
      id <- sid
      output[[paste0("plot_", id)]] <- shiny::renderPlot({
        p <- step_plots[[id]]
        shiny::req(!is.null(p))
        print(p)
      })
    })
  }
}

# ---------------------------------------------------------------------------
# Launch
# ---------------------------------------------------------------------------

shiny::shinyApp(ui = ui, server = server)

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
      "Import CellRanger output (MEX folder or HDF5 file) and create an",
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
      "Infer developmental trajectories and pseudotime ordering with Monocle3.",
      "Optionally estimate RNA velocity (requires CellRanger --include-introns)."
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
    # Run in background so UI stays responsive
    future::future({
      source("R/utils.R")
      check_and_install_dependencies(auto_install = TRUE)
    }, seed = TRUE)
  })

  # ---- Per-step log state ---------------------------------------------------
  step_log <- shiny::reactiveValues(
    input       = "",
    qc          = "",
    normalize   = "",
    features    = "",
    clustering  = "",
    annotation  = "",
    integration = "",
    deg         = "",
    trajectory  = "",
    interactome = "",
    report      = ""
  )

  make_run_handler <- function(id, label) {
    shiny::observeEvent(input[[paste0("run_", id)]], {
      step_log[[id]] <- safe_run(
        {
          log_message(paste0("Step ", label, " triggered (not yet implemented)."))
          paste0(
            format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "),
            "Step '", label, "' is not yet implemented.\n",
            "Module file: R/", label, ".R\n"
          )
        },
        context = paste0("Step ", label)
      )
    })
  }

  make_reset_handler <- function(id) {
    shiny::observeEvent(input[[paste0("reset_", id)]], {
      step_log[[id]] <- ""
    })
  }

  # Register handlers for all module steps
  steps <- list(
    list(id = "input",       label = "00_input"),
    list(id = "qc",          label = "01_qc"),
    list(id = "normalize",   label = "02_normalize"),
    list(id = "features",    label = "03_features"),
    list(id = "clustering",  label = "04_clustering"),
    list(id = "annotation",  label = "05_annotation"),
    list(id = "integration", label = "06_integration"),
    list(id = "deg",         label = "07_deg"),
    list(id = "trajectory",  label = "08_trajectory"),
    list(id = "interactome", label = "09_interactome"),
    list(id = "report",      label = "10_report")
  )

  for (s in steps) {
    local({
      sid   <- s$id
      slabel <- s$label
      make_run_handler(sid, slabel)
      make_reset_handler(sid)
    })
  }

  # ---- Log outputs ----------------------------------------------------------
  for (s in steps) {
    local({
      sid <- s$id
      output[[paste0("log_", sid)]] <- shiny::renderText({
        step_log[[sid]]
      })
    })
  }

  # ---- Plot outputs (stubs) -------------------------------------------------
  plot_steps <- steps[sapply(steps, function(s) s$id != "input" &
                                                 s$id != "normalize" &
                                                 s$id != "report")]
  for (s in plot_steps) {
    local({
      sid <- s$id
      output[[paste0("plot_", sid)]] <- shiny::renderPlot({
        shiny::req(step_log[[sid]] != "")
        # Placeholder: real plots rendered by each module script
        p <- ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = paste0("Plot for step '", sid, "' will appear here."),
                            size = 5, color = "grey50") +
          ggplot2::theme_void()
        print(p)
      })
    })
  }
}

# ---------------------------------------------------------------------------
# Launch
# ---------------------------------------------------------------------------

shiny::shinyApp(ui = ui, server = server)

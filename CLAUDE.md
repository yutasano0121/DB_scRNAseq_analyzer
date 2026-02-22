# CLAUDE.md — scRNA-seq Analysis App

## Session Memory Protocol

Claude Code must follow these steps at the boundaries of every session:

- **Session start:** Read `MEMORY.md` before doing anything else. Use the most recent entry to restore context — what was last worked on, any unfinished TODOs (🔲), and any decisions that override defaults in this file.
- **Session end:** When the user types `/memory`, append a new entry to `MEMORY.md` following the format defined there, then confirm it was saved. This is the only time MEMORY.md should be written to.

If `MEMORY.md` does not exist yet, create it using the template in SKILLS.md before proceeding.

---

## Project Overview

This is an R-based single-cell RNA sequencing (scRNA-seq) analysis application designed for users with little to no coding experience. It provides a guided, end-to-end workflow from raw 10X Genomics output through publication-ready results. The app is built around Seurat, BPCells, and ggplot2, with modular support for trajectory inference (Monocle3, scVelo/RNA velocity) and multi-sample integration.

## Intended Users

Non-coders or researchers with minimal R experience. The app should shield users from implementation details while remaining transparent about what each step does scientifically.

---

## Architecture & Design Principles

### Language & Core Libraries
- **Primary language:** R (exclusively — no user-facing Python)
- **Frontend:** Shiny + bslib (Bootstrap 5 theming); runs locally in RStudio or deployed on a Shiny server
- **Single-cell framework:** Seurat (v5+)
- **On-disk matrix backend:** BPCells (default for all data; avoids RAM overflow on large datasets)
- **Plotting:** ggplot2 (all visualizations; no base R plots)
- **Integration:** Seurat's built-in methods (CCA, RPCA, Harmony); Sketch integration available as optional flag
- **Trajectory inference:** Monocle3 (pure R) and RNA velocity via `velociraptor` + `basilisk` (R interface, automatic Python env management — no user-visible Python)
- **DEG analysis:** Seurat `FindMarkers` / `FindAllMarkers` (Wilcoxon default; support for MAST, DESeq2 pseudobulk)
- **Interactome / cell-cell communication:** CellChat

### Input Formats (10X Genomics)
The app accepts two entry points — users only need one:
1. **CellRanger MEX directory:** folder containing `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` — loaded via `Seurat::Read10X()` then written to BPCells
2. **HDF5 file:** `.h5` file from CellRanger — loaded via `Seurat::Read10X_h5()` with BPCells backend

Read mapping (CellRanger) is a **documented prerequisite**, not an in-app step. See the Prerequisites section below and SKILLS.md §0 for full CellRanger + reference installation instructions.

### User Profiles
Two types of users are expected:

| Profile | Starting point | Extra setup needed |
|---|---|---|
| Starting from FASTQs | Raw sequencing files | CellRanger + transcriptome reference (~10 GB) |
| Starting from CellRanger output | `.h5` or MEX folder | Nothing — just the app |

The majority of users (those whose core facility ran CellRanger) fall into profile 2.

### Prerequisites (Outside the App)

**CellRanger** (only for profile 1 users — starting from FASTQs):
- Install CellRanger from 10X Genomics website (requires free registration)
- Download pre-built transcriptome reference from 10X Genomics: human (GRCh38) or mouse (mm10), ~10 GB each, one-time download per species
- Full step-by-step instructions in SKILLS.md §0

**R package auto-installation:**
- All required R packages are installed automatically on first launch via `check_and_install_dependencies()`
- Uses `BiocManager::install()` for Bioconductor packages and `remotes::install_github()` for GitHub packages
- `basilisk` (used by `velociraptor` for RNA velocity) manages its own internal Python environment — no Python installation required from the user

**Known first-launch friction:**
- `basilisk` environment setup takes ~5–10 minutes on first use
- On HPC systems with restricted internet, `basilisk` may need sysadmin assistance to download its Python environment

### Data Storage
- All raw and processed count matrices stored as BPCells on-disk objects
- Seurat objects use BPCells `dgCMatrix`-like interface transparently
- Intermediate objects saved as `.rds` files in a structured output directory

### Future Modules (not v1)
- Multiomics: CITE-seq (protein/ADT), ATAC-seq (via Signac)
- TCR/BCR clonotype analysis (via scRepertoire)

---

## Analysis Workflow

The pipeline is organized into discrete, numbered modules. Each module saves output before the next begins, allowing restarts.

```
00_input        → Load data, create Seurat object with BPCells backend
01_qc           → QC metrics, filtering (nFeature, nCount, % mito)
02_normalize    → NormalizeData or SCTransform
03_features     → Highly variable features, scaling, PCA
04_clustering   → Nearest-neighbor graph, UMAP/tSNE, clustering (Leiden/Louvain)
05_annotation   → Cell type annotation (manual markers, SingleR optional)
06_integration  → Multi-sample/batch integration (CCA, RPCA, Harmony, Sketch optional)
07_deg          → Differential expression (FindMarkers, pseudobulk)
08_trajectory   → Trajectory inference (Monocle3 and/or RNA velocity)
09_interactome  → Cell-cell communication (CellChat or LIANA)
10_report       → Auto-generate summary report (R Markdown / Quarto)
```

---

## Coding Conventions

- All R scripts use `snake_case` for variables and functions
- Functions must have roxygen2-style comments (even if not building a package)
- Avoid `library()` calls inside functions; use `package::function()` or declare at top of script
- Use `tryCatch` for all file I/O and BPCells operations
- Never hardcode paths; all paths passed as parameters or read from a config file (`config.yaml`)
- Each module script is self-contained and can be run independently
- Log progress to console with timestamps using `message()` or a lightweight logger
- All ggplot2 plots saved as both `.pdf` and `.png` via a shared `save_plot()` helper

## Configuration

A single `config.yaml` (or `config.R`) file controls all parameters:
- Input paths
- Species (human / mouse — affects mitochondrial gene prefix and reference databases)
- QC thresholds (min/max nFeature, nCount, max percent mito)
- Normalization method (LogNormalize vs SCTransform)
- Number of PCs, resolution for clustering
- Integration method and whether Sketch is enabled
- DEG method and thresholds (log2FC, adjusted p-value)
- Trajectory module choice (monocle3, velocity, both)

---

## Testing & Validation

- Use a small public 10X dataset (e.g., PBMC 3k or 10k) as the standard test dataset
- Each module should have a minimal smoke test that runs end-to-end on this dataset
- Plots should be visually inspected as part of QA; save expected outputs for comparison
- Document expected cell counts and cluster counts after QC for the test dataset

---

## Error Handling for Non-Coders

Since target users have minimal coding experience:
- All errors should produce human-readable messages explaining *what went wrong* and *what to check* (e.g., "Could not find barcodes.tsv.gz — make sure you selected the correct CellRanger output folder")
- Avoid exposing raw R stack traces to end users
- Validate all inputs at the start of each module and fail fast with clear messages
- Provide a `check_dependencies()` function that verifies all required packages are installed

---

## Directory Structure

```
project_root/
├── CLAUDE.md               ← This file (architecture, conventions, meta-instructions)
├── SKILLS.md               ← Technical skill reference (R/Seurat patterns)
├── MEMORY.md               ← Living session log (auto-updated via /memory)
├── config.yaml             ← User-editable parameters
├── R/
│   ├── 00_input.R
│   ├── 01_qc.R
│   ├── 02_normalize.R
│   ├── 03_features.R
│   ├── 04_clustering.R
│   ├── 05_annotation.R
│   ├── 06_integration.R
│   ├── 07_deg.R
│   ├── 08_trajectory.R
│   ├── 09_interactome.R
│   ├── 10_report.R
│   └── utils.R             ← Shared helpers (save_plot, logging, etc.)
├── data/
│   ├── raw/                ← User input (CellRanger output or .h5)
│   └── bpcells/            ← On-disk BPCells matrices
├── output/
│   ├── qc/
│   ├── clustering/
│   ├── deg/
│   ├── trajectory/
│   ├── interactome/
│   └── report/
└── tests/
    └── test_pbmc3k.R
```

---

## Key Package Versions to Pin

| Package       | Minimum Version | Source       | Notes                                        |
|---------------|----------------|--------------|----------------------------------------------|
| Seurat        | 5.0.0          | CRAN         | Required for BPCells native support          |
| BPCells       | 0.1.0          | GitHub       | On-disk matrix backend                       |
| SeuratObject  | 5.0.0          | CRAN         | Matches Seurat v5                            |
| SeuratWrappers| latest         | GitHub       | Monocle3 conversion helpers                  |
| bslib         | 0.5.0          | CRAN         | Shiny Bootstrap 5 theming                    |
| monocle3      | 1.3.0          | GitHub/Bioc  | Trajectory inference                         |
| velociraptor  | 1.0.0          | Bioconductor | RNA velocity (R interface to scVelo)         |
| basilisk      | 1.10.0         | Bioconductor | Auto-managed Python env for velociraptor     |
| ggplot2       | 3.4.0          | CRAN         | Plotting                                     |
| harmony       | 1.0.0          | CRAN         | Batch correction (optional)                  |
| SingleR       | 2.0.0          | Bioconductor | Automated cell type annotation               |
| celldex       | 1.0.0          | Bioconductor | Reference datasets for SingleR               |
| CellChat      | 2.0.0          | GitHub       | Cell-cell communication                      |
| DESeq2        | 1.38.0         | Bioconductor | Pseudobulk DEG analysis                      |

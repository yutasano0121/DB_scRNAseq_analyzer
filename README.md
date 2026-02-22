# Democratize Bioinformatics: scRNA-seq Analyzer

A guided, end-to-end single-cell RNA sequencing analysis application built with R and Shiny. Designed for researchers with little or no coding experience, it walks you from precomputed 10X count matrices to publication-ready figures through a point-and-click interface.

---

## Features

- **No coding required** — all steps run through a browser-based UI
- **Memory-efficient** — uses [BPCells](https://github.com/bnprks/BPCells) on-disk matrices; handles datasets that exceed available RAM
- **Full pipeline coverage** — QC → normalization → clustering → annotation → DEG → trajectory → cell-cell communication → report
- **Auto-installer** — missing R packages are detected on startup and can be installed with one click
- **Reproducible** — all parameters live in a single `config.yaml`; every stochastic step uses a fixed seed

---

## Prerequisites

- **R ≥ 4.3** — [Download R](https://cran.r-project.org/)
- **RStudio** (recommended) — [Download RStudio](https://posit.co/download/rstudio-desktop/)

All required R packages are installed automatically on first launch.

---

## Quick Start

```r
# 1. Clone the repository
#    git clone https://github.com/<your-org>/10Xapp.git
#    cd 10Xapp

# 2. Open RStudio and set the working directory to the project root, then:
shiny::runApp("app.R")
```

The app opens in your browser. On the first launch, go to the **Settings** tab to check and install any missing packages, then restart R.

---

## Installation Details

### R packages

All packages are installed automatically. If you prefer to install manually:

```r
# CRAN
install.packages(c("Seurat", "SeuratObject", "ggplot2", "dplyr", "yaml",
                   "harmony", "bslib", "shiny", "promises", "future",
                   "ggsci", "viridis", "patchwork", "ggrepel",
                   "rmarkdown", "remotes"))

# Bioconductor
BiocManager::install(c("SingleR", "celldex", "scuttle", "DESeq2", "BiocParallel"))

# GitHub
remotes::install_github("bnprks/BPCells")
remotes::install_github("satijalab/seurat-wrappers")
remotes::install_github("cole-trapnell-lab/monocle3")
remotes::install_github("sqjin/CellChat")
```

## Workflow

The pipeline is split into numbered modules. Each module saves its output before the next begins, so you can stop and resume at any point.

| Step | Module | What it does |
|------|--------|--------------|
| 00 | Load Data | Import 10X MEX folder or `.h5` file; write BPCells on-disk matrix |
| 01 | QC & Filtering | Compute QC metrics (genes, UMIs, % mito); filter low-quality cells |
| 02 | Normalization | LogNormalize or SCTransform |
| 03 | Features & PCA | Highly variable genes, scaling, PCA, elbow plot |
| 04 | Clustering & UMAP | k-NN graph, Leiden clustering, UMAP/tSNE |
| 05 | Annotation | Marker genes per cluster; optional SingleR automated annotation |
| 06 | Integration | Multi-sample batch correction (Harmony, RPCA, CCA, or Sketch) |
| 07 | DEG Analysis | Differential expression (Wilcoxon or pseudobulk DESeq2); volcano plot |
| 08 | Trajectory | Pseudotime with Monocle3 |
| 09 | Cell-Cell Communication | Ligand-receptor interactions via CellChat |
| 10 | Report | Auto-generated HTML/PDF summary with figures and parameters |

### Input formats

| Format | Path to provide |
|--------|----------------|
| 10X MEX folder | `data/raw/sample1/mex/` |
| 10X HDF5 file | `data/raw/sample1/matrix.h5` |

Set input paths in `config.yaml` before running step 00.

---

## Configuration

All parameters are controlled by `config.yaml` in the project root. Key settings:

```yaml
project_name: "MyProject"
species: "human"          # "human" or "mouse"
seed: 42                  # Global random seed for reproducibility

input:
  type: "mex"             # "mex" or "h5"
  paths:
    sample1: "data/raw/sample1/mex/"

qc:
  min_features: 200
  max_features: 6000
  max_mito: 20            # % mitochondrial reads cutoff

normalization:
  method: "LogNormalize"  # or "SCTransform"

clustering:
  n_pcs: 30
  resolution: 0.5
  algorithm: 4            # 4 = Leiden (preferred), 1 = Louvain

integration:
  enabled: false          # set true for multi-sample datasets
  method: "harmony"
```

See `config.yaml` for the full list of parameters with inline documentation.

---

## Directory Structure

```
10Xapp/
├── app.R               ← Shiny app entry point
├── config.yaml         ← All user-editable parameters
├── R/
│   ├── utils.R         ← Shared helpers (save_plot, logging, installer)
│   ├── 00_input.R
│   ├── 01_qc.R
│   ├── ...
│   └── 10_report.R
├── data/
│   ├── raw/            ← Place input 10X MEX folders or .h5 files here
│   └── bpcells/        ← On-disk BPCells matrices (auto-created)
├── output/             ← All plots and results (auto-created per module)
└── tests/
    └── test_pbmc3k.R   ← Smoke test using the public PBMC 3k dataset
```

---

## Package Versions

| Package | Min. Version | Source |
|---------|-------------|--------|
| Seurat | 5.0.0 | CRAN |
| SeuratObject | 5.0.0 | CRAN |
| BPCells | 0.1.0 | GitHub |
| SeuratWrappers | latest | GitHub |
| bslib | 0.5.0 | CRAN |
| ggplot2 | 3.4.0 | CRAN |
| harmony | 1.0.0 | CRAN |
| monocle3 | 1.3.0 | GitHub |
| SingleR | 2.0.0 | Bioconductor |
| celldex | 1.0.0 | Bioconductor |
| CellChat | 2.0.0 | GitHub |
| DESeq2 | 1.38.0 | Bioconductor |

---

## Troubleshooting

| Problem | Likely cause | Fix |
|---------|-------------|-----|
| App won't start — missing packages | First launch | Go to Settings tab → Install Missing Packages |
| "Could not find barcodes.tsv.gz" | Wrong input folder | Select the MEX folder containing `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` |
| Leiden algorithm not found | `leidenAlg` not installed | App falls back to Louvain automatically; install `leidenAlg` for Leiden |
| BPCells error after moving files | Relative paths stored | Always keep the project folder in the same location; use absolute paths in config |
| UMAP looks different each run | No seed set | Set `seed` in `config.yaml` (top-level key) |
| SCTransform + Harmony fails | Missing `PrepSCTIntegration()` | Handled automatically in module 06; report as a bug if it occurs |

---

## License

GPL-3.0 — see [LICENSE](LICENSE) for details.

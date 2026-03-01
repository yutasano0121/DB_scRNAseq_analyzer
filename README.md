# Democratize Bioinformatics: scRNA-seq Analyzer

A guided, end-to-end single-cell RNA sequencing analysis application built with R and Shiny. Designed for researchers with little or no coding experience, it walks you from precomputed 10X count matrices to publication-ready figures through a point-and-click interface.

---

## Features

- **No coding required** — all steps run through a browser-based UI
- **Memory-efficient** — uses [BPCells](https://github.com/bnprks/BPCells) on-disk matrices; handles datasets that exceed available RAM
- **Full pipeline coverage** — QC → normalization → clustering → annotation → DEG → trajectory → cell-cell communication → report
- **Auto-installer** — missing R packages are installed automatically on first launch, with a confirmation prompt
- **Reproducible** — all parameters live in a single `config.yaml`; every stochastic step uses a fixed seed
- **Desktop launcher** — double-click to open in your browser; no terminal or RStudio required after setup

---

## Prerequisites

- **R ≥ 4.3** — [Download R](https://cran.r-project.org/)

All required R packages are installed automatically on first launch. RStudio is not required.

### System libraries

Some R packages (BPCells, sf, Cairo) include compiled C/C++ code and need system-level libraries to build. Whether these are needed depends on your OS:

| OS | System libraries needed? | How |
|----|--------------------------|-----|
| **Linux** (Ubuntu/Debian, Fedora) | **Yes** — installed automatically by `install_launcher.sh` | `sudo apt` / `sudo dnf` |
| **macOS** | **Sometimes** — Homebrew libraries needed for a few packages | `brew install hdf5 udunits gdal ...` (handled by `install_launcher.sh`) |
| **Windows** | **No** — R packages ship as pre-compiled binaries | Nothing required |

Running `install_launcher.sh` handles this automatically on all platforms.

---

## Quick Start

### Option A — Desktop icon (recommended for non-coders)

Run this **once** from a terminal to create a desktop shortcut:

```bash
bash /path/to/10Xapp/install_launcher.sh
```

After that, double-click **scRNA-seq Analysis** on your Desktop to start the app. The browser opens automatically — no terminal needed.

> **First launch:** a terminal window will appear to install required R packages (10–30 min, one time only). The browser opens automatically when setup is complete.

### Option B — Terminal

```bash
cd /path/to/10Xapp
bash launch.sh
```

### Option C — RStudio

Open the project in RStudio and click **Run App**, or run:

```r
shiny::runApp("app.R")
```

---

## Package Installation

On first launch, the app checks for missing R packages and installs them automatically.

- **Interactive sessions** (RStudio / terminal): a list of missing packages is printed and you are asked `[y/N]` before installation begins.
- **Desktop launcher**: installation runs without a prompt (the terminal window that opens shows live progress).
- **Re-check**: after installation, the **Settings** tab shows any packages that failed to install (e.g. due to network errors) and provides an **Install Missing Packages** button to retry.

---

## Configuration

All parameters are controlled by `config.yaml` in the project root. You can edit it in two ways:

**From within the app (recommended):** open the **Settings** tab → drag-and-drop your `config.yaml` onto the upload zone, or click to browse. Use the **Download current config.yaml** button to get a copy to edit.

**Manually:** open `config.yaml` in any text editor. Key settings:

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

## Directory Structure

```
10Xapp/
├── app.R                  ← Shiny app entry point
├── config.yaml            ← All user-editable parameters
├── launch.sh              ← Desktop launcher script
├── 10Xapp.desktop         ← Linux desktop shortcut definition
├── install_launcher.sh    ← One-time setup: creates the desktop icon
├── launch.log             ← App startup log (auto-created; share when reporting issues)
├── R/
│   ├── utils.R            ← Shared helpers (save_plot, logging, installer)
│   ├── 00_input.R
│   ├── 01_qc.R
│   ├── ...
│   └── 10_report.R
├── data/
│   ├── raw/               ← Place input 10X MEX folders or .h5 files here
│   └── bpcells/           ← On-disk BPCells matrices (auto-created)
├── output/                ← All plots and results (auto-created per module)
└── tests/
    └── test_pbmc3k.R      ← Smoke test using the public PBMC 3k dataset
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
| Desktop icon does nothing | Not marked as trusted | Right-click the icon → **Allow Launching**, or re-run `install_launcher.sh` |
| Browser doesn't open after double-click | Packages still installing | Wait for the setup terminal to finish; browser opens automatically |
| App won't start — missing packages | Install failed | Open Settings tab → **Install Missing Packages**; check `launch.log` for errors |
| "Could not find barcodes.tsv.gz" | Wrong input folder | Select the MEX folder containing `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` |
| Leiden algorithm not found | `leidenAlg` not installed | App falls back to Louvain automatically; install `leidenAlg` for Leiden |
| BPCells error after moving files | Relative paths stored | Keep the project folder in the same location; use absolute paths in config |
| UMAP looks different each run | No seed set | Set `seed` in `config.yaml` (top-level key) |
| SCTransform + Harmony fails | Missing `PrepSCTIntegration()` | Handled automatically in module 06; report as a bug if it occurs |
| App crashed silently | Startup error | Open `launch.log` in the project root and check the last few lines |

---

## License

GPL-3.0 — see [LICENSE](LICENSE) for details.

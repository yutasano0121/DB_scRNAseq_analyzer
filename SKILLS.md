# SKILLS.md — scRNA-seq App Technical Reference

This file documents the key technical patterns, idioms, and gotchas for implementing each module of this app. Consult this before writing any new module.

---

## 0. Prerequisites: Input Matrices

The app starts from precomputed 10X count matrices. Supported formats:
- MEX folder containing `barcodes.tsv.gz`, `features.tsv.gz` (or `genes.tsv.gz`), and `matrix.mtx.gz`
- HDF5 `.h5` count matrix

External read mapping is out of scope for this project.

---

## 1. Loading Data with BPCells Backend

### From 10X MEX directory
```r
# Always use BPCells for on-disk storage
mat <- BPCells::open_matrix_dir(dir = "data/bpcells/sample1")

# If importing fresh from MEX:
raw_mat <- Seurat::Read10X("data/raw/sample1/mex/")
# Write to BPCells on-disk format:
BPCells::write_matrix_dir(mat = raw_mat, dir = "data/bpcells/sample1")
mat <- BPCells::open_matrix_dir("data/bpcells/sample1")
```

### From HDF5 (.h5)
```r
# Read10X_h5 returns a sparse matrix; immediately write to BPCells
raw_mat <- Seurat::Read10X_h5("data/raw/sample.h5")
BPCells::write_matrix_dir(mat = raw_mat$`Gene Expression`, dir = "data/bpcells/sample1")
mat <- BPCells::open_matrix_dir("data/bpcells/sample1")
```

### Creating a Seurat Object with BPCells
```r
seu <- Seurat::CreateSeuratObject(
  counts = mat,       # BPCells IterableMatrix — Seurat v5 handles this natively
  project = "MySample",
  min.cells = 3,
  min.features = 200
)
# Verify the assay uses BPCells:
class(seu[["RNA"]]$counts)  # Should be "BPCells::IterableMatrix" or similar
```

### Multi-sample loading
```r
# Load each sample, then merge — BPCells merge is memory-efficient
samples <- list(
  s1 = BPCells::open_matrix_dir("data/bpcells/sample1"),
  s2 = BPCells::open_matrix_dir("data/bpcells/sample2")
)
seu_list <- lapply(names(samples), function(s) {
  Seurat::CreateSeuratObject(counts = samples[[s]], project = s)
})
seu_merged <- merge(seu_list[[1]], y = seu_list[-1], add.cell.ids = names(samples))
```

---

## 2. Quality Control

### Compute QC metrics
```r
# Mito prefix: "MT-" for human, "mt-" for mouse
seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RP[SL]")

# Visualize — always save with save_plot()
p <- Seurat::VlnPlot(seu, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
save_plot(p, "output/qc/vln_qc.pdf")
```

### Filtering
```r
# Thresholds come from config; do NOT hardcode
seu <- subset(seu,
  subset = nFeature_RNA > config$qc$min_features &
           nFeature_RNA < config$qc$max_features &
           percent.mt   < config$qc$max_mito
)
```

### DoubletFinder (optional but recommended)
```r
# Run after initial clustering; requires PCs and resolution from clustering step
# Use doubletFinder_v3() — parameters pN=0.25, pK from bcmvn optimization
```

---

## 3. Normalization

### LogNormalize (default, faster)
```r
seu <- Seurat::NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4)
seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- Seurat::ScaleData(seu)  # Scales only variable features by default; fine for PCA
```

### SCTransform (preferred for datasets with varying sequencing depth)
```r
# SCTransform replaces NormalizeData + FindVariableFeatures + ScaleData
seu <- Seurat::SCTransform(seu, vars.to.regress = "percent.mt", verbose = FALSE)
# After SCTransform, use assay = "SCT" in downstream steps
```

**Important:** When using SCTransform with integration, use `PrepSCTIntegration()` before `FindIntegrationAnchors()`.

---

## 4. Dimensionality Reduction & Clustering

```r
seu <- Seurat::RunPCA(seu, npcs = 50)

# Choose number of PCs: ElbowPlot + JackStraw (or just use 30 as safe default)
p <- Seurat::ElbowPlot(seu, ndims = 50)
save_plot(p, "output/clustering/elbow_plot.pdf")

seu <- Seurat::FindNeighbors(seu, dims = 1:config$clustering$n_pcs)
seu <- Seurat::FindClusters(seu, resolution = config$clustering$resolution,
                             algorithm = 4)  # algorithm=4 is Leiden (preferred)

seu <- Seurat::RunUMAP(seu, dims = 1:config$clustering$n_pcs)

# Visualize
p <- Seurat::DimPlot(seu, reduction = "umap", label = TRUE) + ggplot2::theme_minimal()
save_plot(p, "output/clustering/umap_clusters.pdf")
```

**Note on Leiden:** Requires the `igraph` package with Leiden support or the `leidenAlg` package. Fall back to `algorithm = 1` (Louvain) if unavailable, with a warning.

---

## 5. Integration (Multi-sample / Batch Correction)

### Seurat CCA / RPCA (default)
```r
# Split by sample identity, then integrate
seu_list <- SplitObject(seu_merged, split.by = "orig.ident")
seu_list <- lapply(seu_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
})
anchors <- Seurat::FindIntegrationAnchors(seu_list, dims = 1:30,
                                           reduction = "rpca")  # rpca faster than cca
seu_int  <- Seurat::IntegrateData(anchors, dims = 1:30)
DefaultAssay(seu_int) <- "integrated"
```

### Harmony (fast, recommended for >10 samples)
```r
# harmony::RunHarmony wraps Seurat integration
seu <- harmony::RunHarmony(seu, group.by.vars = "orig.ident")
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30)
```

### Sketch Integration (optional, for very large datasets)
```r
# Seurat v5 sketch integration — subsample representative cells first
seu <- Seurat::SketchData(seu, ncells = 5000, method = "LeverageScore",
                           sketched.assay = "sketch")
DefaultAssay(seu) <- "sketch"
# Run integration on sketch, then project to full dataset
```

---

## 6. Cell Type Annotation

### Marker-based (manual)
```r
# Find markers for each cluster
markers <- Seurat::FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25,
                                   logfc.threshold = 0.25)
# Visualize top markers
top_markers <- markers |> dplyr::group_by(cluster) |> dplyr::slice_max(avg_log2FC, n = 5)
p <- Seurat::DoHeatmap(seu, features = top_markers$gene) + ggplot2::theme(text = ggplot2::element_text(size = 6))
save_plot(p, "output/annotation/marker_heatmap.pdf", width = 14, height = 10)
```

### SingleR (automated reference-based)
```r
# Requires SingleR and celldex packages
ref <- celldex::HumanPrimaryCellAtlasData()  # or MonacoImmuneData() etc.
sce <- as.SingleCellExperiment(seu)
pred <- SingleR::SingleR(test = sce, ref = ref, labels = ref$label.main)
seu$singler_label <- pred$labels
```

---

## 7. Differential Expression

### Standard (Wilcoxon, default)
```r
Idents(seu) <- "celltype"  # or "seurat_clusters"
deg <- Seurat::FindMarkers(seu,
  ident.1 = "CD4 T",
  ident.2 = "CD8 T",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25
)
```

### Pseudobulk (DESeq2, for multi-sample comparisons)
```r
# Aggregate counts per sample per cell type first
pseudo <- Seurat::AggregateExpression(seu, assays = "RNA",
                                        group.by = c("celltype", "orig.ident"),
                                        return.seurat = FALSE)
# Then run DESeq2 on the pseudobulk matrix
# This is the statistically correct approach when comparing conditions
```

### Volcano plot (ggplot2)
```r
deg$label <- rownames(deg)
p <- ggplot2::ggplot(deg, ggplot2::aes(x = avg_log2FC, y = -log10(p_val_adj),
                                         color = p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)) +
  ggplot2::geom_point(alpha = 0.5, size = 0.8) +
  ggplot2::scale_color_manual(values = c("grey70", "firebrick")) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "DEG Volcano Plot", x = "log2 Fold Change", y = "-log10 adj. p-value")
save_plot(p, "output/deg/volcano.pdf")
```

---

## 8. Trajectory Inference

### Monocle3
```r
library(monocle3)

# Convert Seurat → CellDataSet
cds <- SeuratWrappers::as.cell_data_set(seu)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# User must specify root cells/nodes interactively or via a marker gene
# For automated root selection: use earliest-timepoint cells or stem marker-high cells
cds <- order_cells(cds)  # opens interactive plot if in RStudio

p <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
save_plot(p, "output/trajectory/monocle3_pseudotime.pdf")
```

**Note:** `SeuratWrappers::as.cell_data_set()` transfers UMAP coordinates from the Seurat object. Make sure UMAP is already computed before conversion.

---

## 9. Cell-Cell Communication (Interactome)

### CellChat
```r
library(CellChat)

# Requires normalized data and cell type annotations in Idents(seu)
cellchat <- CellChat::createCellChat(object = seu, group.by = "celltype",
                                      assay = "RNA")
CellChatDB <- CellChat::CellChatDB.human  # or CellChatDB.mouse
cellchat@DB <- CellChatDB
cellchat <- CellChat::subsetData(cellchat)
cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
cellchat <- CellChat::computeCommunProb(cellchat)
cellchat <- CellChat::computeCommunProbPathway(cellchat)
cellchat <- CellChat::aggregateNet(cellchat)

# Visualize
CellChat::netVisual_bubble(cellchat, sources.use = 1:3, targets.use = 4:6)
```

---

## 10. Plotting Conventions (ggplot2)

All plots must use a shared `save_plot()` utility:

```r
# In utils.R
save_plot <- function(plot, path, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = paste0(tools::file_path_sans_ext(path), ".pdf"),
                  plot = plot, width = width, height = height)
  ggplot2::ggsave(filename = paste0(tools::file_path_sans_ext(path), ".png"),
                  plot = plot, width = width, height = height, dpi = dpi)
  invisible(plot)
}
```

**Preferred theme:** `ggplot2::theme_classic()` or `ggplot2::theme_minimal()` — no grey backgrounds in final figures.

**Color palettes:**
- Clusters / categories: `ggsci::pal_d3()` or custom discrete palettes (avoid default ggplot2 hues)
- Continuous (e.g., expression, pseudotime): `viridis::scale_color_viridis_c(option = "magma")`

---

## 11. Config File Pattern

Use `yaml::read_yaml()` to load `config.yaml` at the start of each module:

```r
config <- yaml::read_yaml("config.yaml")
```

Example `config.yaml`:
```yaml
project_name: "MyProject"
species: "human"          # human or mouse
input:
  type: "mex"             # mex or h5
  paths:
    sample1: "data/raw/sample1/mex/"
    sample2: "data/raw/sample2/mex/"
qc:
  min_features: 200
  max_features: 6000
  max_mito: 20
normalization:
  method: "LogNormalize"  # LogNormalize or SCTransform
clustering:
  n_pcs: 30
  resolution: 0.5
  algorithm: 4            # 4=Leiden, 1=Louvain
integration:
  method: "harmony"       # harmony, rpca, cca, sketch
  enabled: true
deg:
  test: "wilcox"
  logfc_threshold: 0.25
  padj_threshold: 0.05
trajectory:
  methods: ["monocle3"]
```

---

## 12. Dependency Checking

Include this function in `utils.R` and call at app startup:

```r
# Package source map — used by auto-installer
PACKAGE_SOURCES <- list(
  cran  = c("Seurat", "SeuratObject", "ggplot2", "yaml", "dplyr", "harmony",
             "bslib", "shiny", "promises", "future", "ggsci", "viridis"),
  bioc  = c("BPCells", "SingleR", "celldex", "scuttle", "DESeq2", "BiocParallel"),
  github = list(
    "satijalab/seurat-wrappers"    = "SeuratWrappers",
    "cole-trapnell-lab/monocle3"   = "monocle3",
    "sqjin/CellChat"               = "CellChat"
  )
)

check_and_install_dependencies <- function(auto_install = FALSE) {
  all_pkgs <- c(PACKAGE_SOURCES$cran, PACKAGE_SOURCES$bioc,
                unlist(PACKAGE_SOURCES$github))
  missing <- all_pkgs[!sapply(all_pkgs, requireNamespace, quietly = TRUE)]

  if (length(missing) == 0) {
    message("All required packages found.")
    return(invisible(TRUE))
  }

  msg <- paste0("Missing packages: ", paste(missing, collapse = ", "))
  if (!auto_install) {
    stop(paste0(msg, "\nRe-run with auto_install = TRUE to install automatically."))
  }

  message(msg, "\nInstalling...")

  # CRAN
  cran_missing <- intersect(missing, PACKAGE_SOURCES$cran)
  if (length(cran_missing)) install.packages(cran_missing)

  # Bioconductor
  bioc_missing <- intersect(missing, PACKAGE_SOURCES$bioc)
  if (length(bioc_missing)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(bioc_missing)
  }

  # GitHub
  for (repo in names(PACKAGE_SOURCES$github)) {
    pkg <- PACKAGE_SOURCES$github[[repo]]
    if (pkg %in% missing) {
      if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
      remotes::install_github(repo)
    }
  }

  message("Installation complete. Please restart R and re-run the app.")
}
```

---

## 13. Common Pitfalls & Gotchas

| Issue | Cause | Fix |
|---|---|---|
| BPCells matrix becomes invalid after subset | BPCells paths are relative | Always use absolute paths for BPCells dirs; call `normalizePath()` |
| SCTransform fails with Harmony | Need `PrepSCTIntegration()` | Run `PrepSCTIntegration()` on list before `FindIntegrationAnchors()` |
| Leiden algorithm not found | `igraph` or `leidenAlg` not installed | Check for `leidenAlg` package; fall back to Louvain with warning |
| Monocle3 root node error | Root not set before `order_cells()` | Provide root via `root_principal_points` or interactive selection |
| Memory spike during `ScaleData` | Scaling all genes | Only scale variable features (default) or use `SCTransform` |
| UMAP non-reproducible | No seed set | Always set `set.seed()` before any stochastic step; document seed in config |

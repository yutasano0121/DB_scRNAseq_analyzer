# SKILLS.md — scRNA-seq App Technical Reference

This file documents the key technical patterns, idioms, and gotchas for implementing each module of this app. Consult this before writing any new module.

---

## 0. Prerequisites: CellRanger Installation & Transcriptome Reference

This section contains the instructions to provide to users who are **starting from raw FASTQ files**. Users who already have a CellRanger output folder or `.h5` file can skip this entirely.

### 0.1 Install CellRanger

CellRanger runs on **Linux only** (macOS and Windows are not supported). Most users will run this on a lab server or HPC cluster.

1. Register for a free account at [10x Genomics support](https://www.10xgenomics.com/support/software/cell-ranger)
2. Download the latest CellRanger tarball (e.g. `cellranger-8.x.x.tar.gz`)
3. Extract and add to PATH:
```bash
tar -xzvf cellranger-8.x.x.tar.gz
export PATH=/path/to/cellranger-8.x.x:$PATH
# Add the export line to ~/.bashrc to make it permanent
```
4. Verify installation:
```bash
cellranger --version
```

**HPC users:** CellRanger is often already available as a module. Check with:
```bash
module avail cellranger
module load cellranger/8.x.x
```

### 0.2 Download Transcriptome Reference

10X Genomics provides pre-built references. Download once per species (~10 GB each).

**Human (GRCh38):**
```bash
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz
```

**Mouse (mm10):**
```bash
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
tar -xzvf refdata-gex-mm10-2020-A.tar.gz
```

Store the extracted folder somewhere permanent (e.g. `/refs/GRCh38/`). You will point CellRanger to this path.

**Building a custom reference** (e.g. for non-human/mouse species or custom genomes):
```bash
cellranger mkref \
  --genome=custom_genome \
  --fasta=genome.fa \
  --genes=genes.gtf
```

### 0.3 Run CellRanger count

Basic usage for a single sample:
```bash
cellranger count \
  --id=sample1 \                          # output folder name
  --transcriptome=/refs/GRCh38/ \         # path to reference
  --fastqs=/path/to/fastqs/ \             # folder with FASTQ files
  --sample=Sample_Name \                  # sample name prefix in FASTQ filenames
  --localcores=16 \                       # CPUs to use
  --localmem=64                            # RAM in GB
```

Output will be in `sample1/outs/`:
- `filtered_feature_bc_matrix/` — MEX format (use this as app input)
- `filtered_feature_bc_matrix.h5` — HDF5 format (alternative app input)
- `possorted_genome_bam.bam` — BAM file (needed if running RNA velocity)

**For RNA velocity:** CellRanger must be run with the `--include-introns` flag to capture unspliced reads:
```bash
cellranger count \
  --id=sample1 \
  --transcriptome=/refs/GRCh38/ \
  --fastqs=/path/to/fastqs/ \
  --sample=Sample_Name \
  --include-introns true \
  --localcores=16 \
  --localmem=64
```

### 0.4 Expected CellRanger Output Structure

```
sample1/outs/
├── filtered_feature_bc_matrix/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── filtered_feature_bc_matrix.h5
├── possorted_genome_bam.bam        ← needed for RNA velocity
├── possorted_genome_bam.bam.bai
├── metrics_summary.csv             ← QC summary — check this first
└── web_summary.html                ← interactive QC report
```

Always check `metrics_summary.csv` or `web_summary.html` before proceeding. Key metrics to verify:
- **Estimated Number of Cells** — should match expected cell count
- **Median Genes per Cell** — typically 1,000–5,000 for good quality
- **Reads Mapped Confidently to Transcriptome** — should be >70%

---

## 1. Loading Data with BPCells Backend

### From CellRanger MEX directory
```r
# Always use BPCells for on-disk storage
mat <- BPCells::open_matrix_dir(dir = "data/bpcells/sample1")

# If importing fresh from CellRanger:
raw_mat <- Seurat::Read10X("data/raw/sample1/filtered_feature_bc_matrix/")
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

### RNA Velocity (velociraptor + basilisk — pure R, no user-visible Python)

`velociraptor` provides an R interface to scVelo via `basilisk`, which manages a self-contained Python environment automatically. Users never install Python or conda.

**First-time setup note:** `basilisk` downloads its Python environment on first use (~5–10 min). On HPC with restricted internet, this may need sysadmin assistance.

**Prerequisite:** CellRanger must have been run with `--include-introns true` (see SKILLS.md §0.3) to capture spliced and unspliced counts.

```r
library(velociraptor)
library(scuttle)

# Convert Seurat to SingleCellExperiment (velociraptor works with SCE)
sce <- Seurat::as.SingleCellExperiment(seu)

# velociraptor expects separate spliced/unspliced assays
# These come from the CellRanger output when --include-introns was used
# Read them in from the CellRanger h5 or loom file:
#   spliced   -> assay(sce, "spliced")
#   unspliced -> assay(sce, "unspliced")

# Run RNA velocity (calls scVelo internally via basilisk)
velo_out <- velociraptor::scvelo(
  x = sce,
  assay.X = "counts",      # use raw counts
  use.dimred = "UMAP"      # use existing UMAP from Seurat
)

# Extract velocity embeddings for plotting
velo_embed <- velociraptor::embedVelocity(
  x = reducedDim(velo_out, "UMAP"),
  ve = velo_out
)

# Plot with ggplot2
velo_df <- as.data.frame(reducedDim(velo_out, "UMAP"))
velo_df$dx <- velo_embed[,1]
velo_df$dy <- velo_embed[,2]
velo_df$celltype <- seu$celltype

p <- ggplot2::ggplot(velo_df, ggplot2::aes(x = V1, y = V2, color = celltype)) +
  ggplot2::geom_point(size = 0.5, alpha = 0.5) +
  ggplot2::geom_segment(ggplot2::aes(xend = V1 + dx * 0.3, yend = V2 + dy * 0.3),
                         arrow = ggplot2::arrow(length = ggplot2::unit(0.05, "inches")),
                         alpha = 0.3) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "RNA Velocity", x = "UMAP 1", y = "UMAP 2")
save_plot(p, "output/trajectory/velocity_umap.pdf", width = 8, height = 6)
```

**Note on spliced/unspliced assays:** When `--include-introns true` is used, CellRanger outputs a single matrix that includes intronic reads. Splitting into spliced/unspliced requires running `velocyto run10x` on the CellRanger BAM file, or using `alevin-fry` / `STARsolo --soloFeatures Velocyto`. Document whichever approach is supported. `velociraptor` can also accept a pre-built loom file:

```r
# If user has a .loom file from velocyto:
loom <- scuttle::readVelocytoLoom("sample1.loom")
# Merge with Seurat metadata by barcode before passing to velociraptor
```

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
  type: "cellranger"      # cellranger or h5
  paths:
    sample1: "data/raw/sample1/filtered_feature_bc_matrix/"
    sample2: "data/raw/sample2/filtered_feature_bc_matrix/"
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
  methods: ["monocle3"]   # monocle3, velocity, or both
```

---

## 12. Dependency Checking

Include this function in `utils.R` and call at app startup:

```r
# Package source map — used by auto-installer
PACKAGE_SOURCES <- list(
  cran  = c("Seurat", "SeuratObject", "ggplot2", "yaml", "dplyr", "harmony",
             "bslib", "shiny", "promises", "future", "ggsci", "viridis"),
  bioc  = c("BPCells", "SingleR", "celldex", "velociraptor", "basilisk",
             "scuttle", "DESeq2", "BiocParallel"),
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
| velociraptor fails on first run | basilisk env not built yet | Expected; takes 5-10 min; show progress message to user |
| basilisk fails on HPC | Restricted internet access | Ask sysadmin to pre-cache basilisk env; see basilisk docs |
| Spliced/unspliced assays missing | CellRanger run without `--include-introns` | Re-run CellRanger with flag, or run velocyto on BAM |
| Memory spike during `ScaleData` | Scaling all genes | Only scale variable features (default) or use `SCTransform` |
| UMAP non-reproducible | No seed set | Always set `set.seed()` before any stochastic step; document seed in config |

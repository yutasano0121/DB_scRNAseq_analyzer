# MEMORY.md — Session Log

## Format
Each entry follows this structure:

```
## YYYY-MM-DD
**Last worked on:** <brief description>
**Status:** <what was completed or left mid-flight>
**Decisions:** <any choices made that override CLAUDE.md defaults>
**TODOs:**
- 🔲 <unfinished item>
- ✅ <completed item>
```

---

## 2026-02-20

**Last worked on:** Initial project skeleton — greenfield setup

**Status:** Completed first-session bootstrap. All foundational files created; no analysis modules implemented yet.

**Decisions:**
- Single-file `app.R` (ui + server together) chosen over `ui.R`/`server.R` split for simplicity at this stage
- bslib theme: `bootswatch = "flatly"` with Inter font
- Each Shiny tab is a stub with Run button + log output + plot placeholder

**TODOs:**
- 🔲 Implement `R/00_input.R` — Load 10X data (MEX or HDF5) into BPCells-backed Seurat object
- 🔲 Implement `R/01_qc.R` — QC metrics, violin plots, cell filtering
- 🔲 Implement `R/02_normalize.R` — LogNormalize or SCTransform
- 🔲 Implement `R/03_features.R` — HVG selection, scaling, PCA, ElbowPlot
- 🔲 Implement `R/04_clustering.R` — FindNeighbors, FindClusters (Leiden), UMAP
- 🔲 Implement `R/05_annotation.R` — FindAllMarkers, optional SingleR
- 🔲 Implement `R/06_integration.R` — Harmony / CCA / RPCA / Sketch
- 🔲 Implement `R/07_deg.R` — FindMarkers, pseudobulk DESeq2, volcano plot
- 🔲 Implement `R/08_trajectory.R` — Monocle3 trajectory and pseudotime
- 🔲 Implement `R/09_interactome.R` — CellChat pipeline
- 🔲 Implement `R/10_report.R` — R Markdown / Quarto summary report
- 🔲 Wire each Shiny tab to its module (replace stubs with real logic)
- 🔲 Add `tests/test_pbmc3k.R` smoke test
- 🔲 Add DoubletFinder optional step in QC tab

---

## 2026-02-22

**Last worked on:** Removed RNA velocity analysis and removed CellRanger/read-mapping assumptions from the project.

**Status:** Completed code, config, UI text, and documentation updates so the app now supports Monocle3-only trajectory and starts from precomputed 10X matrices (`mex`/`h5`) without CellRanger workflow guidance.

**Decisions:**
- RNA velocity was fully removed from runtime logic and docs (no `velociraptor`/`basilisk` references remain).
- Input mode naming was standardized to `mex` and `h5` (replacing `cellranger` + `h5` wording).
- CellRanger setup/read-mapping instructions were removed; project scope now begins at matrix input.
- These edits were implemented by GPT Codex.

**TODOs:**
- ✅ Remove RNA velocity execution path from `R/08_trajectory.R`
- ✅ Remove RNA velocity dependencies from `R/utils.R` and docs
- ✅ Update UI/config/docs for Monocle3-only trajectory
- ✅ Remove CellRanger/read-mapping references across `R/00_input.R`, `config.yaml`, `app.R`, `README.md`, `CLAUDE.md`, and `SKILLS.md`

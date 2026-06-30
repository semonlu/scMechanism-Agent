# Seurat Pipeline Rules

## Runtime

- R 4.3 or newer.
- Core packages: `Seurat`, `SeuratObject`, `Matrix`, `ggplot2`, `patchwork`, `dplyr`, `data.table`.
- Optional packages: `harmony`, `clustree`, `DoubletFinder`, `scDblFinder`, `SingleR`, `SCINA`, `CellChat`, `monocle3`, `clusterProfiler`, `copykat`, `infercnv`, `hdWGCNA`.

## Primary Script

Use the course-adapted Seurat script:

```bash
python scripts/render_template.py --template scripts/course_adapted/01_seurat_v5_core_pipeline.R --out run/01_seurat_v5_core_pipeline.R --define INPUT_PATH=/data --define OUTPUT_DIR=analysis/run1
```

For marker enrichment after a processed Seurat object exists:

```bash
python scripts/render_template.py --template scripts/course_adapted/02_marker_enrichment_from_seurat.R --out run/02_marker_enrichment_from_seurat.R --define SEURAT_RDS=analysis/run1/objects/processed_seurat.rds --define OUTPUT_DIR=analysis/markers
```

## Course-Derived Logic

- Import: `Read10X()`, `Read10X_h5()`, `CreateSeuratObject()`.
- QC: `PercentageFeatureSet()`, `VlnPlot()`, `FeatureScatter()`.
- Normalization: `NormalizeData()`, `FindVariableFeatures()`, `ScaleData()`.
- Embeddings: `RunPCA()`, `RunUMAP()`, `RunTSNE()`.
- Batch: `RunHarmony()` only when metadata supports it.
- Clustering: sweep `FindClusters()` resolution and review with `clustree`.
- Markers: `FindAllMarkers()`, optionally `presto`, `COSG`, `starTracer`.

## Guardrails

- Do not keep interactive working-directory selection.
- Do not hard-code `Type`; use a user-provided `sample_id` or `batch` column.
- Do not hard-code DoubletFinder classification column names.
- Keep raw counts/unintegrated assay for DE, CellChat, CNV, and deconvolution.
- Save `sessionInfo()` and analysis parameters.

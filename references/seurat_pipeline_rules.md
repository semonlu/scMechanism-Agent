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

For Seurat V5 course-derived work, keep the following sequence. Annotation and downstream modules must not be placed before clustering and marker review.

| Order | Step | Course source | Skill script or rule |
|---|---|---|---|
| 1 | Import/object construction | `05_read_10x_standard.R`, `08_read_10x_h5.R`, `09_merge_mixed_inputs.R` | `01_seurat_v5_core_pipeline.R` or `00_multi_sample_merge_harmony.R` |
| 2 | Single-cell QC/filtering | `10_quality_control.R` | `singlecell_qc_rules.md`, `01_seurat_v5_core_pipeline.R` |
| 3 | First normalization/PCA/Harmony | `11_normalization_decontx_harmony.R` | `01_seurat_v5_core_pipeline.R`, `00_multi_sample_merge_harmony.R` |
| 4 | Doublet detection | `12_doublet_finder.R`, `13_scdblfinder.R` | propose when raw counts and sample/loading-batch metadata support it |
| 5 | Doublet removal | `12_doublet_finder.R`, `13_scdblfinder.R` | remove predicted doublets per sample/loading batch |
| 6 | Post-doublet re-normalization | `14_post_doublet_normalization.R` | rerun normalization, variable features, scaling, PCA, neighbors, clustering, UMAP |
| 7 | Clustering/resolution review | `15_clustering_resolution.R` | `01_seurat_v5_core_pipeline.R`; review resolution and marker support |
| 8 | Marker detection | `23_marker_detection_methods.R` | `01_seurat_v5_core_pipeline.R`, `02_marker_enrichment_from_seurat.R` |
| 9 | Cell annotation | `16_manual_cell_annotation.R`, `17_singler_annotation.R`, `18_scina_annotation.R`, `21_transferdata_annotation.R`, `22_scpred_annotation.R` | `05_singler_cell_annotation.R` plus manual marker evidence |
| 10 | Downstream proposal | CellChat/Monocle/CNV/enrichment modules | propose scope and wait for user approval before execution |

Core functions include `Read10X()`, `Read10X_h5()`, `CreateSeuratObject()`, `PercentageFeatureSet()`, `VlnPlot()`, `FeatureScatter()`, `NormalizeData()`, `FindVariableFeatures()`, `ScaleData()`, `RunPCA()`, optional `RunHarmony()`, `FindNeighbors()`, `FindClusters()`, `RunUMAP()`, and `FindAllMarkers()`.

## Guardrails

- Do not keep interactive working-directory selection.
- Do not hard-code `Type`; use a user-provided `sample_id` or `batch` column.
- Do not hard-code DoubletFinder classification column names.
- Keep raw counts/unintegrated assay for DE, CellChat, CNV, and deconvolution.
- Save `sessionInfo()` and analysis parameters.
- If doublet detection is skipped, state why in the plan and report.
- If doublets are removed, downstream clustering, annotation, CellChat, and pseudotime must use the post-doublet reprocessed object.

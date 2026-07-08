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
| 1 | Import/object construction | `05_read_10x_standard.R`, `06_read_10x_nonstandard.R`, `08_read_10x_h5.R`, `09_merge_mixed_inputs.R` | `01_seurat_v5_core_pipeline.R` or `00_multi_sample_merge_harmony.R`; use `INPUT_TYPE=10x_nonstandard` for `count_matrix_*` folders |
| 2 | Single-cell QC/filtering | `10_quality_control.R` | `singlecell_qc_rules.md`, `01_seurat_v5_core_pipeline.R` |
| 3 | First normalization/PCA/Harmony | `11_normalization_decontx_harmony.R` | `01_seurat_v5_core_pipeline.R`, `00_multi_sample_merge_harmony.R` |
| 4 | Doublet detection | `12_doublet_finder.R`, `13_scdblfinder.R` | propose when raw counts and sample/loading-batch metadata support it |
| 5 | Doublet removal | `12_doublet_finder.R`, `13_scdblfinder.R` | remove predicted doublets per sample/loading batch |
| 6 | Post-doublet re-normalization | `14_post_doublet_normalization.R` | rerun normalization, variable features, scaling, PCA, neighbors, clustering, UMAP |
| 7 | Clustering/resolution review | `15_clustering_resolution.R` | `01_seurat_v5_core_pipeline.R`; review resolution and marker support |
| 8 | Marker detection | `23_marker_detection_methods.R` | `01_seurat_v5_core_pipeline.R`, `02_marker_enrichment_from_seurat.R` |
| 9 | Cell annotation | `16_manual_cell_annotation.R`, `17_singler_annotation.R`, `18_scina_annotation.R`, `21_transferdata_annotation.R`, `22_scpred_annotation.R` | `05_singler_cell_annotation.R` plus manual marker evidence |
| 10 | Downstream proposal | CellChat/Monocle/CNV/enrichment modules | propose scope and wait for user approval before execution |

Core functions include `Read10X()`, `Read10X_h5()`, `Matrix::readMM()`, barcode/gene TSV parsing for non-standard 10x, `CreateSeuratObject()`, `PercentageFeatureSet()`, `VlnPlot()`, `FeatureScatter()`, `NormalizeData()`, `FindVariableFeatures()`, `ScaleData()`, `RunPCA()`, optional `RunHarmony()`, `FindNeighbors()`, `FindClusters()`, `RunUMAP()`, and `FindAllMarkers()`.

## Resolution Review

Do not treat `FindClusters(resolution = 0.4)` or any other single value as a universal answer.

For major cell-class annotation, first sweep at least `0.1`, `0.3`, `0.5`, and `0.8`, then review:

1. Cluster count: for typical 10x data with about 5,000 to 20,000 cells, roughly 10 to 20 clusters is often a reasonable major-lineage range. Fewer than about 8 to 10 may miss cell classes; more than about 25 may split one cell type into artificial states.
2. Marker clarity: top markers should show canonical major-lineage genes. If top markers are dominated by ribosomal, heat-shock, mitochondrial, hemoglobin, stress, or cell-cycle genes, review QC, doublets, and over-clustering before accepting the cluster.
3. Biological plausibility: clusters should match tissue context and prior biology, and should not be mainly driven by sample, batch, doublet score, or QC metrics.

Required outputs:

```text
tables/resolution_sweep.tsv
tables/cluster_marker_audit.tsv
figures/umap_clusters.pdf
```

## Subcluster Analysis

Subcluster analysis must not be implemented by only increasing global resolution.

Workflow:

1. Choose a conservative major-lineage resolution and annotate major cell classes.
2. Subset the target cell class.
3. Rerun `FindVariableFeatures()`, `ScaleData()`, `RunPCA()`, `FindNeighbors()`, `FindClusters()`, and `RunUMAP()` on the subset object.
4. Re-check subset markers, doublets, stress, cell cycle, batch, and sample composition.
5. Save subset parameters and object separately.

## Guardrails

- Do not keep interactive working-directory selection.
- Do not hard-code `Type`; use a user-provided `sample_id` or `batch` column.
- Do not hard-code DoubletFinder classification column names.
- Keep raw counts/unintegrated assay for DE, CellChat, CNV, and deconvolution.
- Save `sessionInfo()` and analysis parameters.
- If doublet detection is skipped, state why in the plan and report.
- If doublets are removed, downstream clustering, annotation, CellChat, and pseudotime must use the post-doublet reprocessed object.

# Course Adapted Scripts

These R scripts are the runnable, parameterized versions of the Seurat V5 course workflow. They are adapted from the raw course scripts in `scripts/course_source/` and are intended to be rendered or copied into a project run directory.

| Script | Course modules adapted | Purpose |
|---|---|---|
| `00_multi_sample_merge_harmony.R` | `07_read_multiple_10x_samples.R`, `09_merge_mixed_inputs.R`, `11_normalization_decontx_harmony.R` | Merge multiple samples, preserve sample/batch metadata, and run Harmony when available. |
| `01_seurat_v5_core_pipeline.R` | `05_read_10x_standard.R`, `06_read_10x_nonstandard.R`, `08_read_10x_h5.R`, QC, normalization, clustering, markers | Build and process a Seurat object from 10x MEX, non-standard 10x MEX, H5, RDS, or CSV input. Use `INPUT_TYPE=10x_nonstandard` for `count_matrix_*` folders. |
| `02_marker_enrichment_from_seurat.R` | marker genes, GO/KEGG | Run cluster markers and optional human/mouse enrichment. |
| `03_cellchat_from_seurat.R` | CellChat | Run ligand-receptor communication from a processed Seurat object. |
| `04_monocle3_from_seurat.R` | Monocle3 | Run trajectory and pseudotime analysis from a processed Seurat object. |
| `05_singler_cell_annotation.R` | `17_singler_annotation.R` | Add SingleR cluster/cell labels to a processed Seurat object. |

The expected Seurat V5 course order is: import/object construction -> QC -> first normalization/PCA/Harmony -> doublet detection/removal when applicable -> post-doublet re-normalization -> clustering/resolution review -> marker detection -> cell annotation -> downstream proposal and approval. Downstream CellChat or Monocle3 should use the annotated, post-doublet reprocessed object when doublets were removed.

The Python scripts in the parent `scripts/` folder are helper controllers for GEO diagnosis, planning, template rendering, result checking, and codebase summaries. They do not replace these R analysis scripts.

`templates/` contains generic reusable templates for non-course workflows. `scripts/course_adapted/` contains course-traced runnable modules and should be preferred for Seurat V5 course-derived analyses.

# Course Adapted Scripts

These R scripts are the runnable, parameterized versions of the Seurat V5 course workflow. They are adapted from the raw course scripts in `scripts/course_source/` and are intended to be rendered or copied into a project run directory.

| Script | Course modules adapted | Purpose |
|---|---|---|
| `00_multi_sample_merge_harmony.R` | `07_read_multiple_10x_samples.R`, `09_merge_mixed_inputs.R`, `11_normalization_decontx_harmony.R` | Merge multiple samples, preserve sample/batch metadata, and run Harmony when available. |
| `01_seurat_v5_core_pipeline.R` | import, QC, normalization, clustering, markers | Build and process a Seurat object from 10x MEX/H5/RDS/CSV input. |
| `02_marker_enrichment_from_seurat.R` | marker genes, GO/KEGG | Run cluster markers and optional human/mouse enrichment. |
| `03_cellchat_from_seurat.R` | CellChat | Run ligand-receptor communication from a processed Seurat object. |
| `04_monocle3_from_seurat.R` | Monocle3 | Run trajectory and pseudotime analysis from a processed Seurat object. |
| `05_singler_cell_annotation.R` | `17_singler_annotation.R` | Add SingleR cluster/cell labels to a processed Seurat object. |

The Python scripts in the parent `scripts/` folder are helper controllers for GEO diagnosis, planning, template rendering, result checking, and codebase summaries. They do not replace these R analysis scripts.

`templates/` contains generic reusable templates for non-course workflows or lightweight platform demonstrations. `scripts/course_adapted/` contains course-traced runnable modules and should be preferred for Seurat V5 course-derived analyses.

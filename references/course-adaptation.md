# Seurat V5 Course Adaptation

Use this reference when adapting the local Chinese Seurat V5 course code.

## Source Role

The course archive is reorganized into English-named reference scripts in `scripts/course_source/`. It is not used verbatim as the runnable project code because many original scripts:

- use interactive working-directory selection;
- assume example file names and local objects such as `scRNA.qs`;
- use fixed metadata columns such as `Type`;
- include example thresholds and very low resolution values;
- load many packages regardless of the selected module;
- may call external APIs or package-installation routines that are not appropriate inside an analysis run.

The reference scripts in `scripts/course_source/` are flattened and renamed to English. Each script header records its original course path, and `source_manifest.csv` provides the full mapping.

## New Skill Behavior

Use the preserved course code through parameterized adapted scripts:

| Course module | Adapted skill script/reference |
|---|---|
| 10x/MEX/H5 import | `scripts/course_adapted/01_seurat_v5_core_pipeline.R` |
| QC and filtering | `scripts/course_adapted/01_seurat_v5_core_pipeline.R`, `references/singlecell_qc_rules.md` |
| Normalize/PCA/Harmony/UMAP | `scripts/course_adapted/01_seurat_v5_core_pipeline.R` |
| Doublet detection | `scripts/course_source/12_doublet_finder.R`, `scripts/course_source/13_scdblfinder.R` |
| Doublet removal and post-doublet reprocessing | `scripts/course_source/14_post_doublet_normalization.R`; rerun the adapted core workflow on retained cells |
| Resolution sweep/clustree | `references/seurat_pipeline_rules.md` |
| Clustering before annotation | `scripts/course_adapted/01_seurat_v5_core_pipeline.R`, informed by `scripts/course_source/15_clustering_resolution.R` |
| Marker genes | `scripts/course_adapted/01_seurat_v5_core_pipeline.R`, `scripts/course_adapted/02_marker_enrichment_from_seurat.R` |
| Cell annotation after marker review | `scripts/course_adapted/05_singler_cell_annotation.R`, informed by `scripts/course_source/16_manual_cell_annotation.R`, `17_singler_annotation.R`, `18_scina_annotation.R`, `21_transferdata_annotation.R`, `22_scpred_annotation.R` |
| GO/KEGG | `scripts/course_adapted/02_marker_enrichment_from_seurat.R` |
| CellChat | `scripts/course_adapted/03_cellchat_from_seurat.R` |
| Monocle3 | `scripts/course_adapted/04_monocle3_from_seurat.R` |
| copykat/inferCNV/hdWGCNA/deconvolution | Document as optional project-specific extensions. |

## Required Adaptations

- Replace interactive working directory selection with explicit paths.
- Replace fixed object names with `INPUT_PATH`, `METADATA_PATH`, `OUTPUT_DIR`.
- Replace hard-coded organism resources with `ORGANISM`.
- Replace fixed sample/batch column names with user-provided metadata columns.
- Save logs, session info, major objects, and method-status notes.
- Keep raw counts for DE, CellChat, CNV, and deconvolution.
- Preserve course order: QC -> first normalization/PCA -> doublet detection/removal when applicable -> post-doublet re-normalization -> clustering/resolution review -> marker detection -> cell annotation -> downstream proposal.
- Do not run CellChat, Monocle3, CNV, or deconvolution before annotation evidence and user-approved downstream scope exist.

## Source Traceability

When generating code from course logic, state which source module inspired it. Example:

```text
Source: scripts/course_source/10_quality_control.R
Original course path: see scripts/course_source/source_manifest.csv
Adapted output: scripts/course_adapted/01_seurat_v5_core_pipeline.R
Deviation: thresholds are placeholders and must be chosen from current data distribution.
```

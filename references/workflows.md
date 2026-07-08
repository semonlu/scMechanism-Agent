# Workflow Overview

This skill follows a GEO-to-Mechanism workflow inspired by separate planning, execution, validation, and interpretation layers.

## Standard Flow

1. Parse clinical question with `agents/01_clinical_question_parser.md`.
2. Diagnose accession/file format with `agents/02_geo_dataset_and_format_diagnosis.md` and `scripts/diagnose_geo_inputs.py`.
3. Build a plan with `agents/03_analysis_plan_generator.md` and `scripts/build_analysis_plan.py`.
4. Check or install the environment with `scripts/env_setup/check_environment.ps1` and `scripts/env_setup/install_environment.ps1`.
5. Render local code from `scripts/course_adapted/` or `templates/` using `scripts/render_template.py`.
6. The user runs local R/Python scripts when tools are available.
7. Check returned result folder using `agents/05_result_quality_checker.md` and `scripts/validate_result_bundle.py`.
8. Interpret and report with agents 06-08 and `scripts/write_analysis_report.py`.
9. Before publishing or uploading an example, run `scripts/validate_full_workflow.py` and fix every error.

## Seurat V5 Ordered Course Flow

When the selected route is Seurat V5 course-derived scRNA-seq analysis, use this sequence in the plan and code:

1. Import/object construction: `05_read_10x_standard.R`, `06_read_10x_nonstandard.R`, `08_read_10x_h5.R`, or `09_merge_mixed_inputs.R` -> `01_seurat_v5_core_pipeline.R` or `00_multi_sample_merge_harmony.R`; non-standard `count_matrix_*` folders use `INPUT_TYPE=10x_nonstandard`.
2. Single-cell QC/filtering: `10_quality_control.R` -> `singlecell_qc_rules.md` and `01_seurat_v5_core_pipeline.R`.
3. First normalization/PCA/Harmony: `11_normalization_decontx_harmony.R` -> `01_seurat_v5_core_pipeline.R` or `00_multi_sample_merge_harmony.R`.
4. Doublet detection/removal when applicable: `12_doublet_finder.R` or `13_scdblfinder.R`; document if skipped.
5. Post-doublet normalization/reprocessing: `14_post_doublet_normalization.R`; rerun normalization, PCA, neighbors, UMAP, and clustering on retained cells.
6. Clustering/resolution review: `15_clustering_resolution.R` -> `01_seurat_v5_core_pipeline.R`.
7. Marker detection: `23_marker_detection_methods.R` -> `01_seurat_v5_core_pipeline.R` and `02_marker_enrichment_from_seurat.R`.
8. Cell annotation: `16_manual_cell_annotation.R`, `17_singler_annotation.R`, `18_scina_annotation.R`, `21_transferdata_annotation.R`, `22_scpred_annotation.R` -> `05_singler_cell_annotation.R` plus manual marker review.
9. Downstream proposal: only after annotation evidence is reviewed, propose CellChat, pseudotime, CNV, enrichment, or deconvolution scope and wait for user approval.

## Start Layer By Input Type

- Clinical question only: parse and ask for missing dataset/design.
- Publication URL/PDF/methods text: extract accessions, sample groups, organism, tissue, platform, available file types, and missing metadata before generating code.
- GEO file list: diagnose format and required metadata.
- 10x MEX/H5: generate Seurat or Scanpy template.
- Non-standard 10x MEX: generate Seurat course-adapted code with `INPUT_TYPE=10x_nonstandard`; extract archives first and keep one sample per directory or one row per sample table.
- h5ad: generate Scanpy template and inspect raw/layers.
- Seurat RDS: generate R audit or Seurat downstream template.
- FASTQ/SRA: generate matrix reconstruction advice; do not run downstream directly.
- Result tables: run quality review and interpretation, not reanalysis.

## Output Discipline

Every workflow response should include:

- data readiness;
- missing metadata;
- recommended template/script;
- expected outputs;
- required audit tables: `workflow_step_audit.tsv`, `resolution_sweep.tsv`, and `cluster_marker_audit.tsv`;
- risks and limitations;
- which claims are observations, statistics, or hypotheses.

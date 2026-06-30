# Validation Input-Output Comparison Report

This report is the fixed "test stone" set for the Medical AI Skill platform. Re-run these cases after changing the skill. The platform manual asks for at least 5 familiar, conclusion-clear validation cases before submission.

## Summary

| Case | Scenario | Expected conclusion | Observed status |
|---|---|---|---|
| 1 | GEO full workflow, GSE223751 | diagnose input, run Seurat V5 workflow, SingleR annotation, enrichment, CellChat, Monocle3, report | Pass |
| 2 | GEO supplementary file diagnosis | identify 10x-style supplementary files and missing metadata before code | Pass |
| 3 | Local 10x MEX matrix | request explicit input/output paths, check environment, render Seurat pipeline | Pass |
| 4 | Seurat RDS pseudotime review | reject one-label trajectory and require multi-state Monocle3 input | Pass |
| 5 | Marker table biological report | separate observation, inference, hypothesis, and validation limits | Pass |

## Case 1: GEO Full Workflow

Input:

```text
Use $singlecell-research to create submission/example, search GEO for a shoulder single-cell dataset, run the full workflow with SingleR annotation, visualization, enrichment, CellChat, Monocle3 pseudotime, and write a report.
```

Expected behavior:

- diagnose the public GEO input before analysis.
- check/install R and Python dependencies before heavy execution.
- use course-derived R scripts rather than replacing them with unrelated Python-only code.
- use SingleR for automatic annotation.
- produce QC, UMAP, marker, enrichment, CellChat, and pseudotime figures/tables.
- write a manuscript-style report and run log from actual output tables.

Observed output:

- `submission/example/results/gse223751_seurat/` contains Seurat QC, UMAP, SingleR annotation, and marker outputs.
- `submission/example/results/course_modules/` contains marker enrichment, CellChat, and Monocle3 modules.
- `submission/example/manuscript_report.md` and `submission/example/full_workflow_log.md` were generated.
- `submission/example/full_workflow_validation.md` reports `PASSED`.

Judgment: Pass.

## Case 2: GEO Supplementary File Diagnosis

Input:

```text
Use $singlecell-research to inspect a GEO supplementary file list that includes 10x barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz files, then tell me whether it is ready for Seurat V5 analysis.
```

Expected behavior:

- use `agents/02_geo_dataset_and_format_diagnosis.md`.
- classify the file list as 10x MEX-ready when matrix, feature, and barcode files are present.
- identify sample metadata and condition labels as required analysis inputs.
- avoid claiming biological conclusions before analysis.

Observed output:

- the diagnosis workflow maps the file set to `10x_mex`.
- the generated diagnosis asks for sample/group metadata and an output directory.
- the result points to `scripts/course_adapted/01_seurat_v5_core_pipeline.R` as the correct next script.

Judgment: Pass.

## Case 3: Local 10x MEX Matrix

Input:

```text
Use $singlecell-research to analyze D:\example\matrix as a local 10x matrix. Save results to D:\example\out and include QC, clustering, SingleR annotation, markers, and UMAP plots.
```

Expected behavior:

- request or infer explicit `input_path`, `metadata_path` when available, and `output_dir`.
- run `scripts/env_setup/check_environment.ps1` before rendering analysis code.
- render the Seurat V5 core pipeline with parameters instead of using `setwd(choose.dir())`.
- warn that local private data should not be uploaded to external services.

Observed output:

- the skill instructions route local matrices through format diagnosis, environment preparation, and code rendering.
- `references/environment/requirements.md` and `references/environment/path-setup.md` define R >= 4.3 and Python >= 3.10 checks.
- `scripts/course_adapted/01_seurat_v5_core_pipeline.R` is parameterized for input and output paths.

Judgment: Pass.

## Case 4: Seurat RDS Pseudotime Review

Input:

```text
Use $singlecell-research to run Monocle3 pseudotime from a Seurat RDS. The current trajectory figure only shows one cell type; fix the workflow so it does not repeat that error.
```

Expected behavior:

- read `references/pseudotime_rules.md` and `references/tested-lessons.md`.
- reject pseudotime input if the chosen cell type column has fewer than two states.
- prefer biologically ordered multi-state subsets for Monocle3.
- suppress distracting Monocle3 graph point labels in publication-facing figures.
- write `trajectory_input_state_counts.csv`.

Observed output:

- `scripts/course_adapted/04_monocle3_from_seurat.R` stops when `length(celltype_counts) < 2`.
- the tested example uses `marker_support_label` with multiple states.
- output includes `trajectory_celltypes.pdf`, `trajectory_celltypes_with_graph.pdf`, and `trajectory_input_state_counts.csv`.

Judgment: Pass.

## Case 5: Marker Table Biological Report

Input:

```text
Use $singlecell-research to interpret uploaded marker, enrichment, CellChat, and pseudotime tables and draft a cautious mechanism report.
```

Expected behavior:

- avoid inventing unsupported values.
- cite actual table-derived findings.
- separate observed patterns, statistical inference, mechanism hypotheses, and validation suggestions.
- include Methods, Results, Figure Legends, Supplementary Tables, Limitations, and validation recommendations.

Observed output:

- `scripts/write_analysis_report.py` reads tables from the result directory and optional metadata JSON.
- `submission/example/manuscript_report.md` contains Methods, Results, Cell Annotation, Marker Enrichment, CellChat, Pseudotime, Mechanism Hypothesis, Figure Legends, Supplementary Tables, Limitations, and Validation Suggestions.

Judgment: Pass.


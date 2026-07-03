# scMechanism-Agent

A Multi-Agent Framework for Mechanistic Discovery in Single-Cell Transcriptomics via Tool-Augmented Scientific Reasoning

---

## Overview

scMechanism-Agent is a multi-agent system designed to transform single-cell transcriptomic data into **structured biological mechanism graphs**.

Unlike traditional scRNA-seq analysis pipelines that focus on descriptive outputs (e.g., differential genes, clusters, enrichment tables), this framework aims to generate **mechanistically grounded biological explanations** by integrating:

- Large Language Model (LLM)-based reasoning
- Bioinformatics tools (GO, KEGG, literature databases)
- Structured multi-step scientific reasoning


## Clinical Pain Point

Clinical researchers often find public scRNA-seq datasets in GEO but cannot easily decide:

- whether files are raw reads, 10x matrices, processed Seurat/h5ad objects, or result tables;
- whether disease/control or sample-level metadata is sufficient;
- which analyses are scientifically defensible;
- how to convert example scripts into reproducible local code;
- how far a result can be interpreted as a disease mechanism hypothesis.

## First-Version Capabilities

- Parse clinical research questions into analyzable single-cell tasks.
- Diagnose GEO supplementary file lists and local file formats.
- Generate analysis plans for processed matrices, 10x data, h5ad/RDS objects, and FASTQ/SRA.
- Render or copy Seurat, Scanpy, annotation, marker, and enrichment code.
- Generate CellChat and Monocle3 code only after annotation review and user-approved downstream scope.
- Check result folders for metadata, marker, annotation evidence, DE, enrichment, CellChat, pseudotime, and figure evidence.
- Draft Methods, Results, figure legends, limitations, and validation suggestions.

## Not Promised

- No guarantee of automatic platform-side reanalysis.
- No guarantee that every GEO dataset is complete or analyzable.
- No automatic FASTQ/SRA reconstruction without user-provided reference, chemistry, and compute choices.
- No clinical diagnosis or treatment recommendation.
- No claim that computational cell communication, pseudotime, or CNV inference is experimental proof.

## Runtime Requirements

- Python 3.10 or newer for the helper scripts and Scanpy templates.
- R 4.3 or newer for Seurat-oriented scripts.
- Package needs depend on the selected script/template; see `references/environment/requirements.md`, `references/environment/path-setup.md`, `references/seurat_pipeline_rules.md`, and `references/scanpy_pipeline_rules.md`.
- Before running Seurat V5 course-derived modules, use `scripts/env_setup/check_environment.ps1`; when tools or packages are missing, use `scripts/env_setup/install_environment.ps1`.

## File Structure

```text
SKILL.md
README.md
CODEBASE_SUMMARY.md
environment.yml
agents/
references/
templates/
examples/
scripts/
  env_setup/
  course_source/
  course_adapted/
```

## Main Script Calls

```bash
python scripts/diagnose_geo_inputs.py --file-list supplementary_files.txt --out-json diagnosis.json --out-md diagnosis.md
python scripts/build_analysis_plan.py --diagnosis-json diagnosis.json --question "lung cancer immune microenvironment" --organism human --comparison "tumor vs normal" --out-md analysis_plan.md
python scripts/render_template.py --template scripts/course_adapted/01_seurat_v5_core_pipeline.R --out run/01_seurat_v5_core_pipeline.R --define INPUT_PATH=/data/GSE --define OUTPUT_DIR=analysis/GSE
python scripts/render_template.py --template scripts/course_adapted/05_singler_cell_annotation.R --out run/05_singler_cell_annotation.R --define INPUT_RDS=analysis/GSE/objects/processed_seurat.rds --define OUTPUT_DIR=analysis/GSE/singler
python scripts/render_template.py --template templates/scanpy_batch_annotation_enrichment_template.py --out run/scanpy_optional_modules.py --define INPUT_H5AD=analysis/GSE/objects/processed.h5ad --define OUTPUT_DIR=analysis/GSE/scanpy_optional
python scripts/validate_result_bundle.py --result-dir analysis/GSE --out-md result_quality_check.md
python scripts/propose_downstream_modules.py --result-dir analysis/GSE --out-md analysis/GSE/downstream_proposal.md
python scripts/write_analysis_report.py --result-dir analysis/GSE --metadata-json analysis/GSE/report_metadata.json --out-md analysis/GSE/manuscript_report.md
python scripts/validate_full_workflow.py --project-root . --example-root analysis/GSE --out-md analysis/GSE/full_workflow_validation.md
```

Environment preparation:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\env_setup\check_environment.ps1 -CondaEnv seuratv5-course-py
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\env_setup\install_environment.ps1 -InstallRPackages -InstallPythonEnv
```

Python/Scanpy users can start from the conda environment file:

```powershell
conda env create -f .\environment.yml
conda activate scmechanism-agent
```

The Python scripts are helpers. The Seurat V5 course-derived R analysis code is in:

- `scripts/course_source/`: English-named, lightly adapted course reference R scripts plus `source_manifest.csv`.
- `scripts/course_adapted/`: runnable, parameterized R scripts adapted from the course modules.

## Course Code Reuse

The Seurat V5 course archive was reorganized into English filenames under `scripts/course_source/` and then adapted into runnable scripts under `scripts/course_adapted/`:

- `scripts/course_adapted/01_seurat_v5_core_pipeline.R`: adapted from `05_read_10x_standard.R`, `08_read_10x_h5.R`, `10_quality_control.R`, `11_normalization_decontx_harmony.R`, and `15_clustering_resolution.R`.
- `scripts/course_adapted/00_multi_sample_merge_harmony.R`: adapted from multi-sample import, merge, and Harmony course modules.
- `scripts/course_adapted/02_marker_enrichment_from_seurat.R`: adapted from `23_marker_detection_methods.R` and `24_go_kegg_enrichment.R`.
- `scripts/course_adapted/03_cellchat_from_seurat.R`: adapted from `27_cellchat_analysis.R`.
- `scripts/course_adapted/04_monocle3_from_seurat.R`: adapted from `26_monocle3_pseudotime.R`.
- `scripts/course_adapted/05_singler_cell_annotation.R`: adapted from `17_singler_annotation.R`.

The adapted scripts are parameterized and do not preserve interactive working-directory selection, fixed object names, Chinese output names, or example-specific thresholds. The source reference scripts also have English filenames and a header describing their original source.

`templates/` is for generic reusable templates and lightweight platform demonstrations. `scripts/course_adapted/` is for course-traced runnable Seurat modules. When both are possible, prefer `scripts/course_adapted/` for Seurat V5 course-derived work and `templates/` for Scanpy/report scaffolds.

CellChat and Monocle3 are approval-gated. First run annotation review and generate `downstream_proposal.md`; then ask the user to approve the cell groups, microenvironment, lineage, and root choice before rendering or running those scripts.

## Workflow Layering

The skill keeps data diagnosis, planning, code generation, execution, validation, and interpretation as separate steps. Heavy computation remains external and auditable; the skill coordinates the steps and produces reproducible scripts and review artifacts.

## Validation

Use `scripts/validate_full_workflow.py` against a local analysis folder when you need to check that generated outputs are complete enough for review. Keep large matrices, generated analysis folders, and private metadata out of the repository.

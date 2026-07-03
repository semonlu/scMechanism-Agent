# scMechanism-Agent

A multi-agent framework for mechanistic discovery in single-cell transcriptomics through tool-augmented scientific reasoning.

---

## Overview

**scMechanism-Agent** is a multi-agent system designed to help researchers transform single-cell transcriptomic data into structured biological mechanism hypotheses.

Unlike traditional scRNA-seq analysis pipelines that mainly produce descriptive outputs, such as differential genes, clusters, enrichment tables, and visualization figures, this framework focuses on generating mechanistically grounded biological explanations by integrating:

- Large Language Model (LLM)-based scientific reasoning
- Bioinformatics tools and databases, including GO, KEGG, and literature resources
- Structured multi-step analysis planning
- Reproducible code generation and result validation

The goal of scMechanism-Agent is not to replace bioinformatics analysis platforms, but to help clinical and biomedical researchers diagnose data availability, select defensible analysis strategies, generate reproducible workflows, and organize results into interpretable biological mechanism narratives.

---

## Clinical Research Pain Point

Clinical researchers often identify potentially useful public scRNA-seq datasets from GEO or other repositories, but may have difficulty determining:

- Whether the available files are raw reads, 10x matrices, processed Seurat or h5ad objects, or downstream result tables
- Whether disease/control labels and sample-level metadata are sufficient for analysis
- Which analyses are scientifically appropriate for a given dataset
- How to convert example scripts into reproducible local workflows
- How far computational results can be interpreted as disease mechanism hypotheses

scMechanism-Agent is designed to support these steps through data diagnosis, workflow planning, code generation, result inspection, and mechanism-oriented reporting.

---

## Key Features

- Parse clinical research questions into analyzable single-cell tasks
- Diagnose GEO supplementary file lists and local file formats
- Generate analysis plans for:
  - Processed expression matrices
  - 10x Genomics data
  - h5ad objects
  - Seurat RDS objects
  - FASTQ/SRA-based workflows
- Generate or render reproducible code for:
  - Seurat
  - Scanpy
  - CellChat
  - Monocle3
  - SingleR
- Validate result folders for:
  - Metadata completeness
  - Marker gene evidence
  - Differential expression results
  - Enrichment analysis outputs
  - CellChat results
  - Pseudotime results
  - Figure evidence
- Draft analysis-oriented scientific text, including:
  - Methods
  - Results
  - Figure legends
  - Limitations
  - Validation suggestions

---

## Scope and Limitations

scMechanism-Agent does **not** guarantee:

- Automatic platform-side reanalysis of all datasets
- That every GEO dataset is complete or analyzable
- Automatic reconstruction of FASTQ/SRA workflows without user-provided reference genome, sequencing chemistry, and compute choices
- Clinical diagnosis or treatment recommendations
- That computational cell communication, pseudotime, or CNV inference constitutes experimental proof

The framework is intended for research assistance and hypothesis generation. All biological interpretations should be reviewed by domain experts and, where possible, supported by independent validation.

---

## Runtime Requirements

- Python 3.10 or newer for the helper scripts and Scanpy templates.
- R 4.3 or newer for Seurat-oriented scripts.
- Conda or another Python environment manager.
- Package needs depend on the selected script/template; see `references/environment/cross-platform-setup.md`, `references/environment/requirements.md`, `references/environment/path-setup.md`, `references/seurat_pipeline_rules.md`, and `references/scanpy_pipeline_rules.md`.
- Use `scripts/env_setup/check_environment.py` as the preferred Windows/Linux/macOS environment checker. The older PowerShell scripts remain available for Windows/course compatibility.

Relevant environment and pipeline references are provided in:

```text
references/environment/cross-platform-setup.md
references/environment/requirements.md
references/environment/path-setup.md
references/seurat_pipeline_rules.md
references/scanpy_pipeline_rules.md
```

Cross-platform environment preparation:

```bash
conda env create -f environment.yml
conda activate scmechanism-agent
Rscript scripts/env_setup/install_r_packages.R --profile minimal
python scripts/env_setup/check_environment.py --profile minimal
```

Existing environments can be updated with one command:

```bash
python scripts/env_setup/install_minimal_env.py --env-name scmechanism-agent
```

Use `--profile extended` only when running optional CellChat, Monocle3, or
marker-enrichment course modules.

Legacy Windows/course PowerShell helpers remain available:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\env_setup\check_environment.ps1 -CondaEnv seuratv5-course-py
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\env_setup\install_environment.ps1 -InstallRPackages -InstallPythonEnv
```

---

## Repository Structure

```text
SKILL.md
README.md
CODEBASE_SUMMARY.md
environment.yml
agents/
references/
templates/
scripts/
  env_setup/
  course_source/
  course_adapted/
```

### Main Directories

| Directory                 | Description                                                  |
| ------------------------- | ------------------------------------------------------------ |
| `agents/`                 | Agent definitions and reasoning modules                      |
| `references/`             | Environment notes, pipeline rules, and methodological references |
| `templates/`              | Reusable Scanpy and reporting templates                      |
| `scripts/env_setup/`      | Environment checking and installation scripts                |
| `scripts/course_source/`  | English-named source reference scripts adapted from Seurat V5 course materials |
| `scripts/course_adapted/` | Runnable and parameterized Seurat-oriented analysis scripts  |

---

## Main Script Usage

### 1. Diagnose GEO Input Files

```bash
python scripts/diagnose_geo_inputs.py \
  --file-list supplementary_files.txt \
  --out-json diagnosis.json \
  --out-md diagnosis.md
```

### 2. Build an Analysis Plan

```bash
python scripts/build_analysis_plan.py \
  --diagnosis-json diagnosis.json \
  --question "lung cancer immune microenvironment" \
  --organism human \
  --comparison "tumor vs normal" \
  --out-md analysis_plan.md
```

### 3. Render a Seurat V5 Core Pipeline

```bash
python scripts/render_template.py \
  --template scripts/course_adapted/01_seurat_v5_core_pipeline.R \
  --out run/01_seurat_v5_core_pipeline.R \
  --define INPUT_PATH=/data/GSE \
  --define OUTPUT_DIR=analysis/GSE
```

### 4. Render a SingleR Cell Annotation Pipeline

```bash
python scripts/render_template.py \
  --template scripts/course_adapted/05_singler_cell_annotation.R \
  --out run/05_singler_cell_annotation.R \
  --define INPUT_RDS=analysis/GSE/objects/processed_seurat.rds \
  --define OUTPUT_DIR=analysis/GSE/singler
```

### 5. Render Optional Scanpy Modules

```bash
python scripts/render_template.py \
  --template templates/scanpy_batch_annotation_enrichment_template.py \
  --out run/scanpy_optional_modules.py \
  --define INPUT_H5AD=analysis/GSE/objects/processed.h5ad \
  --define OUTPUT_DIR=analysis/GSE/scanpy_optional
```

### 6. Validate an Analysis Result Folder

```bash
python scripts/validate_result_bundle.py \
  --result-dir analysis/GSE \
  --out-md result_quality_check.md
```

### 7. Generate an Analysis Report

```bash
python scripts/write_analysis_report.py \
  --result-dir analysis/GSE \
  --metadata-json analysis/GSE/report_metadata.json \
  --out-md analysis/GSE/manuscript_report.md
```

### 8. Validate a Full Workflow

```bash
python scripts/validate_full_workflow.py \
  --project-root . \
  --example-root analysis/GSE \
  --out-md analysis/GSE/full_workflow_validation.md
```

---

## Seurat V5 Course Code Adaptation

The Seurat V5 course archive was reorganized into English filenames under `scripts/course_source/` and adapted into runnable scripts under `scripts/course_adapted/`.

The main adapted modules include:

- `scripts/course_adapted/01_seurat_v5_core_pipeline.R`  
  Adapted from 10x import, H5 import, quality control, normalization, DecontX, Harmony, and clustering modules.

- `scripts/course_adapted/00_multi_sample_merge_harmony.R`  
  Adapted from multi-sample import, merge, and Harmony integration modules.

- `scripts/course_adapted/02_marker_enrichment_from_seurat.R`  
  Adapted from marker detection and GO/KEGG enrichment modules.

- `scripts/course_adapted/03_cellchat_from_seurat.R`  
  Adapted from CellChat analysis modules.

- `scripts/course_adapted/04_monocle3_from_seurat.R`  
  Adapted from Monocle3 pseudotime analysis modules.

- `scripts/course_adapted/05_singler_cell_annotation.R`  
  Adapted from SingleR cell annotation modules.

The adapted scripts are parameterized and do not preserve interactive working-directory selection, fixed object names, Chinese output names, or example-specific thresholds.

---

## Templates

The repository contains two types of reusable code resources:

- `scripts/course_adapted/`  
  Runnable Seurat modules derived from the Seurat V5 course materials.

- `templates/`  
  Generic reusable templates for Scanpy workflows, optional downstream modules, and report scaffolds.

When both options are available, use:

- `scripts/course_adapted/` for Seurat V5 course-derived workflows
- `templates/` for Scanpy-based workflows and lightweight reporting templates

---

## Workflow Design

scMechanism-Agent separates the research workflow into six layers:

1. Data diagnosis
2. Analysis planning
3. Code generation
4. External execution
5. Result validation
6. Mechanistic interpretation

Heavy computation remains external and auditable. The framework coordinates the analysis steps, generates reproducible scripts, validates available outputs, and assists with structured scientific interpretation.

For a fuller system view, see `docs/architecture.md`.

---

## Supported Use Cases

The framework is designed for cases covering:

1. Clinical research question only
2. 10x MEX files
3. h5ad/H5 files
4. Seurat RDS files
5. FASTQ/SRA-only datasets
6. Uploaded marker tables
7. Uploaded CellChat results

These use cases demonstrate how the framework diagnoses data structure, selects appropriate analysis strategies, and generates reproducible outputs.

---

## Recommended Workflow

A typical workflow is:

```text
Clinical question
    ↓
GEO or local file diagnosis
    ↓
Analysis feasibility assessment
    ↓
Reproducible analysis plan
    ↓
Seurat or Scanpy code generation
    ↓
External execution
    ↓
Result bundle validation
    ↓
Mechanism-oriented report drafting
```

---

## Security and Privacy

Users should not upload sensitive clinical data unless appropriate de-identification, institutional approval, and data governance procedures have been completed.

The framework is intended for biomedical research support and does not provide clinical diagnosis, treatment decisions, or patient-specific medical advice.

See:

```text
SECURITY_AND_PRIVACY.md
```

for additional safety, privacy, and data handling boundaries.

---

## License

Please specify the license for this repository before public release.

---

## Citation

If this framework is used in research, please cite the repository or associated manuscript once available.

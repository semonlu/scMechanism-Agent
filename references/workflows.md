# Workflow Overview

This skill follows a GEO-to-Mechanism workflow inspired by separate planning, execution, validation, and interpretation layers.

## Standard Flow

1. Parse clinical question with `agents/01_clinical_question_parser.md`.
2. Diagnose accession/file format with `agents/02_geo_dataset_and_format_diagnosis.md` and `scripts/diagnose_geo_inputs.py`.
3. Build a plan with `agents/03_analysis_plan_generator.md` and `scripts/build_analysis_plan.py`.
4. Render local code from `scripts/course_adapted/` or `templates/` using `scripts/render_template.py`.
5. User runs local R/Python scripts.
6. Check returned result folder using `agents/05_result_quality_checker.md` and `scripts/validate_result_bundle.py`.
7. Interpret and report with agents 06-08.

## Start Layer By Input Type

- Clinical question only: parse and ask for missing dataset/design.
- Publication URL/PDF/methods text: extract accessions, sample groups, organism, tissue, platform, available file types, and missing metadata before generating code.
- GEO file list: diagnose format and required metadata.
- 10x MEX/H5: generate Seurat or Scanpy template.
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
- risks and limitations;
- which claims are observations, statistics, or hypotheses.

---
name: scmechanism-agent
description: "Use for GEO/SRA/public or local single-cell research workflows: parse clinical research questions, diagnose supplementary file formats, prepare/check Windows R/Python Seurat V5 environments before running code, plan Seurat/Scanpy analysis, generate local runnable code templates, check uploaded result quality, and draft cautious biological interpretation and manuscript/report text. Also use for local 10x matrices, Seurat RDS, h5ad, loom, marker tables, CellChat, pseudotime, enrichment, CNV, deconvolution, Seurat V5 course-derived workflows, and report/log generation after a full run."
---

# scMechanism Agent: Public Single-cell Research Skill

## When To Use

Use this skill when a user wants to turn public or local single-cell data into a reproducible clinical research workflow, especially when they provide:

- GEO/GSM/SRP/SRR accessions.
- Publication URLs, abstracts, methods text, or PDF-derived notes that mention public single-cell datasets.
- GEO supplementary file lists.
- A disease, tissue, organism, clinical mechanism question, or candidate biomarker question.
- 10x MEX, non-standard 10x MEX, 10x H5, Seurat RDS/RDA/h5Seurat, h5ad, loom, CSV/TSV matrices, FASTQ/SRA files, or result tables.
- Local outputs such as marker tables, DE tables, enrichment tables, cell proportion tables, CellChat results, pseudotime results, UMAPs, dot plots, or report drafts.

## First-Version Boundary

This skill is a decision and code-generation controller. It does not promise that the platform itself will automatically download every GEO dataset or run all heavy single-cell computation. The intended workflow is:

1. Parse the clinical/scientific question.
2. Diagnose GEO/SRA/local file formats and data readiness.
3. Generate an analysis plan with risks and required metadata.
4. Check or install the required R/Python environment before running Seurat V5 course-derived code.
5. Select or render local runnable Seurat/Scanpy/annotation/marker/enrichment code.
6. Review cell annotation evidence and generate a downstream proposal before CellChat or Monocle3.
7. Render or run CellChat/Monocle3 only after the user explicitly approves the proposed cell groups, microenvironment, lineage, and root choice.
8. The user runs approved code locally, in RStudio, Python, Colab, or an institutional server when tools are available.
9. Skill checks result quality, removes known-bad rerun artifacts, writes logs, and writes cautious interpretation, Methods, Results, legends, limitations, and validation suggestions.

Do not claim:

- Fully automatic platform-side reanalysis.
- Guaranteed handling of all GEO datasets.
- FASTQ/SRA to matrix reconstruction without reference/chemistry/compute choices.
- Clinical diagnosis, treatment decision support, or individual patient management.
- Conclusions that do not need human review.

## Required First Step

Start every task by determining which of these modes applies:

1. **Question only**: use `agents/01_clinical_question_parser.md`; ask for disease, tissue, organism, dataset/accession, and comparison design if missing.
2. **Publication or methods text**: extract accessions, organism, tissue, sample groups, platform, matrix/object availability, and missing metadata before writing code.
3. **File list/accession diagnosis**: use `agents/02_geo_dataset_and_format_diagnosis.md` and optionally run `scripts/diagnose_geo_inputs.py`.
4. **Analysis planning**: use `agents/03_analysis_plan_generator.md` and optionally run `scripts/build_analysis_plan.py`.
5. **Environment preparation**: before Seurat V5 course-derived execution, read `references/environment/cross-platform-setup.md`, `references/environment/requirements.md`, and `references/environment/path-setup.md`; run `scripts/env_setup/check_environment.py --profile minimal` for the default local check, or `--profile extended` for CellChat/Monocle3/enrichment modules. Windows PowerShell helpers remain available for legacy course installs.
6. **Code generation**: use `agents/04_code_generator.md`; prefer `scripts/course_adapted/` for Seurat V5 workflows and render placeholders through `scripts/render_template.py`.
7. **Annotation and result review**: use `agents/05_result_quality_checker.md` and optionally run `scripts/validate_result_bundle.py`.
8. **Downstream proposal gate**: before CellChat or Monocle3, run `scripts/propose_downstream_modules.py`, present `downstream_proposal.md`, and wait for user approval.
9. **Approved downstream execution**: only after approval, render `scripts/course_adapted/03_cellchat_from_seurat.R` or `scripts/course_adapted/04_monocle3_from_seurat.R` with the approved scope.
10. **Biological interpretation/report**: use agents 06-08, then run `scripts/write_analysis_report.py` to generate a durable report from the actual result tables.

## Script Layout

The skill script folder has three roles:

```text
scripts/
  course_source/        English-named, lightly adapted reference scripts from the Seurat V5 course
  course_adapted/       runnable R scripts adapted from the course modules
  env_setup/            Windows R/Python/Rtools/JAGS environment install and verification scripts
templates/              generic reusable templates for Scanpy and report scaffolds
  *.py                  helper controllers for diagnosis, planning, rendering, validation, and summaries
```

Use the Python helpers for orchestration:

```bash
python scripts/diagnose_geo_inputs.py --file-list supplementary_files.txt --out-json diagnosis.json --out-md diagnosis.md
python scripts/build_analysis_plan.py --diagnosis-json diagnosis.json --question "disease mechanism" --organism human --comparison "disease vs control" --out-md analysis_plan.md
python scripts/render_template.py --template scripts/course_adapted/01_seurat_v5_core_pipeline.R --out run/01_seurat_v5_core_pipeline.R --define INPUT_PATH=/path/to/data --define OUTPUT_DIR=analysis/run1 --define INPUT_TYPE=10x_mtx
python scripts/validate_result_bundle.py --result-dir analysis/run1 --out-md result_quality_check.md
python scripts/propose_downstream_modules.py --result-dir analysis/run1 --out-md analysis/run1/downstream_proposal.md
python scripts/build_codebase_summary.py --course-root scripts/course_source --out CODEBASE_SUMMARY.md
python scripts/write_analysis_report.py --result-dir analysis/run1 --metadata-json analysis/run1/report_metadata.json --out-md analysis/run1/manuscript_report.md
python scripts/validate_full_workflow.py --project-root . --example-root analysis/run1 --out-md analysis/run1/full_workflow_validation.md
```

Use the cross-platform environment helpers before running Seurat V5 course-derived modules:

```bash
python scripts/env_setup/install_minimal_env.py --env-name scmechanism-agent
python scripts/env_setup/check_environment.py --profile minimal
```

For Windows course-archive compatibility, the PowerShell helpers under `scripts/env_setup/*.ps1` are still available.

Use these course-adapted R scripts for Seurat V5 analysis:

```text
scripts/course_adapted/01_seurat_v5_core_pipeline.R
scripts/course_adapted/02_marker_enrichment_from_seurat.R
scripts/course_adapted/03_cellchat_from_seurat.R
scripts/course_adapted/04_monocle3_from_seurat.R
scripts/course_adapted/05_singler_cell_annotation.R
scripts/course_adapted/00_multi_sample_merge_harmony.R
templates/scanpy_basic_pipeline_template.py
templates/scanpy_batch_annotation_enrichment_template.py
```

The Python scripts do not replace the course code. Heavy Seurat V5 computation belongs in the adapted R scripts, which trace back to the English-named course reference scripts in `scripts/course_source/` and the mapping table `scripts/course_source/source_manifest.csv`.

Use `scripts/course_adapted/` for Seurat V5 course-derived modules. Use `templates/` for generic reusable templates, especially Scanpy workflows and Markdown report scaffolds. Render either family through `scripts/render_template.py` when placeholders need to be filled.

## Required Seurat V5 Course Order

For Seurat V5 course-derived scRNA-seq workflows, do not jump directly from import to annotation or downstream interpretation. The analysis plan and generated code must preserve this order and cite the corresponding course-derived module:

1. Input diagnosis and object construction: `05_read_10x_standard.R`, `06_read_10x_nonstandard.R`, `08_read_10x_h5.R`, or mixed-input logic from `09_merge_mixed_inputs.R` -> `scripts/course_adapted/01_seurat_v5_core_pipeline.R` or `00_multi_sample_merge_harmony.R`.
2. Single-cell QC and filtering: `10_quality_control.R` -> `references/singlecell_qc_rules.md` and `01_seurat_v5_core_pipeline.R`.
3. First normalization, variable features, scaling, PCA, and optional Harmony: `11_normalization_decontx_harmony.R` -> `01_seurat_v5_core_pipeline.R` or `00_multi_sample_merge_harmony.R`.
4. Doublet detection: `12_doublet_finder.R` or `13_scdblfinder.R`; propose it when raw counts and per-sample/loading-batch metadata support it, otherwise explain why it is skipped.
5. Doublet removal: remove predicted doublets per sample/loading batch, without hard-coding DoubletFinder classification column names.
6. Post-doublet re-normalization and reprocessing: `14_post_doublet_normalization.R`; rerun normalization, variable features, scaling, PCA, neighbors, clustering, and UMAP on retained cells.
7. Clustering and resolution review: `15_clustering_resolution.R` -> `01_seurat_v5_core_pipeline.R`; review cluster stability and marker support.
8. Marker detection: `23_marker_detection_methods.R` -> `01_seurat_v5_core_pipeline.R` and `02_marker_enrichment_from_seurat.R`.
9. Cell annotation after clustering and marker review: `16_manual_cell_annotation.R`, `17_singler_annotation.R`, `18_scina_annotation.R`, `21_transferdata_annotation.R`, `22_scpred_annotation.R` -> `05_singler_cell_annotation.R` plus manual marker evidence.
10. Downstream proposal gate: only after annotation evidence is reviewed, propose CellChat, Monocle3/pseudotime, CNV, enrichment, or deconvolution scope and wait for user approval before rendering or running those modules.

If doublet detection is skipped, the report must state the reason. If doublets are removed, all later clustering, annotation, CellChat, and pseudotime steps must use the post-doublet reprocessed object.

`03_cellchat_from_seurat.R` and `04_monocle3_from_seurat.R` are approval-gated scripts. Do not render or run them until `downstream_proposal.md` has been reviewed and the user has explicitly approved the module and cell scope.

## References

Read references selectively:

| Reference | Use when |
|---|---|
| `references/supported_geo_formats.md` | GEO/SRA supplementary file or local input format diagnosis. |
| `references/singlecell_qc_rules.md` | QC threshold planning and result review. |
| `references/seurat_pipeline_rules.md` | Seurat V5 script generation and course-derived R workflow adaptation. |
| `references/scanpy_pipeline_rules.md` | h5ad/Scanpy workflow generation. |
| `references/cell_annotation_rules.md` | Marker-based and reference-assisted annotation. |
| `references/cellchat_rules.md` | CellChat planning, assumptions, and interpretation limits. |
| `references/pseudotime_rules.md` | Monocle3/trajectory planning and review. |
| `references/virtual_knockout_extension.md` | Future perturbation/virtual knockout positioning. |
| `references/clinical_translation_rules.md` | Clinical mechanism interpretation and safety wording. |
| `references/output_file_checklist.md` | Expected output files and quality checks. |
| `references/course-code-index.md` | Exact mapping from Seurat V5 course scripts to skill scripts. |
| `references/course-adaptation.md` | How to adapt the Chinese Seurat V5 course code safely. |
| `references/environment/cross-platform-setup.md` | Windows/Linux/macOS environment profiles and setup commands. |
| `references/environment/requirements.md` | R, Python, R package, and external software inventory. |
| `references/environment/path-setup.md` | PATH, Rscript, Rtools, JAGS, gzip, and Conda setup rules. |
| `references/tested-lessons.md` | Hard-won fixes from real Windows/Seurat V5 test runs. |

## Language Convention

Agent and reference files may mix Chinese and English intentionally:

- user-facing reasoning, clinical framing, and validation guidance can be Chinese.
- file-format routing, code generation, and package/runtime details can be English where it improves exactness.

When editing a file, keep its existing language unless a user explicitly requests translation.

## Course Code Adaptation

The local Seurat V5 course archive is reorganized inside the skill at `scripts/course_source/` as English-named source evidence. Runnable versions live in `scripts/course_adapted/`. The course contributes logic for:

- 10x MEX, non-standard 10x MEX, and H5 import.
- Seurat object construction and QC metrics.
- QC visualization and filtering.
- NormalizeData, variable features, ScaleData, PCA, Harmony, UMAP/tSNE.
- DoubletFinder/scDblFinder options.
- Resolution sweep and clustering.
- Manual/SingleR/SCINA/TransferData/scPred/LLM-assisted annotation.
- Marker detection, GO/KEGG enrichment, CellChat, Monocle2/3, copykat, inferCNV, hdWGCNA, CIBERSORT, and MuSiC.
- Multi-sample merge/Harmony and SingleR annotation are exposed as adapted modules, not hidden assumptions in the core pipeline.
- CellChat and Monocle3 are not automatic next steps; they require annotation review and a user-approved downstream proposal.

Adaptation rules:

- Replace interactive working-directory selection with explicit `input_path`, `metadata_path`, `output_dir`, and config values.
- Replace hard-coded sample column names such as `Type` with user-provided columns.
- For non-standard 10x folders containing `count_matrix_sparse.mtx`, `count_matrix_barcodes.tsv`, and `count_matrix_genes.tsv`, route to `INPUT_TYPE=10x_nonstandard` instead of forcing users to manually rename files to standard 10x names.
- Treat course thresholds and resolution values as examples, not defaults.
- Record English source script, original course path from `source_manifest.csv`, adapted R script/template, deviation, and method status.
- Do not upload private expression matrices or clinical metadata to external LLM/API services.

## Output Standard

A good answer or generated artifact should include:

- Input diagnosis and data readiness.
- Missing metadata/questions.
- Analysis plan with module rationale.
- Cell annotation evidence and uncertainty status before downstream interpretation.
- A downstream module proposal and user approval record before CellChat or pseudotime execution.
- Local runnable code or script names.
- Expected output files.
- Quality-control and interpretation limits.
- Clear separation among observation, statistical inference, and mechanism hypothesis.
- A generated `manuscript_report.md` or equivalent after result review when a full workflow has completed.
- A run log that records environment checks, rendered scripts, executed modules, warnings, and output locations.
- A passing `validate_full_workflow.py` check before sharing a completed local workflow.

## Safety

This skill supports clinical research assistance and hypothesis generation. It must not present public single-cell reanalysis as clinical diagnostic evidence. Cell-cell communication, pseudotime, CNV, and perturbation results require independent validation and should be worded as computational inference unless experimentally proven.

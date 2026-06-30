# Full Workflow Contract

This contract defines what the first version of the skill can produce. It is not a promise that the platform runs all computation automatically.

## Required Modules

| Module | First-version behavior | Script/template |
|---|---|---|
| Question parsing | Structured clinical research question | `agents/01_clinical_question_parser.md` |
| Format diagnosis | Identify GEO/SRA/local formats | `scripts/diagnose_geo_inputs.py` |
| Analysis plan | Produce staged plan with risks | `scripts/build_analysis_plan.py` |
| Seurat code | Render course-adapted local R script | `scripts/course_adapted/01_seurat_v5_core_pipeline.R` |
| Scanpy code | Render local Python template | `templates/scanpy_basic_pipeline_template.py` |
| CellChat | Render course-adapted CellChat script | `scripts/course_adapted/03_cellchat_from_seurat.R` |
| Monocle3 | Render course-adapted trajectory script | `scripts/course_adapted/04_monocle3_from_seurat.R` |
| Result QC | Review uploaded result folder | `scripts/validate_result_bundle.py` |
| Interpretation/report | Draft cautious text | `templates/*.md` and agents 06-08 |

## Optional/Future Modules

| Module | Status |
|---|---|
| GEO automatic download | Future MCP/API extension |
| FASTQ-to-matrix execution | Local/server execution outside first-version platform |
| scVelo/CellRank | Future when spliced/unspliced exists |
| Virtual knockout/perturbation | Future extension; current skill only frames hypotheses |
| Full server-side pipeline orchestration | Future MCP/tooling extension |

## Acceptance Criteria

- No missing function is described as implemented.
- Every generated analysis plan states required metadata and limitations.
- Every code script/template exposes input path, output directory, organism, and key parameters.
- Result interpretation does not invent p values, genes, pathways, or clinical claims.
- Reports separate observation, statistical inference, and mechanism hypothesis.

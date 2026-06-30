# scMechanism Agent

`scMechanism Agent` is a Codex skill project for public and local single-cell mechanism research assistance. It helps users diagnose GEO/SRA/local single-cell data formats, plan reproducible Seurat/Scanpy workflows, generate local runnable code, review uploaded result bundles, and draft cautious biological interpretation and manuscript text.

The skill follows the boundary in `CODEX_BUILD_GEO_SINGLECELL_SKILL.md`: the skill plans, generates code, and checks results; users run expensive computation locally or on their own servers.

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
- Render or copy Seurat, Scanpy, CellChat, and Monocle3 code.
- Check result folders for metadata, marker, DE, enrichment, CellChat, pseudotime, and figure evidence.
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
- Package needs depend on the selected script/template; see `references/seurat_pipeline_rules.md` and `references/scanpy_pipeline_rules.md`.

## File Structure

```text
SKILL.md
README.md
CODEBASE_SUMMARY.md
agents/
references/
templates/
examples/
scripts/
  course_source/
  course_adapted/
```

## Main Script Calls

```bash
python scripts/diagnose_geo_inputs.py --file-list supplementary_files.txt --out-json diagnosis.json --out-md diagnosis.md
python scripts/build_analysis_plan.py --diagnosis-json diagnosis.json --question "lung cancer immune microenvironment" --organism human --comparison "tumor vs normal" --out-md analysis_plan.md
python scripts/render_template.py --template scripts/course_adapted/01_seurat_v5_core_pipeline.R --out run/01_seurat_v5_core_pipeline.R --define INPUT_PATH=/data/GSE --define OUTPUT_DIR=analysis/GSE
python scripts/validate_result_bundle.py --result-dir analysis/GSE --out-md result_quality_check.md
```

The Python scripts are helpers. The Seurat V5 course-derived R analysis code is in:

- `scripts/course_source/`: English-named, lightly adapted course reference R scripts plus `source_manifest.csv`.
- `scripts/course_adapted/`: runnable, parameterized R scripts adapted from the course modules.

## Course Code Reuse

The Seurat V5 course archive was reorganized into English filenames under `scripts/course_source/` and then adapted into runnable scripts under `scripts/course_adapted/`:

- `scripts/course_adapted/01_seurat_v5_core_pipeline.R`: adapted from `05_read_10x_standard.R`, `08_read_10x_h5.R`, `10_quality_control.R`, `11_normalization_decontx_harmony.R`, and `15_clustering_resolution.R`.
- `scripts/course_adapted/02_marker_enrichment_from_seurat.R`: adapted from `23_marker_detection_methods.R` and `24_go_kegg_enrichment.R`.
- `scripts/course_adapted/03_cellchat_from_seurat.R`: adapted from `27_cellchat_analysis.R`.
- `scripts/course_adapted/04_monocle3_from_seurat.R`: adapted from `26_monocle3_pseudotime.R`.

The adapted scripts are parameterized and do not preserve interactive working-directory selection, fixed object names, Chinese output names, or example-specific thresholds. The source reference scripts also have English filenames and a header describing their original source.

## Workflow Layering

The skill keeps data diagnosis, planning, code generation, execution, validation, and interpretation as separate steps. Heavy computation remains external and auditable; the skill coordinates the steps and produces reproducible scripts and review artifacts.

## Validation Cases

See `examples/` for cases covering:

1. Clinical question only.
2. 10x MEX files.
3. h5ad/H5 files.
4. Seurat RDS files.
5. FASTQ/SRA only.
6. Uploaded marker tables.
7. Uploaded CellChat results.

## Packaging

To package the project from its parent directory:

```powershell
Compress-Archive -Path .\scMechanism-Agent -DestinationPath .\scMechanism-Agent.zip -Force
```

The resulting archive contains the skill folder, including `SKILL.md`, `agents/`, `references/`, `templates/`, `examples/`, `scripts/course_source/`, and `scripts/course_adapted/`.

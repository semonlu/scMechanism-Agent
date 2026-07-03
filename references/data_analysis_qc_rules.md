# Data Analysis QC Rules

## Purpose

Before running Seurat, Scanpy, CellChat, Monocle3, or interpretation, confirm that the analyzed data are the same data selected during planning, GEO download, archive extraction, or user upload.

This check is separate from biological QC. It answers: "Are we analyzing the intended dataset?"

## Required Files

```text
data_input_manifest.json
data_analysis_qc.md
```

`data_input_manifest.json` records the project input sources:

- GEO accession or user-provided source label.
- downloaded or registered `workspace/inputs/...` path.
- input type such as `10x_mtx`, `10x_nonstandard`, `10x_h5`, `h5ad`, `rds`, or `csv`.
- representative files, file count, sampled byte size, metadata path, and registration time.

`data_analysis_qc.md` records whether the actual analysis `input_path` matches the manifest.

## MCP Rules

- For GEO data, call `start_geo_download`; the backend registers the download directory automatically.
- If an archive is extracted, call `start_extract_workspace_archive`; the backend registers the extracted directory when it can infer the project from `workspace/inputs/<project>/...`.
- For user-uploaded or manually copied data, call `register_input_dataset(project_id, input_path, input_type, metadata_path)` before `start_seurat_basic` or `start_scanpy_basic`.
- Before running analysis, call `validate_data_analysis_qc(project_id, input_path, input_type, metadata_path)`.
- `start_seurat_basic` and `start_scanpy_basic` also run this QC internally and write `data_analysis_qc.md`.
- If a project has manifest entries and the analysis `input_path` does not match any entry, stop and ask the user to select the correct input.

## Non-MCP Rules

When running locally without MCP, run:

```bash
python scripts/validate_data_sync.py \
  --result-dir analysis/run1 \
  --input-path /path/to/actual/input \
  --input-type 10x_mtx \
  --manifest analysis/run1/data_input_manifest.json \
  --out-md analysis/run1/data_analysis_qc.md
```

If no manifest exists, state that data provenance is incomplete and ask the user to register the intended input or provide the download/extraction record.

## Interpretation Rules

- Do not report a completed workflow unless `data_analysis_qc.md` is present.
- If `data_analysis_qc.md` status is `error`, do not interpret marker, enrichment, CellChat, or pseudotime outputs.
- If status is `warning` because no manifest exists, report the run as locally executed but not fully traceable.
- If the input path matches a registered directory but metadata are missing, continue only for exploratory QC/UMAP; avoid group comparison claims.

## Report Language

Use precise wording:

```text
The analysis input matched the registered/downloaded dataset for this project.
```

or:

```text
The analysis input could not be matched to the project data manifest; results should not be interpreted until the correct input is confirmed.
```

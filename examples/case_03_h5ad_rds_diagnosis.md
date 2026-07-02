# Case 3: h5ad / RDS Object Diagnosis

## Input

```text
杩欎釜GEO鏁版嵁鎻愪緵浜?processed_data.h5ad 鍜?seurat_object.rds锛屽彲浠ユ€庝箞鍒嗘瀽锛?```

## Expected Output

- h5ad: use Scanpy/anndata and inspect `obs`, `var`, `X`, `raw`, and `layers`.
- RDS: use R/Seurat and inspect `meta.data`, assays/layers, reductions, and annotations.
- Confirm raw counts, sample metadata, batch/sample columns, and cell type labels before interpretation or DEG analysis.

## Expected Triggered Files

- `agents/02_geo_dataset_and_format_diagnosis.md`
- `references/supported_geo_formats.md`
- `references/scanpy_pipeline_rules.md`
- `references/seurat_pipeline_rules.md`

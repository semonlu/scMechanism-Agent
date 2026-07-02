# Case 2: 10x MEX File Format Diagnosis

## Input

```text
GEO补充文件包括：
GSM1_matrix.mtx.gz
GSM1_barcodes.tsv.gz
GSM1_features.tsv.gz
GSM2_matrix.mtx.gz
GSM2_barcodes.tsv.gz
GSM2_features.tsv.gz
metadata.csv
```

## Expected Output

- Format: 10x MEX.
- Ready for downstream analysis if matrix/features/barcodes are paired per sample.
- Recommended readers: `Seurat::Read10X()` and `scanpy.read_10x_mtx()`.
- FASTQ reconstruction is not required.
- Next modules: QC, clustering, annotation, marker genes, DEG, enrichment.

## Expected Triggered Files

- `agents/02_geo_dataset_and_format_diagnosis.md`
- `references/supported_geo_formats.md`
- `references/seurat_pipeline_rules.md`
- `references/scanpy_pipeline_rules.md`

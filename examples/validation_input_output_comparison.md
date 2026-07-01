# Validation Input-Output Comparison Report

This is the fixed 5-case validation set for the Medical AI Skill platform. After uploading and enabling the skill on the platform, run these prompts again and paste the platform output summary or screenshot links into the "Platform observed output" fields.

## Summary

| Case | Input type | Expected output | Platform observed output | Status | Revision record |
|---|---|---|---|---|---|
| Case 1 | Clinical research question | Identify disease, tissue, key cells, and analysis route | To be filled after platform run | Pending platform run | None |
| Case 2 | 10x MEX supplementary files | Recognize 10x MEX and recommend Read10X / scanpy.read_10x_mtx | To be filled after platform run | Pending platform run | None |
| Case 3 | h5ad / RDS processed objects | Distinguish Scanpy and Seurat routes | To be filled after platform run | Pending platform run | None |
| Case 4 | FASTQ / SRA raw data | Refuse direct Seurat analysis and route through expression matrix generation | To be filled after platform run | Pending platform run | None |
| Case 5 | Uploaded result tables | Quality-check results and draft cautious Results/mechanism text | To be filled after platform run | Pending platform run | None |

## Case 1: Clinical Question Input

Input:

```text
我想研究膝骨关节炎滑膜组织中巨噬细胞和成纤维细胞的异常变化，想用 GEO 单细胞数据做机制分析。
```

Expected skill behavior:

- Split the question into disease, tissue, organism if available, target cell types, comparison design, and missing metadata.
- Recommend a reproducible analysis route: dataset search/diagnosis, QC, clustering, cell annotation, cell proportion analysis, differential expression, enrichment, cell communication, and pseudotime where justified.
- Avoid claiming that a GEO reanalysis is clinical diagnostic evidence.

Expected output summary:

```text
Disease: knee osteoarthritis
Tissue: synovium
Key cells: macrophages, fibroblasts, T cells, endothelial cells
Recommended modules: cell annotation, cell proportion, DEG, enrichment, cell communication, pseudotime
Missing information: GEO accession or candidate dataset, organism, disease/control design, sample metadata
```

Expected triggered files/rules:

- `agents/01_clinical_question_parser.md`
- `agents/03_analysis_plan_generator.md`
- `references/clinical_translation_rules.md`
- `references/seurat_pipeline_rules.md`

Platform observed output:

```text
To be filled after platform run.
```

Judgment: Pending platform run.

## Case 2: 10x MEX File Format Diagnosis

Input:

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

Expected skill behavior:

- Recognize 10x MEX-style supplementary files.
- State that FASTQ reconstruction is not required for standard downstream analysis.
- Recommend Seurat `Read10X()` and Scanpy `read_10x_mtx()` routes.
- Ask for sample-level metadata and disease/control grouping before differential analysis.

Expected output summary:

```text
Format: 10x MEX
Ready for downstream analysis: yes, if matrix/features/barcodes are paired per sample
Recommended readers: Seurat::Read10X() / scanpy.read_10x_mtx()
FASTQ reconstruction: not required
Next modules: QC, clustering, annotation, marker genes, DEG, enrichment
```

Expected triggered files/rules:

- `agents/02_geo_dataset_and_format_diagnosis.md`
- `references/supported_geo_formats.md`
- `references/seurat_pipeline_rules.md`
- `references/scanpy_pipeline_rules.md`

Platform observed output:

```text
To be filled after platform run.
```

Judgment: Pending platform run.

## Case 3: h5ad / RDS Object Diagnosis

Input:

```text
这个GEO数据提供了 processed_data.h5ad 和 seurat_object.rds，可以怎么分析？
```

Expected skill behavior:

- Distinguish AnnData/h5ad and Seurat RDS objects.
- Recommend Scanpy for h5ad and R/Seurat for RDS.
- Check whether objects contain counts/raw counts, metadata, embeddings, and cell type annotations.
- Explain that existing processed objects may be suitable for interpretation but may limit raw-count differential expression.

Expected output summary:

```text
h5ad: prefer Scanpy/anndata reading and inspect obs, var, X/raw/layers
RDS: prefer R/Seurat reading and inspect meta.data, assays/layers, reductions, annotations
Before analysis: confirm raw counts, sample metadata, batch/sample columns, and cell type labels
```

Expected triggered files/rules:

- `agents/02_geo_dataset_and_format_diagnosis.md`
- `references/supported_geo_formats.md`
- `references/scanpy_pipeline_rules.md`
- `references/seurat_pipeline_rules.md`

Platform observed output:

```text
To be filled after platform run.
```

Judgment: Pending platform run.

## Case 4: FASTQ / SRA Raw Data Boundary

Input:

```text
这个数据只有 SRR 文件，可以直接做 Seurat 分析吗？
```

Expected skill behavior:

- Clearly state that SRR/FASTQ cannot go directly into Seurat/Scanpy downstream analysis.
- Route the user through SRA Toolkit download and expression quantification.
- Ask for organism, reference genome, chemistry/platform, sample metadata, and compute plan.
- Avoid promising fully automatic raw-read processing without those choices.

Expected output summary:

```text
Cannot directly run downstream Seurat/Scanpy analysis from SRR alone.
First generate an expression matrix using SRA Toolkit plus Cell Ranger, kallisto-bustools, or alevin-fry.
Then proceed to QC, clustering, annotation, marker, DEG, and downstream interpretation.
```

Expected triggered files/rules:

- `agents/02_geo_dataset_and_format_diagnosis.md`
- `agents/08_safety_and_limitation_checker.md`
- `references/supported_geo_formats.md`
- `references/full-workflow-contract.md`

Platform observed output:

```text
To be filled after platform run.
```

Judgment: Pending platform run.

## Case 5: Uploaded Result Interpretation

Input:

```text
这是我本地跑出来的 marker_gene.csv、DEG.csv、enrichment.csv 和 cellchat_results.csv，请帮我判断结果是否合理，并生成论文 Results 草稿。
```

Expected skill behavior:

- Switch from code generation to result quality review and biological interpretation.
- Check marker support, differential cell groups, key genes, enrichment terms, and CellChat changes.
- Separate observations, statistical inference, mechanism hypotheses, limitations, and validation suggestions.
- Draft a cautious manuscript Results section without inventing table values.

Expected output summary:

```text
Review: cell annotation support, marker consistency, DEG direction, enrichment plausibility, communication changes
Output: result quality notes, mechanism hypothesis, limitations, validation suggestions, and Results draft
Safety: computational inference only; requires independent validation
```

Expected triggered files/rules:

- `agents/05_result_quality_checker.md`
- `agents/06_biological_interpreter.md`
- `agents/07_report_writer.md`
- `agents/08_safety_and_limitation_checker.md`
- `references/output_file_checklist.md`
- `references/clinical_translation_rules.md`

Platform observed output:

```text
To be filled after platform run.
```

Judgment: Pending platform run.


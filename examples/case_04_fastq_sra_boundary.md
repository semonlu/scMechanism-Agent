# Case 4: FASTQ / SRA Raw Data Boundary

## Input

```text
杩欎釜鏁版嵁鍙湁 SRR 鏂囦欢锛屽彲浠ョ洿鎺ュ仛 Seurat 鍒嗘瀽鍚楋紵
```

## Expected Output

- Do not claim that SRR/FASTQ can directly enter downstream Seurat/Scanpy analysis.
- Explain that SRA Toolkit plus Cell Ranger, kallisto-bustools, or alevin-fry is needed first to generate an expression matrix.
- Ask for organism, reference genome, chemistry/platform, sample metadata, and compute plan.

## Expected Triggered Files

- `agents/02_geo_dataset_and_format_diagnosis.md`
- `agents/08_safety_and_limitation_checker.md`
- `references/supported_geo_formats.md`
- `references/full-workflow-contract.md`

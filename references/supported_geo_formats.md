# Supported GEO And Local Single-cell Formats

| 文件类型 | 典型文件名 | 判断 | Seurat 读取 | Scanpy 读取 | 风险 |
|---|---|---|---|---|---|
| 10x MEX | `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz` | 可直接分析 | `Read10X()` | `sc.read_10x_mtx()` | 需确认每个样本路径和 gene symbol/Ensembl ID |
| 非标准 10x MEX | `count_matrix_sparse.mtx`, `count_matrix_barcodes.tsv`, `count_matrix_genes.tsv`，常见于 GEO tar/tar.gz 解压后的样本目录 | 解压并确认三件套后可直接分析 | `INPUT_TYPE=10x_nonstandard` -> `scripts/course_adapted/01_seurat_v5_core_pipeline.R`；多样本用 `00_multi_sample_merge_harmony.R` 的 `input_type=10x_nonstandard` | 先整理为标准 MEX 或自行用 `scipy.io.mmread` + TSV 读取 | 不能只按文件名判断方向；需检查矩阵维度、barcode 数、gene 数是否匹配 |
| GEO 压缩归档 | `*.tar`, `*.tar.gz`, `*.tgz` | 需先解压再判断 | 解压后按标准 10x、非标准 10x、H5、RDS 或普通矩阵分流 | 同左 | 归档扩展名不能说明是否可直接下游分析；不能误判为 FASTQ |
| 10x H5 | `filtered_feature_bc_matrix.h5`, `raw_feature_bc_matrix.h5` | 可直接分析 | `Read10X_h5()` | `sc.read_10x_h5()` | 需确认是否为 10x 格式而非任意 HDF5 |
| AnnData | `*.h5ad` | 可直接分析 | 需转换或用 Python | `sc.read_h5ad()` | 可能缺 raw counts，DE/CellChat 受限 |
| Seurat | `*.rds`, `*.rda`, `*.h5Seurat` | 可能可用 | `readRDS()`, `load()` | 需转换 | 需确认 assays/layers/metadata/reductions |
| loom | `*.loom` | 可能可用 | 需转换 | `sc.read_loom()` | 常用于 velocity，但需确认 spliced/unspliced |
| 普通矩阵 | `*.csv`, `*.tsv`, `*.txt`, `.gz` | 不确定 | `read.csv/read.delim` | `pandas.read_csv` | 行列方向、raw/normalized、稀疏性不确定 |
| FASTQ/SRA | `SRR*`, `*.fastq.gz`, `*.fq.gz` | 不能直接下游分析 | 先生成矩阵 | 先生成矩阵 | 需要参考基因组、chemistry、计算资源 |
| TCR/BCR | `contig_annotations.csv`, `clonotypes.csv` | 附加分析 | 与表达对象整合 | 与表达对象整合 | 不能替代表达矩阵 |
| 空间转录组 | `spatial/`, `tissue_positions`, `scalefactors` | 需专门流程 | `Load10X_Spatial()` | `squidpy/scanpy` | 坐标、图像和表达矩阵必须匹配 |

## 诊断脚本

```bash
python scripts/diagnose_geo_inputs.py --file-list supplementary_files.txt --out-json diagnosis.json --out-md diagnosis.md
```

## 判断边界

- 只有 FASTQ/SRA 时，不要生成 Seurat 下游分析脚本；先给出矩阵重建方案。
- 有 processed object 时，先审计 raw counts/layers，再决定能否做 DE、CellChat、CNV。
- 缺少 metadata 时，组间比较、比例比较和临床解释都只能作为待补充计划。

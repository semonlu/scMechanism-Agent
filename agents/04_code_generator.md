# 04 Code Generator

## 目标

根据数据格式和分析计划生成本地可运行代码。Seurat V5 相关流程优先使用 `scripts/course_adapted/` 中的课程改造 R 脚本；只有需要替换占位符时才调用 `scripts/render_template.py`。不要临场生成不可解释的大段新代码。

## 调用脚本

```bash
python scripts/render_template.py --template scripts/course_adapted/01_seurat_v5_core_pipeline.R --out run/01_seurat_v5_core_pipeline.R --define INPUT_PATH=/data --define OUTPUT_DIR=analysis/run1
```

## 必须支持

```text
多样本合并/Harmony: scripts/course_adapted/00_multi_sample_merge_harmony.R
Seurat R 基础流程: scripts/course_adapted/01_seurat_v5_core_pipeline.R
SingleR 细胞注释: scripts/course_adapted/05_singler_cell_annotation.R
Marker/GO/KEGG: scripts/course_adapted/02_marker_enrichment_from_seurat.R
CellChat R: scripts/course_adapted/03_cellchat_from_seurat.R
Monocle3 R: scripts/course_adapted/04_monocle3_from_seurat.R
Scanpy Python: templates/scanpy_basic_pipeline_template.py
Scanpy 批次/注释/富集: templates/scanpy_batch_annotation_enrichment_template.py
结果导出与报告模板: templates/*.md
```

## 格式到脚本

- 10x MEX: `scripts/course_adapted/01_seurat_v5_core_pipeline.R` 或 `templates/scanpy_basic_pipeline_template.py`
- 10x H5: Seurat `Read10X_h5` 分支或 Scanpy `read_10x_h5`
- h5ad: `templates/scanpy_basic_pipeline_template.py`
- RDS/RDA: `scripts/course_adapted/01_seurat_v5_core_pipeline.R` 的 RDS 分支或下游脚本
- CSV/TSV: 先判断行列方向和 raw counts/normalized 值
- FASTQ/SRA: 只生成重建矩阵命令建议，不直接进入下游分析
- loom: 建议用 Scanpy `sc.read_loom()` 转为 h5ad 后进入 Scanpy 通用流程

## 决策树

1. 如果用户给出多个样本或多个 GSM 文件，并且每个样本有独立矩阵，优先生成 `00_multi_sample_merge_harmony.R`，保留 `sample_id`、`batch`、`condition`。
2. 如果用户要求 Seurat V5 课程路线，优先选择 `scripts/course_adapted/`，不要转到无溯源模板。
3. 如果用户要求 Python/h5ad/Scanpy，先用 `scanpy_basic_pipeline_template.py`，再按需求追加 `scanpy_batch_annotation_enrichment_template.py`。
4. 如果用户同时有 h5ad 和 RDS，先询问首选生态；没有偏好时，按已有注释和 raw counts 完整度选择，而不是同时生成两套完整流程。
5. 如果输入是 FASTQ/SRA，先生成 raw-to-matrix 方案，不生成 Seurat/Scanpy 下游分析代码。
6. 如果用户上传的是结果表，不生成重跑代码；转到结果质控、解释和报告 agent。
7. 空间转录组数据不要套用普通 scRNA-seq 模板；先要求平台、组织切片信息、spot/cell 坐标和空间图像，再建议 Squidpy/Seurat spatial 路线。

## 规则

- 每段代码必须暴露输入路径、metadata 路径、输出目录和关键参数。
- 默认保存日志、session info、主要中间对象和关键表格。
- 写明用户必须替换的占位符。
- 不要硬编码课程里的 `Type`、固定阈值、固定分辨率或示例对象名。
- 生成代码时标注英文来源脚本和原始课程路径，例如 `scripts/course_source/10_quality_control.R`；原始路径从 `scripts/course_source/source_manifest.csv` 查。

# 04 Code Generator

## 目标

根据数据格式和分析计划生成本地可运行代码。Seurat V5 相关流程优先使用 `scripts/course_adapted/` 中的课程改造 R 脚本；只有需要替换占位符时才调用 `scripts/render_template.py`。不要临场生成不可解释的大段新代码。

## 调用脚本

```bash
python scripts/render_template.py --template scripts/course_adapted/01_seurat_v5_core_pipeline.R --out run/01_seurat_v5_core_pipeline.R --define INPUT_PATH=/data --define OUTPUT_DIR=analysis/run1
```

## 必须支持

```text
Seurat R 基础流程: scripts/course_adapted/01_seurat_v5_core_pipeline.R
Marker/GO/KEGG: scripts/course_adapted/02_marker_enrichment_from_seurat.R
CellChat R: scripts/course_adapted/03_cellchat_from_seurat.R
Monocle3 R: scripts/course_adapted/04_monocle3_from_seurat.R
Scanpy Python: templates/scanpy_basic_pipeline_template.py
结果导出与报告模板: templates/*.md
```

## 格式到脚本

- 10x MEX: `scripts/course_adapted/01_seurat_v5_core_pipeline.R` 或 `templates/scanpy_basic_pipeline_template.py`
- 10x H5: Seurat `Read10X_h5` 分支或 Scanpy `read_10x_h5`
- h5ad: `templates/scanpy_basic_pipeline_template.py`
- RDS/RDA: `scripts/course_adapted/01_seurat_v5_core_pipeline.R` 的 RDS 分支或下游脚本
- CSV/TSV: 先判断行列方向和 raw counts/normalized 值
- FASTQ/SRA: 只生成重建矩阵命令建议，不直接进入下游分析

## 规则

- 每段代码必须暴露输入路径、metadata 路径、输出目录和关键参数。
- 默认保存日志、session info、主要中间对象和关键表格。
- 写明用户必须替换的占位符。
- 不要硬编码课程里的 `Type`、固定阈值、固定分辨率或示例对象名。
- 生成代码时标注英文来源脚本和原始课程路径，例如 `scripts/course_source/10_quality_control.R`；原始路径从 `scripts/course_source/source_manifest.csv` 查。

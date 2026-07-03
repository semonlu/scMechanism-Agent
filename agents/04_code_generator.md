# 04 Code Generator

## 目标

根据数据格式、分析计划、使用者确认状态和已有结果，选择本地可运行脚本。Seurat V5 相关流程优先使用 `scripts/course_adapted/` 中的课程改造 R 脚本；只有需要替换占位符时才调用 `scripts/render_template.py`。

## 可直接生成的脚本

```text
多样本合并/Harmony: scripts/course_adapted/00_multi_sample_merge_harmony.R
Seurat R 基础流程: scripts/course_adapted/01_seurat_v5_core_pipeline.R
SingleR 细胞注释: scripts/course_adapted/05_singler_cell_annotation.R
Marker/GO/KEGG: scripts/course_adapted/02_marker_enrichment_from_seurat.R
Scanpy Python: templates/scanpy_basic_pipeline_template.py
Scanpy 批次/注释/富集: templates/scanpy_batch_annotation_enrichment_template.py
结果检查/报告: scripts/validate_result_bundle.py, scripts/write_analysis_report.py
```

## Seurat V5 脚本顺序

生成 Seurat V5 代码时，必须按课程流程顺序组织，不能把细胞注释放在 QC/标准化/去双/聚类之前：

1. 数据导入和对象构建
   课程代码: `05_read_10x_standard.R`, `06_read_10x_nonstandard.R`, `08_read_10x_h5.R`, `09_merge_mixed_inputs.R`
   调用: `01_seurat_v5_core_pipeline.R` 或 `00_multi_sample_merge_harmony.R`

2. 单细胞 QC 和过滤
   课程代码: `10_quality_control.R`
   调用: `01_seurat_v5_core_pipeline.R` 和 `singlecell_qc_rules.md`

3. 第一次标准化、PCA、Harmony/批次处理
   课程代码: `11_normalization_decontx_harmony.R`
   调用: `01_seurat_v5_core_pipeline.R` 或 `00_multi_sample_merge_harmony.R`

4. 双细胞检测
   课程代码: `12_doublet_finder.R` 或 `13_scdblfinder.R`
   规则: raw counts 和 sample/loading batch 信息足够时要加入；跳过时必须说明原因。

5. 去除双细胞
   课程代码: `12_doublet_finder.R` 或 `13_scdblfinder.R`
   规则: 按样本/批次移除 predicted doublets，不硬编码 DoubletFinder 分类列名。

6. 双细胞后标准化/重新处理
   课程代码: `14_post_doublet_normalization.R`
   规则: 去除双细胞后必须重新 NormalizeData、FindVariableFeatures、ScaleData、PCA、neighbors、clusters、UMAP。

7. 单细胞分群聚类和分辨率检查
   课程代码: `15_clustering_resolution.R`
   调用: `01_seurat_v5_core_pipeline.R`

8. Marker 检测
   课程代码: `23_marker_detection_methods.R`
   调用: `01_seurat_v5_core_pipeline.R` 和 `02_marker_enrichment_from_seurat.R`

9. 细胞注释
   课程代码: `16_manual_cell_annotation.R`, `17_singler_annotation.R`, `18_scina_annotation.R`, `21_transferdata_annotation.R`, `22_scpred_annotation.R`
   调用: `05_singler_cell_annotation.R` 加人工 marker 复核。

## 需要确认后才能生成的脚本

```text
CellChat R: scripts/course_adapted/03_cellchat_from_seurat.R
Monocle3 R: scripts/course_adapted/04_monocle3_from_seurat.R
```

在生成 CellChat 或 Monocle3 脚本前，必须确认已经存在：

- `annotation_evidence.tsv`、`singleR_cluster_labels.csv` 或人工确认的等效注释证据。
- marker 表和 UMAP/DotPlot/FeaturePlot 等注释支撑。
- `downstream_proposal.md`，列出候选细胞群、微环境或拟时序谱系。
- 使用者明确同意运行的模块和范围。

如果缺少这些条件，只能生成“下游模块方案”，不能直接渲染或运行 CellChat/Monocle3。

## 格式到脚本

- 10x MEX: `scripts/course_adapted/01_seurat_v5_core_pipeline.R` 或 `templates/scanpy_basic_pipeline_template.py`
- 非标准 10x MEX: `scripts/course_adapted/01_seurat_v5_core_pipeline.R`，设置 `INPUT_TYPE=10x_nonstandard`；多样本表格用 `00_multi_sample_merge_harmony.R` 并在 `input_type` 列写 `10x_nonstandard`
- 10x H5: Seurat `Read10X_h5` 分支或 Scanpy `read_10x_h5`
- h5ad: `templates/scanpy_basic_pipeline_template.py`
- RDS/RDA: Seurat RDS 分支或下游脚本
- CSV/TSV: 先判断行列方向和 raw counts/normalized 值
- FASTQ/SRA: 只生成 raw-to-matrix 方案，不直接进入下游分析
- loom: 建议用 Scanpy `sc.read_loom()` 转为 h5ad 后进入 Scanpy 通用流程

## 决策规则

1. 多个样本或多个 GSM 文件优先生成 `00_multi_sample_merge_harmony.R`，保留 `sample_id`、`batch`、`condition`。
2. Seurat V5 路线优先选择 `scripts/course_adapted/`。
3. 如果输入和 metadata 支持双细胞检测，必须在计划中加入 DoubletFinder 或 scDblFinder，并说明去双细胞后要重新标准化/聚类。
4. 如果跳过双细胞检测，必须写明原因，例如只有 processed object、细胞数过少、缺少 sample/loading batch 信息、包不可用。
5. Python/h5ad/Scanpy 先用 `scanpy_basic_pipeline_template.py`，再按需求追加批次、注释、富集模板。
6. 同时有 h5ad 和 RDS 时，先询问首选生态；没有偏好时按 raw counts 完整度和已有注释选择。
7. FASTQ/SRA 先生成矩阵重建方案，不生成 Seurat/Scanpy 下游代码。
8. 输入是结果表时，转到结果质控、下游方案、解释和报告，不生成重跑代码。
9. 空间转录组数据不能套用普通 scRNA-seq 模板；先要求平台、组织切片、坐标和图像信息。
10. CellChat/Monocle3 必须先过“下游模块确认门”。没有使用者确认时，不能生成运行命令。

## 输出要求

- 暴露输入路径、metadata 路径、输出目录和关键参数。
- 默认保存日志、session info、主要中间对象和关键表格。
- 写明用户必须替换的占位符。
- 不硬编码课程里的 `Type`、固定阈值、固定分辨率或示例对象名。
- 生成代码时标注英文来源脚本；不要引用本机原始课程目录。

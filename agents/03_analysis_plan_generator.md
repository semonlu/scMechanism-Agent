# 03 Analysis Plan Generator

## 目标

根据研究问题、数据格式、metadata 完整度和已有中间结果，生成可执行但谨慎的单细胞分析计划。计划必须先完成基础 scRNA-seq 流程，再进入细胞注释和下游机制分析。

## 必须输出的层次

```text
A. 数据准备和格式诊断
B. 数据分析质控和输入同步
C. 单细胞 QC 和过滤
D. 第一次标准化、降维和批次处理
E. 双细胞检测与去除方案
F. 去除双细胞后的重新标准化、重新降维和重新聚类
G. 分群聚类、分辨率检查和 marker 检测
H. 细胞注释和注释证据复核
I. 组间比较
J. 功能富集
K. 下游解释性模块候选方案
L. 预期结果文件
M. 风险与限制
N. 需要使用者确认的问题
```

## Seurat V5 课程顺序

当计划选择 Seurat V5 课程来源流程时，必须按下面顺序写明调用依据：

1. 数据导入和对象构建：`05_read_10x_standard.R`、`06_read_10x_nonstandard.R`、`08_read_10x_h5.R`、`09_merge_mixed_inputs.R` -> `01_seurat_v5_core_pipeline.R` 或 `00_multi_sample_merge_harmony.R`。非标准 10x 使用 `INPUT_TYPE=10x_nonstandard`。
2. 数据分析质控：确认实际 `input_path` 与计划下载、解压或用户注册输入一致，输出 `data_input_manifest.json` 和 `data_analysis_qc.md`。
3. 单细胞 QC：`10_quality_control.R` -> `singlecell_qc_rules.md` 和 `01_seurat_v5_core_pipeline.R`。
4. 第一次标准化、PCA、Harmony：`11_normalization_decontx_harmony.R` -> `01_seurat_v5_core_pipeline.R` 或 `00_multi_sample_merge_harmony.R`。
5. 双细胞检测：`12_doublet_finder.R` 或 `13_scdblfinder.R`。如果不跑，必须说明原因。
6. 去除双细胞：按 sample/loading batch 移除 predicted doublets，不能硬编码分类列名。
7. 去除双细胞后重新标准化：`14_post_doublet_normalization.R`；重新跑 NormalizeData、FindVariableFeatures、ScaleData、PCA、neighbors、clusters、UMAP。
8. 分群聚类和分辨率检查：`15_clustering_resolution.R` -> `01_seurat_v5_core_pipeline.R`。
9. Marker 检测：`23_marker_detection_methods.R` -> `01_seurat_v5_core_pipeline.R` 和 `02_marker_enrichment_from_seurat.R`。
10. 细胞注释：`16_manual_cell_annotation.R`、`17_singler_annotation.R`、`18_scina_annotation.R`、`21_transferdata_annotation.R`、`22_scpred_annotation.R` -> `05_singler_cell_annotation.R` 加人工 marker 复核。
11. 下游模块：只有在注释证据复核后，才能提出 CellChat、拟时序、CNV、反卷积等方案，并等待使用者确认。

## 关键规则

- 细胞注释必须先输出 `annotation_evidence.tsv`、`singleR_cluster_labels.csv` 或同等证据表，再进入 CellChat、拟时序、CNV、反卷积等解释性模块。
- 细胞注释不能只依赖 SingleR/CellTypist 等自动标签；必须结合 cluster markers、组织背景、物种、UMAP/DotPlot/FeaturePlot 证据。
- 对证据不足的 cluster 使用 `Unknown`、`Ambiguous` 或上一级粗粒度标签，不要强行命名精细亚型。
- 大群注释先扫 `0.1`、`0.3`、`0.5`、`0.8`，结合 cluster 数量、marker 清晰度和生物学合理性选 resolution。
- 亚群分析必须 subset 目标大群后重新 HVG、scale、PCA、neighbors、clustering 和 UMAP，不能只提高全局 resolution。
- CellChat 和拟时序不是默认必跑模块。必须先根据前序结果写出候选分析方案，并等待使用者确认。
- 下游方案必须说明推荐分析对象、为什么选它、需要哪些前提、预计输出什么、哪些结果不能解释为因果。
- 分析前必须确认 `data_analysis_qc.md` 状态；如果实际 `input_path` 和 `data_input_manifest.json` 不匹配，停止分析并要求用户选择正确输入。

## 下游模块确认门

当计划包含 CellChat 或拟时序时，必须先生成如下内容，而不是直接运行脚本：

```text
候选模块：
1. CellChat / 细胞通讯
   - 推荐微环境/细胞群：
   - 依据：
   - 分组方式：
   - 需要排除的低可信细胞群：
   - 预期输出：
   - 风险：

2. Monocle3 / 拟时序
   - 推荐细胞谱系或连续过程：
   - root 选择建议：
   - 依据：
   - 不建议纳入的细胞：
   - 预期输出：
   - 风险：

请使用者确认：
- 是否运行 CellChat？
- CellChat 使用哪些细胞群/分组？
- 是否运行拟时序？
- 拟时序聚焦哪条谱系，root 是什么？
```

只有收到明确同意后，才进入 `agents/04_code_generator.md` 生成对应脚本。

## 调用脚本

```bash
python scripts/build_analysis_plan.py --diagnosis-json diagnosis.json --question "..." --organism human --comparison "disease vs control" --out-md analysis_plan.md
python scripts/validate_data_sync.py --result-dir analysis/run1 --input-path /path/to/input --input-type 10x_mtx --manifest analysis/run1/data_input_manifest.json --out-md analysis/run1/data_analysis_qc.md
python scripts/propose_downstream_modules.py --result-dir analysis/run1 --out-md analysis/run1/downstream_proposal.md
```

## 必须区分

- 单样本探索分析。
- 多样本疾病组 vs 对照组分析。
- 多批次整合分析。
- 只有 processed data。
- 只有 raw data。
- 缺少分组 metadata。
- 样本量过少或无生物学重复。
- 细胞注释未确认。
- 下游模块已获使用者确认或尚未确认。

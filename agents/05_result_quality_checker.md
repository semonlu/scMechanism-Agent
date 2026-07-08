# 05 Result Quality Checker

## 目标

检查用户上传或本地运行产生的结果是否足以支持继续分析、解释和报告写作。

## 输入可能包括

```text
metadata.csv
data_input_manifest.json
data_analysis_qc.md
qc_summary.tsv
qc_thresholds.tsv
resolution_sweep.tsv
cluster_marker_audit.tsv
workflow_step_audit.tsv
cluster_markers.csv
singleR_cluster_labels.csv
annotation_evidence.tsv
deg_results.csv
celltype_proportion.csv
enrichment_results.csv
cellchat_results.csv
pseudotime_results.csv
umap.png / umap.pdf
dotplot.png / dotplot.pdf
```

## 调用脚本

```bash
python scripts/validate_result_bundle.py --result-dir analysis/run1 --out-md result_quality_check.md
python scripts/validate_data_sync.py --result-dir analysis/run1 --input-path /path/to/input --input-type 10x_mtx --manifest analysis/run1/data_input_manifest.json --out-md analysis/run1/data_analysis_qc.md
python scripts/propose_downstream_modules.py --result-dir analysis/run1 --out-md analysis/run1/downstream_proposal.md
```

## 必须检查

- `data_analysis_qc.md` 是否存在，且实际分析 `input_path` 是否与前一步计划、下载、解压或用户注册的数据一致。
- `data_input_manifest.json` 是否记录了 GEO accession、下载目录、解压目录或用户输入目录。
- 是否有 metadata、分组、样本、批次和 donor/replicate 信息。
- QC 阈值是否有数据分布依据，是否按样本检查。
- 是否考虑环境 RNA、红细胞/血红蛋白污染、线粒体、核糖体、低复杂度和细胞周期。
- 是否做过 doublet 检测；如果跳过，原因是否合理。
- 去除 doublet 后是否重新 NormalizeData、FindVariableFeatures、ScaleData、PCA、neighbors、clusters、UMAP。
- 是否输出 `resolution_sweep.tsv`，并说明为何选择当前 resolution。
- cluster 数量是否异常：常规 10x、5,000 到 20,000 细胞的大群注释通常先期待约 10 到 20 个 cluster；少于约 8 到 10 可能漏掉细胞类型，多于约 25 需要警惕过度拆分。
- marker 是否清楚：top marker 是否能对应经典大群标志，而不是主要由 ribosomal、heat-shock、mitochondrial、hemoglobin、stress 或 cell-cycle 基因驱动。
- 自动注释是否有 `annotation_evidence.tsv`，是否存在低置信度、Ambiguous 或 Unknown cluster。
- 是否存在明显过细、互相矛盾或单一标签覆盖所有细胞的问题。
- DE gene 数量和方向是否异常。
- 富集通路是否与疾病问题和细胞类型相关。
- CellChat 是否建立在可信注释和足够细胞数之上。
- 拟时序是否有合理 root、连续过程和输入细胞子集。

## 进入单细胞流程 10 步审核

每个结果目录必须能回答下面 10 步是否满足。缺证据时要写 `missing_evidence`，不能默认通过。

| 步骤 | 审核内容 | 关键证据 |
|---|---|---|
| 1 | 获取数据，通常 GEO/SRA 或本地可追踪数据 | accession、来源、file manifest、输入诊断 |
| 2 | 读入数据，表达矩阵/对象和 metadata 对齐 | import summary、细胞数、基因数、barcode 对齐 |
| 3 | 质量控制：线粒体、红细胞/血红蛋白、counts、genes、低质量细胞 | QC 图、阈值、过滤前后细胞数 |
| 4 | 标准化和归一化：NormalizeData、ScaleData、FindVariableFeatures、SCT 或 Scanpy 等价步骤 | normalization 参数、HVG、raw counts 保留 |
| 5 | 去除或审核干扰因素：细胞周期、双细胞、不同批次、ambient RNA | doublet、cell-cycle、batch/integration、ambient RNA 状态 |
| 6 | 降维聚类：PCA、UMAP、tSNE、neighbors、resolution | PCA/UMAP/tSNE、resolution sweep、cluster audit |
| 7 | 细胞注释：自动注释和人工注释 | annotation_evidence、marker 图、confidence |
| 8 | 寻找 marker 基因：FindMarkers/FindAllMarkers | marker 表、marker specificity audit |
| 9 | 拟时序分析 | 下游审批、subset/root 依据、trajectory 输出 |
| 10 | 细胞通讯 | 下游审批、细胞数门槛、ligand-receptor 输出 |

步骤 1-8 是进入解释前的核心上游证据；步骤 9-10 是可选下游模块，但一旦运行必须有审批和解释边界。

## 细胞注释判定

注释结果分为：

```text
可用于下游：有 marker 支持、自动注释一致、关键图可复核。
需要人工复核：自动注释和 marker 不一致，或存在 Ambiguous/Unknown cluster。
不能用于下游：缺 marker 表、缺注释证据、所有 cluster 被标成同一类型、或标签明显不符合组织背景。
```

当注释不能用于下游时，不得建议直接运行 CellChat 或拟时序。

当数据分析质控显示输入不匹配时，不得解释任何 marker、DE、富集、CellChat 或拟时序结果；必须先回到数据选择/下载/注册步骤。

## 输出格式

```text
总体质量判断：可靠 / 基本可靠 / 需谨慎 / 不建议解释
数据同步状态：通过 / 缺少 manifest / input_path 不匹配
细胞注释状态：可用于下游 / 需要人工复核 / 不能用于下游
主要支持点：
主要风险：
需要补充的结果：
是否允许进入 CellChat/拟时序候选方案：
是否适合写入论文：
```

# 05 Result Quality Checker

## 目标

检查用户上传或本地运行产生的结果是否足以支持继续分析、解释和报告写作。

## 输入可能包括

```text
metadata.csv
data_input_manifest.json
data_analysis_qc.md
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
- QC 阈值是否过严或过松。
- cluster 数量是否异常。
- marker 是否符合组织背景和细胞类型。
- 自动注释是否有证据表和不确定性标记。
- 是否存在明显过细、互相矛盾或单一标签覆盖所有细胞的问题。
- DE gene 数量和方向是否异常。
- 富集通路是否与疾病问题和细胞类型相关。
- CellChat 是否建立在可信注释和足够细胞数之上。
- 拟时序是否有合理 root、连续过程和输入细胞子集。

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

# 05 Result Quality Checker

## 目标

检查用户上传的结果是否足以支持分析解读。

## 输入可能包括

```text
metadata.csv
cluster_markers.csv
deg_results.csv
celltype_proportion.csv
enrichment_results.csv
cellchat_results.csv
pseudotime_results.csv
umap.png
dotplot.png
```

## 调用脚本

```bash
python scripts/validate_result_bundle.py --result-dir analysis/run1 --out-md result_quality_check.md
```

## 必须检查

- 是否有分组信息。
- 样本量和生物学重复是否足够。
- 是否有批次信息。
- QC 阈值是否过严或过松。
- cluster 数量是否异常。
- marker 是否符合细胞类型。
- DE gene 数量是否异常。
- 富集通路是否与疾病问题相关。
- CellChat 是否被过度解释。
- 拟时序是否有合理 root 和连续生物学过程。

## 输出格式

```text
总体质量判断：可靠/基本可靠/需谨慎/不建议解读
主要支持点：
主要风险：
需要补充的结果：
是否适合写入论文：
```

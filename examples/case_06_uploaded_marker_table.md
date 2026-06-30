# Case 06 Uploaded Marker Table

## 用户输入

上传 `cluster_markers.csv`，希望判断细胞注释是否合理。

## 期望 Skill 判断

进入结果质控和注释复核，不重新跑全流程。

## 期望输出结构

- 检查 marker 表字段：cluster、gene、avg_log2FC、pct、p_val_adj
- 对照 canonical marker
- 输出可信/需谨慎标签

## 验证重点

不能仅凭一个 marker 自动确定所有细胞类型。

## 可能错误输出

把所有 cluster 强行命名。

## 修正规则

对低证据 cluster 标记 ambiguous。

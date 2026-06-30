# Case 03 h5ad Or H5 Format

## 用户输入

GSEXXXX_processed.h5ad 和 filtered_feature_bc_matrix.h5

## 期望 Skill 判断

h5ad 优先 Scanpy 读取；10x H5 可用 Seurat `Read10X_h5()` 或 Scanpy `read_10x_h5()`。必须提醒检查 raw counts。

## 期望输出结构

- 识别格式：AnnData h5ad；10x HDF5
- 风险：processed h5ad 可能缺 raw counts
- 推荐：先审计 `.raw` 和 `.layers`

## 验证重点

不能默认 h5ad 可做所有 DE/CellChat。

## 可能错误输出

直接运行 CellChat 而不确认 raw counts 和 cell labels。

## 修正规则

先用 Scanpy 模板检查对象结构。

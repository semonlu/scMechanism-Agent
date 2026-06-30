# Case 02 10x MEX Format

## 用户输入

GSMxxx_matrix.mtx.gz, GSMxxx_barcodes.tsv.gz, GSMxxx_features.tsv.gz

## 期望 Skill 判断

识别为 10x MEX，可进入 Seurat/Scanpy 下游分析。

## 期望输出结构

- 识别格式：10x MEX
- 推荐读取：Seurat `Read10X()` 或 Scanpy `read_10x_mtx()`
- 下一步：补充 metadata 和样本分组

## 验证重点

必须检查每个样本是否有完整三件套。

## 可能错误输出

把 matrix.mtx 当普通 CSV。

## 修正规则

用 `scripts/diagnose_geo_inputs.py` 验证文件列表。

# Case 04 RDS Format

## 用户输入

作者提供了 `GSEXXXX_seurat_object.rds`。

## 期望 Skill 判断

RDS 可能是 Seurat 对象，需要 `readRDS()` 后检查 class、assays、metadata、reductions、raw counts。

## 期望输出结构

- 识别格式：R/Seurat object
- 推荐读取：Seurat R 模板 RDS 分支
- 风险：可能缺 counts 或样本 metadata

## 验证重点

不能假设任何 RDS 都是可分析 Seurat 对象。

## 可能错误输出

直接说可运行全部流程。

## 修正规则

先生成对象审计步骤。

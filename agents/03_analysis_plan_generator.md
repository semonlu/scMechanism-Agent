# 03 Analysis Plan Generator

## 目标

根据研究问题、数据格式和 metadata 完整度生成单细胞分析计划。

## 必须输出层次

```text
A. 数据准备
B. 质量控制
C. 标准化与降维
D. 聚类与细胞注释
E. 组间比较
F. 功能富集
G. 细胞通讯
H. 拟时序/细胞命运
I. 可选扩展：RNA velocity / 虚拟敲除
J. 预期结果文件
K. 风险与限制
```

## 调用脚本

```bash
python scripts/build_analysis_plan.py --diagnosis-json diagnosis.json --question "..." --organism human --comparison "disease vs control" --out-md analysis_plan.md
```

## 必须区分

- 单样本探索分析。
- 多样本疾病组 vs 对照组分析。
- 多批次整合分析。
- 只有 processed data。
- 只有 raw data。
- 缺失分组 metadata。
- 样本量过小或无生物学重复。

## 规则

- 先说明哪些分析可以做，哪些只能做描述性展示。
- 重分析计划必须保留 raw counts 或说明 raw counts 缺失。
- 组间差异和细胞比例比较必须关注样本级重复。
- CellChat、pseudotime、CNV、velocity 都是有前提的可选模块。

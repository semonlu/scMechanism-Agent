# 03 Analysis Plan Generator

## 目标

根据研究问题、数据格式、metadata 完整度和已产生的中间结果，生成可执行但谨慎的单细胞分析计划。

## 必须输出层次

```text
A. 数据准备
B. 质量控制
C. 标准化、降维和批次处理
D. 聚类与细胞注释
E. 组间比较
F. 功能富集
G. 下游解释性模块候选方案
H. 预期结果文件
I. 风险与限制
J. 需要使用者确认的问题
```

## 关键规则

- 细胞注释必须先输出 `annotation_evidence.tsv` 或同等证据表，再进入 CellChat、拟时序、CNV、反卷积等解释性模块。
- 细胞注释不能只依赖 SingleR/CellTypist 等自动标签；必须结合 cluster markers、组织背景、物种、UMAP/DotPlot/FeaturePlot 证据。
- 对证据不足的 cluster 使用 `Unknown`、`Ambiguous` 或上一级粗粒度标签，不要强行命名精细亚型。
- CellChat 和拟时序不是默认必跑模块。必须先根据前序结果写出候选分析方案，并等待使用者确认。
- 下游方案必须说明：推荐分析对象、为什么选它、需要哪些前提、预计输出什么、哪些结果不应解释为因果。

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

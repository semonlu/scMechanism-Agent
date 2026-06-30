# 07 Report Writer

## 目标

生成论文、课题申报或内部汇报可用文本。

## 必须能生成

```text
Methods
Results
Figure legends
Supplementary table list
Limitations
Future directions
基金/课题依据摘要
```

## 规则

- 不虚构具体结果。
- 用户没有提供结果表时，不写具体 p 值、基因名或通路名。
- Methods 必须写清数据来源、对象构建、QC、降维聚类、注释、差异表达、富集、细胞通讯和拟时序工具。
- Results 要区分“观察到的结果”和“推测的机制”。
- Limitations 必须覆盖样本量、批次、注释质量、数据库推断和外部验证需求。

## 模板

- `templates/methods_template.md`
- `templates/results_template.md`
- `templates/result_interpretation_template.md`
- `templates/validation_report_template.md`

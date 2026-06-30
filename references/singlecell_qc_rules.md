# Single-cell QC Rules

QC 阈值必须来自当前数据分布、组织背景、物种、平台和预期细胞类型。课程脚本中的阈值只能作为示例，不能直接套用。

## 常用指标

| 指标 | 含义 | 风险提示 |
|---|---|---|
| `nFeature_RNA` | 每个细胞检测到的基因数 | 过低可能为空滴/低质量；过高可能双细胞 |
| `nCount_RNA` | 每个细胞 UMI/reads 总数 | 过低低质量；过高可能双细胞或高 RNA content 细胞 |
| `percent.mt` | 线粒体比例 | 高值提示应激/死亡；组织和物种前缀不同 |
| `percent.ribo` | 核糖体比例 | 高值可能反映细胞状态或低质量，需结合组织 |
| `percent.hb` | 血红蛋白比例 | 血液污染或红细胞相关信号 |
| doublet score | 双细胞风险 | 应按样本或加载批次评估 |
| batch | 批次/样本效应 | 不应在未检查前盲目整合 |

## 推荐输出

- QC 前后小提琴图。
- nCount vs percent.mt、nCount vs nFeature 散点图。
- 每个样本过滤前后细胞数。
- 过滤规则和阈值解释。
- doublet 方法、预期双细胞率、每样本移除数。

## Course-Derived Notes

课程脚本使用 `VlnPlot()`、`FeatureScatter()` 和 `subset()`，示例阈值包括 `nCount_RNA <= 25000`、`nFeature_RNA <= 5000`、`percent.mt <= 25`、`percent.rb <= 40`。这些值不是通用默认值。

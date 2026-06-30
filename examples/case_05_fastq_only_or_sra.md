# Case 05 FASTQ Only Or SRA

## 用户输入

只有 SRR12345678 或 `sample_R1.fastq.gz`, `sample_R2.fastq.gz`。

## 期望 Skill 判断

不能直接进入 Seurat/Scanpy 下游分析。需要先用 Cell Ranger、STARsolo、kallisto-bustools 或 alevin-fry 生成表达矩阵。

## 期望输出结构

- 识别格式：FASTQ/SRA raw reads
- 是否直接下游：否
- 下一步：确认物种、参考基因组、chemistry、样本表

## 验证重点

不能说 FASTQ 可以直接 CreateSeuratObject。

## 可能错误输出

跳过比对/定量步骤。

## 修正规则

只生成重建矩阵命令模板和风险说明。

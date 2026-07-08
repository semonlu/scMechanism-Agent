#!/usr/bin/env python3
"""Build a cautious GEO-to-mechanism single-cell analysis plan."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


SECTIONS = [
    ("A. 数据准备", [
        "确认 GEO/SRA 编号、supplementary files、物种、组织、样本数和分组设计。",
        "识别输入格式；FASTQ/SRA 先重建表达矩阵，processed object 先审计 raw counts/layers。",
        "建立 metadata_template.csv，至少包含 sample_id、condition、batch、donor_id、tissue、platform。"
    ]),
    ("B. 数据分析质控和输入同步", [
        "记录计划下载、解压或用户上传的真实输入，生成 data_input_manifest.json。",
        "运行 validate_data_sync.py 或 MCP validate_data_analysis_qc，确认实际 input_path 与 manifest 匹配。",
        "如果 input_path 不匹配，停止 Seurat/Scanpy 分析并要求重新选择正确输入。"
    ]),
    ("C. 质量控制", [
        "生成 QC 分布图和 pre/post cell count 表；阈值来自数据分布而非固定套用课程示例。",
        "计算 mitochondrial/ribosomal/hemoglobin 指标，并按物种调整基因前缀。",
        "根据样本和加载密度选择 DoubletFinder、scDblFinder 或 Scanpy/Scrublet。"
    ]),
    ("D. 标准化与降维", [
        "运行 NormalizeData/SCTransform 或 Scanpy normalize_total/log1p。",
        "选择高变基因、PCA、邻接图、UMAP/tSNE；记录 PC 数和随机种子。",
        "如有批次，先可视化批次结构，再决定 Harmony/CCA/SCVI 等整合方案。"
    ]),
    ("E. 聚类与细胞注释", [
        "大群注释先扫 0.1、0.3、0.5、0.8，输出 resolution_sweep.tsv。",
        "结合 cluster 数量、marker coherence、组织背景和 cluster_marker_audit.tsv 选择粒度。",
        "亚群分析必须 subset 已注释大群后重新 HVG、scale、PCA、neighbors、clustering 和 UMAP。",
        "以 marker 证据为主，SingleR/SCINA/TransferData/scPred/CellTypist 为辅助。",
        "输出 annotation_evidence.tsv，保留 uncertain/ambiguous 标签。"
    ]),
    ("F. 组间比较", [
        "优先在细胞类型内比较 disease vs control；多样本时优先 pseudobulk 或尊重 donor 的模型。",
        "细胞比例比较必须按样本聚合，单细胞级比例仅描述性展示。",
        "保留 raw counts/unintegrated assay 用于表达层面的统计。"
    ]),
    ("G. 功能富集", [
        "对 marker/DE gene list 做 GO/KEGG/GSEA，记录数据库、物种、ID 转换和 gene universe。",
        "富集结果解释为功能线索，不等同通路活性证明。"
    ]),
    ("H. 细胞通讯", [
        "仅在注释可信且每群细胞数足够时运行 CellChat/NicheNet/CellPhoneDB。",
        "CellChat 结果是配体-受体数据库推断，必须标记为机制假说。"
    ]),
    ("I. 拟时序/细胞命运", [
        "仅对存在连续状态假设的细胞子集运行 Monocle3/Slingshot/PAGA。",
        "记录 root choice、subset rationale 和 sensitivity limits。"
    ]),
    ("J. 可选扩展：RNA velocity / 虚拟敲除", [
        "RNA velocity 需要 spliced/unspliced counts；普通表达矩阵不能运行。",
        "虚拟敲除/扰动预测作为后续扩展，第一版仅给出假说生成和外部工具建议。"
    ]),
    ("K. 预期结果文件", [
        "data_input_manifest.json 和 data_analysis_qc.md 必须记录并确认分析输入。",
        "figures/: QC、UMAP、marker、annotation、DE/enrichment、CellChat、pseudotime 图。",
        "tables/: metadata、qc_summary、workflow_step_audit、resolution_sweep、cluster_marker_audit、markers、deg、enrichment、cell_proportion、module_status。",
        "objects/: processed Seurat RDS 或 h5ad；logs/: commands、sessionInfo、package versions。"
    ]),
    ("L. 风险与限制", [
        "缺少数据同步质控时，不能确认结果来自计划下载或用户指定的数据。",
        "缺少 metadata 时无法可靠做 disease/control 比较。",
        "样本量过小或无生物学重复时，统计结论只能作为探索性描述。",
        "公共数据再分析不能替代临床诊断或治疗决策。"
    ]),
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnosis-json", help="Output from diagnose_geo_inputs.py")
    parser.add_argument("--question", default="")
    parser.add_argument("--organism", default="unknown")
    parser.add_argument("--comparison", default="unknown")
    parser.add_argument("--out-md", default="analysis_plan.md")
    args = parser.parse_args()

    diagnosis = {}
    if args.diagnosis_json and Path(args.diagnosis_json).exists():
        diagnosis = json.loads(Path(args.diagnosis_json).read_text(encoding="utf-8"))

    lines = [
        "# Single-cell Analysis Plan",
        "",
        f"- Research question: {args.question or 'not provided'}",
        f"- Organism: {args.organism}",
        f"- Comparison/design: {args.comparison}",
        f"- Detected formats: {', '.join(diagnosis.get('detected_formats', ['not diagnosed']))}",
        f"- Direct downstream analysis: {diagnosis.get('direct_downstream_analysis', 'unknown')}",
        "",
    ]
    for title, bullets in SECTIONS:
        lines.append(f"## {title}")
        lines.extend(f"- {x}" for x in bullets)
        lines.append("")
    Path(args.out_md).write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Check whether uploaded single-cell result files are interpretable."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


EXPECTED = {
    "metadata": ["metadata", "meta"],
    "markers": ["marker", "findallmarkers"],
    "deg": ["deg", "differential"],
    "enrichment": ["enrich", "kegg", "go_"],
    "cellchat": ["cellchat", "communication", "ligand", "receptor"],
    "pseudotime": ["pseudotime", "monocle", "trajectory"],
    "umap": ["umap"],
}


def has_columns(path: Path, columns: list[str]) -> bool:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            sample = handle.read(4096)
            handle.seek(0)
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
            reader = csv.DictReader(handle, dialect=dialect)
            header = [x.lower() for x in (reader.fieldnames or [])]
            return any(col in header for col in columns)
    except Exception:
        return False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", required=True)
    parser.add_argument("--out-md", default="result_quality_check.md")
    args = parser.parse_args()

    root = Path(args.result_dir)
    files = [p for p in root.rglob("*") if p.is_file()]
    lower = {p: p.name.lower() for p in files}
    found = {
        key: [str(p) for p, name in lower.items() if any(token in name for token in tokens)]
        for key, tokens in EXPECTED.items()
    }

    risks = []
    metadata_files = [Path(x) for x in found["metadata"]]
    if not metadata_files:
        risks.append("No metadata table detected; group/batch/sample-level interpretation is limited.")
    elif not any(has_columns(p, ["sample_id", "sample", "condition", "group", "batch", "donor_id"]) for p in metadata_files):
        risks.append("Metadata exists but common sample/group/batch columns were not detected.")
    if not found["markers"]:
        risks.append("No marker table detected; annotation quality cannot be checked.")
    if found["cellchat"] and not found["markers"]:
        risks.append("CellChat results without marker/annotation evidence are prone to over-interpretation.")
    if found["pseudotime"] and not found["umap"]:
        risks.append("Pseudotime results should be reviewed with embedding/trajectory plots.")

    quality = "基本可靠" if not risks else "需谨慎"
    if not files:
        quality = "不建议解读"
        risks.append("No files found in result directory.")

    lines = ["# Result Quality Check", "", f"总体质量判断：{quality}", ""]
    lines.append("## Detected Files")
    for key, paths in found.items():
        lines.append(f"- {key}: {len(paths)}")
    lines.append("")
    lines.append("## 主要风险")
    lines.extend(f"- {x}" for x in risks) if risks else lines.append("- 未发现明显结构性缺失；仍需人工复核生物学合理性。")
    lines.append("")
    lines.append("## 是否适合写入论文")
    lines.append("需要同时具备 metadata、QC、annotation evidence、statistical method 和可追溯图表后，才建议写入正式论文结果。")
    Path(args.out_md).write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()

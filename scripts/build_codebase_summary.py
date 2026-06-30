#!/usr/bin/env python3
"""Summarize course-derived R scripts for this skill."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


MODULES = {
    "GEO download": ["getGEO(", "getGEOSuppFiles", "prefetch", "fasterq-dump"],
    "10x MEX import": ["Read10X", "matrix.mtx"],
    "10x H5 import": ["Read10X_h5"],
    "h5ad import": ["h5ad", "zellkonverter", "anndata"],
    "Seurat object import": ["readRDS", "qread", ".qs"],
    "CSV/TSV import": ["read.csv", "read.delim", "read.table"],
    "QC": ["percent.mt", "VlnPlot", "FeatureScatter"],
    "normalization": ["NormalizeData", "SCTransform"],
    "batch correction": ["RunHarmony", "FindIntegrationAnchors"],
    "clustering": ["FindClusters", "FindNeighbors"],
    "UMAP/tSNE": ["RunUMAP", "RunTSNE"],
    "marker detection": ["FindAllMarkers", "wilcoxauc", "cosg"],
    "cell annotation": ["SingleR", "SCINA", "TransferData", "scPred"],
    "differential expression": ["FindMarkers", "FindAllMarkers"],
    "enrichment": ["enrichGO", "enrichKEGG", "clusterProfiler"],
    "cell proportion": ["prop.table", "table("],
    "CellChat": ["createCellChat", "computeCommunProb"],
    "Monocle": ["monocle3", "learn_graph", "order_cells", "monocle"],
    "CNV": ["copykat", "infercnv"],
    "coexpression": ["hdWGCNA", "WGCNA"],
    "deconvolution": ["CIBERSORT", "MuSiC", "music_prop"],
    "scVelo/CellRank": ["scvelo", "cellrank"],
    "virtual perturbation": ["knockout", "perturb", "CellOracle"],
    "report generation": ["rmarkdown", "quarto", "report"],
}


def read_manifest(root: Path) -> dict[str, dict[str, str]]:
    manifest = root / "source_manifest.csv"
    if not manifest.exists():
        return {}
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        return {row["english_file"]: row for row in csv.DictReader(handle)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--course-root", required=True, help="Directory containing English course-derived R scripts.")
    parser.add_argument("--out", default="CODEBASE_SUMMARY.md")
    args = parser.parse_args()

    root = Path(args.course_root)
    files = sorted([p for p in root.glob("*.R")], key=lambda p: p.name)
    combined = "\n".join(p.read_text(encoding="utf-8", errors="ignore") for p in files)
    packages = sorted(set(re.findall(r"library\(([^),]+)", combined)))
    manifest = read_manifest(root)

    display_root = str(root).replace("\\", "/")
    lines = [
        "# CODEBASE_SUMMARY",
        "",
        f"Course-derived script root: `{display_root}`",
        "",
        "## Main Languages",
        "- R",
        "- Python helper scripts",
        "- Markdown skill instructions",
        "",
        "## Entry Scripts",
    ]

    for p in files:
        text = p.read_text(encoding="utf-8", errors="ignore")
        calls = ", ".join(sorted({token for tokens in MODULES.values() for token in tokens if token in text})[:8])
        row = manifest.get(p.name, {})
        role = row.get("role", "course-derived reference script")
        manifest_note = "; original path listed in `scripts/course_source/source_manifest.csv`" if row else ""
        lines.append(f"- `{p.name}`: {role}; {calls or 'inspect before adaptation'}{manifest_note}")

    lines.append("")
    lines.append("## Dependencies Detected")
    cleaned_packages = [pkg.strip().strip('"').strip("'") for pkg in packages]
    lines.extend(f"- {pkg}" for pkg in cleaned_packages if pkg)

    lines.append("")
    lines.append("## Module Coverage")
    for module, tokens in MODULES.items():
        status = "present" if any(token in combined for token in tokens) else "missing"
        lines.append(f"- {module}: `{status}`")

    lines.append("")
    lines.append("## Notes")
    lines.append("- Course source scripts are renamed to English and flattened under `scripts/course_source/`.")
    lines.append("- The manifest `scripts/course_source/source_manifest.csv` maps English filenames to original course paths.")
    lines.append("- Reference scripts still require review for project-specific object names, thresholds, metadata columns, and filenames.")
    lines.append("- Runnable project workflows should use `scripts/course_adapted/`.")
    lines.append("- Missing modules should be documented as future extensions, not implemented claims.")

    Path(args.out).write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Propose CellChat and pseudotime scopes from existing result tables.

This script does not run CellChat or Monocle3. It writes an approval-gated
proposal that the user must review before downstream scripts are rendered.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


IMMUNE = {"t cell", "cd4 t", "cd8 t", "nk", "b cell", "plasma", "monocyte", "macrophage", "dendritic", "mast"}
STROMAL = {"fibroblast", "endothelial", "pericyte", "smooth muscle", "stromal"}
EPITHELIAL = {"epithelial", "keratinocyte", "basal", "club", "alveolar"}


def find_files(root: Path, tokens: list[str]) -> list[Path]:
    hits = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        name = path.name.lower()
        if any(token in name for token in tokens):
            hits.append(path)
    return sorted(hits)


def read_table(path: Path, limit: int = 20000) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            sample = handle.read(4096)
            handle.seek(0)
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
            reader = csv.DictReader(handle, dialect=dialect)
            rows = []
            for i, row in enumerate(reader):
                if i >= limit:
                    break
                rows.append({(k or "").strip(): (v or "").strip() for k, v in row.items()})
            return rows
    except Exception:
        return []


def norm_label(value: str) -> str:
    return " ".join(value.replace("_", " ").replace("-", " ").lower().split())


def collect_cell_labels(root: Path) -> tuple[Counter[str], list[str]]:
    labels: Counter[str] = Counter()
    warnings: list[str] = []
    evidence_files = find_files(root, ["annotation_evidence"])
    singler_files = find_files(root, ["singler_cluster", "singleR_cluster".lower()])
    proportion_files = find_files(root, ["celltype_proportion", "cell_type_proportion", "proportion"])

    for path in evidence_files:
        for row in read_table(path):
            label = row.get("final_label") or row.get("coarse_label") or row.get("singleR_label") or ""
            confidence = norm_label(row.get("confidence", ""))
            count = row.get("cell_count", "1")
            try:
                n = max(1, int(float(count)))
            except ValueError:
                n = 1
            if label:
                labels[label] += n
            if confidence == "low":
                warnings.append(f"Low-confidence annotation in {path.name}: cluster {row.get('cluster', '?')} {label}")

    if not labels:
        for path in singler_files:
            for row in read_table(path):
                label = row.get("labels") or row.get("pruned.labels") or row.get("singleR_label") or ""
                if label:
                    labels[label] += 1
        if singler_files:
            warnings.append("SingleR labels detected but annotation_evidence.tsv was not found.")

    if not labels:
        for path in proportion_files:
            for row in read_table(path):
                label = row.get("cell_type") or row.get("celltype") or row.get("label") or ""
                if label:
                    labels[label] += 1

    if len(labels) <= 1:
        warnings.append("Only one or zero cell labels detected; CellChat and multi-state pseudotime are usually not appropriate.")
    return labels, warnings


def collect_marker_hints(root: Path) -> dict[str, list[str]]:
    marker_files = find_files(root, ["marker"])
    hints: dict[str, list[str]] = defaultdict(list)
    for path in marker_files[:3]:
        for row in read_table(path, limit=5000):
            cluster = row.get("cluster") or row.get("seurat_clusters") or row.get("group") or "unknown"
            gene = row.get("gene") or row.get("features") or row.get("feature") or row.get("Gene") or ""
            if gene and len(hints[cluster]) < 10:
                hints[cluster].append(gene)
    return dict(hints)


def label_family(label: str) -> str:
    value = norm_label(label)
    if any(token in value for token in IMMUNE):
        return "immune"
    if any(token in value for token in STROMAL):
        return "stromal"
    if any(token in value for token in EPITHELIAL):
        return "epithelial"
    return "other"


def propose_cellchat(labels: Counter[str]) -> list[str]:
    families = {label_family(label) for label in labels}
    selected = [label for label, _ in labels.most_common() if not norm_label(label).startswith(("unknown", "ambiguous"))]
    if len(selected) < 2:
        return ["CellChat is not recommended until at least two credible cell groups are confirmed."]
    if {"immune", "stromal"} <= families:
        scope = "immune-stromal microenvironment"
    elif {"immune", "epithelial"} <= families:
        scope = "immune-epithelial injury microenvironment"
    elif {"stromal", "epithelial"} <= families:
        scope = "epithelial-stromal remodeling microenvironment"
    else:
        scope = "confirmed major cell groups"
    return [
        f"Recommended scope: {scope}.",
        "Candidate groups: " + ", ".join(selected[:8]),
        "Exclude groups labeled Unknown/Ambiguous or groups below the minimum cell threshold.",
        "Ask the user which comparison to run, for example disease vs control or stage-specific networks.",
    ]


def propose_pseudotime(labels: Counter[str]) -> list[str]:
    selected = [label for label, _ in labels.most_common() if not norm_label(label).startswith(("unknown", "ambiguous"))]
    if len(selected) < 2:
        return ["Pseudotime is not recommended until a plausible multi-state lineage is confirmed."]
    families = defaultdict(list)
    for label in selected:
        families[label_family(label)].append(label)
    for family in ("epithelial", "stromal", "immune"):
        if len(families[family]) >= 2:
            root_hint = families[family][0]
            return [
                f"Recommended lineage: {family} states.",
                "Candidate states: " + ", ".join(families[family][:8]),
                f"Proposed root: {root_hint}, pending marker/stage review.",
                "Ask the user to confirm root cells/root state before Monocle3 is rendered.",
            ]
    return [
        "Pseudotime may be inappropriate because detected labels are broad unrelated lineages.",
        "Ask the user whether there is a known differentiation, activation, or time-course hypothesis.",
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", required=True)
    parser.add_argument("--out-md", default="downstream_proposal.md")
    args = parser.parse_args()

    root = Path(args.result_dir)
    labels, warnings = collect_cell_labels(root)
    marker_hints = collect_marker_hints(root)

    lines = [
        "# Downstream Module Proposal",
        "",
        "This proposal is an approval gate. Do not run CellChat or Monocle3 until the user confirms the selected module and scope.",
        "",
        "## Detected Cell Labels",
    ]
    if labels:
        for label, count in labels.most_common():
            lines.append(f"- {label}: {count}")
    else:
        lines.append("- No usable cell labels detected.")

    lines.extend(["", "## Annotation Warnings"])
    lines.extend(f"- {warning}" for warning in warnings) if warnings else lines.append("- No blocking annotation warning detected from available tables.")

    lines.extend(["", "## Marker Hints"])
    if marker_hints:
        for cluster, genes in marker_hints.items():
            lines.append(f"- {cluster}: {', '.join(genes)}")
    else:
        lines.append("- No marker table detected; downstream modules require manual annotation review.")

    lines.extend(["", "## CellChat Candidate Plan"])
    lines.extend(f"- {item}" for item in propose_cellchat(labels))
    lines.extend([
        "- Required user confirmation: run CellChat? Which sender/receiver groups and comparison should be used?",
        "",
        "## Pseudotime Candidate Plan",
    ])
    lines.extend(f"- {item}" for item in propose_pseudotime(labels))
    lines.extend([
        "- Required user confirmation: run pseudotime? Which lineage/subset and root state should be used?",
        "",
        "## Approval Checklist",
        "- [ ] Annotation evidence reviewed.",
        "- [ ] CellChat scope approved or skipped.",
        "- [ ] Pseudotime lineage/root approved or skipped.",
        "- [ ] Low-confidence/ambiguous groups excluded or manually accepted.",
    ])

    Path(args.out_md).write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

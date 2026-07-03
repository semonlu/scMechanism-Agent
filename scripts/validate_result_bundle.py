#!/usr/bin/env python3
"""Check whether uploaded single-cell result files are interpretable."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


EXPECTED = {
    "metadata": ["metadata", "meta"],
    "data_manifest": ["data_input_manifest"],
    "data_qc": ["data_analysis_qc"],
    "markers": ["marker", "findallmarkers"],
    "annotation": ["annotation_evidence", "singler", "celltypist", "cell_labels"],
    "deg": ["deg", "differential"],
    "enrichment": ["enrich", "kegg", "go_"],
    "cellchat": ["cellchat", "communication", "ligand", "receptor"],
    "pseudotime": ["pseudotime", "monocle", "trajectory"],
    "umap": ["umap"],
}


def sniff_rows(path: Path, limit: int = 5000) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            sample = handle.read(4096)
            handle.seek(0)
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
            reader = csv.DictReader(handle, dialect=dialect)
            header = [x or "" for x in (reader.fieldnames or [])]
            rows = []
            for i, row in enumerate(reader):
                if i >= limit:
                    break
                rows.append({(k or "").strip(): (v or "").strip() for k, v in row.items()})
            return header, rows
    except Exception:
        return [], []


def has_columns(path: Path, columns: list[str]) -> bool:
    header, _ = sniff_rows(path, limit=1)
    lower = [x.lower() for x in header]
    return any(col in lower for col in columns)


def annotation_status(annotation_files: list[Path], marker_files: list[Path]) -> tuple[str, list[str], list[str]]:
    risks: list[str] = []
    support: list[str] = []
    if not marker_files:
        risks.append("No marker table detected; annotation cannot be reviewed against marker evidence.")
    if not annotation_files:
        risks.append("No annotation evidence or automatic label table detected.")
        return "not_usable", support, risks

    evidence = [p for p in annotation_files if "annotation_evidence" in p.name.lower()]
    if evidence:
        labels: Counter[str] = Counter()
        low_confidence = 0
        for path in evidence:
            _, rows = sniff_rows(path)
            for row in rows:
                label = row.get("final_label") or row.get("coarse_label") or row.get("singleR_label") or ""
                if label:
                    labels[label] += 1
                if (row.get("confidence") or "").lower() == "low":
                    low_confidence += 1
        if labels:
            support.append(f"Annotation evidence detected with {len(labels)} unique labels.")
        if len(labels) <= 1:
            risks.append("Only one annotation label detected; verify that the dataset is truly single-lineage.")
        if low_confidence:
            risks.append(f"{low_confidence} low-confidence annotation entries require manual review.")
        if marker_files and len(labels) > 1 and not low_confidence:
            return "usable", support, risks
        return "needs_review", support, risks

    risks.append("Automatic annotation table detected, but annotation_evidence.tsv is missing.")
    if marker_files:
        return "needs_review", support, risks
    return "not_usable", support, risks


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", required=True)
    parser.add_argument("--out-md", default="result_quality_check.md")
    args = parser.parse_args()

    root = Path(args.result_dir)
    files = [p for p in root.rglob("*") if p.is_file()]
    lower = {p: p.name.lower() for p in files}
    found = {
        key: [p for p, name in lower.items() if any(token in name for token in tokens)]
        for key, tokens in EXPECTED.items()
    }

    risks: list[str] = []
    support: list[str] = []
    if not found["data_manifest"]:
        risks.append("No data_input_manifest.json detected; analysis input provenance is incomplete.")
    else:
        support.append("data_input_manifest.json detected.")
    if not found["data_qc"]:
        risks.append("No data_analysis_qc.md detected; cannot confirm that analyzed input matches the planned/downloaded/registered dataset.")
    else:
        qc_text = "\n".join(p.read_text(encoding="utf-8", errors="replace")[:4000] for p in found["data_qc"])
        if "Status: error" in qc_text or "does not match" in qc_text:
            risks.append("data_analysis_qc.md indicates an input mismatch; results should not be interpreted.")
        else:
            support.append("data_analysis_qc.md detected without an obvious input mismatch flag.")
    metadata_files = found["metadata"]
    if not metadata_files:
        risks.append("No metadata table detected; group/batch/sample-level interpretation is limited.")
    elif not any(has_columns(p, ["sample_id", "sample", "condition", "group", "batch", "donor_id"]) for p in metadata_files):
        risks.append("Metadata exists but common sample/group/batch columns were not detected.")
    else:
        support.append("Metadata with common sample/group/batch columns was detected.")

    if not found["markers"]:
        risks.append("No marker table detected; annotation quality cannot be checked.")
    else:
        support.append(f"Marker-related files detected: {len(found['markers'])}.")

    ann_status, ann_support, ann_risks = annotation_status(found["annotation"], found["markers"])
    support.extend(ann_support)
    risks.extend(ann_risks)

    if found["cellchat"] and ann_status == "not_usable":
        risks.append("CellChat results exist but annotation is not usable; do not interpret communication results yet.")
    if found["pseudotime"] and ann_status == "not_usable":
        risks.append("Pseudotime results exist but annotation is not usable; root/subset choices require review.")
    if found["pseudotime"] and not found["umap"]:
        risks.append("Pseudotime results should be reviewed with embedding/trajectory plots.")

    if not files:
        quality = "not_interpretable"
        risks.append("No files found in result directory.")
    elif ann_status == "not_usable" or len(risks) >= 4:
        quality = "caution_required"
    elif risks:
        quality = "basically_reliable_with_review"
    else:
        quality = "reliable_for_exploratory_reporting"

    downstream_allowed = ann_status in {"usable", "needs_review"}

    lines = [
        "# Result Quality Check",
        "",
        f"- Overall quality: {quality}",
        f"- Data synchronization status: {'present' if found['data_qc'] else 'missing'}",
        f"- Cell annotation status: {ann_status}",
        f"- Allow CellChat/pseudotime proposal: {'yes' if downstream_allowed else 'no'}",
        "",
        "## Detected Files",
    ]
    for key, paths in found.items():
        lines.append(f"- {key}: {len(paths)}")

    lines.extend(["", "## Main Support"])
    lines.extend(f"- {x}" for x in support) if support else lines.append("- No strong structural support detected yet.")

    lines.extend(["", "## Main Risks"])
    lines.extend(f"- {x}" for x in risks) if risks else lines.append("- No obvious structural blocker detected; biological review is still required.")

    lines.extend([
        "",
        "## Next Step",
        "Before running CellChat or pseudotime, generate `downstream_proposal.md` and ask the user to approve the module and scope.",
        "",
        "## Manuscript Suitability",
        "Use in manuscript text only after metadata, QC, annotation evidence, statistical method, and traceable figures/tables are available.",
    ])
    Path(args.out_md).write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

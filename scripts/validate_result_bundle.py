#!/usr/bin/env python3
"""Check whether uploaded single-cell result files are interpretable."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


EXPECTED = {
    "workflow": ["workflow_step_audit", "pipeline_step_audit"],
    "input": ["input_diagnosis", "data_input_manifest", "file_manifest", "import_summary", "source_traceability"],
    "qc": ["qc"],
    "metadata": ["metadata", "meta"],
    "data_manifest": ["data_input_manifest"],
    "data_qc": ["data_analysis_qc"],
    "normalization": ["normalization", "highly_variable", "hvg", "sctransform", "processed_seurat", "processed.h5ad"],
    "markers": ["marker", "findallmarkers"],
    "resolution": ["resolution_sweep"],
    "cluster_audit": ["cluster_marker_audit"],
    "doublet": ["doublet", "scrublet", "scdblfinder"],
    "cell_cycle": ["cell_cycle", "s_score", "g2m"],
    "batch": ["batch", "harmony", "integration", "combat", "bbknn"],
    "annotation": ["annotation_evidence", "singler", "celltypist", "cell_labels"],
    "deg": ["deg", "differential"],
    "enrichment": ["enrich", "kegg", "go_"],
    "cellchat": ["cellchat", "communication", "ligand", "receptor"],
    "pseudotime": ["pseudotime", "monocle", "trajectory"],
    "umap": ["umap"],
    "pca": ["pca", "elbow"],
    "tsne": ["tsne"],
}


WORKFLOW_STEPS = [
    (1, "data_acquisition", "Get data from GEO/SRA or a traceable local/public source", ["input", "data_manifest"], True),
    (2, "data_import", "Read data and align metadata/barcodes", ["input", "metadata", "data_qc"], True),
    (3, "quality_control", "QC mitochondrial/red-cell/count/gene/low-quality cells", ["qc"], True),
    (4, "normalization_scaling_hvg", "Normalize/scale/find variable features or SCT while retaining raw counts", ["normalization"], True),
    (5, "confounder_review", "Review cell cycle, doublets, batch, ambient RNA, and integration risk", ["doublet", "cell_cycle", "batch"], True),
    (6, "dimensionality_clustering", "Run/audit PCA, UMAP/tSNE, neighbors, clustering, and resolution", ["pca", "umap", "resolution", "cluster_audit"], True),
    (7, "cell_annotation", "Cross-check automatic and manual cell annotation", ["annotation"], True),
    (8, "marker_detection", "Find and review marker genes", ["markers", "cluster_audit"], True),
    (9, "pseudotime", "Optional lineage/root-approved pseudotime analysis", ["pseudotime"], False),
    (10, "cell_communication", "Optional annotation-gated CellChat/cell communication analysis", ["cellchat"], False),
]


def workflow_step_audit(found: dict[str, list[Path]]) -> tuple[list[dict[str, str]], list[str], list[str]]:
    risks: list[str] = []
    support: list[str] = []
    rows: list[dict[str, str]] = []

    for step_id, key, description, evidence_keys, core_required in WORKFLOW_STEPS:
        present_keys = [name for name in evidence_keys if found.get(name)]
        evidence = []
        for name in present_keys:
            evidence.extend(path.name for path in found[name][:3])
        evidence = list(dict.fromkeys(evidence))
        if core_required:
            if present_keys:
                status = "complete"
                note = "Evidence detected for this core upstream step."
                support.append(f"Workflow step {step_id} {key}: evidence detected.")
            else:
                status = "missing_evidence"
                note = "Core upstream step lacks auditable files."
                risks.append(f"Workflow step {step_id} ({key}) lacks auditable evidence.")
        else:
            if present_keys:
                status = "present_requires_review"
                note = "Optional downstream module is present; verify approval, scope, and inference-only wording."
                support.append(f"Workflow step {step_id} {key}: downstream evidence detected.")
            else:
                status = "gated_optional_not_run"
                note = "Optional downstream module not detected; acceptable if not approved or not needed."
        rows.append(
            {
                "step": str(step_id),
                "workflow_step": key,
                "required_before_interpretation": "yes" if core_required else "no_optional_gated",
                "status": status,
                "description": description,
                "evidence_files": ";".join(evidence),
                "review_note": note,
            }
        )
    return rows, support, risks


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
    workflow_rows, workflow_support, workflow_risks = workflow_step_audit(found)
    support.extend(workflow_support)
    risks.extend(workflow_risks)
    workflow_audit_path = root / "tables" / "workflow_step_audit.tsv"
    workflow_audit_path.parent.mkdir(parents=True, exist_ok=True)
    with workflow_audit_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(workflow_rows[0]))
        writer.writeheader()
        writer.writerows(workflow_rows)
    found["workflow"] = list(dict.fromkeys([*found.get("workflow", []), workflow_audit_path]))

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
        "## 10-Step Workflow Audit",
        "",
        "| Step | Workflow step | Required before interpretation | Status | Evidence | Review note |",
        "|---:|---|---|---|---|---|",
    ]
    for row in workflow_rows:
        lines.append(
            f"| {row['step']} | {row['workflow_step']} | {row['required_before_interpretation']} | "
            f"{row['status']} | {row['evidence_files'] or 'none'} | {row['review_note']} |"
        )

    lines.extend([
        "",
        "## Detected Files",
    ])
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

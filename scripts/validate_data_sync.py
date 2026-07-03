#!/usr/bin/env python3
"""Validate that the analysis input matches the planned/downloaded dataset."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def resolve_path(value: str, base: Path) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = base / path
    return path.resolve(strict=False)


def is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def load_manifest(path: Path) -> dict:
    if not path.exists():
        return {"entries": []}
    return json.loads(path.read_text(encoding="utf-8"))


def entry_path(entry: dict, base: Path) -> Path | None:
    rel = (entry.get("input") or {}).get("path")
    if not rel:
        return None
    return resolve_path(rel, base)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", required=True, help="Project result directory.")
    parser.add_argument("--input-path", required=True, help="Actual input used for analysis.")
    parser.add_argument("--input-type", default="")
    parser.add_argument("--metadata-path", default="")
    parser.add_argument("--manifest", default="", help="data_input_manifest.json; defaults to result-dir/data_input_manifest.json")
    parser.add_argument("--workspace-root", default="", help="Workspace root for relative manifest paths.")
    parser.add_argument("--out-md", default="data_analysis_qc.md")
    parser.add_argument("--out-json", default="")
    args = parser.parse_args()

    result_dir = Path(args.result_dir).resolve()
    workspace_root = Path(args.workspace_root).resolve() if args.workspace_root else result_dir.parent.parent
    manifest_path = Path(args.manifest).resolve() if args.manifest else result_dir / "data_input_manifest.json"
    input_path = resolve_path(args.input_path, workspace_root)
    metadata_path = resolve_path(args.metadata_path, workspace_root) if args.metadata_path else None
    manifest = load_manifest(manifest_path)
    entries = manifest.get("entries", [])

    matched = []
    for entry in entries:
        registered = entry_path(entry, workspace_root)
        if registered is None:
            continue
        if registered.is_dir() and is_relative_to(input_path, registered):
            matched.append(entry)
        elif registered == input_path:
            matched.append(entry)

    risks: list[str] = []
    support: list[str] = []
    if not input_path.exists():
        risks.append("Analysis input_path does not exist.")
    if metadata_path and not metadata_path.exists():
        risks.append("metadata_path was provided but does not exist.")
    if not entries:
        risks.append("No data_input_manifest.json entries found.")
    elif not matched:
        risks.append("Analysis input_path does not match any registered/downloaded input entry.")
    else:
        support.append(f"Matched {len(matched)} registered/downloaded input entry.")

    status = "ok" if matched and input_path.exists() and not any("does not match" in r for r in risks) else ("warning" if not entries else "error")
    payload = {
        "status": status,
        "input_path": str(input_path),
        "input_type": args.input_type,
        "metadata_path": str(metadata_path) if metadata_path else "",
        "manifest": str(manifest_path),
        "matched_entries": matched,
        "risks": risks,
        "support": support,
    }

    lines = [
        "# Data Analysis QC",
        "",
        f"- Status: {status}",
        f"- Input path: {input_path}",
        f"- Input type: {args.input_type or 'not specified'}",
        f"- Metadata path: {metadata_path if metadata_path else 'not provided'}",
        f"- Manifest: {manifest_path}",
        f"- Manifest entries: {len(entries)}",
        f"- Matched entries: {len(matched)}",
        "",
        "## Support",
        *(f"- {item}" for item in support),
        *([] if support else ["- No manifest-backed input match detected."]),
        "",
        "## Risks",
        *(f"- {item}" for item in risks),
        *([] if risks else ["- No structural input synchronization risk detected."]),
        "",
        "## Rule",
        "The analyzed input must match the dataset selected during planning, GEO download, archive extraction, or manual input registration.",
    ]
    Path(args.out_md).write_text("\n".join(lines) + "\n", encoding="utf-8")
    if args.out_json:
        Path(args.out_json).write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

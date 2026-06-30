#!/usr/bin/env python3
"""Validate full scMechanism workflow deliverables and known regression traps."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def fail(message: str, errors: list[str]) -> None:
    errors.append(message)


def count_files(root: Path, pattern: str) -> int:
    return sum(1 for _ in root.rglob(pattern)) if root.exists() else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--example-root", type=Path)
    parser.add_argument("--out-md", type=Path)
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    example_root = (args.example_root or project_root / "submission" / "example").resolve()
    result_dir = example_root / "results"
    errors: list[str] = []
    warnings: list[str] = []

    required_files = [
        project_root / "scripts" / "env_setup" / "install_environment.ps1",
        project_root / "scripts" / "env_setup" / "check_environment.ps1",
        project_root / "references" / "environment" / "requirements.md",
        project_root / "references" / "environment" / "path-setup.md",
        project_root / "references" / "tested-lessons.md",
        project_root / "scripts" / "write_analysis_report.py",
        result_dir / "course_modules" / "monocle3" / "figures" / "trajectory_celltypes.pdf",
        result_dir / "course_modules" / "monocle3" / "figures" / "trajectory_celltypes_with_graph.pdf",
        result_dir / "course_modules" / "monocle3" / "figures" / "pseudotime_umap.pdf",
        result_dir / "course_modules" / "monocle3" / "tables" / "pseudotime.csv",
        result_dir / "course_modules" / "monocle3" / "tables" / "trajectory_input_state_counts.csv",
        example_root / "manuscript_report.md",
        example_root / "full_workflow_log.md",
        example_root / "logs" / "environment_check_20260630_203000" / "environment_check_summary.md",
        example_root / "logs" / "environment_check_20260630_203000" / "r_packages.tsv",
        example_root / "logs" / "environment_check_20260630_203000" / "python_packages.tsv",
    ]
    for path in required_files:
        if not path.exists():
            fail(f"Missing required artifact: {path}", errors)

    pseudotime = read_csv(result_dir / "course_modules" / "monocle3" / "tables" / "pseudotime.csv")
    if pseudotime:
        state_field = "marker_support_label" if "marker_support_label" in pseudotime[0] else "singleR_label"
        states = {row.get(state_field, "") for row in pseudotime if row.get(state_field, "")}
        if len(states) < 2:
            fail(f"Monocle3 pseudotime has fewer than 2 states in {state_field}: {sorted(states)}", errors)
        if len(pseudotime) < 100:
            fail(f"Monocle3 pseudotime has suspiciously few cells: {len(pseudotime)}", errors)
        finite = [row for row in pseudotime if re.fullmatch(r"-?\d+(\.\d+)?([eE][+-]?\d+)?", row.get("pseudotime", ""))]
        if len(finite) != len(pseudotime):
            fail(f"Monocle3 pseudotime contains non-finite or missing values: {len(finite)}/{len(pseudotime)} finite", errors)
    else:
        fail("Missing or unreadable Monocle3 pseudotime table", errors)

    report_text = (example_root / "manuscript_report.md").read_text(encoding="utf-8", errors="ignore") if (example_root / "manuscript_report.md").exists() else ""
    for heading in ["## Methods", "## Results", "## Figure Legends", "## Supplementary Tables", "## Limitations"]:
        if heading not in report_text:
            fail(f"Report is missing heading: {heading}", errors)

    rendered_r = example_root / "rendered_course_modules" / "04_monocle3_from_seurat.R"
    if rendered_r.exists():
        rendered_text = rendered_r.read_text(encoding="utf-8", errors="ignore")
        if "\\单" in rendered_text or re.search(r"[A-Z]:\\", rendered_text):
            fail("Rendered Monocle3 R script contains Windows backslash paths", errors)
        for guard in ["label_principal_points = FALSE", "length(celltype_counts) < 2"]:
            if guard not in rendered_text:
                fail(f"Rendered Monocle3 script is missing guard: {guard}", errors)

    bad_artifacts = [
        result_dir / "course_modules" / "monocle3" / "root_cells.txt",
        result_dir / "course_modules" / "monocle3" / "sessionInfo_monocle3.txt",
        example_root / "pdf_inspection",
    ]
    for path in bad_artifacts:
        if path.exists():
            fail(f"Known-bad or scratch artifact still exists: {path}", errors)

    if count_files(result_dir, "*.rds"):
        warnings.append("RDS objects exist in results; confirm they are gitignored before publishing.")

    lines = ["# Full Workflow Validation", ""]
    lines.append("## Status")
    lines.append("")
    lines.append("FAILED" if errors else "PASSED")
    lines.append("")
    lines.append("## Errors")
    lines.append("")
    lines.extend(f"- {item}" for item in errors) if errors else lines.append("- None")
    lines.append("")
    lines.append("## Warnings")
    lines.append("")
    lines.extend(f"- {item}" for item in warnings) if warnings else lines.append("- None")
    lines.append("")

    if args.out_md:
        args.out_md.parent.mkdir(parents=True, exist_ok=True)
        args.out_md.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())

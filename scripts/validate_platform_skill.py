#!/usr/bin/env python3
"""Validate the skill against the Medical AI Skill platform upload checklist."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


REQUIRED_PATHS = [
    "SKILL.md",
    "README.md",
    "SECURITY_AND_PRIVACY.md",
    "agents",
    "references",
    "examples",
    "examples/validation_input_output_comparison.md",
    "submission/platform",
    "submission/platform/skill_brief.md",
    "submission/platform/clinical_pain_point.md",
    "submission/platform/implementation_plan.md",
    "submission/platform/business_value.md",
    "submission/platform/demo_script.md",
    "submission/platform/validation_report.md",
    "submission/platform/github_link.md",
]

FORBIDDEN_TOP_LEVEL = {
    ".git",
    "seuratv5_environment_check",
    "__pycache__",
}

FORBIDDEN_PACKAGE_PARTS = {
    ".git",
    "__pycache__",
    "data",
    "objects",
    "seuratv5_environment_check",
}

FORBIDDEN_SUFFIXES = {
    ".rds",
    ".rda",
    ".h5ad",
    ".loom",
    ".mtx",
    ".qs",
}


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def parse_frontmatter(skill_md: Path) -> tuple[dict[str, str], list[str]]:
    errors: list[str] = []
    text = read_text(skill_md)
    if not text.startswith("---\n"):
        return {}, ["SKILL.md must start with YAML frontmatter"]
    end = text.find("\n---", 4)
    if end == -1:
        return {}, ["SKILL.md frontmatter is not closed"]
    block = text[4:end].strip().splitlines()
    meta: dict[str, str] = {}
    for line in block:
        if not line.strip():
            continue
        if ":" not in line:
            errors.append(f"Invalid frontmatter line: {line}")
            continue
        key, value = line.split(":", 1)
        meta[key.strip()] = value.strip().strip('"').strip("'")
    extra = sorted(set(meta) - {"name", "description"})
    if extra:
        errors.append(f"SKILL.md frontmatter may only contain name and description, found: {', '.join(extra)}")
    for key in ("name", "description"):
        if not meta.get(key):
            errors.append(f"SKILL.md frontmatter missing {key}")
    if meta.get("name") and not re.fullmatch(r"[a-z0-9-]{1,63}", meta["name"]):
        errors.append("Skill name must be lowercase letters, digits, and hyphens only, max 63 chars")
    return meta, errors


def count_validation_cases(report: Path) -> int:
    text = read_text(report)
    return len(re.findall(r"^## Case\s+\d+:", text, flags=re.MULTILINE))


def validate(skill_root: Path) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []

    for rel in REQUIRED_PATHS:
        if not (skill_root / rel).exists():
            errors.append(f"Missing required path: {rel}")

    if (skill_root / "SKILL.md").exists():
        _, fm_errors = parse_frontmatter(skill_root / "SKILL.md")
        errors.extend(fm_errors)

    agent_files = sorted((skill_root / "agents").glob("*.md")) if (skill_root / "agents").exists() else []
    if not agent_files:
        errors.append("agents/ must contain at least one .md role file")

    if not (skill_root / "agents" / "openai.yaml").exists():
        warnings.append("agents/openai.yaml is recommended for platform UI metadata")

    case_report = skill_root / "examples" / "validation_input_output_comparison.md"
    if case_report.exists():
        n_cases = count_validation_cases(case_report)
        if n_cases < 5:
            errors.append(f"Validation report must contain at least 5 cases, found {n_cases}")

    for name in FORBIDDEN_TOP_LEVEL:
        if (skill_root / name).exists() and name == "__pycache__":
            warnings.append(f"Top-level generated cache exists and should not be packaged: {name}")

    excluded_path_seen = False
    excluded_object_seen = False
    for path in skill_root.rglob("*"):
        if ".git" in path.parts:
            continue
        rel_parts = set(path.relative_to(skill_root).parts)
        if rel_parts & FORBIDDEN_PACKAGE_PARTS:
            excluded_path_seen = True
        if path.is_file() and path.suffix.lower() in FORBIDDEN_SUFFIXES:
            excluded_object_seen = True

    if excluded_path_seen:
        warnings.append("Local/generated folders such as data/, objects/, __pycache__/, and seuratv5_environment_check/ exist locally and must be excluded from the upload package")
    if excluded_object_seen:
        warnings.append("Large/private analysis object files exist locally and must be excluded from the upload package")

    package_script = skill_root / "scripts" / "package_platform_skill.py"
    if not package_script.exists():
        errors.append("Missing scripts/package_platform_skill.py for reproducible upload package generation")

    return sorted(set(errors)), sorted(set(warnings))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skill-root", default=".", help="Path to the skill root")
    args = parser.parse_args()

    skill_root = Path(args.skill_root).resolve()
    errors, warnings = validate(skill_root)

    print("# Platform Skill Validation")
    print()
    print(f"Skill root: {skill_root}")
    print(f"Status: {'PASSED' if not errors else 'FAILED'}")
    print()
    print("## Errors")
    if errors:
        for item in errors:
            print(f"- {item}")
    else:
        print("- None")
    print()
    print("## Warnings")
    if warnings:
        for item in warnings:
            print(f"- {item}")
    else:
        print("- None")

    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())

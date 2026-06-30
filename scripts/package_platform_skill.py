#!/usr/bin/env python3
"""Create a platform upload zip using an explicit skill-file allowlist."""

from __future__ import annotations

import argparse
import sys
import zipfile
from pathlib import Path


ALLOWLIST = [
    "SKILL.md",
    "agents",
    "references",
    "examples",
    "scripts",
    "templates",
]

EXCLUDE_DIRS = {
    ".git",
    "__pycache__",
    ".pytest_cache",
    ".venv",
    ".conda",
    "renv",
    "seuratv5_environment_check",
    "data",
    "objects",
}

EXCLUDE_SUFFIXES = {
    ".pyc",
    ".rds",
    ".rda",
    ".h5ad",
    ".loom",
    ".mtx",
    ".qs",
    ".log",
}

EXCLUDE_FILENAMES = {
    "README.md",
}


def should_include(path: Path, root: Path) -> bool:
    rel = path.relative_to(root)
    if any(part in EXCLUDE_DIRS for part in rel.parts):
        return False
    if path.is_file() and path.suffix.lower() in EXCLUDE_SUFFIXES:
        return False
    if path.is_file() and path.name in EXCLUDE_FILENAMES:
        return False
    if path.name.startswith(".") and path.name not in {".gitkeep"}:
        return False
    return True


def iter_package_files(root: Path):
    for item in ALLOWLIST:
        path = root / item
        if not path.exists():
            continue
        if path.is_file():
            if should_include(path, root):
                yield path
            continue
        for child in sorted(path.rglob("*")):
            if child.is_file() and should_include(child, root):
                yield child


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skill-root", default=".", help="Path to the skill root")
    parser.add_argument("--out", required=True, help="Output zip path")
    parser.add_argument("--package-name", default=None, help="Top-level folder name inside the zip")
    args = parser.parse_args()

    root = Path(args.skill_root).resolve()
    out = Path(args.out).resolve()
    package_name = args.package_name or root.name

    if not (root / "SKILL.md").exists():
        print(f"Missing SKILL.md under {root}", file=sys.stderr)
        return 1

    files = list(iter_package_files(root))
    if not files:
        print("No files selected for packaging", file=sys.stderr)
        return 1

    out.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for path in files:
            rel = path.relative_to(root).as_posix()
            zf.write(path, f"{package_name}/{rel}")

    print(f"Wrote {out}")
    print(f"Included files: {len(files)}")
    print("Excluded local analysis outputs, data folders, object files, caches, and Git metadata.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

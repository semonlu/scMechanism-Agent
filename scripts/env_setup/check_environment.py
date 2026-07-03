#!/usr/bin/env python3
"""Cross-platform environment checker for scMechanism-Agent."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import platform
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class EnvironmentProfile:
    name: str
    python_imports: list[str]
    r_packages: list[str]
    required_tools: list[str]


PYTHON_IMPORTS = [
    "scanpy",
    "anndata",
    "pandas",
    "numpy",
    "scipy",
    "sklearn",
    "matplotlib",
    "seaborn",
    "igraph",
    "leidenalg",
    "celltypist",
    "gseapy",
    "harmonypy",
    "bbknn",
    "scanorama",
    "liana",
]

MINIMAL_R_PACKAGES = [
    "Seurat",
    "SeuratObject",
    "Matrix",
    "ggplot2",
    "patchwork",
    "dplyr",
    "data.table",
    "tibble",
    "ggrepel",
    "SingleR",
    "celldex",
    "SingleCellExperiment",
    "SummarizedExperiment",
]

EXTENDED_R_PACKAGES = [
    "clusterProfiler",
    "org.Mm.eg.db",
    "org.Hs.eg.db",
    "CellChat",
    "monocle3",
    "BiocParallel",
]


def platform_name() -> str:
    system = platform.system().lower()
    if system == "darwin":
        return "macos"
    if system == "windows":
        return "windows"
    if system == "linux":
        return "linux"
    return system or "unknown"


def get_profile(name: str) -> EnvironmentProfile:
    normalized = name.lower().strip()
    if normalized == "python":
        return EnvironmentProfile("python", PYTHON_IMPORTS, [], [])
    if normalized == "minimal":
        return EnvironmentProfile("minimal", PYTHON_IMPORTS, MINIMAL_R_PACKAGES, ["Rscript"])
    if normalized == "extended":
        return EnvironmentProfile(
            "extended",
            PYTHON_IMPORTS,
            MINIMAL_R_PACKAGES + EXTENDED_R_PACKAGES,
            ["Rscript"],
        )
    raise ValueError(f"Unknown profile: {name}. Choose python, minimal, or extended.")


def check_python_import(name: str) -> dict[str, Any]:
    return {"name": name, "available": importlib.util.find_spec(name) is not None}


def check_tools(names: list[str]) -> list[dict[str, Any]]:
    return [{"name": name, "path": shutil.which(name), "available": shutil.which(name) is not None} for name in names]


def check_r_packages(packages: list[str], rscript: str = "Rscript") -> list[dict[str, Any]]:
    if not packages:
        return []
    if shutil.which(rscript) is None and not Path(rscript).exists():
        return [{"name": pkg, "available": False, "version": "", "error": f"{rscript} not found"} for pkg in packages]

    r_expr = """
packages <- strsplit(Sys.getenv("SCMECH_R_PACKAGES"), "\n", fixed = TRUE)[[1]]
packages <- packages[nzchar(packages)]
for (pkg in packages) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  version <- if (ok) as.character(utils::packageVersion(pkg)) else ""
  cat(pkg, ok, version, sep = "\t")
  cat("\n")
}
"""
    env = {"SCMECH_R_PACKAGES": "\n".join(packages)}
    try:
        proc = subprocess.run(
            [rscript, "-e", r_expr],
            check=False,
            capture_output=True,
            text=True,
            env={**os.environ, **env},
        )
    except OSError as exc:
        return [{"name": pkg, "available": False, "version": "", "error": str(exc)} for pkg in packages]

    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        name, available = parts[0], parts[1] == "TRUE"
        seen.add(name)
        rows.append({"name": name, "available": available, "version": parts[2] if len(parts) > 2 else ""})
    for pkg in packages:
        if pkg not in seen:
            rows.append({"name": pkg, "available": False, "version": "", "error": proc.stderr.strip()})
    return rows


def build_report(profile: EnvironmentProfile, json_ready: bool = False, rscript: str = "Rscript") -> dict[str, Any]:
    tools = check_tools(profile.required_tools)
    python_rows = [check_python_import(pkg) for pkg in profile.python_imports]
    r_rows = check_r_packages(profile.r_packages, rscript=rscript)
    all_rows = tools + python_rows + r_rows
    status = "PASSED" if all(item.get("available") for item in all_rows) else "FAILED"
    return {
        "status": status,
        "profile": profile.name,
        "platform": platform_name(),
        "python": python_rows,
        "tools": tools,
        "r": r_rows,
        "json_ready": json_ready,
    }


def print_markdown(report: dict[str, Any]) -> None:
    print("# scMechanism Environment Check")
    print()
    print(f"Status: {report['status']}")
    print(f"Profile: {report['profile']}")
    print(f"Platform: {report['platform']}")
    for section in ("tools", "python", "r"):
        print()
        print(f"## {section.title()}")
        rows = report[section]
        if not rows:
            print("- Not checked")
            continue
        for row in rows:
            marker = "OK" if row["available"] else "MISSING"
            version = f" {row.get('version')}" if row.get("version") else ""
            path = f" ({row.get('path')})" if row.get("path") else ""
            print(f"- {marker}: {row['name']}{version}{path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=["python", "minimal", "extended"], default="minimal")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable or absolute path")
    parser.add_argument("--json", action="store_true", help="Print JSON instead of Markdown")
    parser.add_argument("--out-json", default="", help="Optional path to write the JSON report")
    args = parser.parse_args()

    report = build_report(get_profile(args.profile), json_ready=args.json, rscript=args.rscript)
    if args.out_json:
        Path(args.out_json).write_text(json.dumps(report, indent=2), encoding="utf-8")
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print_markdown(report)
    return 0 if report["status"] == "PASSED" else 1


if __name__ == "__main__":
    sys.exit(main())

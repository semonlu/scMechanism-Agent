#!/usr/bin/env python3
"""Standalone synthetic demo pipeline for local MCP testing."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    out = Path(args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    tables = out / "tables"

    write_csv(tables / "marker_genes.csv", [{"cluster": "0", "gene": "LYZ", "cell_type_hint": "Macrophage"}])
    write_csv(tables / "deg_results.csv", [{"gene": "IL1B", "cell_type": "Macrophage", "log2FC": 1.5}])
    write_csv(tables / "enrichment_results.csv", [{"term": "inflammatory response", "p_adj": 0.01}])
    write_csv(tables / "celltype_proportion.csv", [{"sample_id": "S1", "cell_type": "Macrophage", "proportion": 0.32}])
    (out / "report_skeleton.md").write_text("# Demo Report\n\nSynthetic demo output only.\n", encoding="utf-8")
    (out / "run_log.txt").write_text(f"{datetime.now(timezone.utc).isoformat()} demo completed\n", encoding="utf-8")
    print(f"Demo files written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


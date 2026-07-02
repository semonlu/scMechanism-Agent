#!/usr/bin/env python3
"""Guarded placeholder for a future Scanpy local pipeline."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-path", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--species", default="human")
    args = parser.parse_args()

    result = {
        "status": "error",
        "message": "scanpy_basic execution is disabled in the first local MCP version.",
        "input_exists": Path(args.input_path).exists(),
        "scanpy_installed": importlib.util.find_spec("scanpy") is not None,
        "species": args.species,
    }
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    (Path(args.output_dir) / "scanpy_basic_status.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())


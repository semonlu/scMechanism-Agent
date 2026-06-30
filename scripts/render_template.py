#!/usr/bin/env python3
"""Render a skill template with simple {{PLACEHOLDER}} replacements."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_define(values: list[str]) -> dict[str, str]:
    out = {}
    for item in values:
        if "=" not in item:
            raise SystemExit(f"--define must use KEY=VALUE, got: {item}")
        key, value = item.split("=", 1)
        out[key.strip()] = value
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--define", action="append", default=[], help="Placeholder replacement, e.g. PROJECT_ID=GSE123.")
    args = parser.parse_args()

    text = Path(args.template).read_text(encoding="utf-8")
    for key, value in parse_define(args.define).items():
        text = text.replace("{{" + key + "}}", value)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()

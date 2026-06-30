#!/usr/bin/env python3
"""Organize the flat GSE223751 GEO tar into one 10x MEX folder per sample."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
import tarfile
from pathlib import Path


NAME_RE = re.compile(r"^(GSM\d+)_([^.]+(?:\.\d+)?)\.(barcodes|features|matrix)\.(tsv|mtx)\.gz$")


def sample_id_from_label(label: str) -> str:
    return label.replace(".", "_")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tar", required=True, help="Path to GSE223751_RAW.tar")
    parser.add_argument("--out-dir", required=True, help="Output directory for organized 10x folders")
    parser.add_argument("--metadata", required=True, help="CSV metadata file to write or update")
    args = parser.parse_args()

    tar_path = Path(args.tar)
    out_dir = Path(args.out_dir)
    metadata_path = Path(args.metadata)
    out_dir.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)

    rows: dict[str, dict[str, str]] = {}
    extracted = 0
    with tarfile.open(tar_path, "r") as tar:
        for member in tar.getmembers():
            name = Path(member.name).name
            match = NAME_RE.match(name)
            if not match:
                continue
            gsm, label, kind, _suffix = match.groups()
            sample_id = sample_id_from_label(label)
            sample_dir = out_dir / sample_id
            sample_dir.mkdir(parents=True, exist_ok=True)
            target_name = {
                "barcodes": "barcodes.tsv.gz",
                "features": "features.tsv.gz",
                "matrix": "matrix.mtx.gz",
            }[kind]
            source = tar.extractfile(member)
            if source is None:
                raise RuntimeError(f"Cannot extract {member.name}")
            with (sample_dir / target_name).open("wb") as handle:
                shutil.copyfileobj(source, handle)
            rows.setdefault(sample_id, {
                "sample_id": sample_id,
                "gsm": gsm,
                "stage": label.split(".")[0].upper(),
                "replicate": label.split(".")[1] if "." in label else "1",
                "condition": "developmental_timecourse",
                "tissue": "rotator_cuff_enthesis",
                "organism": "Mus musculus",
                "platform": "10x Genomics Chromium v3",
            })
            extracted += 1

    if not rows:
        raise RuntimeError(f"No GSE223751 10x files were detected in {tar_path}")

    fieldnames = ["sample_id", "gsm", "stage", "replicate", "condition", "tissue", "organism", "platform"]
    with metadata_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for sample_id in sorted(rows):
            writer.writerow(rows[sample_id])

    for sample_dir in sorted(out_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        required = ["matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"]
        missing = [x for x in required if not (sample_dir / x).exists()]
        if missing:
            raise RuntimeError(f"{sample_dir} is missing {missing}")
        with gzip.open(sample_dir / "barcodes.tsv.gz", "rt") as handle:
            n_barcodes = sum(1 for _ in handle)
        print(f"{sample_dir.name}: {n_barcodes} barcodes")

    print(f"Organized {len(rows)} samples and {extracted} files under {out_dir}")
    print(f"Wrote metadata: {metadata_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Diagnose GEO/SRA/single-cell input file lists for the skill workflow."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


ACCESSION_PATTERNS = {
    "GEO Series": re.compile(r"\bGSE\d+\b", re.I),
    "GEO Sample": re.compile(r"\bGSM\d+\b", re.I),
    "SRA Run": re.compile(r"\bSRR\d+\b", re.I),
    "SRA Project": re.compile(r"\bSRP\d+\b", re.I),
}


def normalize_name(value: str) -> str:
    return Path(value.strip().strip('"').strip("'")).name.lower()


def classify_files(items: list[str]) -> dict:
    names = [normalize_name(x) for x in items if x.strip()]
    joined = "\n".join(items)
    accessions = {
        label: sorted(set(match.group(0).upper() for match in pattern.finditer(joined)))
        for label, pattern in ACCESSION_PATTERNS.items()
    }

    has_matrix = any("matrix.mtx" in x for x in names)
    has_barcodes = any("barcodes.tsv" in x for x in names)
    has_features = any("features.tsv" in x or "genes.tsv" in x for x in names)
    has_nonstandard_matrix = any("count_matrix_sparse.mtx" in x for x in names)
    has_nonstandard_barcodes = any("count_matrix_barcodes.tsv" in x for x in names)
    has_nonstandard_genes = any("count_matrix_genes.tsv" in x for x in names)
    has_archive = any(x.endswith((".tar", ".tar.gz", ".tgz")) for x in names)
    has_h5 = any(x.endswith(".h5") or x.endswith(".h5.gz") for x in names)
    has_h5ad = any(x.endswith(".h5ad") for x in names)
    has_r_object = any(x.endswith((".rds", ".rda", ".h5seurat")) for x in names)
    has_loom = any(x.endswith(".loom") for x in names)
    has_table = any(x.endswith((".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz")) for x in names)
    has_fastq = any(x.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")) for x in names)
    has_tcr_bcr = any("contig_annotations" in x or "clonotypes" in x for x in names)
    has_spatial = any("spatial" in x or "tissue_positions" in x or "scalefactors" in x for x in names)

    formats = []
    if has_matrix and has_barcodes and has_features:
        formats.append("10x MEX")
    if has_nonstandard_matrix and has_nonstandard_barcodes and has_nonstandard_genes:
        formats.append("10x non-standard MEX")
    elif has_archive:
        formats.append("compressed GEO archive requiring extraction")
    if has_h5:
        formats.append("10x HDF5 or generic HDF5")
    if has_h5ad:
        formats.append("AnnData h5ad")
    if has_r_object:
        formats.append("R/Seurat object")
    if has_loom:
        formats.append("loom")
    if has_fastq or accessions["SRA Run"] or accessions["SRA Project"]:
        formats.append("FASTQ/SRA raw reads")
    has_standard_mex = has_matrix and has_barcodes and has_features
    has_nonstandard_mex = has_nonstandard_matrix and has_nonstandard_barcodes and has_nonstandard_genes
    if has_table and not (has_standard_mex or has_nonstandard_mex):
        formats.append("flat expression/result table")
    if has_tcr_bcr:
        formats.append("TCR/BCR contig/clonotype files")
    if has_spatial:
        formats.append("spatial transcriptomics files")

    direct = "uncertain"
    needs_rebuild = False
    if any(x in formats for x in ["10x MEX", "10x non-standard MEX", "10x HDF5 or generic HDF5", "AnnData h5ad", "R/Seurat object", "loom"]):
        direct = "yes"
    if "FASTQ/SRA raw reads" in formats and len(formats) == 1:
        direct = "no"
        needs_rebuild = True
    elif "FASTQ/SRA raw reads" in formats:
        needs_rebuild = True

    recommendations = []
    if "10x MEX" in formats:
        recommendations.append("Use Seurat::Read10X or scanpy.read_10x_mtx.")
    if "10x non-standard MEX" in formats:
        recommendations.append("Use input_type=10x_nonstandard with scripts/course_adapted/01_seurat_v5_core_pipeline.R or 00_multi_sample_merge_harmony.R; the reader uses Matrix::readMM plus barcode/gene TSV files.")
    if "compressed GEO archive requiring extraction" in formats:
        recommendations.append("Extract archives into a project input folder, then inspect for standard 10x MEX, non-standard count_matrix_* MEX, H5, RDS, h5ad, or plain matrices before choosing the reader.")
    if "10x HDF5 or generic HDF5" in formats:
        recommendations.append("Use Seurat::Read10X_h5 for 10x H5 or scanpy.read_10x_h5 after confirming structure.")
    if "AnnData h5ad" in formats:
        recommendations.append("Use scanpy.read_h5ad; verify raw counts/layers before differential expression.")
    if "R/Seurat object" in formats:
        recommendations.append("Use readRDS/load and inspect class, assays, layers, reductions, and metadata.")
    if "FASTQ/SRA raw reads" in formats:
        recommendations.append("Do not run downstream Seurat/Scanpy directly; first rebuild matrices with Cell Ranger, STARsolo, kallisto-bustools, or alevin-fry.")
    if "flat expression/result table" in formats:
        recommendations.append("Inspect orientation, gene identifiers, delimiter, and whether values are raw counts, normalized expression, or results.")

    risks = []
    if not any("metadata" in x or "phenotype" in x or "sample" in x for x in names):
        risks.append("Sample metadata/group design not detected.")
    if needs_rebuild:
        risks.append("Raw-read reconstruction requires reference genome, chemistry, and compute resources.")
    if has_archive:
        risks.append("Archive contents must be inspected after extraction; file extension alone is not enough to choose Seurat/Scanpy reader.")
    if has_h5ad or has_r_object:
        risks.append("Processed objects may not retain raw counts needed for DE, CellChat, CNV, or deconvolution.")
    if not formats:
        risks.append("No supported single-cell file pattern was detected.")

    return {
        "detected_formats": formats or ["unknown"],
        "accessions": accessions,
        "direct_downstream_analysis": direct,
        "needs_fastq_rebuild": needs_rebuild,
        "recommended_reading": recommendations,
        "major_risks": risks,
        "next_step": "Generate an analysis plan only after organism, metadata columns, and research comparison are known.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--file-list", help="Text file containing GEO supplementary file names or local paths.")
    parser.add_argument("--items", nargs="*", default=[], help="Inline file names/accessions.")
    parser.add_argument("--out-json", default="geo_input_diagnosis.json")
    parser.add_argument("--out-md", default="geo_input_diagnosis.md")
    args = parser.parse_args()

    items = list(args.items)
    if args.file_list:
        items.extend(Path(args.file_list).read_text(encoding="utf-8-sig").splitlines())

    result = classify_files(items)
    Path(args.out_json).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    md = ["# GEO/SRA Input Diagnosis", ""]
    for key, value in result.items():
        md.append(f"## {key.replace('_', ' ').title()}")
        if isinstance(value, list):
            md.extend([f"- {x}" for x in value] or ["- none"])
        else:
            md.append(f"`{value}`")
        md.append("")
    Path(args.out_md).write_text("\n".join(md), encoding="utf-8")


if __name__ == "__main__":
    main()

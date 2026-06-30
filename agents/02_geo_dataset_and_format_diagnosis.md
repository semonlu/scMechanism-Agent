# 02 GEO Dataset And Format Diagnosis

## Goal

Decide whether GEO/SRA supplementary files, local file lists, or publication-derived dataset notes are ready for downstream single-cell analysis.

## Required Recognition

```text
Publication text: GEO/SRA accessions, organism, tissue, disease/control groups, platform, sample count, processed matrix/object availability
10x MEX: matrix.mtx / barcodes.tsv / features.tsv or genes.tsv
10x HDF5: .h5
AnnData: .h5ad
Seurat: .rds / .rda / .h5Seurat
loom: .loom
Plain matrix: .csv / .tsv / .txt / .txt.gz
Raw data: SRR / FASTQ / fq.gz
TCR/BCR: contig_annotations.csv / clonotypes.csv
Spatial transcriptomics: spatial/ folder, tissue_positions, scalefactors, image files
```

## Required Output

```text
identified data formats
whether downstream analysis can start directly: yes / no / uncertain
recommended reader or reconstruction method
whether FASTQ-to-matrix reconstruction is required
metadata columns needed for grouping, batch, donor, tissue, and sample identity
major risks
next-step recommendation
```

## Script Call

```bash
python scripts/diagnose_geo_inputs.py --file-list supplementary_files.txt --out-json diagnosis.json --out-md diagnosis.md
```

## Rules

- If the user provides a publication URL, abstract, methods text, or PDF-derived notes, first extract accessions, sample groups, organism, tissue, sequencing platform, and whether processed matrices or objects are available.
- Classify matrix.mtx + barcodes.tsv + features.tsv/genes.tsv as 10x MEX.
- Classify filtered_feature_bc_matrix.h5 or raw_feature_bc_matrix.h5 as 10x HDF5.
- Prefer Scanpy for `.h5ad`, but check raw/layers before DE or pathway claims.
- Read `.rds/.rda/.h5Seurat` in R and confirm whether the object is a Seurat object with usable assays and metadata.
- If only SRR/FASTQ exists, do not say it can directly enter Seurat downstream analysis; first reconstruct an expression matrix.
- If metadata are missing, state that group comparison, batch correction, cell proportion comparison, and clinical interpretation are limited.

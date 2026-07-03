# 02 GEO Dataset And Format Diagnosis

## Goal

Decide whether GEO/SRA supplementary files, local file lists, or publication-derived dataset notes are ready for downstream single-cell analysis.

## Required Recognition

```text
Publication text: GEO/SRA accessions, organism, tissue, disease/control groups, platform, sample count, processed matrix/object availability
10x MEX: matrix.mtx / barcodes.tsv / features.tsv or genes.tsv
10x non-standard MEX: count_matrix_sparse.mtx / count_matrix_barcodes.tsv / count_matrix_genes.tsv, often after extracting GEO tar/tar.gz archives
10x HDF5: .h5
AnnData: .h5ad
Seurat: .rds / .rda / .h5Seurat
loom: .loom
Plain matrix: .csv / .tsv / .txt / .txt.gz
Raw data: SRR / FASTQ / fq.gz
Compressed GEO archive: .tar / .tar.gz / .tgz; extract first, then classify the contents
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
- Classify count_matrix_sparse.mtx + count_matrix_barcodes.tsv + count_matrix_genes.tsv as 10x non-standard MEX; route Seurat runs with `INPUT_TYPE=10x_nonstandard` to `scripts/course_adapted/01_seurat_v5_core_pipeline.R`, or use `00_multi_sample_merge_harmony.R` with `input_type=10x_nonstandard` in the sample table.
- If GEO provides `.tar`, `.tar.gz`, or `.tgz` archives, do not assume FASTQ and do not choose a reader only from the archive name. Ask to extract into a project input folder and inspect whether each sample directory contains standard 10x, non-standard 10x, H5, RDS/h5ad, plain matrix, or raw FASTQ files.
- Classify filtered_feature_bc_matrix.h5 or raw_feature_bc_matrix.h5 as 10x HDF5.
- Prefer Scanpy for `.h5ad`, but check raw/layers before DE or pathway claims.
- Read `.rds/.rda/.h5Seurat` in R and confirm whether the object is a Seurat object with usable assays and metadata.
- If only SRR/FASTQ exists, do not say it can directly enter Seurat downstream analysis; first reconstruct an expression matrix.
- If metadata are missing, state that group comparison, batch correction, cell proportion comparison, and clinical interpretation are limited.

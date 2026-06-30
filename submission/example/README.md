# GSE223751 Shoulder-Related Single-Cell Example

This example demonstrates a reproducible scMechanism Agent workflow on a public GEO single-cell dataset related to shoulder biology.

## Dataset

- GEO accession: GSE223751
- Title: Single-cell RNA sequencing reveals cellular and molecular heterogeneity in fibrocartilaginous enthesis formation
- Organism: Mus musculus
- Tissue/context: rotator cuff fibrocartilaginous enthesis, a shoulder tendon-bone junction model
- Platform: 10x Genomics Chromium v3
- Supplementary data: GSE223751_RAW.tar with 7 10x MEX samples

Important scope note: precise GEO searches for `adhesive capsulitis`, `frozen shoulder`, and `shoulder periarthritis` found public shoulder/frozen-shoulder transcriptomic datasets, but not a directly reusable single-cell shoulder periarthritis matrix. GSE223751 is used here as the closest shoulder-region single-cell example and should not be described as a frozen-shoulder case-control cohort.

## Files

```text
scripts/01_download_gse223751.ps1
scripts/02_prepare_gse223751_mex.py
scripts/03_run_gse223751_seurat.R
metadata/gse223751_sample_metadata.csv
geo_supplementary_files.txt
geo_input_diagnosis.md
analysis_plan.md
```

Generated raw data and organized 10x folders are written under `data/` and are intentionally gitignored. Result figures and tables are written under `results/gse223751_seurat/`.

## Run

From this directory:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\01_download_gse223751.ps1
python .\scripts\02_prepare_gse223751_mex.py --tar .\data\raw\GSE223751_RAW.tar --out-dir .\data\processed\10x_by_sample --metadata .\metadata\gse223751_sample_metadata.csv
Rscript .\scripts\03_run_gse223751_seurat.R --input-dir .\data\processed\10x_by_sample --metadata .\metadata\gse223751_sample_metadata.csv --output-dir .\results\gse223751_seurat
```

The first R script performs multi-sample import, merge, QC, normalization, PCA, UMAP, clustering, SingleR automated annotation, marker detection, marker-support scoring, and visualization. The course-derived downstream modules are then rendered and run with:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\04_run_course_modules.ps1 -Rscript Rscript
```

The course-module runner renders:

- marker enrichment from `scripts/course_adapted/02_marker_enrichment_from_seurat.R`
- CellChat from `scripts/course_adapted/03_cellchat_from_seurat.R`
- Monocle3 pseudotime from `scripts/course_adapted/04_monocle3_from_seurat.R`

For this example, Monocle3 uses a marker-supported multi-state trajectory input (`tenocyte_fibroblast`, `chondrocyte`, `osteoblast`, and `cycling`), downsampled to 5,000 cells balanced by `marker_support_label`. It imports Seurat UMAP coordinates, learns the graph, and orders cells with E15 cells as the root. The verified pseudotime table contains 5,000 cells, all with finite pseudotime.

The integrated environment checker was run from this project after the environment scripts were merged into the skill. A portable summary is stored in `logs/environment_check_20260630_203000` and confirms R 4.4.2, 83/83 required R imports, and all configured Python imports in `seuratv5-course-py`.

See `FULL_WORKFLOW_STATUS.md` and `course_source_traceability.tsv` for the exact course-source mapping.

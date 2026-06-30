# GEO Input Diagnosis

## Dataset

- Accession: GSE223751
- Organism: Mus musculus
- Context: rotator cuff fibrocartilaginous enthesis, shoulder tendon-bone junction
- Platform: 10x Genomics Chromium v3
- Supplementary file: `GSE223751_RAW.tar`

## Detected Formats

- 10x MEX
- Flat GEO tar containing per-sample `matrix.mtx.gz`, `barcodes.tsv.gz`, and `features.tsv.gz`

## Direct Downstream Analysis

Yes, after reorganizing the flat tar into one 10x directory per sample.

## Recommended Reader

- Use `Seurat::Read10X(data.dir = sample_dir)` for each organized sample folder.
- Merge sample-level Seurat objects with sample-prefixed cell IDs.
- Use SingleR for automated cluster-level annotation, with marker tables and dot plots as supporting evidence.

## Major Risks

- The dataset is shoulder-region biology but not adhesive capsulitis/frozen shoulder.
- Metadata represent developmental time points, not disease/control groups.
- Group comparisons should be descriptive unless additional validated metadata are added.

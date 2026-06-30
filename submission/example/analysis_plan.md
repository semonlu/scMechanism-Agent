# Analysis Plan

## A. Data Readiness

- Use GEO accession GSE223751 as a public shoulder-region single-cell example.
- Process `GSE223751_RAW.tar`, which contains 7 sample-level 10x MEX file triplets.
- Organize the flat GEO tar into one 10x directory per sample before calling Seurat.
- Keep sample metadata columns: `sample_id`, `gsm`, `stage`, `replicate`, `condition`, `tissue`, `organism`, and `platform`.

## B. QC

- Calculate `nFeature_RNA`, `nCount_RNA`, mitochondrial percentage, and ribosomal percentage.
- Save pre-filter and post-filter violin plots plus sample-level QC summaries.
- Default filters are exploratory: `nFeature_RNA > 200`, `nFeature_RNA < 8000`, and `percent.mt < 20`.

## C. Normalization And Clustering

- Merge all samples with sample-prefixed cell barcodes.
- Run `NormalizeData`, `FindVariableFeatures`, `ScaleData`, PCA, neighbor graph, clustering, and UMAP.
- Visualize UMAP by cluster, sample, and developmental stage.

## D. SingleR Annotation And Marker Evidence

- Run cluster-level SingleR annotation with a mouse reference from `celldex`.
- Save the SingleR label table, SingleR score heatmap, and UMAP grouped by SingleR labels.
- Run cluster marker detection after Seurat v5 layers are joined.
- Score musculoskeletal and immune marker sets only as auxiliary evidence.
- Export marker tables, `annotation_evidence.tsv`, cell metadata, and SingleR-based cell-type proportions.

## E. Interpretation Limits

- This is not a frozen-shoulder case-control dataset.
- Stage comparisons are descriptive unless sample-level biological replication supports a formal model.
- Cell-type labels are marker-based hypotheses that require domain review.

# GSE223751 Example Run Summary

## Environment

- R: 4.4.2
- Seurat: 5.5.1
- SingleR: 2.8.0
- celldex: 1.16.0
- JAGS: 4.3.1
- SeuratV5 environment check: `seuratv5_environment_check/20260630_154614`
- Integrated scMechanism environment check: `logs/environment_check_20260630_203000`
- R package check: 83/83 required imports installed
- Python environment: `seuratv5-course-py`
- Python package check: all required modules imported successfully

## Dataset

- GEO accession: GSE223751
- Data: 7 sample-level 10x MEX matrices from `GSE223751_RAW.tar`
- Scope: shoulder-region rotator cuff enthesis single-cell example, not a frozen-shoulder case-control cohort

## Completed Workflow

1. Organized the flat GEO tar into per-sample 10x MEX folders.
2. Imported 7 samples with Seurat.
3. Merged samples and calculated QC metrics.
4. Filtered cells with exploratory thresholds.
5. Ran normalization, variable feature selection, scaling, PCA, clustering, and UMAP.
6. Ran SingleR cluster-level automated annotation with `celldex::MouseRNAseqData()`.
7. Ran marker detection and marker-support scoring.
8. Rendered and ran course-derived marker enrichment, CellChat, and Monocle3 modules from `scripts/course_adapted/`.
9. Re-ran Monocle3 pseudotime with the course-style trajectory logic: multi-state marker-supported trajectory input, Seurat UMAP import, graph learning, and E15 root selection.
10. Exported figures, tables, logs, and processed objects.

## Key Run Metrics

- Cells retained after QC: 46,955
- Genes in processed object: 21,187
- Seurat clusters: 25
- Course-derived downstream modules completed: marker/GO/KEGG, CellChat, Monocle3 pseudotime
- Monocle3 trajectory cells: 5,000
- Monocle3 cell states: chondrocyte 820; cycling 762; osteoblast 292; tenocyte_fibroblast 3,126
- Finite pseudotime cells: 5,000
- Unique finite pseudotime values: 4,537
- Monocle3 root rule: E15 cells within the sampled trajectory object

## Main Outputs

- `results/gse223751_seurat/figures/umap_clusters.pdf`
- `results/gse223751_seurat/figures/umap_singleR_labels.pdf`
- `results/gse223751_seurat/figures/singleR_score_heatmap.pdf`
- `results/gse223751_seurat/figures/marker_dotplot.pdf`
- `results/gse223751_seurat/tables/singleR_cluster_labels.csv`
- `results/gse223751_seurat/tables/annotation_evidence.tsv`
- `results/gse223751_seurat/tables/cluster_markers.csv`
- `results/gse223751_seurat/tables/celltype_proportions.tsv`
- `results/gse223751_seurat/objects/gse223751_processed_seurat.rds`
- `results/course_modules/marker_enrichment/tables/GO_BP_top_markers.csv`
- `results/course_modules/marker_enrichment/tables/KEGG_top_markers.csv`
- `results/course_modules/cellchat/tables/cellchat_ligand_receptor.csv`
- `results/course_modules/monocle3/tables/pseudotime.csv`
- `results/course_modules/monocle3/figures/trajectory_celltypes.pdf`
- `results/course_modules/monocle3/figures/trajectory_celltypes_with_graph.pdf`
- `results/course_modules/monocle3/figures/seurat_imported_umap.pdf`
- `results/course_modules/monocle3/figures/pseudotime_umap.pdf`
- `results/course_modules/monocle3/figures/top_trajectory_genes_pseudotime.pdf`
- `course_source_traceability.tsv`
- `FULL_WORKFLOW_STATUS.md`

## Caveats

- GSE223751 is a shoulder enthesis development dataset; it should not be interpreted as shoulder periarthritis or frozen-shoulder disease-control evidence.
- SingleR labels are automated reference-based labels and should be reviewed against marker genes and tissue context.
- Stage-level trends are descriptive unless the analysis design is extended with appropriate replication-aware statistical models.

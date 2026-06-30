# Full Workflow Status

This example now separates the full workflow into completed modules and explicitly records which course-derived scripts were used or referenced.

## Completed

| Stage | Status | Main outputs |
|---|---|---|
| GEO diagnosis and 10x preparation | completed | `data/processed/10x_by_sample`, `geo_input_diagnosis.md` |
| Seurat import, merge, QC, PCA, UMAP, clustering | completed | `results/gse223751_seurat/figures`, `results/gse223751_seurat/objects/gse223751_processed_seurat.rds` |
| SingleR automated annotation | completed | `singleR_cluster_labels.csv`, `umap_singleR_labels.pdf`, `singleR_score_heatmap.pdf` |
| Marker detection and GO/KEGG enrichment | completed | `results/course_modules/marker_enrichment/tables` |
| CellChat ligand-receptor inference | completed | `results/course_modules/cellchat/tables/cellchat_ligand_receptor.csv`, CellChat network PDFs |
| Monocle3 pseudotime | completed | `results/course_modules/monocle3/tables/pseudotime.csv`, `trajectory_genes.csv`, native/Seurat-imported UMAP trajectory PDFs |

## Course-Derived Calls

The downstream modules were rendered from:

- `scripts/course_adapted/02_marker_enrichment_from_seurat.R`
- `scripts/course_adapted/03_cellchat_from_seurat.R`
- `scripts/course_adapted/04_monocle3_from_seurat.R`

Their corresponding source references are listed in `course_source_traceability.tsv`. The rendered runnable scripts are stored in `rendered_course_modules/`.

## Monocle3 Pseudotime Check

The current pseudotime example follows the course-style Monocle3 workflow with a biologically related multi-state mesenchymal trajectory input, balanced downsampling, graph learning, explicit root selection, and Seurat UMAP import for visual consistency with the upstream Seurat object.

- Subset: `marker_support_label %in% c('tenocyte_fibroblast', 'chondrocyte', 'osteoblast', 'cycling')`
- Downsampling: 5,000 cells, balanced by `marker_support_label`
- Root selection: `grepl('^E15$', stage)`
- Cells in trajectory: 5,000
- Finite pseudotime cells: 5,000
- Unique finite pseudotime values: 4,537
- Cell states: chondrocyte 820; cycling 762; osteoblast 292; tenocyte_fibroblast 3,126
- Stage distribution: E15 1,813; P1 1,388; P7 915; P14 301; P28 583

Main figures:

- `results/course_modules/monocle3/figures/monocle3_native_umap.pdf`
- `results/course_modules/monocle3/figures/seurat_imported_umap.pdf`
- `results/course_modules/monocle3/figures/monocle3_vs_seurat_umap.pdf`
- `results/course_modules/monocle3/figures/trajectory_celltypes.pdf`
- `results/course_modules/monocle3/figures/trajectory_celltypes_with_graph.pdf`
- `results/course_modules/monocle3/figures/pseudotime_umap.pdf`
- `results/course_modules/monocle3/figures/trajectory_stage.pdf`
- `results/course_modules/monocle3/figures/top_trajectory_genes_pseudotime.pdf`

## Not Applicable Or Not Run

- DoubletFinder/scDblFinder: not run for this public processed GEO example because raw loading and chemistry metadata were not provided.
- SCINA/scPred/TransferData: optional annotation alternatives; SingleR was used as the automatic annotation method.
- copykat/inferCNV: not biologically appropriate for this non-tumor mouse developmental enthesis dataset.
- CIBERSORT/MuSiC: requires matched bulk expression input and is not applicable to the scRNA-only example.

These modules remain represented in the skill's course source scripts and can be activated for suitable datasets.

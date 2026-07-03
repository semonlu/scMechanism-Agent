# Output File Checklist

## Core Outputs

```text
figures/qc_pre_filter.pdf
figures/qc_post_filter.pdf
figures/umap_clusters.pdf
figures/umap_celltypes.pdf
tables/metadata_checked.tsv
tables/qc_summary.tsv
tables/cluster_markers.tsv
tables/annotation_evidence.tsv
objects/processed_seurat.rds or objects/processed.h5ad
logs/session_info.txt
logs/commands.txt
data_input_manifest.json
data_analysis_qc.md
```

## Optional Module Outputs

```text
tables/deg_results.tsv
tables/enrichment_results.tsv
tables/celltype_proportions.tsv
tables/cellchat_ligand_receptor.tsv
tables/cellchat_pathway_summary.tsv
tables/pseudotime.tsv
tables/trajectory_genes.tsv
figures/cellchat_network.pdf
figures/pseudotime_umap.pdf
figures/enrichment_dotplot.pdf
```

## Review Questions

- Are sample-level metadata and group labels present?
- Does `data_analysis_qc.md` confirm that the analyzed `input_path` matches the planned/downloaded/extracted/registered dataset?
- Are raw counts retained?
- Are QC thresholds justified?
- Are batch effects assessed?
- Are annotation labels backed by markers?
- Are group comparisons sample-aware?
- Are optional modules marked as inference/hypothesis?
- Are package versions and commands recorded?

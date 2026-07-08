# Single-cell QC Rules

QC thresholds must come from the current data distribution, tissue background, organism, platform, and expected cell types. Course-script values are examples only and must not be copied as universal defaults.

## Core QC Metrics

| Metric | Meaning | Review risk |
|---|---|---|
| `nFeature_RNA` | Detected genes per cell | Very low suggests empty droplets or low-quality cells; very high can indicate doublets or high-RNA-content cells. |
| `nCount_RNA` | UMI/read count per cell | Very low suggests poor quality; very high can indicate doublets or large/high-RNA cells. |
| `percent.mt` | Mitochondrial fraction | High values can indicate dying, damaged, stressed, or metabolically active cells. Do not apply a fixed 5% or 10% cutoff blindly. |
| `percent.ribo` | Ribosomal fraction | Can reflect technical noise, low quality, or biological state; review by tissue and cluster. |
| `percent.hb` | Hemoglobin fraction | Suggests blood contamination, red-cell carryover, or ambient RNA. |
| doublet score | Doublet risk | Evaluate by sample or loading batch; high-scoring clusters need marker review. |
| batch/sample | Technical or biological source | Do not integrate blindly before checking whether group and batch are confounded. |

## Required QC Outputs

- Pre-filter and post-filter violin plots for `nFeature_RNA`, `nCount_RNA`, `percent.mt`, and when available `percent.ribo`/`percent.hb`.
- Scatter plots: `nCount_RNA` vs `nFeature_RNA`, and `nCount_RNA` vs `percent.mt`.
- Per-sample pre/post cell retention table.
- A threshold table explaining each filter and why it was chosen.
- Doublet method/status table if doublet detection is possible.
- Ambient RNA, hemoglobin, ribosomal, mitochondrial, and cell-cycle review notes when suspicious clusters appear.

## Threshold Selection

Use common values only as starting points:

- `min.features`: often 200 to 300, but lower may be needed for low-RNA cell types or small validation data.
- `max.features`: often 2,500 to 6,000 or higher depending on platform, tissue, and cell size.
- `percent.mt`: often 5% to 25%, but tissue and preservation method matter.
- `min.cells` per gene: often 3.
- `min.umi`: often 500 or 1,000 when sequencing depth supports it.

Choose thresholds after looking at each sample separately. If uncertain, use permissive filtering first, then re-check whether low-quality clusters remain after clustering.

## Ambient RNA And Contamination

Raise an ambient RNA warning when:

- Red-cell, epithelial, immune, or other high-expression markers appear weakly across many unrelated clusters.
- Hemoglobin genes such as `HBA1`, `HBA2`, `HBB`, or `HBD` are widespread outside erythroid-like cells.
- A high-expression marker appears everywhere rather than in a coherent cell population.

Possible actions include SoupX/DecontX, removing obvious contamination-driven genes from interpretation, or marking the analysis as limited when raw droplets are unavailable.

## Doublet Review

Doublets should be evaluated per sample or loading batch whenever raw counts and metadata support it.

Review:

- Doublet score distribution.
- Predicted doublets on UMAP.
- Clusters co-expressing mutually exclusive lineage markers.
- Per-sample doublet rates.
- Whether the workflow was rerun after doublet removal.

After removing predicted doublets, rerun normalization, HVG selection, scaling, PCA, neighbors, clustering, and UMAP.

## Course-Derived Notes

Course scripts use `VlnPlot()`, `FeatureScatter()`, and `subset()`. Example thresholds such as `nCount_RNA <= 25000`, `nFeature_RNA <= 5000`, `percent.mt <= 25`, or `percent.rb <= 40` are not universal defaults. They must be justified against the current dataset.

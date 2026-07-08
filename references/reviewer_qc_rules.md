# Reviewer-Style scRNA-seq QC Rules

The skill should behave like a single-cell reviewer, not like a script runner. At every step ask:

- Could this be technical noise rather than biology?
- Could a wrong annotation be propagated into DE, CellChat, pseudotime, or mechanism text?
- Is the figure attractive but unsupported by QC, markers, metadata, or statistics?

## Mandatory 10-Step Workflow Audit

Before interpreting any single-cell result, audit whether the workflow has evidence for these 10 steps:

| Step | Required audit question | Minimum evidence |
|---|---|---|
| 1. Data acquisition | Was the dataset obtained from GEO/SRA or another traceable local/public source? | accession, source URL/file manifest, `input_diagnosis.md`, or `data_input_manifest.json` |
| 2. Data import | Was the matrix/object read correctly with metadata aligned to barcodes? | import summary, object construction log, metadata check, cell/gene counts |
| 3. Quality control | Were mitochondrial, hemoglobin/red-cell, counts, genes, low-quality cells, and sample retention reviewed? | QC violin/scatter plots, threshold table, pre/post cell counts |
| 4. Normalization and scaling | Were NormalizeData/ScaleData/FindVariableFeatures or SCT/Scanpy equivalents run, and were raw counts retained? | normalization parameters, HVG table, raw-count layer/assay note |
| 5. Confounder removal/review | Were cell cycle, doublets, batch, ambient RNA, and over-integration risks reviewed? | doublet table, cell-cycle score/status, batch/integration plots, ambient RNA note |
| 6. Dimensionality reduction and clustering | Were PCA, UMAP/tSNE, neighbors, and resolution choices audited? | elbow/PCA outputs, UMAP/tSNE, `resolution_sweep.tsv`, `cluster_marker_audit.tsv` |
| 7. Cell annotation | Were automatic and manual annotations cross-checked? | `annotation_evidence.tsv`, marker dot/feature plots, confidence/unknown labels |
| 8. Marker detection | Were marker genes found and reviewed for specificity? | `cluster_markers.*`, top marker tables, marker audit |
| 9. Pseudotime | If used, was it restricted to a plausible lineage with a justified root? | pseudotime proposal/approval, trajectory table/plot, root rationale |
| 10. Cell communication | If used, did it pass annotation/cell-count gates and remain a computational inference? | downstream proposal/approval, cell counts per group, ligand-receptor table |

The result reviewer should export or infer `tables/workflow_step_audit.tsv` and report every missing or gated step. Steps 9 and 10 are not mandatory for every project, but they must not run before steps 1-8 are reviewable.

## Resolution Decision Rules

Resolution is a scientific decision, not a fixed parameter.

For major cell-class annotation, first sweep at least `0.1`, `0.3`, `0.5`, and `0.8`, then review:

- Cluster count: for typical 10x datasets with about 5,000 to 20,000 cells, roughly 10 to 20 clusters is often a reasonable major-lineage range. Fewer than about 8 to 10 can miss cell types; more than about 25 can split one biological class into artificial states.
- Marker clarity: top markers should show stable canonical major-lineage genes, for example `CD3D` for T cells or `MS4A1` for B cells. If top markers are dominated by ribosomal, heat-shock, mitochondrial, hemoglobin, or generic stress genes, the clustering is not ready for fine annotation.
- Biological plausibility: clusters should match tissue context, known cell classes, sample composition, and marker plots. Clusters driven mainly by sample, batch, QC metrics, or doublet score require review.

The core Seurat workflow should export `tables/resolution_sweep.tsv` and `tables/cluster_marker_audit.tsv`. Do not treat the selected resolution as final until these tables and UMAP/marker plots have been reviewed.

## Subcluster Rules

Subcluster analysis is not simply increasing global resolution.

To study heterogeneity inside one cell class:

1. Annotate major cell classes first with a conservative major-lineage resolution.
2. Subset the target class, for example T cells, fibroblasts, endothelial cells, or macrophages.
3. On the subset, rerun variable feature selection, scaling, PCA, neighbor graph construction, clustering, and UMAP.
4. Re-check markers, doublets, stress, cell cycle, sample dominance, and batch effects inside the subset.
5. Report the subcluster object and parameters separately from the global object.

Rerunning HVG selection on the subset is required because global HVGs often miss within-lineage programs.

## Annotation Evidence Logic

Use this chain:

1. Biological prior knowledge.
2. Marker expression patterns.
3. Cell-type hypothesis.
4. Cross-check with automated references and tissue context.
5. Final label with confidence and review note.

Do not force every cluster into a named fine subtype. Use `Unknown`, `Ambiguous`, or a coarse lineage when evidence is weak.

## Unclear Cluster Review

For any unclear cluster, audit these possibilities before naming it:

- Residual doublet: mixed mutually exclusive lineage markers.
- Stress or dissociation damage: heat-shock, ribosomal, mitochondrial, immediate early, or low-complexity signatures.
- Cell cycle: proliferating clusters often express `MKI67`, `TOP2A`, `UBE2C`, `CCNB1`, `CCNB2`, and `CENPF`. If S/G2M scores are high relative to other clusters, annotate the state rather than inventing a new type.
- Transitional state: use trajectory only when a plausible continuum and root exist.
- True novel population: requires internally consistent markers, enough cells, at least about 1% of total cells, recurrence across multiple samples, and biologically meaningful pathway enrichment. Generic stress, ribosomal, mitochondrial, or hemoglobin signals do not support novelty.

For unclear clusters that remain biologically plausible after QC review, check cell-type marker databases and primary literature before naming them. At minimum, consider CellMarker, PanglaoDB, CellTaxonomy, and tissue-specific papers; record the evidence source in `annotation_evidence.tsv`.

When a targeted contrast against a neighboring or biologically related cluster is available:

- More than about 300 DE genes can indicate a meaningful functional state difference worth deeper review.
- Fewer than about 100 DE genes usually supports `Unknown`, `Ambiguous`, or merging with a nearby cluster unless strong marker evidence exists.

These numbers are review heuristics, not statistical laws; always check sample balance, cell count, effect size, and adjusted p values.

## Downstream Interpretation Boundaries

- CellChat depends on correct annotation and enough cells per group. It generates ligand-receptor hypotheses, not proof of physical interaction.
- Pseudotime requires a biologically plausible subset, continuum, and root. It suggests possible state ordering, not experimental time.
- DE should be sample-aware where possible. Treating thousands of cells as independent biological replicates inflates significance.
- Enrichment is functional context, not mechanism proof.

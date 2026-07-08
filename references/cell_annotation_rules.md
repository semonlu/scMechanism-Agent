# Cell Annotation Rules

Cell annotation must combine canonical markers, tissue context, organism, cluster-level expression patterns, and automated reference labels. Automated labels are supporting evidence, not final authority.

## Minimum Evidence

Every completed annotation should export an evidence table equivalent to `annotation_evidence.tsv` with:

```text
cluster
final_label
coarse_label
singleR_label
singleR_pruned_label
singleR_delta_next
top_markers
canonical_marker_support
conflicting_markers
cell_count
confidence
review_note
```

When available, add these reviewer fields:

```text
marker_count_vs_all
relative_de_gene_count
cell_cycle_warning
doublet_warning
stress_warning
sample_recurrence
novelty_support
```

Confidence should use:

- `high`: marker evidence and reference label agree.
- `medium`: broad lineage is clear but subtype is uncertain.
- `low`: marker evidence is weak, conflicting, or reference label is missing.

Clusters with low confidence must keep `Unknown`, `Ambiguous`, or a coarse lineage label. Do not force fine subtype labels.

## Common Markers

| Cell type | Markers |
|---|---|
| T cells | CD3D, CD3E, TRAC |
| CD4 T cells | CD4, IL7R, CCR7 |
| CD8 T cells | CD8A, CD8B, GZMB, NKG7 |
| NK cells | NKG7, GNLY, KLRD1 |
| B cells | MS4A1, CD79A, CD79B |
| Plasma cells | MZB1, JCHAIN, XBP1 |
| Monocytes/macrophages | LST1, S100A8, S100A9, CD68, C1QA |
| Dendritic cells | FCER1A, CLEC10A, LILRA4 |
| Endothelial cells | PECAM1, VWF, KDR |
| Fibroblasts | COL1A1, COL1A2, DCN, LUM |
| Epithelial cells | EPCAM, KRT8, KRT18, KRT19 |
| Smooth muscle/pericytes | ACTA2, RGS5, PDGFRB |
| Mast cells | TPSAB1, CPA3 |

## Review Rules

- If all clusters receive the same label, treat annotation as failed until marker plots prove that the dataset truly contains one lineage.
- If SingleR/CellTypist label conflicts with top markers, prefer marker evidence and mark the cluster for review.
- If a cluster has mixed lineage markers, mark `Ambiguous` and consider doublet, low-quality cells, or over-clustering.
- If the study tissue lacks a claimed cell type biologically, mark the label as suspicious and require manual review.
- If top markers are mostly ribosomal, mitochondrial, hemoglobin, heat-shock, stress, or cell-cycle genes, do not assign a new cell type label without additional evidence.
- If a cluster is unclear, compare it to neighboring or biologically related clusters when possible. More than about 300 DE genes can justify deeper functional review; fewer than about 100 DE genes usually supports `Unknown`, `Ambiguous`, or merge unless canonical markers are strong.
- For unclear but plausible clusters, review marker databases and literature before final naming. Suggested sources include CellMarker, PanglaoDB, CellTaxonomy, and tissue-specific papers; record the evidence source and confidence in `annotation_evidence.tsv`.
- Proliferating clusters should be audited with cell-cycle markers such as `MKI67`, `TOP2A`, `UBE2C`, `CCNB1`, `CCNB2`, and `CENPF`.
- A claim of a novel cell type requires internally consistent markers, enough cells, at least about 1% of total cells, recurrence across multiple samples, and meaningful functional enrichment. Generic stress or ribosomal signals are not novelty evidence.
- LLM-assisted annotation may use marker summaries only. Do not send raw expression matrices or private metadata to external APIs.

## Downstream Gate

CellChat, pseudotime, CNV, deconvolution, and microenvironment interpretation require annotation status of `usable` or `needs_review_with_approved_scope`.

If annotation status is `not_usable`, the skill must stop and ask for corrected annotation or manual confirmation before generating downstream scripts.

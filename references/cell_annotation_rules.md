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
- LLM-assisted annotation may use marker summaries only. Do not send raw expression matrices or private metadata to external APIs.

## Downstream Gate

CellChat, pseudotime, CNV, deconvolution, and microenvironment interpretation require annotation status of `usable` or `needs_review_with_approved_scope`.

If annotation status is `not_usable`, the skill must stop and ask for corrected annotation or manual confirmation before generating downstream scripts.

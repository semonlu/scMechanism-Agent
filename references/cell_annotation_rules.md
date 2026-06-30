# Cell Annotation Rules

Marker annotation must combine canonical markers, tissue context, organism, and cluster-level expression patterns. Automated tools are supporting evidence, not final authority.

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

## Expected Evidence

- DotPlot/FeaturePlot for marker genes.
- Top marker table per cluster.
- Automated labels from SingleR/SCINA/TransferData/scPred/CellTypist when appropriate.
- `annotation_evidence.tsv` with label, marker support, automated support, and uncertainty.

## Guardrails

- Do not force precise subtype labels when evidence is weak.
- Do not upload raw matrices to external LLM/API annotation services without explicit approval.
- Prefer marker summaries over expression matrices for LLM-assisted annotation.

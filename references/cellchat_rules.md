# CellChat Rules

## When To Propose

Propose CellChat only when:

- Cell-type labels are credible and supported by an annotation evidence table.
- Each sender/receiver group has enough cells after filtering.
- Disease/control or other comparison groups are clear.
- The biological question involves microenvironment, ligand-receptor signaling, stromal-immune interaction, epithelial-immune interaction, or niche remodeling.

## Confirmation Gate

Do not run or render CellChat scripts immediately after clustering/annotation. First generate a proposal based on existing results:

```text
Recommended microenvironment:
Recommended sender/receiver groups:
Groups to exclude:
Comparison design:
Rationale from markers/proportions/DE/enrichment:
Minimum cell count filter:
Expected outputs:
Interpretation limits:
Question for user approval:
```

Only after the user agrees may the skill render:

```bash
python scripts/render_template.py --template scripts/course_adapted/03_cellchat_from_seurat.R --out run/03_cellchat_from_seurat.R --define SEURAT_RDS=objects/processed_seurat.rds --define OUTPUT_DIR=analysis/cellchat
```

## Recommended Scope Examples

- Immune-stromal microenvironment: macrophage, T cell, fibroblast, endothelial, pericyte.
- Inflammatory niche: monocyte/macrophage, dendritic cell, T/NK, fibroblast.
- Fibrotic remodeling: fibroblast, smooth muscle/pericyte, endothelial, macrophage.
- Epithelial injury niche: epithelial, immune, fibroblast, endothelial.

The selected scope must be driven by observed cell types and disease question, not by a generic all-cell default.

## Interpretation

CellChat results are computational ligand-receptor hypotheses. They do not prove physical interactions or causality. For disease comparisons, run equivalent preprocessing, annotation, and filtering in all groups.

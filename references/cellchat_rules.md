# CellChat Rules

## When To Run

Run CellChat only when:

- Cell-type labels are credible.
- Each sender/receiver group has enough cells.
- Normalized expression or raw counts are available in a compatible Seurat object.
- Organism-specific database can be selected.

## Primary Script

```bash
python scripts/render_template.py --template scripts/course_adapted/03_cellchat_from_seurat.R --out run/03_cellchat_from_seurat.R --define SEURAT_RDS=objects/processed_seurat.rds --define OUTPUT_DIR=analysis/cellchat
```

## Course-Derived Steps

- `createCellChat()`
- `addMeta()`
- `setIdent()`
- `CellChatDB.human` or `CellChatDB.mouse`
- `subsetDB()`
- `subsetData()`
- `identifyOverExpressedGenes()`
- `identifyOverExpressedInteractions()`
- `computeCommunProb()`
- `filterCommunication()`
- `computeCommunProbPathway()`
- `aggregateNet()`
- export ligand-receptor and pathway tables.

## Interpretation

CellChat results are computational ligand-receptor hypotheses. They do not prove physical interactions or causality. For disease comparisons, run equivalent preprocessing, annotation, and filtering in all groups.

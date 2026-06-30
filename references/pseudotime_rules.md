# Pseudotime Rules

## When To Run

Run pseudotime only for plausible continua such as differentiation, activation, disease progression, treatment response, or developmental transitions.

## Primary Script

```bash
python scripts/render_template.py --template scripts/course_adapted/04_monocle3_from_seurat.R --out run/04_monocle3_from_seurat.R --define SEURAT_RDS=objects/processed_seurat.rds --define OUTPUT_DIR=analysis/monocle3
```

## Course-Derived Steps

- Extract counts and metadata from Seurat.
- Build `new_cell_data_set()`.
- `preprocess_cds()`.
- `reduce_dimension()`.
- Optionally import Seurat UMAP coordinates for consistency.
- `cluster_cells()`.
- `learn_graph()`.
- `order_cells()`.
- Export pseudotime and trajectory genes.

## Guardrails

- Root choice must be biologically justified.
- Pseudotime is relative ordering, not measured time.
- Cross-sectional data do not prove causal transitions.
- Do not run trajectory on unrelated broad cell types just because UMAP has clusters.

# Pseudotime Rules

## When To Propose

Propose pseudotime only for plausible continua:

- differentiation or maturation.
- activation or exhaustion.
- developmental time course.
- epithelial or stromal state transition.
- disease progression within a related lineage.

Do not run pseudotime on unrelated broad cell types just because UMAP has clusters.

## Confirmation Gate

Do not run or render Monocle3 scripts immediately after clustering/annotation. First generate a proposal based on existing results:

```text
Recommended lineage/subset:
States included:
States excluded:
Root choice and biological reason:
Evidence from markers/stage metadata/UMAP/proportions:
Maximum cell count/downsampling:
Expected outputs:
Interpretation limits:
Question for user approval:
```

Only after the user agrees may the skill render:

```bash
python scripts/render_template.py --template scripts/course_adapted/04_monocle3_from_seurat.R --out run/04_monocle3_from_seurat.R --define SEURAT_RDS=objects/processed_seurat.rds --define OUTPUT_DIR=analysis/monocle3
```

## Guardrails

- Root choice must be biologically justified.
- Pseudotime is relative ordering, not measured time.
- Cross-sectional data do not prove causal transitions.
- Input should be a biologically related subset, not all cells by default.
- If only one state/cell type is present, report that trajectory analysis is not appropriate unless the user explicitly asks for a single-lineage state audit.

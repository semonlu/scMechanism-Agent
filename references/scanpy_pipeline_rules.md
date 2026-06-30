# Scanpy Pipeline Rules

## Runtime

- Python 3.10 or newer.
- Packages: `scanpy`, `anndata`, `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`, `leidenalg`, `igraph`.
- Optional: `scrublet`, `celltypist`, `gprofiler-official`, `scvelo`, `cellrank`, `squidpy`.

## Template

```bash
python scripts/render_template.py --template templates/scanpy_basic_pipeline_template.py --out run/scanpy_pipeline.py --define INPUT_PATH=/data/file.h5ad --define OUTPUT_DIR=analysis/run1
```

## Rules

- h5ad: `sc.read_h5ad()`, then inspect `.raw` and `.layers`.
- 10x MEX: `sc.read_10x_mtx()`.
- 10x H5: `sc.read_10x_h5()`.
- CSV/TSV: read with pandas, then confirm gene/cell orientation before creating AnnData.
- FASTQ/SRA: do not load directly; rebuild a matrix first.

## Outputs

- `objects/processed.h5ad`
- `tables/qc_summary.tsv`
- `tables/cluster_markers.tsv`
- `figures/umap_clusters.png`
- `logs/session_info.txt`

## Guardrails

- Do not run velocity unless spliced/unspliced layers exist.
- Do not use CellRank without velocity/transition information.
- Record raw-count availability before DE or enrichment.

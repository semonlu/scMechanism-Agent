# Scanpy Pipeline Rules

## Runtime

- Python 3.10 or newer.
- Packages: `scanpy`, `anndata`, `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`, `leidenalg`, `igraph`.
- Optional: `bbknn`, `scanorama`, `harmonypy`, `celltypist`, `gseapy`, `decoupler`, `liana`, `scvelo`, `cellrank`, `palantir`, `squidpy`.

## Template

```bash
python scripts/render_template.py --template templates/scanpy_basic_pipeline_template.py --out run/scanpy_pipeline.py --define INPUT_PATH=/data/file.h5ad --define OUTPUT_DIR=analysis/run1
python scripts/render_template.py --template templates/scanpy_batch_annotation_enrichment_template.py --out run/scanpy_optional_modules.py --define INPUT_H5AD=analysis/run1/objects/processed.h5ad --define OUTPUT_DIR=analysis/run1/optional
```

## Rules

- h5ad: `sc.read_h5ad()`, then inspect `.raw` and `.layers`.
- 10x MEX: `sc.read_10x_mtx()`.
- 10x H5: `sc.read_10x_h5()`.
- CSV/TSV: read with pandas, then confirm gene/cell orientation before creating AnnData.
- FASTQ/SRA: do not load directly; rebuild a matrix first.
- Loom: read with `sc.read_loom()`, then write h5ad and continue through the generic Scanpy workflow.

## Optional Modules

- Batch correction: prefer BBKNN when `bbknn` is installed and a valid batch column exists; otherwise use Scanpy ComBat as a conservative fallback. Mention Scanorama/Harmony Python as install-dependent alternatives.
- Annotation: use CellTypist when the model is appropriate; otherwise write marker tables for manual/reference-assisted annotation.
- Enrichment: use GSEApy/Enrichr when internet/package access is available; otherwise export ranked markers for later enrichment.
- Communication: mention LIANA/CellPhoneDB as optional Python routes; do not claim they ran unless result tables exist.
- Pseudotime/RNA velocity: only run scVelo/CellRank when spliced/unspliced layers or a valid transition model exists. Palantir can be suggested for diffusion-map pseudotime when inputs are appropriate.

## Outputs

- `objects/processed.h5ad`
- `tables/qc_summary.tsv`
- `tables/cluster_markers.tsv`
- `figures/umap_clusters.png`
- `logs/scanpy_versions.json`
- optional module logs under `logs/scanpy_optional_modules.json`

## Guardrails

- Do not run velocity unless spliced/unspliced layers exist.
- Do not use CellRank without velocity/transition information.
- Record raw-count availability before DE or enrichment.
- Ensure every figure call writes to `OUTPUT_DIR/figures`; `show=False` alone is not a saved output.

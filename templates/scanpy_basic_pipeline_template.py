#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import scanpy as sc


INPUT_PATH = "{{INPUT_PATH}}"
METADATA_PATH = "{{METADATA_PATH}}"
OUTPUT_DIR = Path("{{OUTPUT_DIR}}")
INPUT_TYPE = "{{INPUT_TYPE}}"  # h5ad, 10x_mtx, 10x_h5, csv
ORGANISM = "{{ORGANISM}}"

for sub in ["figures", "tables", "objects", "logs"]:
    (OUTPUT_DIR / sub).mkdir(parents=True, exist_ok=True)

if INPUT_TYPE == "h5ad":
    adata = sc.read_h5ad(INPUT_PATH)
elif INPUT_TYPE == "10x_mtx":
    adata = sc.read_10x_mtx(INPUT_PATH, var_names="gene_symbols", cache=False)
elif INPUT_TYPE == "10x_h5":
    adata = sc.read_10x_h5(INPUT_PATH)
elif INPUT_TYPE == "csv":
    matrix = pd.read_csv(INPUT_PATH, index_col=0)
    adata = sc.AnnData(matrix.T)
else:
    raise ValueError(f"Unsupported INPUT_TYPE: {INPUT_TYPE}")

if METADATA_PATH and Path(METADATA_PATH).exists():
    meta = pd.read_csv(METADATA_PATH, index_col=0)
    common = adata.obs_names.intersection(meta.index)
    adata = adata[common].copy()
    adata.obs = adata.obs.join(meta.loc[common])

adata.var_names_make_unique()
mt_prefix = "mt-" if ORGANISM.lower() == "mouse" else "MT-"
adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].to_csv(OUTPUT_DIR / "tables" / "qc_metrics.tsv", sep="\t")

sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], multi_panel=True, save=None, show=False)
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", show=False)

# Replace thresholds after reviewing current distributions.
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs["pct_counts_mt"] < 25].copy()

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata, resolution=0.5)
sc.tl.umap(adata)

sc.pl.umap(adata, color=["leiden"], save="_clusters.png", show=False)
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv(OUTPUT_DIR / "tables" / "cluster_markers.tsv", sep="\t", index=False)

adata.write_h5ad(OUTPUT_DIR / "objects" / "processed.h5ad")
(OUTPUT_DIR / "logs" / "scanpy_versions.json").write_text(json.dumps({"scanpy": sc.__version__}, indent=2), encoding="utf-8")

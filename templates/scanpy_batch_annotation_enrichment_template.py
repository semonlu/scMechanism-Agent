#!/usr/bin/env python3
"""Optional Scanpy modules: batch correction, CellTypist annotation, and enrichment."""

from __future__ import annotations

import json
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


INPUT_H5AD = "{{INPUT_H5AD}}"
OUTPUT_DIR = Path("{{OUTPUT_DIR}}")
BATCH_COL = "{{BATCH_COL}}"
ANNOTATION_METHOD = "{{ANNOTATION_METHOD}}"  # celltypist, marker_summary, none
CELLTYPIST_MODEL = "{{CELLTYPIST_MODEL}}"  # e.g. Immune_All_Low.pkl
ENRICHMENT_ORGANISM = "{{ENRICHMENT_ORGANISM}}"  # human, mouse, none

for sub in ["figures", "tables", "objects", "logs"]:
    (OUTPUT_DIR / sub).mkdir(parents=True, exist_ok=True)

sc.settings.figdir = str(OUTPUT_DIR / "figures")
adata = sc.read_h5ad(INPUT_H5AD)
adata.var_names_make_unique()

if "X_pca" not in adata.obsm:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")

batch_method = "none"
if BATCH_COL and BATCH_COL in adata.obs:
    if importlib.util.find_spec("bbknn") is not None:
        import bbknn

        bbknn.bbknn(adata, batch_key=BATCH_COL)
        batch_method = "bbknn"
    else:
        sc.pp.combat(adata, key=BATCH_COL)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
        batch_method = "combat"
else:
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

sc.tl.leiden(adata, resolution=0.5)
sc.tl.umap(adata)
sc.pl.umap(adata, color=["leiden"], save="_scanpy_batch_corrected_clusters.png", show=False)
plt.close("all")

annotation_status = {"method": ANNOTATION_METHOD, "status": "not_run"}
if ANNOTATION_METHOD == "celltypist":
    try:
        import celltypist

        predictions = celltypist.annotate(adata, model=CELLTYPIST_MODEL, majority_voting=True)
        adata.obs["celltypist_label"] = predictions.predicted_labels["majority_voting"].astype(str).values
        adata.obs[["celltypist_label"]].to_csv(OUTPUT_DIR / "tables" / "celltypist_labels.tsv", sep="\t")
        sc.pl.umap(adata, color=["celltypist_label"], save="_celltypist_labels.png", show=False)
        plt.close("all")
        annotation_status["status"] = "ok"
    except Exception as exc:
        annotation_status = {"method": "celltypist", "status": "failed", "error": str(exc)}
elif ANNOTATION_METHOD == "marker_summary":
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(OUTPUT_DIR / "tables" / "cluster_markers_for_annotation.tsv", sep="\t", index=False)
    annotation_status["status"] = "marker_table_written"

sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv(OUTPUT_DIR / "tables" / "scanpy_cluster_markers.tsv", sep="\t", index=False)

enrichment_status = {"organism": ENRICHMENT_ORGANISM, "status": "not_run"}
if ENRICHMENT_ORGANISM in {"human", "mouse"}:
    try:
        import gseapy as gp

        gene_sets = "GO_Biological_Process_2023" if ENRICHMENT_ORGANISM == "human" else "GO_Biological_Process_2023"
        rows = []
        for cluster, group in markers.groupby("group"):
            genes = group.sort_values("scores", ascending=False)["names"].head(100).dropna().astype(str).tolist()
            if len(genes) < 5:
                continue
            enr = gp.enrichr(gene_list=genes, gene_sets=gene_sets, organism=ENRICHMENT_ORGANISM)
            table = enr.results.copy()
            table.insert(0, "cluster", cluster)
            rows.append(table)
        if rows:
            pd.concat(rows, ignore_index=True).to_csv(OUTPUT_DIR / "tables" / "gseapy_enrichment.tsv", sep="\t", index=False)
        enrichment_status["status"] = "ok"
    except Exception as exc:
        enrichment_status = {"organism": ENRICHMENT_ORGANISM, "status": "failed", "error": str(exc)}

adata.write_h5ad(OUTPUT_DIR / "objects" / "scanpy_batch_annotation_enrichment.h5ad")
(OUTPUT_DIR / "logs" / "scanpy_optional_modules.json").write_text(
    json.dumps(
        {
            "scanpy": sc.__version__,
            "batch_method": batch_method,
            "annotation": annotation_status,
            "enrichment": enrichment_status,
        },
        indent=2,
    ),
    encoding="utf-8",
)

#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


INPUT_PATH = "{{INPUT_PATH}}"
METADATA_PATH = "{{METADATA_PATH}}"
OUTPUT_DIR = Path("{{OUTPUT_DIR}}")
INPUT_TYPE = "{{INPUT_TYPE}}"  # h5ad, 10x_mtx, 10x_h5, csv
ORGANISM = "{{ORGANISM}}"
RESOLUTION_VALUES = [0.1, 0.3, 0.5, 0.8]

for sub in ["figures", "tables", "objects", "logs"]:
    (OUTPUT_DIR / sub).mkdir(parents=True, exist_ok=True)

sc.settings.figdir = str(OUTPUT_DIR / "figures")


def resolution_note(total_cells: int, cluster_count: int) -> str:
    if 5000 <= total_cells <= 20000:
        if cluster_count < 8:
            return "possible_under_clustering_for_major_cell_types"
        if cluster_count > 25:
            return "possible_over_clustering_or_state_splitting"
        if 10 <= cluster_count <= 20:
            return "typical_major_lineage_range_for_10x_5k_20k_cells"
    if cluster_count <= 1:
        return "single_cluster_or_no_clear_structure"
    if cluster_count > 25:
        return "many_clusters_review_marker_specificity"
    return "review_marker_clarity_and_biology"


def select_resolution(summary: pd.DataFrame) -> float:
    preference = [0.3, 0.5, 0.1, 0.8]
    total_cells = int(summary["total_cells"].iloc[0])
    if 5000 <= total_cells <= 20000:
        candidates = summary[(summary["cluster_count"] >= 10) & (summary["cluster_count"] <= 20)]
        for res in preference:
            if res in set(candidates["resolution"]):
                return res
    candidates = summary[(summary["cluster_count"] > 1) & (summary["cluster_count"] <= 25)]
    for res in preference:
        if res in set(candidates["resolution"]):
            return res
    return float(summary["resolution"].iloc[0])


def run_resolution_sweep(adata) -> float:
    rows = []
    for res in RESOLUTION_VALUES:
        key = f"leiden_res_{str(res).replace('.', '_')}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        counts = adata.obs[key].value_counts()
        rows.append(
            {
                "resolution": res,
                "cluster_count": int(counts.shape[0]),
                "total_cells": int(adata.n_obs),
                "min_cluster_cells": int(counts.min()),
                "median_cluster_cells": float(counts.median()),
                "max_cluster_cells": int(counts.max()),
                "review_note": resolution_note(int(adata.n_obs), int(counts.shape[0])),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(OUTPUT_DIR / "tables" / "resolution_sweep.tsv", sep="\t", index=False)
    selected = select_resolution(summary)
    selected_key = f"leiden_res_{str(selected).replace('.', '_')}"
    adata.obs["leiden"] = adata.obs[selected_key].astype("category")
    pd.DataFrame(
        [{"parameter": "resolution_values", "value": ",".join(map(str, RESOLUTION_VALUES))}, {"parameter": "selected_resolution", "value": selected}]
    ).to_csv(OUTPUT_DIR / "tables" / "clustering_parameters.tsv", sep="\t", index=False)
    return selected


def write_cluster_marker_audit(adata, markers: pd.DataFrame) -> None:
    generic_prefixes = ("MT-", "mt-", "RPS", "RPL", "Rps", "Rpl", "HBA", "HBB", "HBD", "HSP", "Hsp")
    cell_cycle_genes = {"MKI67", "TOP2A", "UBE2C", "CCNB1", "CCNB2", "CENPF", "PCNA", "MCM2", "MCM5"}
    stress_genes = {"FOS", "JUN", "JUNB", "DUSP1", "HSPA1A", "HSPA1B", "HSP90AA1", "IER2", "ATF3"}
    rows = []
    for cluster in sorted(adata.obs["leiden"].astype(str).unique()):
        group = markers[markers["group"].astype(str) == str(cluster)].copy()
        if "scores" in group:
            group = group.sort_values("scores", ascending=False)
        top_genes = group["names"].dropna().astype(str).head(10).tolist() if "names" in group else []
        generic_count = sum(gene.startswith(generic_prefixes) for gene in top_genes)
        cell_cycle_count = sum(gene.upper() in cell_cycle_genes for gene in top_genes)
        stress_count = sum(gene.upper() in stress_genes for gene in top_genes)
        cell_count = int((adata.obs["leiden"].astype(str) == str(cluster)).sum())
        pct_total = round(100 * cell_count / max(1, adata.n_obs), 3)
        notes = []
        if len(group) < 100:
            notes.append("few_positive_markers_consider_unknown_or_merge_if_unclear")
        if len(group) > 300:
            notes.append("many_positive_markers_possible_functional_state_review")
        if generic_count >= 5:
            notes.append("top_markers_generic_ribo_mt_hb_or_heat_shock_review_qc")
        if cell_cycle_count >= 2:
            notes.append("cell_cycle_signature_review_s_g2m_scores")
        if stress_count >= 2:
            notes.append("stress_signature_review_dissociation_or_low_quality")
        if pct_total < 1:
            notes.append("rare_cluster_under_1_percent_requires_recurrence_before_novelty_claim")
        rows.append(
            {
                "cluster": cluster,
                "cell_count": cell_count,
                "pct_total": pct_total,
                "marker_count_vs_all": int(len(group)),
                "top_markers": ";".join(top_genes),
                "generic_top_marker_count": generic_count,
                "cell_cycle_top_marker_count": cell_cycle_count,
                "stress_top_marker_count": stress_count,
                "review_note": ";".join(dict.fromkeys(notes)) or "marker_review_required",
            }
        )
    pd.DataFrame(rows).to_csv(OUTPUT_DIR / "tables" / "cluster_marker_audit.tsv", sep="\t", index=False)

if INPUT_TYPE == "h5ad":
    adata = sc.read_h5ad(INPUT_PATH)
elif INPUT_TYPE == "10x_mtx":
    adata = sc.read_10x_mtx(INPUT_PATH, var_names="gene_symbols", cache=False)
elif INPUT_TYPE == "10x_h5":
    adata = sc.read_10x_h5(INPUT_PATH)
elif INPUT_TYPE == "csv":
    matrix = pd.read_csv(INPUT_PATH, index_col=0, encoding="utf-8-sig")
    adata = sc.AnnData(matrix.T)
else:
    raise ValueError(f"Unsupported INPUT_TYPE: {INPUT_TYPE}")

if METADATA_PATH and Path(METADATA_PATH).exists():
    meta = pd.read_csv(METADATA_PATH, index_col=0, encoding="utf-8-sig")
    common = adata.obs_names.intersection(meta.index)
    adata = adata[common].copy()
    adata.obs = adata.obs.join(meta.loc[common])

adata.var_names_make_unique()
mt_prefix = "mt-" if ORGANISM.lower() == "mouse" else "MT-"
hb_pattern = r"^Hb[abq]" if ORGANISM.lower() == "mouse" else r"^HB[ABDEGQMZ]"
adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
adata.var["hb"] = adata.var_names.str.match(hb_pattern)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "hb"], inplace=True)

qc_cols = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_hb"]
adata.obs[qc_cols].to_csv(OUTPUT_DIR / "tables" / "qc_metrics.tsv", sep="\t")

sc.pl.violin(
    adata,
    qc_cols,
    multi_panel=True,
    save="_qc_pre_filter.png",
    show=False,
)
plt.close("all")
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", save="_qc_counts_mt.png", show=False)
plt.close("all")

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
run_resolution_sweep(adata)
sc.tl.umap(adata)

sc.pl.umap(adata, color=["leiden"], save="_clusters.png", show=False)
plt.close("all")
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv(OUTPUT_DIR / "tables" / "cluster_markers.tsv", sep="\t", index=False)
write_cluster_marker_audit(adata, markers)

adata.write_h5ad(OUTPUT_DIR / "objects" / "processed.h5ad")
(OUTPUT_DIR / "logs" / "scanpy_versions.json").write_text(json.dumps({"scanpy": sc.__version__}, indent=2), encoding="utf-8")

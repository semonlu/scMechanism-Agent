#!/usr/bin/env python3
"""Whitelisted Scanpy basic workflow for the local MCP backend.

This script is launched only by local_mcp_server.py after input/output paths
have been constrained to mcp_server/workspace. It does not accept shell
commands and only writes to the selected project output directory.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

RESOLUTION_VALUES = [0.1, 0.3, 0.5, 0.8]


def require_scanpy():
    if importlib.util.find_spec("scanpy") is None:
        raise RuntimeError("scanpy is not installed. Create the environment from environment.yml first.")
    import matplotlib.pyplot as plt
    import pandas as pd
    import scanpy as sc

    return sc, pd, plt


def read_input(sc, input_path: Path, input_type: str):
    if input_type == "h5ad":
        return sc.read_h5ad(input_path)
    if input_type == "10x_mtx":
        return sc.read_10x_mtx(input_path, var_names="gene_symbols", cache=False)
    if input_type == "10x_h5":
        return sc.read_10x_h5(input_path)
    raise ValueError("scanpy_basic supports h5ad, 10x_mtx, and 10x_h5 inputs.")


def merge_metadata(pd, adata, metadata_path: str, sample_id: str, group_col: str) -> tuple[str, str, str]:
    sample_id = sample_id or Path(metadata_path).stem if metadata_path else sample_id
    sample_id = sample_id or "sample_1"
    if metadata_path:
        meta = pd.read_csv(metadata_path, encoding="utf-8-sig")
        if "cell" in meta.columns:
            meta = meta.set_index("cell")
        elif "barcode" in meta.columns:
            meta = meta.set_index("barcode")
        else:
            meta = meta.set_index(meta.columns[0])
        common = adata.obs_names.intersection(meta.index.astype(str))
        if len(common) == 0 and len(meta) == 1:
            row = meta.iloc[0]
            for col, value in row.items():
                adata.obs[col] = value
            if not sample_id and "sample_id" in row.index:
                sample_id = str(row["sample_id"])
        elif len(common) == 0:
            raise ValueError("metadata_path was provided, but no metadata row names match cell barcodes.")
        else:
            meta = meta.loc[common].copy()
            adata._inplace_subset_obs(common)
            for col in meta.columns:
                adata.obs[col] = meta[col].values
    if "sample_id" not in adata.obs:
        adata.obs["sample_id"] = sample_id
    if "group" not in adata.obs:
        if group_col and group_col in adata.obs:
            adata.obs["group"] = adata.obs[group_col].astype(str)
        else:
            adata.obs["group"] = adata.obs["sample_id"].astype(str)
    if "batch" not in adata.obs:
        adata.obs["batch"] = adata.obs["sample_id"].astype(str)
    return sample_id, group_col if group_col in adata.obs else "group", "metadata_merged" if metadata_path else "default_metadata_created"


def write_status(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


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


def select_resolution(pd, summary) -> float:
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


def run_resolution_sweep(sc, pd, adata, output_dir: Path) -> float:
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
    summary.to_csv(output_dir / "tables" / "resolution_sweep.tsv", sep="\t", index=False)
    selected = select_resolution(pd, summary)
    selected_key = f"leiden_res_{str(selected).replace('.', '_')}"
    adata.obs["leiden"] = adata.obs[selected_key].astype("category")
    pd.DataFrame(
        [{"parameter": "resolution_values", "value": ",".join(map(str, RESOLUTION_VALUES))}, {"parameter": "selected_resolution", "value": selected}]
    ).to_csv(output_dir / "tables" / "clustering_parameters.tsv", sep="\t", index=False)
    return selected


def write_cluster_marker_audit(pd, adata, markers, output_dir: Path) -> None:
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
    pd.DataFrame(rows).to_csv(output_dir / "tables" / "cluster_marker_audit.tsv", sep="\t", index=False)


CANONICAL_MARKERS = {
    "T_NK_cell": {"CD3D", "CD3E", "CD2", "TRAC", "NKG7", "GNLY", "PRF1", "GZMB", "KLRD1"},
    "B_cell": {"MS4A1", "CD79A", "CD79B", "BANK1", "CD74", "HLA-DRA"},
    "Plasma_cell": {"MZB1", "JCHAIN", "XBP1", "SDC1", "IGHG1", "IGHG3", "IGKC"},
    "Myeloid_macrophage": {"LYZ", "LST1", "AIF1", "CD68", "C1QA", "C1QB", "FCGR3A", "S100A8", "S100A9"},
    "Dendritic_cell": {"FCER1A", "CLEC10A", "CD1C", "LILRA4", "IRF8"},
    "Endothelial_cell": {"PECAM1", "VWF", "KDR", "CLDN5", "ENG", "ESAM", "RAMP2"},
    "Fibroblast_CAF": {"COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "FAP", "THY1"},
    "Osteoblastic_tumor_like": {"RUNX2", "ALPL", "BGLAP", "IBSP", "SPP1", "POSTN", "SPARC", "MMP13", "COL10A1"},
    "Pericyte_smooth_muscle": {"RGS5", "PDGFRB", "MCAM", "ACTA2", "TAGLN", "MYH11", "NOTCH3"},
    "Mast_cell": {"TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"},
    "Cycling_cell": {"MKI67", "TOP2A", "PCNA", "STMN1", "UBE2C", "HMGB2"},
    "Epithelial_like": {"EPCAM", "KRT8", "KRT18", "KRT19", "KRT7"},
}


def infer_cluster_label(top_genes: list[str]) -> dict:
    upper_genes = [str(gene).upper() for gene in top_genes]
    top_set = set(upper_genes[:50])
    scores = {
        label: len(top_set.intersection(markers))
        for label, markers in CANONICAL_MARKERS.items()
    }
    ranked = sorted(scores.items(), key=lambda item: item[1], reverse=True)
    best_label, best_score = ranked[0]
    second_label, second_score = ranked[1] if len(ranked) > 1 else ("", 0)
    if best_score == 0:
        return {
            "label": "",
            "coarse": "unassigned",
            "confidence": "low",
            "support": "no canonical marker match in top markers",
            "conflict": "",
        }
    support = sorted(top_set.intersection(CANONICAL_MARKERS[best_label]))
    conflict = ""
    if second_score > 0 and second_score >= max(2, best_score - 1):
        conflict_markers = sorted(top_set.intersection(CANONICAL_MARKERS[second_label]))
        conflict = f"{second_label}:{','.join(conflict_markers)}"
    confidence = "high" if best_score >= 4 and second_score <= best_score - 2 else "medium" if best_score >= 2 else "low"
    return {
        "label": best_label,
        "coarse": best_label,
        "confidence": confidence,
        "support": ",".join(support),
        "conflict": conflict,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-path", required=True)
    parser.add_argument("--input-type", required=True, choices=["h5ad", "10x_mtx", "10x_h5"])
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--species", default="human")
    parser.add_argument("--batch-col", default="")
    parser.add_argument("--metadata-path", default="")
    parser.add_argument("--sample-id", default="")
    parser.add_argument("--group-col", default="")
    parser.add_argument("--annotation-method", default="marker_summary", choices=["marker_summary", "celltypist", "none"])
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    for sub in ["figures", "tables", "objects", "logs"]:
        (output_dir / sub).mkdir(parents=True, exist_ok=True)
    status_file = output_dir / "logs" / "scanpy_basic_status.json"

    try:
        sc, pd, plt = require_scanpy()
        adata = read_input(sc, Path(args.input_path), args.input_type)
        adata.var_names_make_unique()
        sample_id, effective_group_col, metadata_status = merge_metadata(pd, adata, args.metadata_path, args.sample_id, args.group_col)
        adata.obs[["sample_id", "group", "batch"]].to_csv(output_dir / "tables" / "cell_metadata_minimal.tsv", sep="\t")

        mt_prefix = "mt-" if args.species.lower() == "mouse" else "MT-"
        hb_pattern = r"^Hb[abq]" if args.species.lower() == "mouse" else r"^HB[ABDEGQMZ]"
        adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
        adata.var["hb"] = adata.var_names.str.match(hb_pattern)
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "hb"], percent_top=None, log1p=False, inplace=True)
        qc_cols = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_hb"]
        adata.obs[qc_cols].to_csv(output_dir / "tables" / "qc_metrics.tsv", sep="\t")

        sc.settings.figdir = str(output_dir / "figures")
        sc.pl.violin(adata, qc_cols, jitter=0.4, multi_panel=True, save="_qc_prefilter.png", show=False)
        sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", save="_qc_counts_mt.png", show=False)
        plt.close("all")

        min_genes = min(200, max(1, int(adata.n_vars * 0.05)))
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=3)
        adata = adata[adata.obs["pct_counts_mt"] < 25].copy()
        if adata.n_obs == 0 or adata.n_vars == 0:
            raise RuntimeError(
                f"QC filtering removed all data. Try lower thresholds or inspect input quality. "
                f"Adaptive min_genes={min_genes}, cells={adata.n_obs}, genes={adata.n_vars}."
            )

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata
        n_top_genes = min(2000, max(1, adata.n_vars))
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
        adata = adata[:, adata.var["highly_variable"]].copy()
        if adata.n_obs < 3 or adata.n_vars < 3:
            raise RuntimeError(f"Too few cells or genes after filtering for PCA/clustering: cells={adata.n_obs}, genes={adata.n_vars}.")
        sc.pp.scale(adata, max_value=10)
        n_comps = min(50, adata.n_obs - 1, adata.n_vars - 1)
        sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")
        n_pcs = min(30, n_comps)

        batch_method = "none"
        batch_col = args.batch_col if args.batch_col and args.batch_col in adata.obs else ""
        if batch_col:
            if importlib.util.find_spec("bbknn") is not None:
                import bbknn

                bbknn.bbknn(adata, batch_key=batch_col, n_pcs=n_pcs)
                batch_method = "bbknn"
            else:
                sc.pp.combat(adata, key=batch_col)
                sc.pp.neighbors(adata, n_neighbors=min(15, max(2, adata.n_obs - 1)), n_pcs=n_pcs)
                batch_method = "combat"
        else:
            sc.pp.neighbors(adata, n_neighbors=min(15, max(2, adata.n_obs - 1)), n_pcs=n_pcs)

        selected_resolution = run_resolution_sweep(sc, pd, adata, output_dir)
        sc.tl.umap(adata)
        sc.pl.umap(adata, color=["leiden"], save="_clusters.png", show=False)
        plt.close("all")

        annotation_status = {"method": args.annotation_method, "status": "not_run"}
        if args.annotation_method == "celltypist":
            try:
                import celltypist

                predictions = celltypist.annotate(adata, majority_voting=True)
                adata.obs["celltypist_label"] = predictions.predicted_labels["majority_voting"].astype(str).values
                adata.obs[["celltypist_label"]].to_csv(output_dir / "tables" / "celltypist_labels.tsv", sep="\t")
                sc.pl.umap(adata, color=["celltypist_label"], save="_celltypist_labels.png", show=False)
                plt.close("all")
                annotation_status = {"method": "celltypist", "status": "ok"}
            except Exception as exc:
                annotation_status = {"method": "celltypist", "status": "failed", "error": str(exc)}

        sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
        markers = sc.get.rank_genes_groups_df(adata, group=None)
        if "group" not in markers.columns:
            clusters = sorted(adata.obs["leiden"].astype(str).unique().tolist())
            markers.insert(0, "group", clusters[0] if clusters else "0")
        markers.to_csv(output_dir / "tables" / "scanpy_cluster_markers.tsv", sep="\t", index=False)
        write_cluster_marker_audit(pd, adata, markers, output_dir)

        evidence_rows = []
        cluster_to_label = {}
        for cluster, group in markers.groupby("group"):
            top = group.sort_values("scores", ascending=False).head(10)
            top_genes = group.sort_values("scores", ascending=False)["names"].head(50).astype(str).tolist()
            inferred = infer_cluster_label(top_genes)
            final_label = inferred["label"] or f"Cluster_{cluster}"
            cluster_to_label[str(cluster)] = final_label
            evidence_rows.append(
                {
                    "cluster": cluster,
                    "final_label": final_label,
                    "coarse_label": inferred["coarse"],
                    "singleR_label": "",
                    "singleR_pruned_label": "",
                    "singleR_delta_next": "",
                    "top_markers": ";".join(top["names"].astype(str).tolist()),
                    "canonical_marker_support": inferred["support"],
                    "conflicting_markers": inferred["conflict"],
                    "cell_count": int((adata.obs["leiden"].astype(str) == str(cluster)).sum()),
                    "confidence": inferred["confidence"],
                    "review_note": "Broad marker-rule annotation from Scanpy cluster markers; review manually before manuscript or downstream modules.",
                }
            )
        pd.DataFrame(evidence_rows).to_csv(output_dir / "tables" / "annotation_evidence.tsv", sep="\t", index=False)
        adata.obs["marker_rule_label"] = adata.obs["leiden"].astype(str).map(cluster_to_label).astype(str)
        adata.obs[["marker_rule_label"]].to_csv(output_dir / "tables" / "marker_rule_labels.tsv", sep="\t")
        sc.pl.umap(adata, color=["marker_rule_label"], save="_marker_rule_labels.png", show=False)
        plt.close("all")

        enrichment_status = {"status": "not_run"}
        if importlib.util.find_spec("gseapy") is not None:
            try:
                import gseapy as gp

                rows = []
                organism = "mouse" if args.species.lower() == "mouse" else "human"
                for cluster, group in markers.groupby("group"):
                    genes = group.sort_values("scores", ascending=False)["names"].head(100).dropna().astype(str).tolist()
                    if len(genes) < 5:
                        continue
                    enr = gp.enrichr(gene_list=genes, gene_sets="GO_Biological_Process_2023", organism=organism)
                    table = enr.results.copy()
                    table.insert(0, "cluster", cluster)
                    rows.append(table)
                if rows:
                    pd.concat(rows, ignore_index=True).to_csv(output_dir / "tables" / "gseapy_enrichment.tsv", sep="\t", index=False)
                enrichment_status = {"status": "ok"}
            except Exception as exc:
                enrichment_status = {"status": "failed", "error": str(exc)}

        adata.write_h5ad(output_dir / "objects" / "scanpy_basic_processed.h5ad")
        status = {
            "status": "ok",
            "scanpy": sc.__version__,
            "input_type": args.input_type,
            "species": args.species,
            "sample_id": sample_id,
            "effective_group_col": effective_group_col,
            "metadata_status": metadata_status,
            "batch_method": batch_method,
            "adaptive_min_genes": min_genes,
            "n_top_genes": n_top_genes,
            "n_pcs": n_pcs,
            "selected_resolution": selected_resolution,
            "annotation": annotation_status,
            "enrichment": enrichment_status,
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
        }
        write_status(status_file, status)
        print(json.dumps(status, ensure_ascii=False, indent=2))
        return 0
    except Exception as exc:
        status = {"status": "error", "message": str(exc), "input_path": args.input_path, "input_type": args.input_type}
        write_status(status_file, status)
        print(json.dumps(status, ensure_ascii=False, indent=2))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

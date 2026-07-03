#!/usr/bin/env python3
"""Build a table-backed manuscript-style report from a scMechanism result bundle."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median


def read_table(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def as_float(value: str, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def format_float(value: float, digits: int = 3) -> str:
    if math.isnan(value):
        return "NA"
    return f"{value:.{digits}g}"


def top_rows(rows: list[dict[str, str]], key: str, n: int = 5, reverse: bool = False) -> list[dict[str, str]]:
    return sorted(rows, key=lambda row: as_float(row.get(key, "")), reverse=reverse)[:n]


def rel(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


def portable_display_path(path: Path) -> str:
    parts = list(path.parts)
    for anchor in ("submission", "results", "analysis"):
        if anchor in parts:
            return Path(*parts[parts.index(anchor):]).as_posix()
    return path.name


def load_metadata(path: Path | None) -> dict[str, str]:
    defaults = {
        "dataset": "unspecified dataset",
        "accession": "unspecified accession",
        "organism": "unspecified organism",
        "tissue": "unspecified tissue",
        "scope_note": "Public single-cell reanalysis; interpret as exploratory mechanism evidence.",
        "pipeline": "Seurat V5, SingleR, cluster marker scoring, GO/KEGG enrichment, CellChat, and Monocle3",
    }
    if path and path.exists():
        with path.open("r", encoding="utf-8") as handle:
            defaults.update(json.load(handle))
    return defaults


def summarize_qc(result_dir: Path) -> tuple[list[dict[str, str]], int]:
    rows = read_table(result_dir / "gse223751_seurat" / "tables" / "qc_summary_post_filter.tsv")
    cells = sum(int(as_float(row.get("cells", "0"), 0)) for row in rows)
    return rows, cells


def summarize_celltypes(result_dir: Path) -> tuple[Counter, list[dict[str, str]]]:
    rows = read_table(result_dir / "gse223751_seurat" / "tables" / "celltype_proportions.tsv")
    counts: Counter[str] = Counter()
    for row in rows:
        label = row.get("singleR_label", "unlabeled")
        counts[label] += int(as_float(row.get("cells", "0"), 0))
    return counts, rows


def summarize_pseudotime(result_dir: Path) -> tuple[Counter, dict[str, float], int, int]:
    rows = read_table(result_dir / "course_modules" / "monocle3" / "tables" / "pseudotime.csv")
    counts: Counter[str] = Counter()
    values: defaultdict[str, list[float]] = defaultdict(list)
    finite_values: list[float] = []
    for row in rows:
        state = row.get("marker_support_label") or row.get("singleR_label") or "unlabeled"
        pt = as_float(row.get("pseudotime", ""))
        counts[state] += 1
        if math.isfinite(pt):
            values[state].append(pt)
            finite_values.append(pt)
    medians = {state: median(vals) for state, vals in values.items() if vals}
    unique = len(set(round(value, 9) for value in finite_values))
    return counts, medians, len(rows), unique


def make_report(result_dir: Path, metadata: dict[str, str]) -> str:
    qc_rows, retained_cells = summarize_qc(result_dir)
    celltype_counts, _ = summarize_celltypes(result_dir)
    pseudotime_counts, pseudotime_medians, trajectory_cells, unique_pseudotime = summarize_pseudotime(result_dir)
    go_rows = top_rows(read_table(result_dir / "course_modules" / "marker_enrichment" / "tables" / "GO_BP_top_markers.csv"), "p.adjust", 5)
    kegg_rows = top_rows(read_table(result_dir / "course_modules" / "marker_enrichment" / "tables" / "KEGG_top_markers.csv"), "p.adjust", 5)
    lr_rows = top_rows(read_table(result_dir / "course_modules" / "cellchat" / "tables" / "cellchat_ligand_receptor.csv"), "prob", 5, reverse=True)

    figures = [
        "gse223751_seurat/figures/umap_clusters.pdf",
        "gse223751_seurat/figures/umap_singleR_labels.pdf",
        "gse223751_seurat/figures/umap_marker_support_labels.pdf",
        "gse223751_seurat/figures/marker_dotplot.pdf",
        "course_modules/cellchat/figures/cellchat_interaction_count.pdf",
        "course_modules/cellchat/figures/cellchat_interaction_weight.pdf",
        "course_modules/monocle3/figures/trajectory_celltypes.pdf",
        "course_modules/monocle3/figures/trajectory_celltypes_with_graph.pdf",
        "course_modules/monocle3/figures/pseudotime_umap.pdf",
    ]
    existing_figures = [figure for figure in figures if (result_dir / figure).exists()]

    lines: list[str] = []
    lines.append(f"# Manuscript-Style Report: {metadata['accession']} Single-Cell Reanalysis")
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append(f"- Dataset: {metadata['dataset']}")
    lines.append(f"- Accession: {metadata['accession']}")
    lines.append(f"- Organism: {metadata['organism']}")
    lines.append(f"- Tissue/context: {metadata['tissue']}")
    lines.append(f"- Scope note: {metadata['scope_note']}")
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append(
        f"Public processed single-cell RNA-seq matrices were analyzed locally using {metadata['pipeline']}. "
        "Sample-level 10x matrices were imported into Seurat, merged with sample metadata, filtered after QC inspection, normalized, scaled, reduced by PCA/UMAP, and clustered. "
        "Reference-assisted annotation used SingleR, and marker-support labels were computed from curated marker sets. "
        "Cluster markers were used for GO/KEGG enrichment. Cell-cell communication was inferred with CellChat. "
        "Trajectory analysis used Monocle3 with reviewed Seurat UMAP coordinates, a marker-supported multi-state trajectory input, graph learning, and an E15 root rule."
    )
    lines.append("")
    lines.append("## Results")
    lines.append("")
    lines.append("### Data Processing And QC")
    lines.append("")
    lines.append(f"The processed Seurat object retained {retained_cells:,} cells across {len(qc_rows)} sample entries.")
    if qc_rows:
        stage_counts = Counter()
        for row in qc_rows:
            stage_counts[row.get("stage", "NA")] += int(as_float(row.get("cells", "0"), 0))
        stage_summary = "; ".join(f"{stage} {count:,}" for stage, count in sorted(stage_counts.items()))
        lines.append(f"Retained cells by stage were: {stage_summary}.")
    lines.append("")
    lines.append("### Cell Annotation")
    lines.append("")
    if celltype_counts:
        total = sum(celltype_counts.values())
        top_summary = "; ".join(
            f"{label} {count:,} ({count / total:.1%})" for label, count in celltype_counts.most_common(6)
        )
        lines.append(f"SingleR annotation identified the following major labels: {top_summary}.")
        lines.append("These labels should be reviewed with marker evidence and tissue context before biological claims are finalized.")
    else:
        lines.append("Cell annotation tables were not found.")
    lines.append("")
    lines.append("### Marker Enrichment")
    lines.append("")
    if go_rows:
        lines.append("Top GO biological-process terms among selected marker genes included:")
        for row in go_rows:
            lines.append(f"- {row.get('Description', 'NA')} (adjusted p={row.get('p.adjust', 'NA')}, genes={row.get('Count', 'NA')})")
    if kegg_rows:
        lines.append("")
        lines.append("Top KEGG terms included:")
        for row in kegg_rows:
            lines.append(f"- {row.get('Description', 'NA')} (adjusted p={row.get('p.adjust', 'NA')}, genes={row.get('Count', 'NA')})")
    if not go_rows and not kegg_rows:
        lines.append("No enrichment result table was detected.")
    lines.append("")
    lines.append("### Cell-Cell Communication")
    lines.append("")
    if lr_rows:
        lines.append("CellChat inferred ligand-receptor interactions, including the highest-probability exported interactions:")
        for row in lr_rows:
            lines.append(
                f"- {row.get('source', 'NA')} -> {row.get('target', 'NA')}: "
                f"{row.get('ligand', 'NA')} / {row.get('receptor', 'NA')} "
                f"({row.get('pathway_name', 'NA')}, prob={row.get('prob', 'NA')})"
            )
        lines.append("These interactions are database-based computational hypotheses, not experimental proof of signaling.")
    else:
        lines.append("No CellChat ligand-receptor table was detected.")
    lines.append("")
    lines.append("### Pseudotime And Trajectory")
    lines.append("")
    if trajectory_cells:
        state_summary = "; ".join(f"{state} {count:,}" for state, count in pseudotime_counts.most_common())
        median_summary = "; ".join(
            f"{state} median {format_float(value)}" for state, value in sorted(pseudotime_medians.items(), key=lambda item: item[1])
        )
        lines.append(
            f"Monocle3 trajectory analysis used {trajectory_cells:,} cells with {unique_pseudotime:,} unique finite pseudotime values. "
            f"Trajectory input states were: {state_summary}."
        )
        lines.append(f"Median pseudotime by state was: {median_summary}.")
        lines.append(
            "The trajectory should be interpreted as relative transcriptional ordering across marker-supported states, not measured chronological time or causal proof."
        )
    else:
        lines.append("No pseudotime table was detected.")
    lines.append("")
    lines.append("## Mechanism Hypothesis")
    lines.append("")
    lines.append(
        "Together, marker enrichment, inferred communication, and trajectory ordering support an exploratory hypothesis that mesenchymal, chondrogenic/osteogenic, and proliferative states are organized along a developmental enthesis-associated continuum. "
        "This hypothesis is dataset- and annotation-dependent and requires validation in independent shoulder/enthesis or disease-relevant cohorts."
    )
    lines.append("")
    lines.append("## Figure Legends")
    lines.append("")
    for i, figure in enumerate(existing_figures, start=1):
        if "trajectory_celltypes.pdf" in figure and "with_graph" not in figure:
            legend = "UMAP of marker-supported trajectory input states used for Monocle3 analysis."
        elif "pseudotime_umap" in figure:
            legend = "Monocle3 pseudotime projected onto the imported Seurat UMAP with principal graph overlay."
        elif "cellchat" in figure:
            legend = "CellChat network summary from exported ligand-receptor inference."
        elif "marker" in figure:
            legend = "Marker-gene visualization supporting cluster or state annotation."
        else:
            legend = "Single-cell visualization generated by the workflow."
        lines.append(f"- Figure {i}. `{figure}`: {legend}")
    lines.append("")
    lines.append("## Supplementary Tables")
    lines.append("")
    table_paths = [
        "gse223751_seurat/tables/qc_summary_post_filter.tsv",
        "gse223751_seurat/tables/singleR_cluster_labels.csv",
        "gse223751_seurat/tables/annotation_evidence.tsv",
        "gse223751_seurat/tables/cluster_markers.csv",
        "course_modules/marker_enrichment/tables/GO_BP_top_markers.csv",
        "course_modules/marker_enrichment/tables/KEGG_top_markers.csv",
        "course_modules/cellchat/tables/cellchat_ligand_receptor.csv",
        "course_modules/monocle3/tables/pseudotime.csv",
        "course_modules/monocle3/tables/trajectory_genes.csv",
    ]
    for table in table_paths:
        if (result_dir / table).exists():
            lines.append(f"- `{table}`")
    lines.append("")
    lines.append("## Limitations")
    lines.append("")
    lines.append("- This is a public processed-data reanalysis and does not replace a prospectively designed clinical study.")
    lines.append("- The example dataset is shoulder-region rotator cuff enthesis biology, not a direct frozen-shoulder case-control cohort.")
    lines.append("- SingleR and marker-support annotations require manual review and independent validation.")
    lines.append("- Enrichment and CellChat analyses depend on database coverage and exported marker/interaction tables.")
    lines.append("- Pseudotime is a relative ordering inferred from expression similarity and graph structure, not a direct time measurement.")
    lines.append("- Mechanism hypotheses should be validated in independent cohorts and targeted experiments.")
    lines.append("")
    lines.append("## Validation Suggestions")
    lines.append("")
    lines.append("- Review marker expression for each trajectory state and adjust state labels if tissue expertise suggests alternatives.")
    lines.append("- Repeat key marker, enrichment, CellChat, and pseudotime analyses in an independent disease- or injury-relevant dataset.")
    lines.append("- Validate candidate pathways and ligand-receptor interactions experimentally or with orthogonal spatial/protein evidence.")
    return "\n".join(lines) + "\n"


def write_log(path: Path, result_dir: Path, report_path: Path, environment_report: str = "") -> None:
    lines = [
        "# Full Workflow Log",
        "",
        f"- Result directory: `{portable_display_path(result_dir)}`",
        f"- Report: `{portable_display_path(report_path)}`",
        "- Environment setup/check scripts: `scripts/env_setup/*.ps1`",
        "- Seurat core workflow: `scripts/course_adapted/01_seurat_v5_core_pipeline.R`",
        "- Course modules: `scripts/course_adapted/02_marker_enrichment_from_seurat.R`, `03_cellchat_from_seurat.R`, `04_monocle3_from_seurat.R`, `05_singler_cell_annotation.R`",
        "- Report generator: `scripts/write_analysis_report.py`",
        "- Known-bad older pseudotime artifacts removed: old root cell file, duplicate old sessionInfo file, and PDF inspection scratch images.",
        "- Current Monocle3 report figure uses multi-state `marker_support_label` input rather than a single Fibroblasts-only subset.",
        "",
    ]
    if environment_report:
        lines.insert(5, f"- Latest integrated environment check: `{environment_report}`")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", required=True, type=Path)
    parser.add_argument("--metadata-json", type=Path)
    parser.add_argument("--out-md", required=True, type=Path)
    parser.add_argument("--out-log", type=Path)
    parser.add_argument("--environment-report", default="")
    args = parser.parse_args()

    result_dir = args.result_dir.resolve()
    metadata = load_metadata(args.metadata_json)
    report = make_report(result_dir, metadata)
    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    args.out_md.write_text(report, encoding="utf-8")
    if args.out_log:
        write_log(args.out_log, result_dir, args.out_md.resolve(), args.environment_report)
    print(f"Wrote report: {args.out_md}")
    if args.out_log:
        print(f"Wrote log: {args.out_log}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

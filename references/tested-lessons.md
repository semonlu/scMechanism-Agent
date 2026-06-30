# Tested Lessons From Real Runs

Use this reference when a Seurat V5 course-derived run fails, produces suspicious figures, or needs a reproducible full-run example.

## Environment Lessons

- Check `Rscript` before any R workflow. On Windows, `Rscript` may be installed but absent from `PATH`; prefer passing `-RHome` or adding `<R_HOME>\bin`.
- Require R 4.3 or newer, not one fixed R version. Rtools must match the R minor version.
- Require Python 3.10 or newer for Scanpy/reticulate-compatible helper environments.
- Add JAGS to `PATH` before installing packages that link to JAGS.
- Keep `R_LIBS_USER` stable when repairing package libraries; package checks and analysis runs must use the same library path.
- Treat PowerShell profile execution-policy warnings as shell noise unless they stop the command.

## Course-Code Adaptation Lessons

- Replace Seurat v5-defunct `slot=` arguments with `layer=` in `GetAssayData()`.
- Avoid Windows backslash paths inside rendered R strings, especially with non-ASCII directories. Render paths with `/`.
- Do not preserve interactive `choose.dir()` or fixed object names from course scripts.
- Keep `scripts/course_source/` as traceability evidence and run heavy analysis from `scripts/course_adapted/`.

## Pseudotime Lessons

- Do not run Monocle3 on one coarse label just because it is abundant. A trajectory figure meant to show cell states should use a biologically related continuum with multiple states.
- If reference annotation is too coarse, use marker-supported states such as `marker_support_label` for trajectory coloring, while documenting the caveat.
- Import reviewed Seurat UMAP coordinates into Monocle3 when the upstream Seurat layout is the audited visualization.
- Choose a biologically justified root, such as an early stage or early cell state. Document the root rule.
- In recent Monocle3 versions, disable all node labels when making publication-facing figures: `label_leaves = FALSE`, `label_branch_points = FALSE`, `label_roots = FALSE`, and `label_principal_points = FALSE`.
- Keep a separate graph overlay figure for method checking and a clean cell-state figure for reporting.

## CellChat And Enrichment Lessons

- CellChat and enrichment are database-based inference. Reports must describe them as computational hypotheses.
- CellChat may require Seurat v5 layer-compatible expression extraction.
- KEGG REST calls can require network access; save status logs and do not silently omit failures.

## Reporting Lessons

- A full workflow is incomplete until `validate_result_bundle.py` and `write_analysis_report.py` have run.
- Reports must be table-backed. Do not invent p values, pathways, ligands, receptors, or gene names that are absent from output tables.
- Public processed GEO examples should clearly state scope limits, especially when the dataset is adjacent to but not identical to the clinical disease requested.

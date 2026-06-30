# Course Source Scripts

This folder contains course-derived R reference scripts from the Seurat V5 course codebase.

The scripts are not copied verbatim as Chinese nested folders. They have been:

- renamed to English filenames;
- flattened into one script directory;
- given a header that records the original course path;
- lightly adapted so interactive working-directory selection is replaced with `SC_WORK_DIR`/`getwd()` when present.

Use `source_manifest.csv` to map each English filename back to its original course path.

These files are reference material for code generation and review. For runnable project workflows, use `scripts/course_adapted/`.

## Main References

| English file | Use |
|---|---|
| `05_read_10x_standard.R` | Standard 10x MEX import reference. |
| `08_read_10x_h5.R` | 10x H5 import reference. |
| `10_quality_control.R` | QC plots and filtering reference. |
| `11_normalization_decontx_harmony.R` | Normalization, decontamination notes, Harmony, UMAP/tSNE reference. |
| `15_clustering_resolution.R` | Neighbor graph, clustering, resolution sweep reference. |
| `23_marker_detection_methods.R` | Multiple marker detection method reference. |
| `24_go_kegg_enrichment.R` | GO/KEGG enrichment reference. |
| `26_monocle3_pseudotime.R` | Monocle3 trajectory reference. |
| `27_cellchat_analysis.R` | CellChat ligand-receptor analysis reference. |
| `28_copykat_cnv.R` to `31_infercnv_score.R` | CNV inference references. |
| `32_hdwgcna_coexpression.R` | hdWGCNA coexpression reference. |
| `34_cibersort_deconvolution.R`, `35_music_deconvolution.R` | Deconvolution references. |

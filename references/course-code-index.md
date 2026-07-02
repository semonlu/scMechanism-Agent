# Seurat V5 Course Code Index

The current skill does not keep the original Chinese nested filenames as runnable package paths. Course-derived scripts are flattened, renamed to English, lightly adapted, and stored under `scripts/course_source/`. The source-module label for every English file is recorded in `scripts/course_source/source_manifest.csv`; no absolute local course directory is required.

| Stage | English reference script | Adapted use |
|---|---|---|
| Package setup | `01_install_core_r_packages.R`, `02_install_github_r_packages.R` | Dependency inventory only. |
| Single matrix import | `03_read_single_matrix.R` | Seurat object construction pattern. |
| Multiple matrix import | `04_read_multiple_matrices.R` | Merge pattern with sample metadata. |
| Standard 10x import | `05_read_10x_standard.R` | `scripts/course_adapted/01_seurat_v5_core_pipeline.R`. |
| Non-standard 10x import | `06_read_10x_nonstandard.R` | File normalization guidance. |
| Multiple 10x samples | `07_read_multiple_10x_samples.R` | Multi-sample import pattern. |
| 10x H5 import | `08_read_10x_h5.R` | `scripts/course_adapted/01_seurat_v5_core_pipeline.R`. |
| Mixed input merge | `09_merge_mixed_inputs.R` | Metadata harmonization guidance. |
| QC | `10_quality_control.R` | `scripts/course_adapted/01_seurat_v5_core_pipeline.R`. |
| Normalization/Harmony | `11_normalization_decontx_harmony.R` | `scripts/course_adapted/01_seurat_v5_core_pipeline.R`. |
| DoubletFinder | `12_doublet_finder.R` | Optional extension; dynamic column discovery. |
| scDblFinder | `13_scdblfinder.R` | Optional extension. |
| Post-doublet normalization | `14_post_doublet_normalization.R` | Reprocessing after doublet removal. |
| Clustering | `15_clustering_resolution.R` | `scripts/course_adapted/01_seurat_v5_core_pipeline.R` and resolution guidance. |
| Manual annotation | `16_manual_cell_annotation.R` | Marker-based annotation rules. |
| SingleR annotation | `17_singler_annotation.R` | Automated annotation support. |
| SCINA annotation | `18_scina_annotation.R` | Automated annotation support. |
| LLM annotation | `19_llm_annotation_kimi.R`, `20_llm_annotation_deepseek.R` | Privacy-bounded marker summary only; do not upload private matrices. |
| TransferData/scPred | `21_transferdata_annotation.R`, `22_scpred_annotation.R` | Reference mapping when compatible. |
| Marker detection | `23_marker_detection_methods.R` | `scripts/course_adapted/02_marker_enrichment_from_seurat.R`. |
| GO/KEGG | `24_go_kegg_enrichment.R` | `scripts/course_adapted/02_marker_enrichment_from_seurat.R`. |
| Pseudotime | `25_monocle2_pseudotime.R`, `26_monocle3_pseudotime.R` | `scripts/course_adapted/04_monocle3_from_seurat.R`. |
| CellChat | `27_cellchat_analysis.R` | `scripts/course_adapted/03_cellchat_from_seurat.R`. |
| copykat | `28_copykat_cnv.R` | Tumor CNV optional extension. |
| inferCNV | `29_infercnv_basic.R`, `30_infercnv_with_normal_reference.R`, `31_infercnv_score.R` | Tumor CNV optional extension. |
| hdWGCNA | `32_hdwgcna_coexpression.R` | Co-expression optional extension. |
| CIBERSORT/MuSiC | `33_cibersort_immune_helper.R`, `34_cibersort_deconvolution.R`, `35_music_deconvolution.R` | Bulk deconvolution optional extension. |

Use `scripts/build_codebase_summary.py` to regenerate a local summary from `scripts/course_source/` after the source references change.

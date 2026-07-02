# Course-derived reference script
# English filename: 20_llm_annotation_deepseek.R
# Original course path: see source_manifest.csv
# Role: LLM annotation reference
# Adaptations applied for the public skill package:
# - Removed hard-coded API keys and private helper-package dependency.
# - Exports de-identified marker summaries only.
# - Use SingleR/manual annotation for default runs; use an approved LLM client only when privacy rules allow it.

WORK_DIR <- Sys.getenv("SC_WORK_DIR", unset = getwd())
setwd(WORK_DIR)

required_packages <- c("Seurat", "qs", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

input_qs <- Sys.getenv("SC_SEURAT_QS", unset = "scRNA_clustered.qs")
if (!file.exists(input_qs)) {
  stop("Seurat QS object not found: ", input_qs)
}

output_dir <- Sys.getenv("SC_OUTPUT_DIR", unset = "llm_annotation_safe_export")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

scRNA <- qs::qread(input_qs, nthreads = parallel::detectCores())
markers <- Seurat::FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

marker_summary <- markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)

utils::write.csv(marker_summary, file.path(output_dir, "marker_summary_for_llm_review.csv"), row.names = FALSE)
writeLines(c(
  "LLM annotation safety note:",
  "- Do not upload raw count matrices, cell barcodes, patient identifiers, or protected clinical metadata.",
  "- Review only de-identified cluster-level marker summaries.",
  "- Record the model, prompt, date, and human review decision before accepting labels."
), file.path(output_dir, "llm_annotation_safety_note.txt"))

message("Wrote safe marker summary to: ", normalizePath(output_dir, winslash = "/"))

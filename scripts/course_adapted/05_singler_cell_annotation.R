#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(SummarizedExperiment)
  library(ggplot2)
})

input_rds <- "{{INPUT_RDS}}"
output_dir <- "{{OUTPUT_DIR}}"
organism <- "{{ORGANISM}}" # human or mouse
reference_name <- "{{REFERENCE_NAME}}" # hpca, blueprint_encode, mouse_rnaseq
cluster_col <- "{{CLUSTER_COL}}"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "objects"), showWarnings = FALSE)
dir.create(file.path(output_dir, "logs"), showWarnings = FALSE)

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("INPUT_RDS must contain a Seurat object")
}

if (!requireNamespace("celldex", quietly = TRUE)) {
  stop("Package celldex is required for built-in SingleR references")
}

ref <- switch(
  reference_name,
  hpca = celldex::HumanPrimaryCellAtlasData(),
  blueprint_encode = celldex::BlueprintEncodeData(),
  mouse_rnaseq = celldex::MouseRNAseqData(),
  {
    if (tolower(organism) == "mouse") celldex::MouseRNAseqData() else celldex::HumanPrimaryCellAtlasData()
  }
)

expr <- GetAssayData(obj, assay = DefaultAssay(obj), layer = "data")
if (ncol(expr) == 0) {
  expr <- GetAssayData(obj, assay = DefaultAssay(obj), layer = "counts")
}

clusters <- NULL
if (nzchar(cluster_col) && cluster_col %in% colnames(obj@meta.data)) {
  clusters <- obj@meta.data[[cluster_col]]
} else if ("seurat_clusters" %in% colnames(obj@meta.data)) {
  clusters <- obj$seurat_clusters
}

pred <- SingleR(test = expr, ref = ref, labels = ref$label.main, clusters = clusters)
pred_df <- as.data.frame(pred)
pred_df$cluster <- rownames(pred_df)
write.csv(pred_df, file.path(output_dir, "tables", "singleR_cluster_labels.csv"), row.names = FALSE)

if (!is.null(clusters)) {
  pruned <- if (!is.null(pred$pruned.labels)) pred$pruned.labels else pred$labels
  final_cluster_label <- ifelse(is.na(pruned) | !nzchar(pruned), paste0("Ambiguous_", pred$labels), pruned)
  cluster_to_label <- setNames(final_cluster_label, rownames(pred))
  obj$singleR_label <- unname(cluster_to_label[as.character(clusters)])
} else {
  pruned <- if (!is.null(pred$pruned.labels)) pred$pruned.labels else pred$labels
  obj$singleR_label <- ifelse(is.na(pruned) | !nzchar(pruned), paste0("Ambiguous_", pred$labels), pruned)
}

write.csv(
  data.frame(cell = colnames(obj), singleR_label = obj$singleR_label),
  file.path(output_dir, "tables", "singleR_cell_labels.csv"),
  row.names = FALSE
)

marker_rows <- NULL
if (!is.null(clusters)) {
  Idents(obj) <- as.factor(clusters)
  marker_rows <- tryCatch(
    FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25),
    error = function(e) NULL
  )
}

if (!is.null(clusters)) {
  cluster_ids <- rownames(pred)
  top_markers <- vapply(cluster_ids, function(cl) {
    if (is.null(marker_rows)) return("")
    hits <- marker_rows[as.character(marker_rows$cluster) == as.character(cl), , drop = FALSE]
    if (!nrow(hits)) return("")
    hits <- hits[order(hits$avg_log2FC, decreasing = TRUE), , drop = FALSE]
    paste(head(hits$gene, 10), collapse = ";")
  }, character(1))
  cell_counts <- as.integer(table(as.factor(clusters))[as.character(cluster_ids)])
  pruned <- if (!is.null(pred$pruned.labels)) pred$pruned.labels else pred$labels
  delta_next <- if ("delta.next" %in% colnames(pred_df)) pred_df$delta.next else NA_real_
  confidence <- ifelse(is.na(pruned) | !nzchar(pruned), "low", ifelse(is.na(delta_next) | delta_next < 0.05, "medium", "high"))
  final_label <- ifelse(confidence == "low", paste0("Ambiguous_", pred$labels), pruned)
  evidence <- data.frame(
    cluster = cluster_ids,
    final_label = final_label,
    coarse_label = pred$labels,
    singleR_label = pred$labels,
    singleR_pruned_label = pruned,
    singleR_delta_next = delta_next,
    top_markers = top_markers,
    canonical_marker_support = "manual_review_required",
    conflicting_markers = "",
    cell_count = cell_counts,
    confidence = confidence,
    review_note = ifelse(confidence == "low", "SingleR pruning failed or label is weak; keep ambiguous until marker review.", "Review marker plots before downstream modules."),
    stringsAsFactors = FALSE
  )
  write.table(
    evidence,
    file.path(output_dir, "tables", "annotation_evidence.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

if ("umap" %in% names(obj@reductions)) {
  pdf(file.path(output_dir, "figures", "umap_singleR_labels.pdf"), width = 8, height = 6)
  print(DimPlot(obj, reduction = "umap", group.by = "singleR_label", label = TRUE, repel = TRUE) + ggtitle("SingleR labels"))
  dev.off()
}

saveRDS(obj, file.path(output_dir, "objects", "seurat_singleR_annotated.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "logs", "sessionInfo_singler_annotation.txt"))

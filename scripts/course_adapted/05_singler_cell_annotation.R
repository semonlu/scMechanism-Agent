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
  cluster_to_label <- setNames(pred$labels, rownames(pred))
  obj$singleR_label <- unname(cluster_to_label[as.character(clusters)])
} else {
  obj$singleR_label <- pred$labels
}

write.csv(
  data.frame(cell = colnames(obj), singleR_label = obj$singleR_label),
  file.path(output_dir, "tables", "singleR_cell_labels.csv"),
  row.names = FALSE
)

if ("umap" %in% names(obj@reductions)) {
  pdf(file.path(output_dir, "figures", "umap_singleR_labels.pdf"), width = 8, height = 6)
  print(DimPlot(obj, reduction = "umap", group.by = "singleR_label", label = TRUE, repel = TRUE) + ggtitle("SingleR labels"))
  dev.off()
}

saveRDS(obj, file.path(output_dir, "objects", "seurat_singleR_annotated.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "logs", "sessionInfo_singler_annotation.txt"))

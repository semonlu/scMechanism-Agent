#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

input_path <- "{{INPUT_PATH}}"
metadata_path <- "{{METADATA_PATH}}"
output_dir <- "{{OUTPUT_DIR}}"
organism <- "{{ORGANISM}}" # human or mouse
sample_id <- "{{SAMPLE_ID}}"
batch_col <- "{{BATCH_COL}}"
condition_col <- "{{CONDITION_COL}}"
input_type <- "{{INPUT_TYPE}}" # 10x_mtx, 10x_h5, rds, csv
resolution_values <- c(0.1, 0.3, 0.5, 0.8)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "objects"), showWarnings = FALSE)
dir.create(file.path(output_dir, "logs"), showWarnings = FALSE)

mt_pattern <- if (tolower(organism) == "mouse") "^mt-" else "^MT-"
ribo_pattern <- if (tolower(organism) == "mouse") "^Rp[sl]" else "^RP[SL]"
hb_pattern <- if (tolower(organism) == "mouse") "^Hb[abq]" else "^HB[ABDEGQMZ]"

resolution_col <- function(res) {
  paste0("resolution_", gsub("\\.", "_", as.character(res)))
}

resolution_note <- function(total_cells, cluster_count) {
  if (total_cells >= 5000 && total_cells <= 20000) {
    if (cluster_count < 8) return("possible_under_clustering_for_major_cell_types")
    if (cluster_count > 25) return("possible_over_clustering_or_state_splitting")
    if (cluster_count >= 10 && cluster_count <= 20) return("typical_major_lineage_range_for_10x_5k_20k_cells")
  }
  if (cluster_count <= 1) return("single_cluster_or_no_clear_structure")
  if (cluster_count > 25) return("many_clusters_review_marker_specificity")
  "review_marker_clarity_and_biology"
}

select_resolution <- function(summary_df) {
  preference <- c(0.3, 0.5, 0.1, 0.8)
  total_cells <- summary_df$total_cells[1]
  if (total_cells >= 5000 && total_cells <= 20000) {
    candidates <- summary_df[summary_df$cluster_count >= 10 & summary_df$cluster_count <= 20, , drop = FALSE]
    if (nrow(candidates) > 0) {
      for (res in preference) if (res %in% candidates$resolution) return(res)
    }
  }
  candidates <- summary_df[summary_df$cluster_count > 1 & summary_df$cluster_count <= 25, , drop = FALSE]
  if (nrow(candidates) > 0) {
    for (res in preference) if (res %in% candidates$resolution) return(res)
  }
  preference[preference %in% summary_df$resolution][1]
}

write_resolution_sweep <- function(obj, values) {
  summary_rows <- list()
  for (res in values) {
    obj <- FindClusters(obj, resolution = res, verbose = FALSE)
    col <- resolution_col(res)
    assignments <- as.character(Idents(obj))
    obj[[col]] <- assignments
    sizes <- table(assignments)
    summary_rows[[as.character(res)]] <- data.frame(
      resolution = res,
      cluster_count = length(sizes),
      total_cells = ncol(obj),
      min_cluster_cells = min(as.integer(sizes)),
      median_cluster_cells = median(as.integer(sizes)),
      max_cluster_cells = max(as.integer(sizes)),
      review_note = resolution_note(ncol(obj), length(sizes)),
      stringsAsFactors = FALSE
    )
  }
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  write.table(summary_df, file.path(output_dir, "tables", "resolution_sweep.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  list(obj = obj, summary = summary_df, selected_resolution = select_resolution(summary_df))
}

message("Reading input: ", input_path)
if (input_type == "10x_mtx") {
  counts <- Read10X(data.dir = input_path)
  obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200, project = sample_id)
} else if (input_type == "10x_h5") {
  counts <- Read10X_h5(filename = input_path)
  obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200, project = sample_id)
} else if (input_type == "rds") {
  obj <- readRDS(input_path)
} else if (input_type == "csv") {
  counts <- read.csv(input_path, row.names = 1, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  obj <- CreateSeuratObject(counts = as.matrix(counts), min.cells = 3, min.features = 200, project = sample_id)
} else {
  stop("Unsupported input_type: ", input_type)
}

if (nzchar(metadata_path) && file.exists(metadata_path)) {
  meta <- read.csv(metadata_path, row.names = 1, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  common <- intersect(colnames(obj), rownames(meta))
  obj <- obj[, common]
  obj <- AddMetaData(obj, meta[common, , drop = FALSE])
}

if (!"sample_id" %in% colnames(obj@meta.data)) {
  obj$sample_id <- sample_id
}

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = ribo_pattern)
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = hb_pattern)

qc_before <- data.frame(cells = ncol(obj), genes = nrow(obj), stage = "before_filter")
write.table(qc_before, file.path(output_dir, "tables", "qc_cell_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
qc_cols <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb")
write.table(obj@meta.data[, qc_cols, drop = FALSE], file.path(output_dir, "tables", "cell_qc_metrics_pre_filter.tsv"), sep = "\t", quote = FALSE, col.names = NA)

pdf(file.path(output_dir, "figures", "qc_pre_filter.pdf"), width = 14, height = 6)
print(VlnPlot(obj, features = qc_cols, ncol = 5))
dev.off()

pdf(file.path(output_dir, "figures", "qc_scatter.pdf"), width = 15, height = 5)
print(FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.hb") + FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()

# Replace these thresholds after reviewing current data distributions.
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)

qc_after <- data.frame(cells = ncol(obj), genes = nrow(obj), stage = "after_filter")
write.table(qc_after, file.path(output_dir, "tables", "qc_cell_counts_after.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(obj@meta.data[, qc_cols, drop = FALSE], file.path(output_dir, "tables", "cell_qc_metrics_post_filter.tsv"), sep = "\t", quote = FALSE, col.names = NA)

pdf(file.path(output_dir, "figures", "qc_post_filter.pdf"), width = 14, height = 6)
print(VlnPlot(obj, features = qc_cols, ncol = 5))
dev.off()

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj))

pdf(file.path(output_dir, "figures", "elbow_plot.pdf"), width = 6, height = 5)
print(ElbowPlot(obj, ndims = 40))
dev.off()

dims_use <- 1:30
if (batch_col %in% colnames(obj@meta.data) && requireNamespace("harmony", quietly = TRUE)) {
  obj <- harmony::RunHarmony(obj, group.by.vars = batch_col)
  reduction_use <- "harmony"
} else {
  reduction_use <- "pca"
}

obj <- FindNeighbors(obj, reduction = reduction_use, dims = dims_use)
sweep <- write_resolution_sweep(obj, resolution_values)
obj <- sweep$obj
selected_resolution <- sweep$selected_resolution
selected_col <- resolution_col(selected_resolution)
obj$seurat_clusters <- factor(obj@meta.data[[selected_col]])
Idents(obj) <- "seurat_clusters"
write.table(
  data.frame(parameter = c("resolution_values", "selected_resolution"), value = c(paste(resolution_values, collapse = ","), selected_resolution)),
  file.path(output_dir, "tables", "clustering_parameters.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
obj <- RunUMAP(obj, reduction = reduction_use, dims = dims_use)

pdf(file.path(output_dir, "figures", "umap_clusters.pdf"), width = 7, height = 6)
print(DimPlot(obj, reduction = "umap", label = TRUE) + ggtitle("Clusters"))
dev.off()

if (condition_col %in% colnames(obj@meta.data)) {
  pdf(file.path(output_dir, "figures", "umap_condition.pdf"), width = 7, height = 6)
  print(DimPlot(obj, reduction = "umap", group.by = condition_col) + ggtitle(condition_col))
  dev.off()
}

markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(output_dir, "tables", "cluster_markers.csv"), row.names = FALSE)

saveRDS(obj, file.path(output_dir, "objects", "processed_seurat.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "logs", "sessionInfo.txt"))

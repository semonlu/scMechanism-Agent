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

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "objects"), showWarnings = FALSE)
dir.create(file.path(output_dir, "logs"), showWarnings = FALSE)

mt_pattern <- if (tolower(organism) == "mouse") "^mt-" else "^MT-"
ribo_pattern <- if (tolower(organism) == "mouse") "^Rp[sl]" else "^RP[SL]"

message("Reading input: ", input_path)
if (input_type == "10x_mtx") {
  counts <- Read10X(data.dir = input_path)
  obj <- CreateSeuratObject(counts = counts, min.cells = 1, min.features = 1, project = sample_id)
} else if (input_type == "10x_h5") {
  counts <- Read10X_h5(filename = input_path)
  obj <- CreateSeuratObject(counts = counts, min.cells = 1, min.features = 1, project = sample_id)
} else if (input_type == "rds") {
  obj <- readRDS(input_path)
} else if (input_type == "csv") {
  counts <- read.csv(input_path, row.names = 1, check.names = FALSE)
  obj <- CreateSeuratObject(counts = as.matrix(counts), min.cells = 1, min.features = 1, project = sample_id)
} else {
  stop("Unsupported input_type: ", input_type)
}

if (nzchar(metadata_path) && file.exists(metadata_path)) {
  meta <- read.csv(metadata_path, row.names = 1, check.names = FALSE)
  common <- intersect(colnames(obj), rownames(meta))
  obj <- obj[, common]
  obj <- AddMetaData(obj, meta[common, , drop = FALSE])
}

if (!"sample_id" %in% colnames(obj@meta.data)) {
  obj$sample_id <- sample_id
}

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = ribo_pattern)

qc_before <- data.frame(cells = ncol(obj), genes = nrow(obj), stage = "before_filter")
write.table(qc_before, file.path(output_dir, "tables", "qc_cell_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(output_dir, "figures", "qc_pre_filter.pdf"), width = 12, height = 6)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4))
dev.off()

pdf(file.path(output_dir, "figures", "qc_scatter.pdf"), width = 10, height = 5)
print(FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()

# Adaptive defaults keep small validation datasets runnable while preserving
# conservative upper mitochondrial filtering for ordinary datasets.
min_features <- if (nrow(obj) < 200) max(1, floor(nrow(obj) * 0.1)) else 200
max_features <- max(6000, nrow(obj) + 1)
obj <- subset(obj, subset = nFeature_RNA >= min_features & nFeature_RNA < max_features & percent.mt < 25)
if (ncol(obj) < 3 || nrow(obj) < 3) {
  stop("Too few cells or genes remain after QC filtering; inspect input or lower QC thresholds.")
}

qc_after <- data.frame(cells = ncol(obj), genes = nrow(obj), stage = "after_filter")
write.table(qc_after, file.path(output_dir, "tables", "qc_cell_counts_after.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(output_dir, "figures", "qc_post_filter.pdf"), width = 12, height = 6)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4))
dev.off()

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = rownames(obj))
npcs <- min(30, ncol(obj) - 1, length(VariableFeatures(obj)) - 1)
if (npcs < 2) {
  stop("Too few cells or variable genes for PCA.")
}
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs)

pdf(file.path(output_dir, "figures", "elbow_plot.pdf"), width = 6, height = 5)
print(ElbowPlot(obj, ndims = min(40, npcs)))
dev.off()

dims_use <- 1:npcs
has_batch_levels <- batch_col %in% colnames(obj@meta.data) && length(unique(obj@meta.data[[batch_col]])) >= 2
if (has_batch_levels && requireNamespace("harmony", quietly = TRUE)) {
  obj <- harmony::RunHarmony(obj, group.by.vars = batch_col)
  reduction_use <- "harmony"
} else {
  reduction_use <- "pca"
}

obj <- FindNeighbors(obj, reduction = reduction_use, dims = dims_use)
obj <- FindClusters(obj, resolution = 0.4)
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

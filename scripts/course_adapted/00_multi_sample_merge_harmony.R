#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
})

sample_table_path <- "{{SAMPLE_TABLE}}"
output_dir <- "{{OUTPUT_DIR}}"
organism <- "{{ORGANISM}}" # human or mouse
batch_col <- "{{BATCH_COL}}"
condition_col <- "{{CONDITION_COL}}"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "objects"), showWarnings = FALSE)
dir.create(file.path(output_dir, "logs"), showWarnings = FALSE)

samples <- read.csv(sample_table_path, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("sample_id", "input_path", "input_type")
missing_cols <- setdiff(required_cols, colnames(samples))
if (length(missing_cols) > 0) {
  stop("sample table missing required columns: ", paste(missing_cols, collapse = ", "))
}

read_one <- function(sample_id, input_path, input_type) {
  message("Reading sample ", sample_id, ": ", input_path)
  if (input_type == "10x_mtx") {
    counts <- Read10X(data.dir = input_path)
  } else if (input_type == "10x_h5") {
    counts <- Read10X_h5(filename = input_path)
  } else if (input_type == "csv") {
    counts <- read.csv(input_path, row.names = 1, check.names = FALSE)
    counts <- as.matrix(counts)
  } else {
    stop("Unsupported input_type for merge: ", input_type)
  }
  obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200, project = sample_id)
  obj$sample_id <- sample_id
  obj
}

objects <- vector("list", nrow(samples))
for (i in seq_len(nrow(samples))) {
  objects[[i]] <- read_one(samples$sample_id[i], samples$input_path[i], samples$input_type[i])
  sample_meta <- samples[i, setdiff(colnames(samples), c("input_path", "input_type")), drop = FALSE]
  rownames(sample_meta) <- samples$sample_id[i]
  for (col in colnames(sample_meta)) {
    objects[[i]][[col]] <- sample_meta[[col]][1]
  }
}

if (length(objects) == 1) {
  obj <- objects[[1]]
} else {
  obj <- merge(objects[[1]], y = objects[-1], add.cell.ids = samples$sample_id, project = "merged_scRNA")
}

mt_pattern <- if (tolower(organism) == "mouse") "^mt-" else "^MT-"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)

qc <- obj@meta.data |>
  group_by(sample_id) |>
  summarise(cells = dplyr::n(), median_features = median(nFeature_RNA), median_counts = median(nCount_RNA), .groups = "drop")
write.csv(qc, file.path(output_dir, "tables", "multi_sample_qc_by_sample.csv"), row.names = FALSE)

pdf(file.path(output_dir, "figures", "multi_sample_qc_pre_filter.pdf"), width = 12, height = 6)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample_id", ncol = 3))
dev.off()

obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj))

if (nzchar(batch_col) && batch_col %in% colnames(obj@meta.data) && requireNamespace("harmony", quietly = TRUE)) {
  obj <- harmony::RunHarmony(obj, group.by.vars = batch_col)
  reduction_use <- "harmony"
} else {
  reduction_use <- "pca"
}

dims_use <- 1:30
obj <- FindNeighbors(obj, reduction = reduction_use, dims = dims_use)
obj <- FindClusters(obj, resolution = 0.4)
obj <- RunUMAP(obj, reduction = reduction_use, dims = dims_use)

pdf(file.path(output_dir, "figures", "umap_by_sample.pdf"), width = 7, height = 6)
print(DimPlot(obj, reduction = "umap", group.by = "sample_id"))
dev.off()

if (nzchar(condition_col) && condition_col %in% colnames(obj@meta.data)) {
  pdf(file.path(output_dir, "figures", "umap_by_condition.pdf"), width = 7, height = 6)
  print(DimPlot(obj, reduction = "umap", group.by = condition_col))
  dev.off()
}

markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(output_dir, "tables", "multi_sample_cluster_markers.csv"), row.names = FALSE)
saveRDS(obj, file.path(output_dir, "objects", "multi_sample_harmony_seurat.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "logs", "sessionInfo_multi_sample_harmony.txt"))

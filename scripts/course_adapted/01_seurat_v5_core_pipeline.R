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
input_type <- "{{INPUT_TYPE}}" # 10x_mtx, 10x_nonstandard, 10x_h5, rds, csv

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "objects"), showWarnings = FALSE)
dir.create(file.path(output_dir, "logs"), showWarnings = FALSE)

mt_pattern <- if (tolower(organism) == "mouse") "^mt-" else "^MT-"
ribo_pattern <- if (tolower(organism) == "mouse") "^Rp[sl]" else "^RP[SL]"

read_maybe_gz_table <- function(path) {
  read.delim(path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
}

read_maybe_gz_mtx <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    Matrix::readMM(con)
  } else {
    Matrix::readMM(path)
  }
}

first_existing <- function(paths) {
  hits <- paths[file.exists(paths)]
  if (length(hits) == 0) {
    return(NA_character_)
  }
  hits[[1]]
}

read_nonstandard_10x <- function(path, fallback_sample_id) {
  candidate_dirs <- unique(c(path, list.dirs(path, recursive = FALSE, full.names = TRUE)))
  count_list <- list()
  sample_by_cell <- character()

  for (sample_dir in candidate_dirs) {
    matrix_file <- first_existing(file.path(sample_dir, c("count_matrix_sparse.mtx", "count_matrix_sparse.mtx.gz", "matrix.mtx", "matrix.mtx.gz")))
    barcode_file <- first_existing(file.path(sample_dir, c("count_matrix_barcodes.tsv", "count_matrix_barcodes.tsv.gz", "barcodes.tsv", "barcodes.tsv.gz")))
    gene_file <- first_existing(file.path(sample_dir, c("count_matrix_genes.tsv", "count_matrix_genes.tsv.gz", "features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz")))

    if (is.na(matrix_file) || is.na(barcode_file) || is.na(gene_file)) {
      next
    }

    this_sample <- basename(normalizePath(sample_dir, winslash = "/", mustWork = FALSE))
    if (!nzchar(this_sample) || this_sample == "." || identical(normalizePath(sample_dir, winslash = "/", mustWork = FALSE), normalizePath(path, winslash = "/", mustWork = FALSE))) {
      this_sample <- fallback_sample_id
    }

    counts <- read_maybe_gz_mtx(matrix_file)
    barcodes <- read_maybe_gz_table(barcode_file)
    genes <- read_maybe_gz_table(gene_file)
    gene_names <- if (ncol(genes) >= 2) genes[[2]] else genes[[1]]

    if (ncol(counts) != nrow(barcodes) || nrow(counts) != length(gene_names)) {
      stop("Non-standard 10x dimensions do not match in ", sample_dir)
    }

    rownames(counts) <- make.unique(as.character(gene_names))
    colnames(counts) <- paste(this_sample, as.character(barcodes[[1]]), sep = "_")
    count_list[[this_sample]] <- counts
    sample_by_cell[colnames(counts)] <- this_sample
  }

  if (length(count_list) == 0) {
    stop("No non-standard 10x matrix set found. Expected count_matrix_sparse.mtx, count_matrix_barcodes.tsv, and count_matrix_genes.tsv under INPUT_PATH or its immediate sample subdirectories.")
  }

  list(counts = do.call(cbind, count_list), sample_by_cell = sample_by_cell)
}

message("Reading input: ", input_path)
if (input_type == "10x_mtx") {
  counts <- Read10X(data.dir = input_path)
  obj <- CreateSeuratObject(counts = counts, min.cells = 1, min.features = 1, project = sample_id)
} else if (input_type == "10x_nonstandard") {
  nonstandard <- read_nonstandard_10x(input_path, sample_id)
  obj <- CreateSeuratObject(counts = nonstandard$counts, min.cells = 1, min.features = 1, project = sample_id)
  obj$sample_id <- unname(nonstandard$sample_by_cell[colnames(obj)])
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

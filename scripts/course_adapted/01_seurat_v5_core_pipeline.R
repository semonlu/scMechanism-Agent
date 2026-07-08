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
resolution_values <- c(0.1, 0.3, 0.5, 0.8)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "objects"), showWarnings = FALSE)
dir.create(file.path(output_dir, "logs"), showWarnings = FALSE)

mt_pattern <- if (tolower(organism) == "mouse") "^mt-" else "^MT-"
ribo_pattern <- if (tolower(organism) == "mouse") "^Rp[sl]" else "^RP[SL]"
hb_pattern <- if (tolower(organism) == "mouse") "^Hb[abq]" else "^HB[ABDEGQMZ]"

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
      for (res in preference) {
        if (res %in% candidates$resolution) return(res)
      }
    }
  }
  candidates <- summary_df[summary_df$cluster_count > 1 & summary_df$cluster_count <= 25, , drop = FALSE]
  if (nrow(candidates) > 0) {
    for (res in preference) {
      if (res %in% candidates$resolution) return(res)
    }
  }
  preference[preference %in% summary_df$resolution][1]
}

write_resolution_sweep <- function(obj, values) {
  summary_rows <- list()
  size_rows <- list()
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
    size_rows[[as.character(res)]] <- data.frame(
      resolution = res,
      cluster = names(sizes),
      cells = as.integer(sizes),
      stringsAsFactors = FALSE
    )
  }
  summary_df <- do.call(rbind, summary_rows)
  size_df <- do.call(rbind, size_rows)
  rownames(summary_df) <- NULL
  rownames(size_df) <- NULL
  write.table(summary_df, file.path(output_dir, "tables", "resolution_sweep.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(size_df, file.path(output_dir, "tables", "cluster_sizes_by_resolution.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  list(obj = obj, summary = summary_df, selected_resolution = select_resolution(summary_df))
}

write_cluster_marker_audit <- function(obj, markers) {
  clusters <- levels(Idents(obj))
  generic_patterns <- c("^MT-", "^mt-", "^RP[SL]", "^Rp[sl]", "^HB[ABDEGMQZ]", "^HSP", "^Hsp")
  cell_cycle_genes <- c("MKI67", "TOP2A", "UBE2C", "CCNB1", "CCNB2", "CENPF", "PCNA", "MCM2", "MCM5")
  stress_genes <- c("FOS", "JUN", "JUNB", "DUSP1", "HSPA1A", "HSPA1B", "HSP90AA1", "IER2", "ATF3")
  order_col <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else if ("avg_logFC" %in% colnames(markers)) "avg_logFC" else ""
  audit_rows <- lapply(clusters, function(cl) {
    hits <- markers[as.character(markers$cluster) == as.character(cl), , drop = FALSE]
    if (nrow(hits) > 0 && nzchar(order_col)) {
      hits <- hits[order(hits[[order_col]], decreasing = TRUE), , drop = FALSE]
    }
    top_genes <- head(as.character(hits$gene), 10)
    generic_count <- sum(vapply(top_genes, function(g) any(grepl(paste(generic_patterns, collapse = "|"), g)), logical(1)))
    cell_cycle_count <- sum(toupper(top_genes) %in% cell_cycle_genes)
    stress_count <- sum(toupper(top_genes) %in% stress_genes)
    cell_count <- sum(as.character(Idents(obj)) == as.character(cl))
    pct_total <- round(100 * cell_count / ncol(obj), 3)
    notes <- character(0)
    if (nrow(hits) < 100) notes <- c(notes, "few_positive_markers_consider_unknown_or_merge_if_unclear")
    if (nrow(hits) > 300) notes <- c(notes, "many_positive_markers_possible_functional_state_review")
    if (generic_count >= 5) notes <- c(notes, "top_markers_generic_ribo_mt_hb_or_heat_shock_review_qc")
    if (cell_cycle_count >= 2) notes <- c(notes, "cell_cycle_signature_review_s_g2m_scores")
    if (stress_count >= 2) notes <- c(notes, "stress_signature_review_dissociation_or_low_quality")
    if (pct_total < 1) notes <- c(notes, "rare_cluster_under_1_percent_requires_recurrence_before_novelty_claim")
    if (!length(notes)) notes <- "marker_review_required"
    data.frame(
      cluster = cl,
      cell_count = cell_count,
      pct_total = pct_total,
      marker_count_vs_all = nrow(hits),
      top_markers = paste(top_genes, collapse = ";"),
      generic_top_marker_count = generic_count,
      cell_cycle_top_marker_count = cell_cycle_count,
      stress_top_marker_count = stress_count,
      review_note = paste(unique(notes), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  audit <- do.call(rbind, audit_rows)
  write.table(audit, file.path(output_dir, "tables", "cluster_marker_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

write_workflow_step_audit <- function(output_dir, input_path, input_type, metadata_path, reduction_use, has_batch_levels) {
  batch_note <- if (has_batch_levels) {
    paste0("batch metadata detected; reduction used: ", reduction_use)
  } else {
    "batch metadata with at least two levels not detected"
  }
  audit <- data.frame(
    step = as.character(1:10),
    workflow_step = c("data_acquisition", "data_import", "quality_control", "normalization_scaling_hvg", "confounder_review", "dimensionality_clustering", "cell_annotation", "marker_detection", "pseudotime", "cell_communication"),
    required_before_interpretation = c(rep("yes", 8), "no_optional_gated", "no_optional_gated"),
    status = c("external_input_provided", "complete", "complete", "complete", "partial_review", "complete", "pending_annotation_module", "complete", "gated_optional_not_run", "gated_optional_not_run"),
    evidence_files = c(
      input_path,
      paste(c(input_type, if (nzchar(metadata_path)) metadata_path else "metadata_not_provided"), collapse = ";"),
      "qc_cell_counts.tsv;cell_qc_metrics_pre_filter.tsv;cell_qc_metrics_post_filter.tsv;qc_pre_filter.pdf;qc_scatter.pdf;qc_post_filter.pdf",
      "processed_seurat.rds;sessionInfo.txt",
      batch_note,
      "elbow_plot.pdf;umap_clusters.pdf;resolution_sweep.tsv;cluster_marker_audit.tsv",
      "run scripts/course_adapted/05_singler_cell_annotation.R or equivalent manual annotation next",
      "cluster_markers.csv;cluster_marker_audit.tsv",
      "requires downstream_proposal.md and approved lineage/root",
      "requires downstream_proposal.md and approved sender/receiver groups"
    ),
    review_note = c(
      "Verify GEO/SRA/local provenance before interpretation.",
      "Matrix/object was read into Seurat; metadata alignment should be checked in logs/tables.",
      "QC was performed, but thresholds still require reviewer inspection.",
      "LogNormalize, variable features, scaling, and PCA were run; raw counts should remain in the Seurat object.",
      "Doublet, cell-cycle, and ambient RNA review are not completed by this core script; document or run dedicated modules before strong interpretation.",
      "PCA, neighbors, resolution sweep, selected clustering, and UMAP were generated.",
      "Cell annotation is intentionally not completed in the core pipeline.",
      "Cluster marker detection was run and marker specificity audit was exported.",
      "Pseudotime is optional and must not run before annotation and root approval.",
      "Cell communication is optional and must not run before annotation and scope approval."
    ),
    stringsAsFactors = FALSE
  )
  write.table(audit, file.path(output_dir, "tables", "workflow_step_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
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
  counts <- read.csv(input_path, row.names = 1, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  obj <- CreateSeuratObject(counts = as.matrix(counts), min.cells = 1, min.features = 1, project = sample_id)
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
write.table(obj@meta.data[, qc_cols, drop = FALSE], file.path(output_dir, "tables", "cell_qc_metrics_post_filter.tsv"), sep = "\t", quote = FALSE, col.names = NA)

pdf(file.path(output_dir, "figures", "qc_post_filter.pdf"), width = 14, height = 6)
print(VlnPlot(obj, features = qc_cols, ncol = 5))
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
sweep <- write_resolution_sweep(obj, resolution_values)
obj <- sweep$obj
selected_resolution <- sweep$selected_resolution
selected_col <- resolution_col(selected_resolution)
obj$seurat_clusters <- factor(obj@meta.data[[selected_col]])
Idents(obj) <- "seurat_clusters"
write.table(
  data.frame(
    parameter = c("resolution_values", "selected_resolution", "reduction", "dims"),
    value = c(paste(resolution_values, collapse = ","), selected_resolution, reduction_use, paste(range(dims_use), collapse = ":")),
    stringsAsFactors = FALSE
  ),
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
write_cluster_marker_audit(obj, markers)
write_workflow_step_audit(output_dir, input_path, input_type, metadata_path, reduction_use, has_batch_levels)

saveRDS(obj, file.path(output_dir, "objects", "processed_seurat.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "logs", "sessionInfo.txt"))

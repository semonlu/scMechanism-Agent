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
resolution_values <- c(0.1, 0.3, 0.5, 0.8)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "objects"), showWarnings = FALSE)
dir.create(file.path(output_dir, "logs"), showWarnings = FALSE)

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

write_workflow_step_audit <- function(output_dir, sample_table_path, reduction_use, has_batch_levels) {
  batch_note <- if (has_batch_levels) {
    paste0("batch metadata detected; reduction used: ", reduction_use)
  } else {
    "batch metadata with at least two levels not detected"
  }
  audit <- data.frame(
    step = as.character(1:10),
    workflow_step = c(
      "data_acquisition",
      "data_import",
      "quality_control",
      "normalization_scaling_hvg",
      "confounder_review",
      "dimensionality_clustering",
      "cell_annotation",
      "marker_detection",
      "pseudotime",
      "cell_communication"
    ),
    required_before_interpretation = c(rep("yes", 8), "no_optional_gated", "no_optional_gated"),
    status = c(
      "external_sample_table_provided",
      "complete",
      "complete",
      "complete",
      "partial_review",
      "complete",
      "pending_annotation_module",
      "complete",
      "gated_optional_not_run",
      "gated_optional_not_run"
    ),
    evidence_files = c(
      sample_table_path,
      "multi-sample Seurat merge object construction",
      "multi_sample_qc_by_sample.csv;multi_sample_cell_qc_metrics_pre_filter.csv;multi_sample_cell_qc_metrics_post_filter.csv;multi_sample_qc_pre_filter.pdf",
      "multi_sample_harmony_seurat.rds;sessionInfo_multi_sample_harmony.txt",
      batch_note,
      "umap_by_sample.pdf;resolution_sweep.tsv;cluster_marker_audit.tsv",
      "run scripts/course_adapted/05_singler_cell_annotation.R or equivalent manual annotation next",
      "multi_sample_cluster_markers.csv;cluster_marker_audit.tsv",
      "requires downstream_proposal.md and approved lineage/root",
      "requires downstream_proposal.md and approved sender/receiver groups"
    ),
    review_note = c(
      "Verify GEO/SRA/local provenance for every sample before interpretation.",
      "Samples were read, metadata columns copied, and objects merged.",
      "QC was performed by sample, but thresholds still require reviewer inspection.",
      "LogNormalize, variable features, scaling, and PCA were run; raw counts should remain in the Seurat object.",
      "Doublet, cell-cycle, and ambient RNA review are not completed by this merge script; document or run dedicated modules before strong interpretation.",
      "PCA, neighbors, resolution sweep, selected clustering, and UMAP were generated.",
      "Cell annotation is intentionally not completed in this merge pipeline.",
      "Cluster marker detection was run and marker specificity audit was exported.",
      "Pseudotime is optional and must not run before annotation and root approval.",
      "Cell communication is optional and must not run before annotation and scope approval."
    ),
    stringsAsFactors = FALSE
  )
  write.table(audit, file.path(output_dir, "tables", "workflow_step_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

samples <- read.csv(sample_table_path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
required_cols <- c("sample_id", "input_path", "input_type")
missing_cols <- setdiff(required_cols, colnames(samples))
if (length(missing_cols) > 0) {
  stop("sample table missing required columns: ", paste(missing_cols, collapse = ", "))
}

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

read_nonstandard_10x_counts <- function(path, sample_id) {
  matrix_file <- first_existing(file.path(path, c("count_matrix_sparse.mtx", "count_matrix_sparse.mtx.gz", "matrix.mtx", "matrix.mtx.gz")))
  barcode_file <- first_existing(file.path(path, c("count_matrix_barcodes.tsv", "count_matrix_barcodes.tsv.gz", "barcodes.tsv", "barcodes.tsv.gz")))
  gene_file <- first_existing(file.path(path, c("count_matrix_genes.tsv", "count_matrix_genes.tsv.gz", "features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz")))

  if (is.na(matrix_file) || is.na(barcode_file) || is.na(gene_file)) {
    stop("No non-standard 10x matrix set found for sample ", sample_id, ". Expected count_matrix_sparse.mtx, count_matrix_barcodes.tsv, and count_matrix_genes.tsv in the sample input_path.")
  }

  counts <- read_maybe_gz_mtx(matrix_file)
  barcodes <- read_maybe_gz_table(barcode_file)
  genes <- read_maybe_gz_table(gene_file)
  gene_names <- if (ncol(genes) >= 2) genes[[2]] else genes[[1]]

  if (ncol(counts) != nrow(barcodes) || nrow(counts) != length(gene_names)) {
    stop("Non-standard 10x dimensions do not match for sample ", sample_id)
  }

  rownames(counts) <- make.unique(as.character(gene_names))
  colnames(counts) <- as.character(barcodes[[1]])
  counts
}

read_one <- function(sample_id, input_path, input_type) {
  message("Reading sample ", sample_id, ": ", input_path)
  if (input_type == "10x_mtx") {
    counts <- Read10X(data.dir = input_path)
  } else if (input_type == "10x_nonstandard") {
    counts <- read_nonstandard_10x_counts(input_path, sample_id)
  } else if (input_type == "10x_h5") {
    counts <- Read10X_h5(filename = input_path)
  } else if (input_type == "csv") {
    counts <- read.csv(input_path, row.names = 1, check.names = FALSE, fileEncoding = "UTF-8-BOM")
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
ribo_pattern <- if (tolower(organism) == "mouse") "^Rp[sl]" else "^RP[SL]"
hb_pattern <- if (tolower(organism) == "mouse") "^Hb[abq]" else "^HB[ABDEGQMZ]"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = ribo_pattern)
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = hb_pattern)

qc <- obj@meta.data |>
  group_by(sample_id) |>
  summarise(cells = dplyr::n(), median_features = median(nFeature_RNA), median_counts = median(nCount_RNA), median_mt = median(percent.mt), median_ribo = median(percent.ribo), median_hb = median(percent.hb), .groups = "drop")
write.csv(qc, file.path(output_dir, "tables", "multi_sample_qc_by_sample.csv"), row.names = FALSE)
qc_cols <- c("sample_id", "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb")
write.csv(obj@meta.data[, qc_cols, drop = FALSE], file.path(output_dir, "tables", "multi_sample_cell_qc_metrics_pre_filter.csv"))

pdf(file.path(output_dir, "figures", "multi_sample_qc_pre_filter.pdf"), width = 14, height = 6)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), group.by = "sample_id", ncol = 5))
dev.off()

obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
write.csv(obj@meta.data[, qc_cols, drop = FALSE], file.path(output_dir, "tables", "multi_sample_cell_qc_metrics_post_filter.csv"))
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj))

has_batch_levels <- nzchar(batch_col) && batch_col %in% colnames(obj@meta.data) && length(unique(obj@meta.data[[batch_col]])) >= 2
if (has_batch_levels && requireNamespace("harmony", quietly = TRUE)) {
  obj <- harmony::RunHarmony(obj, group.by.vars = batch_col)
  reduction_use <- "harmony"
} else {
  reduction_use <- "pca"
}

dims_use <- seq_len(min(30, ncol(Embeddings(obj, "pca"))))
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
write_cluster_marker_audit(obj, markers)
write_workflow_step_audit(output_dir, sample_table_path, reduction_use, has_batch_levels)
saveRDS(obj, file.path(output_dir, "objects", "multi_sample_harmony_seurat.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "logs", "sessionInfo_multi_sample_harmony.txt"))

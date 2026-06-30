#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(data.table)
  library(SingleR)
  library(celldex)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  values <- list(
    input_dir = "data/processed/10x_by_sample",
    metadata = "metadata/gse223751_sample_metadata.csv",
    output_dir = "results/gse223751_seurat",
    min_features = 200,
    max_features = 8000,
    max_mt = 20,
    resolution = 0.4,
    dims = 30,
    max_cells = 0
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Arguments must be provided as --key value pairs")
    }
    name <- sub("^--", "", key)
    name <- gsub("-", "_", name)
    values[[name]] <- args[[i + 1]]
    i <- i + 2
  }
  values$min_features <- as.numeric(values$min_features)
  values$max_features <- as.numeric(values$max_features)
  values$max_mt <- as.numeric(values$max_mt)
  values$resolution <- as.numeric(values$resolution)
  values$dims <- as.integer(values$dims)
  values$max_cells <- as.integer(values$max_cells)
  values
}

opts <- parse_args()

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
for (subdir in c("figures", "tables", "objects", "logs")) {
  dir.create(file.path(opts$output_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}

sink(file.path(opts$output_dir, "logs", "run_messages.txt"), split = TRUE)
on.exit({
  writeLines(capture.output(sessionInfo()), file.path(opts$output_dir, "logs", "session_info.txt"))
  sink()
}, add = TRUE)

writeLines(commandArgs(), file.path(opts$output_dir, "logs", "commands.txt"))
message("Input directory: ", opts$input_dir)
message("Metadata: ", opts$metadata)
message("Output directory: ", opts$output_dir)

metadata <- read.csv(opts$metadata, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("sample_id", "gsm", "stage", "replicate", "condition", "tissue")
missing_cols <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0) {
  stop("Metadata missing required columns: ", paste(missing_cols, collapse = ", "))
}

objects <- list()
sample_qc <- list()
for (row_i in seq_len(nrow(metadata))) {
  row <- metadata[row_i, , drop = FALSE]
  sample_id <- row$sample_id[[1]]
  sample_dir <- file.path(opts$input_dir, sample_id)
  if (!dir.exists(sample_dir)) {
    stop("Missing organized 10x directory: ", sample_dir)
  }
  message("Reading sample ", sample_id, " from ", sample_dir)
  counts <- Read10X(data.dir = sample_dir)
  obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200, project = sample_id)
  obj$sample_id <- sample_id
  obj$gsm <- row$gsm[[1]]
  obj$stage <- row$stage[[1]]
  obj$replicate <- row$replicate[[1]]
  obj$condition <- row$condition[[1]]
  obj$tissue <- row$tissue[[1]]
  objects[[sample_id]] <- obj
  sample_qc[[sample_id]] <- data.frame(
    sample_id = sample_id,
    cells_raw = ncol(obj),
    genes_raw = nrow(obj),
    stringsAsFactors = FALSE
  )
}

merged <- Reduce(function(x, y) merge(x, y), objects)
DefaultAssay(merged) <- "RNA"

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^mt-")
merged[["percent.ribo"]] <- PercentageFeatureSet(merged, pattern = "^Rp[sl]")

qc_pre <- merged@meta.data |>
  tibble::rownames_to_column("cell_id") |>
  group_by(sample_id, stage) |>
  summarise(
    cells = n(),
    median_features = median(nFeature_RNA),
    median_counts = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    .groups = "drop"
  )
write.table(qc_pre, file.path(opts$output_dir, "tables", "qc_summary_pre_filter.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(opts$output_dir, "figures", "qc_pre_filter.pdf"), width = 12, height = 7)
print(VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), group.by = "sample_id", ncol = 2, pt.size = 0))
dev.off()

pdf(file.path(opts$output_dir, "figures", "qc_scatter.pdf"), width = 12, height = 5)
print(FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "sample_id") +
  FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample_id"))
dev.off()

merged <- subset(merged, subset = nFeature_RNA > opts$min_features & nFeature_RNA < opts$max_features & percent.mt < opts$max_mt)

if (opts$max_cells > 0 && ncol(merged) > opts$max_cells) {
  set.seed(20260630)
  cells <- unlist(lapply(split(colnames(merged), merged$sample_id), function(x) {
    n <- ceiling(opts$max_cells * length(x) / ncol(merged))
    sample(x, min(length(x), n))
  }), use.names = FALSE)
  merged <- subset(merged, cells = cells)
  message("Downsampled to ", ncol(merged), " cells for an example run")
}

qc_post <- merged@meta.data |>
  tibble::rownames_to_column("cell_id") |>
  group_by(sample_id, stage) |>
  summarise(
    cells = n(),
    median_features = median(nFeature_RNA),
    median_counts = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    .groups = "drop"
  )
write.table(qc_post, file.path(opts$output_dir, "tables", "qc_summary_post_filter.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(opts$output_dir, "figures", "qc_post_filter.pdf"), width = 12, height = 7)
print(VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), group.by = "sample_id", ncol = 2, pt.size = 0))
dev.off()

merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged, features = VariableFeatures(merged))
merged <- RunPCA(merged, features = VariableFeatures(merged))

pdf(file.path(opts$output_dir, "figures", "elbow_plot.pdf"), width = 7, height = 5)
print(ElbowPlot(merged, ndims = min(50, opts$dims + 10)))
dev.off()

dims_use <- seq_len(min(opts$dims, ncol(Embeddings(merged, "pca"))))
merged <- FindNeighbors(merged, reduction = "pca", dims = dims_use)
merged <- FindClusters(merged, resolution = opts$resolution)
merged <- RunUMAP(merged, reduction = "pca", dims = dims_use)
merged <- JoinLayers(merged)

pdf(file.path(opts$output_dir, "figures", "umap_clusters.pdf"), width = 7, height = 6)
print(DimPlot(merged, reduction = "umap", label = TRUE) + ggtitle("GSE223751 Seurat clusters"))
dev.off()

pdf(file.path(opts$output_dir, "figures", "umap_by_sample.pdf"), width = 8, height = 6)
print(DimPlot(merged, reduction = "umap", group.by = "sample_id") + ggtitle("Samples"))
dev.off()

pdf(file.path(opts$output_dir, "figures", "umap_by_stage.pdf"), width = 8, height = 6)
print(DimPlot(merged, reduction = "umap", group.by = "stage") + ggtitle("Developmental stage"))
dev.off()

message("Running SingleR cluster-level annotation")
ref <- celldex::MouseRNAseqData()
test_data <- GetAssayData(merged, assay = "RNA", layer = "data")
single_r <- SingleR::SingleR(
  test = test_data,
  ref = ref,
  labels = ref$label.main,
  clusters = merged$seurat_clusters
)
single_r_table <- data.frame(
  seurat_clusters = rownames(single_r),
  singleR_label = single_r$labels,
  singleR_pruned_label = if ("pruned.labels" %in% colnames(single_r)) single_r$pruned.labels else single_r$labels,
  stringsAsFactors = FALSE
)
single_r_table$singleR_pruned_label[is.na(single_r_table$singleR_pruned_label)] <- single_r_table$singleR_label[is.na(single_r_table$singleR_pruned_label)]
write.csv(single_r_table, file.path(opts$output_dir, "tables", "singleR_cluster_labels.csv"), row.names = FALSE)

pdf(file.path(opts$output_dir, "figures", "singleR_score_heatmap.pdf"), width = 9, height = 7)
print(SingleR::plotScoreHeatmap(single_r, clusters = rownames(single_r), order.by = "cluster"))
dev.off()

single_r_map <- setNames(single_r_table$singleR_pruned_label, single_r_table$seurat_clusters)
merged$singleR_label <- unname(single_r_map[as.character(merged$seurat_clusters)])

pdf(file.path(opts$output_dir, "figures", "umap_singleR_labels.pdf"), width = 8, height = 6)
print(DimPlot(merged, reduction = "umap", group.by = "singleR_label", label = TRUE) + ggtitle("SingleR annotation"))
dev.off()

markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, max.cells.per.ident = 500)
write.csv(markers, file.path(opts$output_dir, "tables", "cluster_markers.csv"), row.names = FALSE)

marker_sets <- list(
  chondrocyte = c("Col2a1", "Acan", "Sox9", "Matn4", "Comp"),
  tenocyte_fibroblast = c("Scx", "Tnmd", "Col1a1", "Col1a2", "Dcn", "Lum", "Pdgfra"),
  osteoblast = c("Sp7", "Runx2", "Bglap", "Ibsp", "Spp1"),
  endothelial = c("Pecam1", "Kdr", "Cdh5", "Vwf"),
  myeloid = c("Lyz2", "Csf1r", "Adgre1", "Cd68"),
  t_cell = c("Cd3d", "Cd3e", "Trac", "Lck"),
  b_cell = c("Ms4a1", "Cd79a", "Cd79b", "Bank1"),
  cycling = c("Mki67", "Top2a", "Pclaf", "Stmn1"),
  erythroid = c("Hbb-bs", "Hba-a1", "Hba-a2", "Alas2")
)
marker_sets <- lapply(marker_sets, function(x) intersect(x, rownames(merged)))
marker_sets <- marker_sets[lengths(marker_sets) > 0]

score_cols <- character()
if (length(marker_sets) > 0) {
  for (set_name in names(marker_sets)) {
    merged <- AddModuleScore(merged, features = list(marker_sets[[set_name]]), name = paste0("score_", set_name))
    score_cols <- c(score_cols, paste0("score_", set_name, "1"))
  }
  names(score_cols) <- names(marker_sets)
  score_df <- merged@meta.data |>
    tibble::rownames_to_column("cell_id") |>
    group_by(seurat_clusters) |>
    summarise(across(all_of(unname(score_cols)), mean), .groups = "drop")
  annotation <- data.frame(
    seurat_clusters = score_df$seurat_clusters,
    coarse_annotation = names(score_cols)[max.col(as.matrix(score_df[, unname(score_cols), drop = FALSE]))],
    stringsAsFactors = FALSE
  )
  score_out <- merge(score_df, annotation, by = "seurat_clusters")
  write.table(score_out, file.path(opts$output_dir, "tables", "cluster_annotation_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cluster_map <- setNames(annotation$coarse_annotation, annotation$seurat_clusters)
  merged$marker_support_label <- unname(cluster_map[as.character(merged$seurat_clusters)])

  annotation_evidence <- merge(single_r_table, annotation, by = "seurat_clusters", all.x = TRUE)
  colnames(annotation_evidence)[colnames(annotation_evidence) == "coarse_annotation"] <- "marker_support_label"
  annotation_evidence$final_example_label <- annotation_evidence$singleR_pruned_label
  annotation_evidence$note <- "SingleR is used as the automated annotation; marker scores are auxiliary evidence for review."
  write.table(annotation_evidence, file.path(opts$output_dir, "tables", "annotation_evidence.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  pdf(file.path(opts$output_dir, "figures", "umap_marker_support_labels.pdf"), width = 8, height = 6)
  print(DimPlot(merged, reduction = "umap", group.by = "marker_support_label", label = TRUE) + ggtitle("Marker support labels"))
  dev.off()

  dot_features <- unique(unlist(marker_sets))
  dot_features <- dot_features[dot_features %in% rownames(merged)]
  dot_features <- dot_features[seq_len(min(length(dot_features), 35))]
  pdf(file.path(opts$output_dir, "figures", "marker_dotplot.pdf"), width = 12, height = 7)
  print(DotPlot(merged, features = dot_features, group.by = "seurat_clusters") + RotatedAxis())
  dev.off()
}

cell_meta <- merged@meta.data |>
  tibble::rownames_to_column("cell_id")
write.table(cell_meta, file.path(opts$output_dir, "tables", "cell_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

prop_group <- if ("coarse_annotation" %in% colnames(cell_meta)) "coarse_annotation" else "seurat_clusters"
if ("singleR_label" %in% colnames(cell_meta)) {
  prop_group <- "singleR_label"
}
proportions <- cell_meta |>
  group_by(sample_id, stage, .data[[prop_group]]) |>
  summarise(cells = n(), .groups = "drop") |>
  group_by(sample_id) |>
  mutate(fraction = cells / sum(cells)) |>
  ungroup()
write.table(proportions, file.path(opts$output_dir, "tables", "celltype_proportions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(opts$output_dir, "figures", "celltype_proportions_by_sample.pdf"), width = 9, height = 5)
print(ggplot(proportions, aes(x = sample_id, y = fraction, fill = .data[[prop_group]])) +
  geom_col(width = 0.8) +
  theme_bw() +
  labs(x = "Sample", y = "Fraction", fill = prop_group))
dev.off()

saveRDS(merged, file.path(opts$output_dir, "objects", "gse223751_processed_seurat.rds"))
message("Done. Cells retained: ", ncol(merged), "; genes: ", nrow(merged))

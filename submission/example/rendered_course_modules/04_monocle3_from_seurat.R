#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SingleCellExperiment)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

seurat_rds <- "submission/example/results/gse223751_seurat/objects/gse223751_processed_seurat.rds"
output_dir <- "submission/example/results/course_modules/monocle3"
celltype_col <- "marker_support_label"
subset_query <- "marker_support_label %in% c('tenocyte_fibroblast', 'chondrocyte', 'osteoblast', 'cycling')"
root_cells_file <- ""
root_query <- "grepl('^E15$', stage)"
max_cells <- suppressWarnings(as.integer("5000"))
balance_col <- "marker_support_label"
use_seurat_umap <- tolower("true") %in% c("true", "t", "1", "yes", "y")
random_seed <- suppressWarnings(as.integer("20260630"))

if (is.na(max_cells)) max_cells <- 0
if (is.na(random_seed)) random_seed <- 20260630

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "objects"), showWarnings = FALSE)
dir.create(file.path(output_dir, "logs"), showWarnings = FALSE)

obj <- readRDS(seurat_rds)
stopifnot(celltype_col %in% colnames(obj@meta.data))

if (nzchar(subset_query)) {
  keep <- eval(parse(text = subset_query), envir = obj@meta.data)
  if (!is.logical(keep) || length(keep) != nrow(obj@meta.data)) {
    stop("SUBSET_QUERY must evaluate to one logical value per cell in obj@meta.data")
  }
  obj <- subset(obj, cells = rownames(obj@meta.data)[keep])
}

if (max_cells > 0 && ncol(obj) > max_cells) {
  set.seed(random_seed)
  if (nzchar(balance_col) && balance_col %in% colnames(obj@meta.data)) {
    split_cells <- split(colnames(obj), obj[[balance_col]][, 1])
    sampled <- unlist(lapply(split_cells, function(cells) {
      n <- ceiling(max_cells * length(cells) / ncol(obj))
      sample(cells, min(length(cells), n))
    }), use.names = FALSE)
    if (length(sampled) > max_cells) sampled <- sample(sampled, max_cells)
  } else {
    sampled <- sample(colnames(obj), max_cells)
  }
  obj <- subset(obj, cells = sampled)
}

celltype_counts <- sort(table(obj@meta.data[[celltype_col]]), decreasing = TRUE)
if (length(celltype_counts) < 2) {
  stop(
    "Monocle3 trajectory input has fewer than 2 states in ",
    celltype_col,
    ". Use a biologically related multi-state continuum or explicitly document this as a single-state trajectory audit."
  )
}
write.csv(
  data.frame(state = names(celltype_counts), cells = as.integer(celltype_counts)),
  file.path(output_dir, "tables", "trajectory_input_state_counts.csv"),
  row.names = FALSE
)

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
cell_metadata <- obj@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))

cds <- new_cell_data_set(counts, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = "UMAP")

p_native <- plot_cells(cds, color_cells_by = celltype_col, label_cell_groups = FALSE, show_trajectory_graph = FALSE) +
  ggtitle("Monocle3 native UMAP")
ggsave(file.path(output_dir, "figures", "monocle3_native_umap.pdf"), p_native, width = 7, height = 6)

if (use_seurat_umap && "umap" %in% Reductions(obj)) {
  seurat_umap <- Embeddings(obj, reduction = "umap")
  seurat_umap <- seurat_umap[colnames(cds), , drop = FALSE]
  colnames(seurat_umap) <- c("UMAP_1", "UMAP_2")
  SingleCellExperiment::reducedDims(cds)$UMAP <- seurat_umap
  p_seurat <- plot_cells(cds, color_cells_by = celltype_col, label_cell_groups = FALSE, show_trajectory_graph = FALSE) +
    ggtitle("Seurat UMAP imported into Monocle3")
  ggsave(file.path(output_dir, "figures", "seurat_imported_umap.pdf"), p_seurat, width = 7, height = 6)
  ggsave(file.path(output_dir, "figures", "monocle3_vs_seurat_umap.pdf"), p_native + p_seurat, width = 12, height = 6)
}

cds <- cluster_cells(cds, cluster_method = "louvain")
cds <- learn_graph(cds, use_partition = FALSE)

if (file.exists(root_cells_file)) {
  root_cells <- readLines(root_cells_file)
  cds <- order_cells(cds, root_cells = intersect(root_cells, colnames(cds)))
} else if (nzchar(root_query) && !grepl("^\\{\\{", root_query)) {
  root_keep <- eval(parse(text = root_query), envir = as.data.frame(colData(cds)))
  if (!is.logical(root_keep) || length(root_keep) != ncol(cds)) {
    stop("ROOT_QUERY must evaluate to one logical value per cell in colData(cds)")
  }
  root_cells <- colnames(cds)[root_keep]
  if (length(root_cells) == 0) stop("ROOT_QUERY selected no root cells")
  cds <- order_cells(cds, root_cells = root_cells)
} else {
  message("No root cell file found; order_cells() will open an interactive selector in some environments.")
  cds <- order_cells(cds)
}

p <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_principal_points = FALSE
)
ggsave(file.path(output_dir, "figures", "pseudotime_umap.pdf"), p, width = 7, height = 6)

p2_data <- data.frame(
  UMAP_1 = SingleCellExperiment::reducedDim(cds, "UMAP")[, 1],
  UMAP_2 = SingleCellExperiment::reducedDim(cds, "UMAP")[, 2],
  cell_state = colData(cds)[[celltype_col]]
)
p2_labels <- p2_data %>%
  group_by(cell_state) %>%
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
p2 <- ggplot(p2_data, aes(UMAP_1, UMAP_2, color = cell_state)) +
  geom_point(size = 0.45, alpha = 0.85) +
  ggrepel::geom_text_repel(data = p2_labels, aes(label = cell_state), color = "black", size = 3, show.legend = FALSE) +
  coord_equal() +
  labs(title = "Trajectory input cell states", color = celltype_col, x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 12)
ggsave(file.path(output_dir, "figures", "trajectory_celltypes.pdf"), p2, width = 7, height = 6)

p2_graph <- plot_cells(
  cds,
  color_cells_by = celltype_col,
  label_cell_groups = TRUE,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_principal_points = FALSE,
  cell_size = 0.45
)
ggsave(file.path(output_dir, "figures", "trajectory_celltypes_with_graph.pdf"), p2_graph, width = 7, height = 6)

p3 <- plot_cells(
  cds,
  color_cells_by = "stage",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_principal_points = FALSE
)
ggsave(file.path(output_dir, "figures", "trajectory_stage.pdf"), p3, width = 7, height = 6)

pseudotime_table <- data.frame(cell = colnames(cds), pseudotime = pseudotime(cds), colData(cds))
write.csv(pseudotime_table, file.path(output_dir, "tables", "pseudotime.csv"), row.names = FALSE)

trajectory_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1)
write.csv(trajectory_genes, file.path(output_dir, "tables", "trajectory_genes.csv"))

top_genes <- trajectory_genes %>%
  filter(!is.na(q_value), q_value < 1e-3) %>%
  arrange(desc(morans_I)) %>%
  head(10) %>%
  pull(gene_short_name) %>%
  as.character()
if (length(top_genes) > 0) {
  p4 <- plot_genes_in_pseudotime(cds[top_genes, ], color_cells_by = celltype_col, min_expr = 0.5, ncol = 2)
  ggsave(file.path(output_dir, "figures", "top_trajectory_genes_pseudotime.pdf"), p4, width = 8, height = 8)
}

summary_lines <- c(
  paste("Seurat object:", seurat_rds),
  paste("Cells in trajectory:", ncol(cds)),
  paste("Genes:", nrow(cds)),
  paste("Subset query:", subset_query),
  paste("Max cells:", max_cells),
  paste("Balance column:", balance_col),
  paste("Use Seurat UMAP:", use_seurat_umap),
  paste("Finite pseudotime cells:", sum(is.finite(pseudotime(cds)))),
  paste("Unique finite pseudotime values:", length(unique(pseudotime(cds)[is.finite(pseudotime(cds))])))
)
writeLines(summary_lines, file.path(output_dir, "logs", "monocle3_status.txt"))

saveRDS(cds, file.path(output_dir, "objects", "monocle3_cds.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "logs", "sessionInfo_monocle3.txt"))

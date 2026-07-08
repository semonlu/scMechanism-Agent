#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(ggplot2)
})

seurat_rds <- "{{SEURAT_RDS}}"
output_dir <- "{{OUTPUT_DIR}}"
celltype_col <- "{{CELLTYPE_COL}}"
subset_query <- "{{SUBSET_QUERY}}"
root_cells_file <- "{{ROOT_CELLS_FILE}}"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "objects"), showWarnings = FALSE)

obj <- readRDS(seurat_rds)
stopifnot(celltype_col %in% colnames(obj@meta.data))

if (nzchar(subset_query)) {
  obj <- subset(obj, subset = eval(parse(text = subset_query)))
}

counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
cell_metadata <- obj@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))

cds <- new_cell_data_set(counts, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds, cluster_method = "louvain")
cds <- learn_graph(cds)

if (file.exists(root_cells_file)) {
  root_cells <- readLines(root_cells_file, encoding = "UTF-8")
  cds <- order_cells(cds, root_cells = intersect(root_cells, colnames(cds)))
} else {
  message("No root cell file found; order_cells() will open an interactive selector in some environments.")
  cds <- order_cells(cds)
}

p <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
ggsave(file.path(output_dir, "figures", "pseudotime_umap.pdf"), p, width = 7, height = 6)

p2 <- plot_cells(cds, color_cells_by = celltype_col, label_cell_groups = FALSE)
ggsave(file.path(output_dir, "figures", "trajectory_celltypes.pdf"), p2, width = 7, height = 6)

pseudotime_table <- data.frame(cell = colnames(cds), pseudotime = pseudotime(cds), colData(cds))
write.csv(pseudotime_table, file.path(output_dir, "tables", "pseudotime.csv"), row.names = FALSE)

trajectory_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1)
write.csv(trajectory_genes, file.path(output_dir, "tables", "trajectory_genes.csv"))

saveRDS(cds, file.path(output_dir, "objects", "monocle3_cds.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_monocle3.txt"))

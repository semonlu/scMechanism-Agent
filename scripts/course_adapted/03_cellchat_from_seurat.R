#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
})

seurat_rds <- "{{SEURAT_RDS}}"
output_dir <- "{{OUTPUT_DIR}}"
celltype_col <- "{{CELLTYPE_COL}}"
organism <- "{{ORGANISM}}" # human or mouse
min_cells <- as.integer("{{MIN_CELLS}}")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "objects"), showWarnings = FALSE)

obj <- readRDS(seurat_rds)
stopifnot(celltype_col %in% colnames(obj@meta.data))

tab <- table(obj[[celltype_col]][, 1])
keep_groups <- names(tab)[tab >= min_cells]
obj <- subset(obj, cells = rownames(obj@meta.data)[obj[[celltype_col]][, 1] %in% keep_groups])

data.input <- GetAssayData(obj, assay = "RNA", slot = "data")
meta <- data.frame(cell_type = obj[[celltype_col]][, 1], row.names = colnames(obj))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")
CellChatDB <- if (tolower(organism) == "mouse") CellChatDB.mouse else CellChatDB.human
cellchat@DB <- subsetDB(CellChatDB, search = "Secreted Signaling")

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = min_cells)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

lr <- subsetCommunication(cellchat)
write.csv(lr, file.path(output_dir, "tables", "cellchat_ligand_receptor.csv"), row.names = FALSE)

groupSize <- as.numeric(table(cellchat@idents))
pdf(file.path(output_dir, "figures", "cellchat_interaction_count.pdf"))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE)
dev.off()

pdf(file.path(output_dir, "figures", "cellchat_interaction_weight.pdf"))
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE)
dev.off()

saveRDS(cellchat, file.path(output_dir, "objects", "cellchat.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_cellchat.txt"))

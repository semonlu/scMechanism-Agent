#!/usr/bin/env Rscript

# Adapted from course modules:
# - 17.多种算法计算marker基因/多种算法计算marker基因.R
# - 18.GO_KEGG/GO_KEGG.r

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

seurat_rds <- "D:/单细胞测试/scMechanism-Agent/submission/example/results/gse223751_seurat/objects/gse223751_processed_seurat.rds"
output_dir <- "D:/单细胞测试/scMechanism-Agent/submission/example/results/course_modules/marker_enrichment"
cluster_col <- "singleR_label"
organism <- "mouse" # human or mouse
top_n <- suppressWarnings(as.integer("20"))
if (is.na(top_n) || top_n <= 0) {
  top_n <- 20
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "logs"), showWarnings = FALSE)

obj <- readRDS(seurat_rds)
if (nzchar(cluster_col) && cluster_col %in% colnames(obj@meta.data)) {
  Idents(obj) <- obj[[cluster_col]][, 1]
} else {
  cluster_col <- "active.ident"
}

markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(markers, file.path(output_dir, "tables", "cluster_markers_all.csv"), row.names = FALSE)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) %>%
  ungroup()
write.csv(top_markers, file.path(output_dir, "tables", "cluster_markers_top.csv"), row.names = FALSE)

status <- c(
  paste("Seurat object:", seurat_rds),
  paste("Identity column:", cluster_col),
  paste("Organism:", organism),
  paste("Top markers per cluster:", top_n)
)

can_enrich <- requireNamespace("clusterProfiler", quietly = TRUE)
is_human <- tolower(organism) %in% c("human", "homo_sapiens", "homo sapiens", "hs")
is_mouse <- tolower(organism) %in% c("mouse", "mus_musculus", "mus musculus", "mm")
org_pkg <- if (is_mouse) "org.Mm.eg.db" else "org.Hs.eg.db"
can_org <- requireNamespace(org_pkg, quietly = TRUE)

if (can_enrich && can_org) {
  OrgDb <- getExportedValue(org_pkg, org_pkg)
  genes <- unique(top_markers$gene)
  gene_map <- clusterProfiler::bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = OrgDb
  )

  if (nrow(gene_map) > 0) {
    go <- clusterProfiler::enrichGO(
      gene = unique(gene_map$ENTREZID),
      OrgDb = OrgDb,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      readable = TRUE
    )
    write.csv(as.data.frame(go), file.path(output_dir, "tables", "GO_BP_top_markers.csv"), row.names = FALSE)

    if (is_human || is_mouse) {
      kegg_code <- if (is_mouse) "mmu" else "hsa"
      kegg <- clusterProfiler::enrichKEGG(
        gene = unique(gene_map$ENTREZID),
        organism = kegg_code,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05
      )
      write.csv(as.data.frame(kegg), file.path(output_dir, "tables", "KEGG_top_markers.csv"), row.names = FALSE)
    }
    status <- c(status, "Enrichment: completed for top marker set.")
  } else {
    status <- c(status, "Enrichment: skipped because no SYMBOL to ENTREZID mapping was found.")
  }
} else {
  status <- c(
    status,
    paste("Enrichment: skipped; requires clusterProfiler and", org_pkg)
  )
}

writeLines(status, file.path(output_dir, "logs", "marker_enrichment_status.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "logs", "sessionInfo_marker_enrichment.txt"))

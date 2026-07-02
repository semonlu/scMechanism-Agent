# Course-derived reference script
# English filename: 26_monocle3_pseudotime.R
# Original course path: see source_manifest.csv
# Role: trajectory reference
# Adaptations applied for the skill package:
# - Filename flattened and translated to English.
# - Interactive working-directory selection replaced with SC_WORK_DIR/getwd() when present.
# - Example-specific object names, thresholds, metadata columns, and local filenames still require review.
# - Prefer scripts/course_adapted/ for runnable project workflows.
####################################################
####################################################

#设置工作路径
WORK_DIR <- Sys.getenv("SC_WORK_DIR", unset = getwd())
setwd(WORK_DIR)
#加载R包
library(Seurat)
library(SeuratObject)
library(celldex)
library(assertthat)
library(monocle3)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(tricycle)   
library(scattermore)
library(scater)
library(patchwork)
library(ggplot2)
library(CCA)
library(clustree)
library(cowplot)
library(monocle3)
library(tidyverse)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)
library(future)
library(ggplot2)
library(ggforce)
library(ggsci)
library(parallel)
library(doParallel)
library(data.table)
library(qs)
#生成随机颜色
randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}

# 生成100个随机颜色
randomColors <- replicate(100,randomColor())
scRNA=qread("deepseek注释.qs",nthreads = detectCores())#读取数据
###拟时序分析
##随机抽取1%个细胞，仅作为演示，实际操作或高性能工作站直接运行scRNA_tpm<-scRNA
N=length(colnames(scRNA))/100
N=round(N)
scRNA_tpm<-scRNA[,sample(x=colnames(scRNA),size = N,replace=F)]
data <-GetAssayData(scRNA_tpm, assay = 'RNA', slot = 'counts')
cell_metadata <- scRNA_tpm@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data),row.names = row.names(data))
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#标准化
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
#umap降维
#umap.fast_sgd=T进行随机降维,cores为多线程运算数量,这两个参数都能加快处理速度,但是windows系统下无效
cds <- reduce_dimension(cds,reduction_method = "UMAP",cores = 1,umap.fast_sgd=F)
plot_cells(cds)
colnames(colData(cds))
#以之前的Seurat坐标来添加颜色
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="cellType") + ggtitle('cds.umap')
p1
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNA_tpm, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="cellType") + ggtitle('int.umap')
p2
#绘制与原来UMAP的对比
p = p1|p2
p
ggsave("1.Monocle3与seurat对比.pdf", plot = p, width = 10, height = 5)
#展示指定基因
genes <- c("TP53","FAM41C","FAM87B")
pdf("2.特定基因表达.pdf")
plot_cells(cds,
           genes=genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()
## 识别轨迹
cds <- cluster_cells(cds,cluster_method = "louvain")#如果出现Error in leidenbase::leiden_find_partition(graph_result[["g"]], partition_type = partition_type,  : REAL() can only be applied to a 'numeric', not a 'NULL'报错，将cluster_method设置为“louvain“，或者将 igraph 软件包降级到 1.4.3 版本
plot_cells(cds, color_cells_by = "cellType")
cds <- learn_graph(cds)
p = plot_cells(cds, color_cells_by = "cellType", label_groups_by_cluster=FALSE,
               label_leaves=FALSE, label_branch_points=FALSE)
ggsave("3.细胞轨迹.pdf", plot = p, width = 8, height = 6)
p=plot_cells(cds,
             color_cells_by = "cellType",
             label_cell_groups=FALSE,
             label_leaves=TRUE,
             label_branch_points=TRUE,
             graph_label_size=1.5)
ggsave("4.发育节点.pdf", plot = p, width = 8, height = 6)

#选择起始点
cds <- order_cells(cds)
p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
               label_leaves = FALSE,  label_branch_points = FALSE)
ggsave("5.手动节点发育时间.pdf", plot = p, width = 8, height = 6)
#筛选时序差异基因
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
write.csv(Track_genes, "6.时序差异基.csv", row.names = F)
#挑选前10个基因进行绘图
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
#绘制这10个基因的趋势图
p <- plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="cellType", 
                              min_expr=0.5, ncol = 2)
ggsave("7.基因趋势图.pdf", plot = p, width = 8, height = 6)
####################################################
####################################################
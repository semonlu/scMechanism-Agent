# Course-derived reference script
# English filename: 15_clustering_resolution.R
# Original course path: see source_manifest.csv
# Role: clustering reference
# Adaptations applied for the skill package:
# - Filename flattened and translated to English.
# - Interactive working-directory selection replaced with SC_WORK_DIR/getwd() when present.
# - Example-specific object names, thresholds, metadata columns, and local filenames still require review.
# - Prefer scripts/course_adapted/ for runnable project workflows.
####################################################
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
####################################################

#设置工作路径
WORK_DIR <- Sys.getenv("SC_WORK_DIR", unset = getwd())
setwd(WORK_DIR)

#加载R包
library(Seurat)
library(SeuratObject)
library(fastSave)
library(celldex)
library(assertthat)
library(monocle)
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
library(ggforce)
library(ggsci)
library(parallel)
library(doParallel)
library(data.table)
library(qs)
scRNA=qread("scRNA常规标准化后.qs",nthreads = detectCores())#读取数据
#选择PC
scRNAPC=30
##对细胞聚类
scRNA=FindNeighbors(scRNA, dims = 1:scRNAPC, reduction = "harmony")
#挑选分辨率
#首先查看0-0.1节点
for (res in seq(0,0.1,by=0.01)) { 
  scRNA=FindClusters(scRNA, graph.name = "RNA_snn", resolution = res, algorithm = 1)}
apply(scRNA@meta.data[,grep("RNA_snn_res",colnames(scRNA@meta.data))],2,table)
p2_tree=clustree(scRNA@meta.data, prefix = "RNA_snn_res.")
pdf(file = "1. 0至0.1节点分辨率.pdf",width =12,height =10)
p2_tree
dev.off()
#先选择较低分辨率进行粗注释
scRNA=FindNeighbors(scRNA, dims = 1:scRNAPC, reduction = "harmony")
scRNA <- FindClusters(scRNA, resolution = 0.01)#选择分辨率进行降维

#选定PC后降维
scRNA <- RunUMAP(scRNA, dims = 1:scRNAPC, reduction = "harmony")
scRNA <- RunTSNE(scRNA, dims = 1:scRNAPC, reduction = "harmony")

pdf(file = "2.选定PC降维UMAP图.pdf",width =6.5,height = 5.5)
DimPlot(scRNA, reduction = "umap", label = T, label.size = 3.5,pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "3.选定PC降维tsne图.pdf",width =6.5,height = 5.5)
DimPlot(scRNA, reduction = "tsne", label = T, label.size = 3.5,pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()

qsave(scRNA,"scRNA分群聚类后.qs",nthreads = detectCores())
####################################################
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
####################################################
# Course-derived reference script
# English filename: 14_post_doublet_normalization.R
# Original course path: see source_manifest.csv
# Role: normalization reference
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
library(monocle)
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
library(ggforce)
library(ggsci)
library(parallel)
library(doParallel)
library(data.table)
library(qs)
library(decontX)
#生成随机颜色
randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}

# 生成100个随机颜色
randomColors <- replicate(100,randomColor())
scRNA=qread("scRNA剔除双细胞后.qs",nthreads = detectCores())#读取数据
##标准化,使用LogNormalize方法
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
## 鉴定高变基因
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
# 提取前10的高变基因
top10 <- head(VariableFeatures(scRNA), 10)
# 选择全部基因归一化
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
# PCA降维
scRNA <- Seurat::RunPCA(scRNA, features = VariableFeatures(object = scRNA))
scRNA <- Seurat::RunTSNE(scRNA,dims = 1:20)
scRNA <- Seurat::RunUMAP(scRNA,dims = 1:20)
#harmony 去批次
scRNA <- RunHarmony(scRNA, group.by.vars = "Type")
##重新鉴定高变基因
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA <- Seurat::RunTSNE(scRNA,dims = 1:20,reduction ='harmony')
scRNA <- Seurat::RunUMAP(scRNA,dims = 1:20,reduction ='harmony')
scRNA <- Seurat::RunPCA(scRNA, features = VariableFeatures(object = scRNA))
qsave(scRNA,"scRNA常规标准化后.qs",nthreads = detectCores())
####################################################
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
####################################################

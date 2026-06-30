# Course-derived reference script
# English filename: 13_scdblfinder.R
# Original course path: see source_manifest.csv
# Role: doublet reference
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
library(scDblFinder)
library(tidyverse)
library(patchwork)
library(glmGamPoi)
library(BiocParallel)
options(future.globals.maxSize = 1024^1000)
#生成随机颜色
randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}

# 生成100个随机颜色
randomColors <- replicate(100,randomColor())
scRNA=qread("scRNA标准化后.qs",nthreads = detectCores())#读取数据
#转换数据结构
scRNAdbin<-as.SingleCellExperiment(scRNA)
scRNAdb <- scDblFinder(scRNAdbin, samples ="Type")
#转换回Seurat对象
scRNA <- as.Seurat(scRNAdb)
#查看双细胞分布
pdf("双细胞分布图.pdf")
DimPlot(scRNA, group.by = "scDblFinder.class", raster = FALSE)
dev.off()
#剔除双细胞
scRNA <- subset(scRNA, subset = (scDblFinder.class== "singlet"))
qsave(scRNA,"scRNA剔除双细胞后.qs",nthreads = detectCores())
####################################################
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
####################################################
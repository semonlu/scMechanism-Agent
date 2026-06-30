# Course-derived reference script
# English filename: 12_doublet_finder.R
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
library(DoubletFinder)
library(tidyverse)
library(patchwork)
library(glmGamPoi)
options(future.globals.maxSize = 1024^1000)
#生成随机颜色
randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}

# 生成100个随机颜色
randomColors <- replicate(100,randomColor())
scRNA=qread("scRNA标准化后.qs",nthreads = detectCores())#读取数据
#进行pN-pK参数扫描
sweeplist <- paramSweep(scRNA, PCs = 1:20, sct = T)
#计算双细胞检测指标
sweepstats <- summarizeSweep(sweeplist, GT = FALSE) 
sweepstats[order(sweepstats$BCreal),]
#找到最优pK参数
bcmvn <- find.pK(sweepstats) 
#提取最优pK值
pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)]) 
#估计的同源双细胞    
homotypic.prop <- modelHomotypic(scRNA$Type) 
#计算总的双细胞数量（假设双细胞形成率为 5%）
nExp_poi <- round(0.05 *nrow(scRNA@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # 计算异源双细胞数量
#鉴定doublets
scRNA <- doubletFinder(scRNA, PCs = 1:20, pN = 0.25,  pK = pK_bcmvn,
                                  nExp = nExp_poi.adj, reuse.pANN = F, sct = T)

#去除双细胞
scRNA <- subset(scRNA, subset = (DF.classifications_0.25_19_417== "Singlet"))
qsave(scRNA,"scRNA剔除双细胞后.qs",nthreads = detectCores())
####################################################
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
####################################################
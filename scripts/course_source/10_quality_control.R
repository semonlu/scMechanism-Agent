# Course-derived reference script
# English filename: 10_quality_control.R
# Original course path: see source_manifest.csv
# Role: QC reference
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
#生成随机颜色
randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}

# 生成100个随机颜色
randomColors <- replicate(100,randomColor())
scRNA=qread("scRNA.qs",nthreads = detectCores())

# QC指标使用小提琴图可视化,ncol为图片排列列数
pdf(file = "1.QC前小提琴图.pdf",width = 13,height = 8)
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4)+scale_fill_manual(values =randomColors)
dev.off()

# 指标之间的相关性
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.rb")+ RotatedAxis()
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ RotatedAxis()
#组图
pdf(file = "2.QC指标相关性.pdf",width =20,height = 8)
plot2+plot3+plot_layout(ncol = 2)      
dev.off()
#调整QC
mask1 <- scRNA$nCount_RNA <= 25000 #细胞内检测到的分子总数
mask2 <- scRNA$nFeature_RNA >= 0 & scRNA$nFeature_RNA <= 5000#每个细胞中检测到的基因数量
mask3 <- scRNA$percent.mt <= 25#结合线粒体基因（percent.mt）去除异常值
mask4<-scRNA$percent.rb<= 40#核糖体基因（percent.rb）除去异常值
mask <-mask1 & mask2 & mask3 & mask4
scRNA_out <- scRNA[, mask]
pdf(file = "3.QC后小提琴图.pdf",width = 13,height = 8)
VlnPlot(scRNA_out, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4)+scale_fill_manual(values =randomColors)
dev.off()
plot2 <- FeatureScatter(scRNA_out, feature1 = "nCount_RNA", feature2 = "percent.rb")+ RotatedAxis()
plot3 <- FeatureScatter(scRNA_out, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ RotatedAxis()
pdf(file = "4.QC后指标相关性.pdf",width =9,height = 8)
plot2+plot3+plot_layout(ncol = 2)      
dev.off()

#QC后保存数据
scRNA<-scRNA_out
qsave(scRNA,"scRNA数据质控后.qs",nthreads = detectCores())
####################################################
####################################################
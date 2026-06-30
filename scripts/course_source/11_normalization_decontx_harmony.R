# Course-derived reference script
# English filename: 11_normalization_decontx_harmony.R
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
scRNA=qread("scRNA数据质控后.qs",nthreads = detectCores())#读取数据

##标准化,使用LogNormalize方法
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
## 鉴定高变基因
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
# 提取前10的高变基因
top10 <- head(VariableFeatures(scRNA), 10)
# 展示高变基因
plot1 <- VariableFeaturePlot(scRNA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(file = "1.高变基因.pdf",width =7,height = 6)
plot2                   
dev.off()
# 选择全部基因归一化
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
# PCA降维
scRNA <- Seurat::RunPCA(scRNA, features = VariableFeatures(object = scRNA))
scRNA <- Seurat::RunTSNE(scRNA,dims = 1:20)
scRNA <- Seurat::RunUMAP(scRNA,dims = 1:20)
pdf(file = "2.tsne图.pdf",width =7.5,height = 5.5)
DimPlot(scRNA, reduction = "tsne",pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right") #top为图列位置最上方，除此之外还有right、left、bottom(意思同英文)
dev.off()
pdf(file = "3.PCA图.pdf",width =7.5,height = 5.5)
DimPlot(scRNA, reduction = "pca",pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "4.UMAP图.pdf",width =7.5,height = 5.5)
DimPlot(scRNA, reduction = "umap",pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "5.每个样本tsne图.pdf",width =12,height = 4)
do_DimPlot(sample = scRNA,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1,2),split.by = "Type",pt.size = 0.01
) 
dev.off()

#去除游离RNA污染
counts<-scRNA@assays$RNA@layers$counts
decontxcounts<-decontX(counts)
scRNA$Contamination=decontxcounts$contamination
scRNAde=scRNA[,scRNA$Contamination<0.3]
#harmony 去批次
scRNAde <- RunHarmony(scRNAde, group.by.vars = "Type")
##重新鉴定高变基因
scRNAde <- FindVariableFeatures(scRNAde, selection.method = "vst", nfeatures = 2000)
#提取前10的高变基因
top10 <- head(VariableFeatures(scRNAde), 10)
# 展示高变基因
plot1 <- VariableFeaturePlot(scRNAde)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(file = "6.去批后高变基因.pdf",width =7,height = 6)
plot2                  
dev.off()
#降维可视化
pdf(file = "7.Harmony去批次图.pdf",width =7.5,height = 5.5)
DimPlot(scRNAde, reduction = "harmony",pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
scRNAde <- Seurat::RunTSNE(scRNAde,dims = 1:20,reduction ='harmony')
scRNAde <- Seurat::RunUMAP(scRNAde,dims = 1:20,reduction ='harmony')
pdf(file = "8.去批次后tsne图.pdf",width =7.5,height = 5.5)
DimPlot(scRNAde, reduction = "tsne",pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "9.去批次后umap图.pdf",width =7.5,height = 5.5)
DimPlot(scRNAde, reduction = "tsne",pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "10.去批次后每个样本tsne图.pdf",width =12,height = 7.5)
do_DimPlot(sample = scRNAde,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1,2),split.by = "Type",pt.size =0.01
) 
dev.off()

#热图可视化前15个PC
pdf(file = "11.前15个PC可视化.pdf",width =7.5,height = 9)
DimHeatmap(scRNAde, dims = 1:15, cells = 1000, balanced = TRUE)
dev.off()
##确定使用PC个数
scRNAde <- JackStraw(scRNAde, num.replicate = 100)
scRNAde <- ScoreJackStraw(scRNAde, dims = 1:20)
pdf(file = "12.jackstrawplot.pdf",width =7.5,height = 5.5)
JackStrawPlot(scRNAde, dims = 1:20)
dev.off()
pdf(file = "13.ElbowPlot.pdf",width =5,height = 4)
ElbowPlot(scRNAde,ndims = 30)
dev.off()
qsave(scRNAde,"scRNA标准化后.qs",nthreads = detectCores())
####################################################
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
####################################################

# Course-derived reference script
# English filename: 22_scpred_annotation.R
# Original course path: see source_manifest.csv
# Role: annotation reference
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
library(Seurat)
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
library(ggplot2)
library(ggforce)
library(ggsci)
library(parallel)
library(doParallel)
library(data.table)
library(scPred)
library(qs)
options(future.globals.maxSize = 1024^1000)
#读取已注释数据(参考数据集)
scRNA=qread("deepseek注释.qs",nthreads = detectCores())#读取数据
#读取未注释数据
unscRNA=qread("正常细胞.qs",nthreads = detectCores())#读取数据
#进行模型训练
scRNAmodle<-getFeatureSpace(scRNA,"cellType")
cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)
scRNAmodle<-trainModel(scRNAmodle,allowParallel = T)
stopCluster(cl)
#进行自动注释
unscRNA<-scPredict(unscRNA,scRNAmodle)
#注释后可视化
pdf(file = "1.scPredict注释UMAP.pdf",width =7.5,height = 5.5)
DimPlot(unscRNA,group.by = "scpred_prediction", reduction = "umap", label = T,label.box = T, label.size = 3.5,pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "2.scPredict注释TSEN.pdf",width =7.5,height = 5.5)
DimPlot(unscRNA,group.by = "scpred_prediction", reduction = "tsne", , label = T,label.box = T, label.size = 3.5,pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
#基因分布图
seuratgenes=c("CCL5","GZMB")      #写入你想展示的基因
pdf(file = "3.基因分布UMAP.pdf",width =12,height = 6)
FeaturePlot(unscRNA, features = seuratgenes, reduction ="umap",  cols = c("grey", "red"),min.cutoff = 0, max.cutoff = 1,ncol=2,pt.size = 0.01)
dev.off()
pdf(file = "4.基因分布TSNE.pdf",width =12,height = 6)
FeaturePlot(unscRNA, features = seuratgenes, reduction ="tsne",  cols = c("grey", "red"),min.cutoff = 0, max.cutoff = 1,ncol=2,pt.size = 0.01)
dev.off()
#高亮度基因分布
pdf(file = "5.高亮基因分布TSNE.pdf",width =12,height = 6)
ggrastr::rasterize(Nebulosa::plot_density(unscRNA, seuratgenes,size=0.01,reduction = "tsne",shape = 2),dpi=300)
dev.off()
pdf(file = "6.高亮基因分布UMAP.pdf",width =12,height = 6)
ggrastr::rasterize(Nebulosa::plot_density(unscRNA, seuratgenes,size=0.01,reduction = "umap",shape = 2),dpi=300)
dev.off()
qsave(unscRNA,"scPredict注释.qs",,nthreads = detectCores())
####################################################
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
####################################################
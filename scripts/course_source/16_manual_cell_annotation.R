# Course-derived reference script
# English filename: 16_manual_cell_annotation.R
# Original course path: see source_manifest.csv
# Role: annotation reference
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
library(qs)
scRNA=qread("scRNA分群聚类后.qs",nthreads = detectCores())#读取数据
#初步注释
#将正常组肿瘤组分开
scRNAf<-SplitObject(scRNA,split.by = "Type")#将所有细胞类型分开
Normaldata<-scRNAf$GSM5580155
Tumordata<-scRNAf$GSM5580154
#独立分开后常规标准流程
##标准化,使用LogNormalize方法
Normaldata <- NormalizeData(Normaldata, normalization.method = "LogNormalize", scale.factor = 10000)
## 鉴定高变基因
Normaldata <- FindVariableFeatures(Normaldata, selection.method = "vst", nfeatures = 2000)
# 选择全部基因归一化
Normaldataall.genes <- rownames(Normaldata)
Normaldata <- ScaleData(Normaldata, features = Normaldataall.genes)
# PCA降维
Normaldata <- Seurat::RunPCA(Normaldata, features = VariableFeatures(object = Normaldata))
Normaldata <- Seurat::RunTSNE(Normaldata,dims = 1:20)
Normaldata <- Seurat::RunUMAP(Normaldata,dims = 1:20)
#将细胞注释为normal
Normaldataidens=mapvalues(Idents(Normaldata), from = levels(Idents(Normaldata)), to = rep("Normal", times = length(levels(Idents(Normaldata)))))
Idents(Normaldata)=Normaldataidens
Normaldata$cellType=Idents(Normaldata)
#导出备用
qsave(Normaldata,"正常细胞.qs",nthreads = detectCores())
#独立分开后常规标准流程
##标准化,使用LogNormalize方法
Tumordata <- NormalizeData(Tumordata, normalization.method = "LogNormalize", scale.factor = 10000)
## 鉴定高变基因
Tumordata <- FindVariableFeatures(Tumordata, selection.method = "vst", nfeatures = 2000)
# 选择全部基因归一化
Tumordataall.genes <- rownames(Tumordata)
Tumordata <- ScaleData(Tumordata, features = Tumordataall.genes)
# PCA降维
Tumordata <- Seurat::RunPCA(Tumordata, features = VariableFeatures(object = Tumordata))
Tumordata <- Seurat::RunTSNE(Tumordata,dims = 1:20)
Tumordata <- Seurat::RunUMAP(Tumordata,dims = 1:20)
scRNA<-Tumordata
genes <- list("immune"=c("PTPRC"),
              "epithelial"=c("EPCAM"),
              "stromal"=c("MME","PECAM1")
)
pdf(file = "1.初步注释标识.pdf",width = 20,height = 15)
do_DotPlot(sample = scRNA,features = genes,dot.scale = 12,
           legend.framewidth = 2, font.size =10)
dev.off()
ann.ids <- c("immune",#0
             "immune",#1
             "immune",#2
             "epithelial",#3
             "unknow",#4
             "unknow",#5
             "stromal",#6
             "immune"#7
)
scRNAidens=mapvalues(Idents(scRNA), from = levels(Idents(scRNA)), to = ann.ids)
Idents(scRNA)=scRNAidens
scRNA$cellType=Idents(scRNA)
#提取免疫细胞进行二次注释
scRNAf<-SplitObject(scRNA,split.by = "cellType")#将所有细胞类型分开
immune<-scRNAf$immune
#独立分开后常规标准流程
##标准化,使用LogNormalize方法
immune <- NormalizeData(immune, normalization.method = "LogNormalize", scale.factor = 10000)
## 鉴定高变基因
immune <- FindVariableFeatures(immune, selection.method = "vst", nfeatures = 2000)
# 选择全部基因归一化
immuneall.genes <- rownames(immune)
immune <- ScaleData(immune, features = immuneall.genes)
# PCA降维
immune <- Seurat::RunPCA(immune, features = VariableFeatures(object = immune))
immune <- Seurat::RunTSNE(immune,dims = 1:20)
immune <- Seurat::RunUMAP(immune,dims = 1:20)
immune=FindNeighbors(immune, dims = 1:20, reduction = "harmony")
#挑选分辨率
for (res in seq(0.1,1.0,by=0.1)) { 
  immune=FindClusters(immune, graph.name = "RNA_snn", resolution = res, algorithm = 1)}
apply(immune@meta.data[,grep("RNA_snn_res",colnames(immune@meta.data))],2,table)
p2_tree=clustree(immune@meta.data, prefix = "RNA_snn_res.")
pdf(file = "2.immune分辨率.pdf",width =12,height =10)
p2_tree
dev.off()
#第二次注释
immune <- FindClusters(immune, resolution =0.2, reduction = "harmony")
genes <- list("T Cells"= c("CD3D", "CD3E", "CD8A","PTPRC","CD4","CD2"),
              "B cells "=c("CD19", "CD79A", "MS4A1"),
              "Plasma cells"= c("IGHG1", "MZB1", "SDC1", "CD79A"),
              "Monocytes"=c("CD68", "CD163", "CD14"),
              "NK Cells"= c("GNLY","NKG7","KLRD1"),
              "Fibroblasts"= c("FGF7", "MME","GSN","LUM","DCN"),
              "Endothelial cells"= c("PECAM1", "VWF"),
             "Mast"=c("CPA3", "CST3", "KIT", "TPSAB1", "TPSB2", "MS4A2"))
pdf(file = "3.immune第二次注释.pdf",width = 13,height = 15)
do_DotPlot(sample = immune,features = genes,dot.scale = 12,
           legend.framewidth = 2, font.size =10)
dev.off()
ann.ids <- c("T cells",#0
             "T cells",#1
             "T cells",#2
             "Mast",#3
             "B cells",#4
             "Monocytes",#5
             "T cells",#6
             "T cells"#7
)
immuneidens=mapvalues(Idents(immune), from = levels(Idents(immune)), to = ann.ids)
Idents(immune)=immuneidens
immune$cellType=Idents(immune)
pdf(file = "4.immune第二次注释后.pdf",width = 15,height = 8)
do_DotPlot(sample = immune,features = genes,dot.scale = 12,
           legend.framewidth = 2, font.size =10)
dev.off()
#########人工注释后结果可视化
# 可视化UMAP/tSNE
pdf(file = "5.immune第二次注释UMAP.pdf",width =7.5,height = 5.5)
DimPlot(immune,group.by = "cellType", reduction = "umap", label = T,label.box = T, label.size = 3.5,pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "6.immune第二次注释TSEN.pdf",width =7.5,height = 5.5)
DimPlot(immune,group.by = "cellType", reduction = "tsne", label = T,label.box = T, label.size = 3.5,pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
#基因分布图
seuratgenes=c("CCL5","GZMB")      #写入你想展示的基因
pdf(file = "7.基因分布UMAP.pdf",width =12,height = 6)
FeaturePlot(immune, features = seuratgenes, reduction ="umap",  cols = c("grey", "red"),min.cutoff = 0, max.cutoff = 1,ncol=2,pt.size = 0.01)
dev.off()
pdf(file = "8.基因分布TSNE.pdf",width =12,height = 6)
FeaturePlot(immune, features = seuratgenes, reduction ="tsne",  cols = c("grey", "red"),min.cutoff = 0, max.cutoff = 1,ncol=2,pt.size = 0.01)
dev.off()
#高亮度基因分布
pdf(file = "9.高亮基因分布TSNE.pdf",width =12,height = 6)
ggrastr::rasterize(Nebulosa::plot_density(immune, seuratgenes,size=0.01,reduction = "tsne",shape = 2),dpi=300)
dev.off()
pdf(file = "10.高亮基因分布UMAP.pdf",width =12,height = 6)
ggrastr::rasterize(Nebulosa::plot_density(immune, seuratgenes,size=0.01,reduction = "umap",shape = 2),dpi=300)
dev.off()
qsave(immune,"手动注释.qs",nthreads = detectCores())

####################################################
####################################################
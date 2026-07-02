# Course-derived reference script
# English filename: 18_scina_annotation.R
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
library(SCINA)
scRNA=qread("scRNA分群聚类后.qs",nthreads = detectCores())#读取数据
#构建marker基因list
markerlist<-list("T Cells"= c("CD3D", "CD3E", "CD8A","PTPRC","CD4","CD2"),
                 "B cells "=c("CD19", "CD79A", "MS4A1"),
                 "Plasma cells"= c("IGHG1", "MZB1", "SDC1", "CD79A"),
                 "Monocytes"=c("CD68", "CD163", "CD14"),
                 "NK Cells"= c("GNLY","NKG7","KLRD1"),
                 "Fibroblasts"= c("FGF7", "MME","GSN","LUM","DCN"),
                 "Endothelial cells"= c("PECAM1", "VWF"),
                 "Mast"=c("CPA3", "CST3", "KIT", "TPSAB1", "TPSB2", "MS4A2"))
#获取单细胞表达矩阵
scRNAcounts<-scRNA@assays[["RNA"]]@counts
scRNAout=SCINA(scRNAcounts,markerlist,rm_overlap=TRUE,allow_unknown=F,log_file='SCINA.log')
#结果可视化
seurat<-scRNA
newLabels=scRNAout[["cell_labels"]]
seuratidens=mapvalues(Idents(seurat), from = Idents(seurat), to = newLabels)
Idents(seurat)=seuratidens
seurat$cellType=Idents(seurat)
# 可视化UMAP/tSNE
pdf(file = "1.SCINA注释后UMAP.pdf",width =7.5,height = 5.5)
DimPlot(seurat,group.by = "cellType" ,reduction = "umap", label = T,label.box = T, label.size = 3.5,pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "2.SCINA注释后TSEN.pdf",width =7.5,height = 5.5)
DimPlot(seurat,group.by = "cellType" , reduction = "tsne", label = T,label.box = T, label.size = 3.5,pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
#基因分布图
seuratgenes=c("CCL5","GZMB")      #写入你想展示的基因
pdf(file = "3.基因分布UMAP.pdf",width =12,height = 6)
FeaturePlot(seurat, features = seuratgenes, reduction ="umap",  cols = c("grey", "red"),min.cutoff = 0, max.cutoff = 1,ncol=2,pt.size = 0.01)
dev.off()
pdf(file = "4.基因分布TSNE.pdf",width =12,height = 6)
FeaturePlot(seurat, features = seuratgenes, reduction ="tsne",  cols = c("grey", "red"),min.cutoff = 0, max.cutoff = 1,ncol=2,pt.size = 0.01)
dev.off()
#高亮度基因分布
pdf(file = "5.高亮基因分布TSNE.pdf",width =12,height = 6)
ggrastr::rasterize(Nebulosa::plot_density(seurat, seuratgenes,size=0.01,reduction = "tsne",shape = 2),dpi=300)
dev.off()
pdf(file = "6.高亮基因分布UMAP.pdf",width =12,height = 6)
ggrastr::rasterize(Nebulosa::plot_density(seurat, seuratgenes,size=0.01,reduction = "umap",shape = 2),dpi=300)
dev.off()
qsave(seurat,"SCINA注释.qs",nthreads = detectCores())
####################################################
####################################################
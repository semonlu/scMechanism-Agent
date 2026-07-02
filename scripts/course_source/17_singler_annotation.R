# Course-derived reference script
# English filename: 17_singler_annotation.R
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
library(SingleR)
library(qs)
seurat=qread("scRNA分群聚类后.qs",nthreads = detectCores())#读取数据
refdata=celldex::HumanPrimaryCellAtlasData()#获取人类示例
# 使用的数据为标化后的数据
testdata <- GetAssayData(seurat)
#进行自动注释
cellpred <- SingleR(test = testdata,  
                    ref = refdata, 
                    labels = refdata$label.main)
#保存注释结果表格
celltype = data.frame(ClusterID = rownames(cellpred), 
                      celltype = cellpred$labels, 
                      stringsAsFactors = F)
write.csv(celltype, "1.SingleR细胞注释.csv")

pdf(file = "2.singleR热图.pdf",width =7.5,height = 5.5)
plotScoreHeatmap(cellpred, clusters = rownames(cellpred), order.by = "cluster")
dev.off()
#sigleR注释后结果可视化
newLabels=cellpred$labels
seuratidens=mapvalues(Idents(seurat), from = Idents(seurat), to = newLabels)
Idents(seurat)=seuratidens
seurat$cellType=Idents(seurat)
# 可视化UMAP/tSNE
pdf(file = "3.singleR注释后UMAP.pdf",width =7.5,height = 5.5)
DimPlot(seurat,group.by = "cellType" ,reduction = "umap", label = T,label.box = T, label.size = 3.5,pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "4.singleR注释后TSEN.pdf",width =7.5,height = 5.5)
DimPlot(seurat,group.by = "cellType" , reduction = "tsne", label = T,label.box = T, label.size = 3.5,pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
#基因分布图
seuratgenes=c("CCL5","GZMB")      #写入你想展示的基因
pdf(file = "5.基因分布UMAP.pdf",width =12,height = 6)
FeaturePlot(seurat, features = seuratgenes, reduction ="umap",  cols = c("grey", "red"),min.cutoff = 0, max.cutoff = 1,ncol=2,pt.size = 0.01)
dev.off()
pdf(file = "6.基因分布TSNE.pdf",width =12,height = 6)
FeaturePlot(seurat, features = seuratgenes, reduction ="tsne",  cols = c("grey", "red"),min.cutoff = 0, max.cutoff = 1,ncol=2,pt.size = 0.01)
dev.off()
#高亮度基因分布
pdf(file = "7.高亮基因分布TSNE.pdf",width =12,height = 6)
ggrastr::rasterize(Nebulosa::plot_density(seurat, seuratgenes,size=0.01,reduction = "tsne",shape = 2),dpi=300)
dev.off()
pdf(file = "8.高亮基因分布UMAP.pdf",width =12,height = 6)
ggrastr::rasterize(Nebulosa::plot_density(seurat, seuratgenes,size=0.01,reduction = "umap",shape = 2),dpi=300)
dev.off()
qsave(scRNA,"SingleR注释.qs",nthreads = detectCores())
####################################################
####################################################
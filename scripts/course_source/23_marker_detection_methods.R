# Course-derived reference script
# English filename: 23_marker_detection_methods.R
# Original course path: see source_manifest.csv
# Role: marker reference
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
library(presto)
library(COSG)
library(starTracer)

scRNA=qread("deepseek注释.qs",nthreads = detectCores())#读取数据
#FindAllMarkers算法
scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)  
write.csv(scRNA.markers,file = "1.FindAllMarkers每个细胞Maeker基因.csv")
top5scRNA.markers <- scRNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
col <- c(ggsci::pal_npg()(9),ggsci::pal_jco()(9),ggsci::pal_jama()(7),ggsci::pal_nejm()(8))
pdf(file = "1.FindAllMarkers每个细胞前5基因热图.pdf",width =22,height = 16)
DoHeatmap(scRNA,features = top5scRNA.markers$gene,
          group.colors = col) +
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
dev.off()

#presto-wilcoxauc算法
scRNA.markers <- wilcoxauc(scRNA,group_by = "cellType")  
write.csv(scRNA.markers,file = "2.presto-wilcoxauc每个细胞Maeker基因.csv")
top5scRNA.markers <- scRNA.markers %>%
  group_by(group) %>%
  top_n(n = 5, wt = logFC)
col <- c(ggsci::pal_npg()(9),ggsci::pal_jco()(9),ggsci::pal_jama()(7),ggsci::pal_nejm()(8))
pdf(file = "2.presto-wilcoxauc每个细胞前5基因热图.pdf",width =22,height = 16)
DoHeatmap(scRNA,features = top5scRNA.markers$feature,
          group.colors = col) +
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
dev.off()

#COSG算法
scRNA.markers<-cosg(
  scRNA,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  remove_lowly_expressed=TRUE,
  expressed_pct=0.1,
  n_genes_user=5
)

write.csv(scRNA.markers,file = "3.COSG每个细胞Maeker基因.csv")
top5scRNA.markers <- lapply(1:length(scRNA.markers[["names"]]), function(i){
  genes<-data.frame(gene=scRNA.markers[["names"]][[i]])
  return(genes)
})
top5scRNA.markers<-do.call(rbind,top5scRNA.markers)
col <- c(ggsci::pal_npg()(9),ggsci::pal_jco()(9),ggsci::pal_jama()(7),ggsci::pal_nejm()(8))
pdf(file = "3.COSG每个细胞前5基因热图.pdf",width =22,height = 16)
DoHeatmap(scRNA,features = top5scRNA.markers$gene,
          group.colors = col) +
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
dev.off()

#starTracer算法
starTracer<-searchMarker(
  x=scRNA,
  thresh.1=0.5,
  thresh.2=0.3,
  method="pos",
  num=5,
  gene.use=NULL,
  meta.data=NULL,
  ident.use=NULL
)
top5scRNA.markers<-starTracer$para_frame[starTracer$genes.markers,c("max.X","gene")]
write.csv(top5scRNA.markers,file = "4.starTracer每个细胞Maeker基因.csv")
col <- c(ggsci::pal_npg()(9),ggsci::pal_jco()(9),ggsci::pal_jama()(7),ggsci::pal_nejm()(8))
pdf(file = "4.starTracer每个细胞前5基因热图.pdf",width =22,height = 16)
DoHeatmap(scRNA,features = top5scRNA.markers$gene,
          group.colors = col) +
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
dev.off()
####################################################
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
####################################################
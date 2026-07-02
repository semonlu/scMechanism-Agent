# Course-derived reference script
# English filename: 07_read_multiple_10x_samples.R
# Original course path: see source_manifest.csv
# Role: 10x import reference
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
library(dplyr)
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
library(parallel)
library(doParallel)
library(data.table)
library(qs)
#数据读取##
dataid<-c("GSE178858",
          "GSE184198")
scRNAs<-lapply(dataid, function(ids){
data_dir <- paste0(getwd(),"/",ids)       
samples=list.files(data_dir)#读取文件夹下所有文件
dataGSE<-gsub("(.*?)\\_.*", "\\1", samples)#读取GSM样本号
#对每个GSM样本构建一个文件夹
lapply(levels(factor(dataGSE)), function(idGSM){
  dir.create(paste0(data_dir,"/",idGSM))
})
#自动识别并重命名每个GSM样本的10X数据
lapply(levels(factor(dataGSE)), function(idGSE){
  #将每个GSM样品移入对应文件夹
  files0 <- list.files(path =paste0(data_dir,"/") ,pattern = paste0(idGSE,".*barcodes.*"), full.names = TRUE)
  file.rename(from = files0,to =paste0(data_dir,"/",idGSE,"/barcodes.tsv.gz"))
  files1 <- list.files(path =paste0(data_dir,"/") ,pattern = paste0(idGSE,".*features.*"), full.names = TRUE)
  file.rename(from = files1,to =paste0(data_dir,"/",idGSE,"/features.tsv.gz"))
  files2 <- list.files(path =paste0(data_dir,"/") ,pattern = paste0(idGSE,".*matrix.*"), full.names = TRUE)
  file.rename(from = files2,to =paste0(data_dir,"/",idGSE,"/matrix.mtx.gz"))
})
samples=list.files(data_dir)#重新读取文件夹下所有文件
dir=file.path(data_dir,samples)
scRNAdata <- Read10X(data.dir = dir)
scRNA <- CreateSeuratObject(counts = scRNAdata,
                            min.cells = 3,
                            min.features = 200)
scRNAidens=mapvalues(Idents(scRNA), from = levels(Idents(scRNA)), to = samples)
Idents(scRNA)=scRNAidens
scRNA$Type=Idents(scRNA)
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP")
return(scRNA)
})
scRNA<-Reduce(function(x, y) merge(x, y, all = TRUE), scRNAs)
qsave(scRNA, "scRNA.qs", nthreads = detectCores())
####################################################
####################################################
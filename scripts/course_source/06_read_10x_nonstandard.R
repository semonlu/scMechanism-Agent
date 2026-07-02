# Course-derived reference script
# English filename: 06_read_10x_nonstandard.R
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
library(R.utils)
#数据读取##
data_dir <- paste0(getwd(),"/GSE176078")    
samples=list.files(data_dir)#读取文件夹下所有文件
lapply(samples, function(id){
gunzip(paste0(data_dir,"/",id))
})
samples=list.files(data_dir)#重新读取文件夹下所有文件
lapply(samples, function(id){
untar(paste0(data_dir,"/",id),exdir =paste0(data_dir,"/"))
})
samples=list.files(data_dir)#读取文件夹下所有文件夹
scRNAdatas<-lapply(samples, function(dirs){
scdata_dir <- paste0(data_dir,"/",dirs)       
scRNAdata<-readMM(paste0(scdata_dir,"/count_matrix_sparse.mtx"))#表达数据
metadata=fread(paste0(scdata_dir,"/count_matrix_barcodes.tsv"),header = F)#细胞标签
colnames(scRNAdata)<-metadata$V1
gene<-fread(paste0(scdata_dir,"/count_matrix_genes.tsv"),header = F)#基因ID
rownames(scRNAdata)<-gene$V2
return(scRNAdata)
})
scRNAdata<-do.call(cbind,scRNAdatas)
scRNA <- CreateSeuratObject(counts = scRNAdata, 
                         min.cells = 3,
                         min.features = 200)
scRNAidens=mapvalues(Idents(scRNA), from = levels(Idents(scRNA)), to = samples)
Idents(scRNA)=scRNAidens
scRNA$Type=Idents(scRNA)
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP")
qsave(scRNA, "scRNA.qs", nthreads = detectCores())
####################################################
####################################################
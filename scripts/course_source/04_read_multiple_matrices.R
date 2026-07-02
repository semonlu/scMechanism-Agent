# Course-derived reference script
# English filename: 04_read_multiple_matrices.R
# Original course path: see source_manifest.csv
# Role: input import reference
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
library(limma)
library(dplyr)
library(patchwork)
library(ggplot2)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(parallel)
library(doParallel)
library(data.table)
library(qs)
setDTthreads(threads =detectCores())
dataid<-c("GSE118389_counts_rsem.txt",
          "GSE158631_count.csv")
scRNAs = lapply(dataid,function(id){ 
  data=fread(id,check.names = F,data.table = T)
  dataGSE<-gsub("(.*?)\\_.*", "\\1", id)
  rownames(data)<-data$V1
  data<-data[,2:ncol(data)]
  scRNA=CreateSeuratObject(counts =data, 
                           project = dataGSE, 
                           min.cells = 3,
                           min.features = 200)
  scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")#线粒体基因
  scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP")#核糖体基因
  scRNA$Type=Idents(scRNA)
  return(scRNA)
})
scRNA<-Reduce(function(x, y) merge(x, y, all = TRUE), scRNAs)
qsave(scRNA, "scRNA.qs", nthreads = detectCores())
####################################################
####################################################

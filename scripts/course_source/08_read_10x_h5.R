# Course-derived reference script
# English filename: 08_read_10x_h5.R
# Original course path: see source_manifest.csv
# Role: h5 import reference
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
library(tidyverse)
library(Matrix)
library(stringr)
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
library(hdf5r)   
library(parallel)
library(doParallel)
library(data.table)
library(qs)
ids <-list.files(path =getwd() ,pattern = ".h5.*", full.names = F)
scRNAs<-lapply(ids, function(h5id){
data1 = Read10X_h5(h5id, use.names = T)
scRNAcol=gsub("(.*?)\\_.*", "\\1",h5id)
scRNA<-CreateSeuratObject(counts =data1, 
                       project = scRNAcol, 
                       min.cells = 3,
                       min.features = 200)
scRNA$Type=Idents(scRNA)
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP")
return(scRNA)
})
scRNA<-Reduce(function(x, y) merge(x, y, all = TRUE), scRNAs)
qsave(scRNA, "scRNA.qs", nthreads = detectCores())
####################################################
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
####################################################
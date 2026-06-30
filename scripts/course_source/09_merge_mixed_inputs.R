# Course-derived reference script
# English filename: 09_merge_mixed_inputs.R
# Original course path: see source_manifest.csv
# Role: merge reference
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
library(limma)
library(dplyr)
library(patchwork)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(Matrix)
library(stringr)
library(ggplot2)
library(stringr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(hdf5r)
library(harmony)
library(plyr)
library(parallel)
library(doParallel)
library(data.table)
library(qs)
ids <-list.files(path =getwd() ,pattern = ".qs.*", full.names = F)
scRNAs<-lapply(ids, function(scid){
  scRNA=qread(scid, nthreads = detectCores())
  return(scRNA)
})
scRNA<-Reduce(function(x, y) merge(x,y),scRNAs)
qsave(scRNA, "scRNA.qs", nthreads = detectCores())
####################################################
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
####################################################
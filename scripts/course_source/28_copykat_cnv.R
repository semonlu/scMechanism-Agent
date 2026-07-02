# Course-derived reference script
# English filename: 28_copykat_cnv.R
# Original course path: see source_manifest.csv
# Role: CNV reference
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
library(copykat)
library(parallel)
library(doParallel)
library(data.table)
library(qs)
scRNA=qread("deepseek注释.qs",nthreads = detectCores())#读取数据

##随机抽取2%个细胞，仅作为演示，实际操作或高性能工作站直接运行scRNA_tpm<-scRNA
#5000个细胞预计8小时，4000个细胞预计5小时
N=length(colnames(scRNA))/50
N=round(N)
scRNA_tpm<-scRNA[,sample(x=colnames(scRNA),size = N,replace=F)]
#提取细胞表达量
expr <- as.matrix(scRNA_tpm@assays$RNA@counts)
expr<-na.omit(expr)
copykat.test <- copykat(rawmat=expr, 
                        id.type="S", 
                        cell.line="no", 
                        ngene.chr=3, 
                        LOW.DR = 0.01,
                        win.size=25, 
                        KS.cut=0.01, 
                        sam.name="st", 
                        distance="euclidean", 
                        n.cores=1)
qsave(copykat.test,"copykat.qs")

#记录坐标信息
pre <-copykat.test$prediction
pre <-as.data.frame(pre)
rownames(pre)<-pre$cell.names
pre <-pre[rownames(scRNA@meta.data),]
scRNA$copykat.pred <-pre$copykat.pred
scRNA$copykat.tumor.pred<-
pdf(file = "copykat_tsne图.pdf",width =8.5,height = 7.5)
DimPlot(scRNA,group.by = "copykat.pred",reduction = "tsne",pt.size = 0.01)
dev.off()
####################################################
####################################################

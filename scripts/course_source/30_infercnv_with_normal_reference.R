# Course-derived reference script
# English filename: 30_infercnv_with_normal_reference.R
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
library(phylogram)
library(gridExtra)
library(grid)
require(dendextend)
library(ape)
require(ggthemes)
library(tidyverse)
library(infercnv)
library(miscTools)
library(AnnoProbe)
library(remotes)
library(infercnv)
library(future)
library(future.apply)
library(DoubletFinder)
library(parallel)
library(doParallel)
library(data.table)
library(qs)
scRNA=qread("手动注释.qs",nthreads = detectCores())#读取数据
normal<-qread("正常细胞.qs",nthreads = detectCores())
#合并正常
scRNA<-merge(scRNA,normal)
#合并后常规流程
##标准化,使用LogNormalize方法
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
## 鉴定高变基因
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
# 选择全部基因归一化
scRNAall.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = scRNAall.genes)
# PCA降维
scRNA <- Seurat::RunPCA(scRNA, features = VariableFeatures(object = scRNA))
scRNA <- Seurat::RunTSNE(scRNA,dims = 1:20)
scRNA <- Seurat::RunUMAP(scRNA,dims = 1:20)
##随机抽取5%个细胞，仅作为演示，实际操作或高性能工作站直接运行
N=length(colnames(scRNA))/20
N=round(N)
scRNA<-scRNA[,sample(x=colnames(scRNA),size = N,replace=F)]
#记录细胞信息
scRNA$cellType=scRNA@active.ident
table(scRNA$cellType)
yourcell=names(table(scRNA$cellType))
cells.use=colnames(scRNA)[which(scRNA$cellType%in% yourcell)]
#构建细胞表达表格
dat=as.data.frame(GetAssayData(subset(scRNA, cells=cells.use)))
groupinfo=data.frame(v1=colnames(dat),
                     v2=scRNA@active.ident[cells.use])
#构建染色体及分组信息
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
dat=dat[rownames(dat) %in% geneInfor[,1],] ##去除性染色体
dat=dat[match(geneInfor[,1], rownames(dat) ),] #给染色体排序
write.table(groupinfo,file = 'groupFiles.txt',sep = '\t',quote = F,col.names = F,row.names = F)

dat=na.omit(dat)
write.table(dat,file ='expFile.txt',sep = '\t',quote = F)
write.table(geneInfor,file = 'geneFile.txt',sep = '\t',quote = F,col.names = F,row.names = F)
#构建Infer对象
infercnv_obj = CreateInfercnvObject(delim = '\t',
                                    raw_counts_matrix = 'expFile.txt',
                                    annotations_file = 'groupFiles.txt',
                                    gene_order_file = 'geneFile.txt',
                                    ref_group_names = "Normal")
#ref_group_names根据细胞注释进行填写，将哪些细胞定义为参考，即输入的细胞应是作为一个参考值，其余细胞根据该参考值进行判断是否为恶性
#简而言之，如果分群中有正常细胞，则填入。如果无正常细胞，则填NULL
#如果填写NULL则会选用所有细胞的平均表达作为参照，这样就需要确保细胞中有足够的差异

#运行InferCNV
#10x数据cutoff推荐使用0.1
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, 
                             out_dir='inferCNV/', 
                             cluster_by_groups=T, #先区分细胞源再进行聚类
                             denoise=TRUE,
                             plot_steps =T,
                             analysis_mode = "subclusters",#默认是“samples"
                             HMM=F,#如果设置HMM会运行得比较久
                             num_threads=6)#设置线程数

##绘制矢量图
plot_cnv(infercnv_obj,
         obs_title="Observations (Cells)",
         ref_title="References (Cells)",
         cluster_by_groups=TRUE,
         plot_chr_scale=T,#绘制染色体全长
         x.center=1,
         x.range="auto",
         hclust_method='ward.D',
         custom_color_pal = color.palette(c("#0071B2", "white", "#C3250A"), c(2, 2)),
         color_safe_pal=FALSE,
         output_filename="infercnv",
         output_format="pdf",
         png_res=600,
         dynamic_resize=0)
####################################################
####################################################
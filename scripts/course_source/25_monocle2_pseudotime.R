# Course-derived reference script
# English filename: 25_monocle2_pseudotime.R
# Original course path: see source_manifest.csv
# Role: trajectory reference
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
#生成随机颜色
randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}

# 生成100个随机颜色
randomColors <- replicate(100,randomColor())
scRNA=qread("deepseek注释.qs",nthreads = detectCores())#读取数据
###拟时序分析
##随机抽取1%个细胞，仅作为演示，实际操作或高性能工作站直接运行scRNA_tpm<-scRNA
N=length(colnames(scRNA))/100
N=round(N)
scRNA_tpm<-scRNA[,sample(x=colnames(scRNA),size = N,replace=F)]
#计算marker基因
scRNA_tpm.markers <- FindAllMarkers(scRNA_tpm, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
logFCfilter=1           
adjPvalFilter=0.05
#按条件筛选marker基因
scRNA_tpm.markers=scRNA_tpm.markers[(abs(as.numeric(as.vector(scRNA_tpm.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(scRNA_tpm.markers$p_val_adj))<adjPvalFilter),]
#保存为表格
monocle.matrix=as.matrix(scRNA_tpm@assays$RNA@counts, 'sparseMatrix')
monocle.sample=scRNA_tpm@meta.data[,8,drop=F]
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
marker=scRNA_tpm.markers

#将Seurat结果转换为monocle需要的细胞矩阵，细胞注释表和基因注释表表
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
#添加细胞聚类数据
pData(cds)$Type=scRNA_tpm$Type
pData(cds)$cellType=scRNA_tpm$cellType
#伪时间分析流程
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
#monocle选择高变基因
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')
cds <- orderCells(cds)
cds <- setOrderingFilter(cds, marker$gene)
pdf(file="1.细胞拟时序.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cellType")
dev.off()
pdf(file="2.模拟发育时间.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "Pseudotime")
dev.off()
#细胞轨迹差异分析
groups=subset(pData(cds),select='State')
pbmc=AddMetaData(object=scRNA, metadata=groups, col.name="group")
geneList=list()
for(i in levels(factor(groups$State))){
  pbmc.markers=FindMarkers(pbmc, ident.1 = i, group.by = 'group')
  sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
  sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
  write.table(sig.markers,file=paste0("05.monocleDiff.", i, ".txt"),sep="\t",row.names=F,quote=F)
  geneList[[i]]=row.names(sig.markers)
}
#保存交集基因
unionGenes=Reduce(union,geneList)
write.table(file="3.细胞轨迹差异基因.txt",unionGenes,sep="\t",quote=F,col.names=F,row.names=F)
#特定基因时序
pdf(file = "4.特定基因拟时序.pdf",width =8,height = 7)
pData(cds)[,'RUVBL1'] = scRNA_tpm@assays$RNA@scale.data['RUVBL1',]
plot_cell_trajectory(cds, color_by = 'RUVBL1') + scale_color_gsea()
dev.off()

####################################################
####################################################
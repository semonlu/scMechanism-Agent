# Course-derived reference script
# English filename: 34_cibersort_deconvolution.R
# Original course path: see source_manifest.csv
# Role: deconvolution reference
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
library(tcltk)
library(limma)
library(reshape2)
library(ggpubr)
library(vioplot)
library(ggExtra)
library(data.table)
library(tidyverse)
library(doParallel)
setDTthreads(threads =detectCores())
scRNA=qread("手动注释.qs",nthreads = detectCores())#读取数据
#构建伪表达矩阵
exps<-scRNA@assays[["RNA"]]@counts
#构建细胞参考矩阵
cellexp<-lapply(levels(scRNA$cellType),function(cells){
  identifiers <- colnames(subset(scRNA,subset =cellType %in% cells))
  exp<- exps[, identifiers]
  avexp <- apply(exp, 1, mean)
  avexp<-as.data.frame(avexp)
  colnames(avexp)<-cells
  return(avexp)
})
cellexps<-do.call(cbind,cellexp)
cellexps<-as.data.frame(cellexps)
#数据过小则进行放大
normalized_x <- (cellexps - min(cellexps)) / (max(cellexps) - min(cellexps))
normalized_scaled_x <- normalized_x * 1000000
normalized_scaled_x<-round(normalized_scaled_x)
normalized_scaled_x<-normalized_scaled_x[normalized_scaled_x>0,]
normalized_scaled_x<-na.omit(normalized_scaled_x)
rownames(normalized_scaled_x)<-make.names(rownames(normalized_scaled_x))
LMgenes<-fread("LMgenes.txt",sep = "\t",header = T)
sameLMgenes<-intersect(LMgenes$ID,rownames(normalized_scaled_x))
normalized_scaled_x<-normalized_scaled_x[sameLMgenes,]
#导出为参考矩阵
fwrite(normalized_scaled_x,file = "scimmune.txt",quote = F,sep = "\t",row.names = T)


#读取TCGA表达文件
data=fread("TCGAexp.txt",header = T,check.names = F,data.table = F)
rownames(data)=data$V1
data=data[,2:ncol(data)] 
Tumordata=data%>% dplyr::select(str_which(colnames(.), "-01A")) 
colnames(Tumordata )=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(Tumordata))
out=Tumordata
out=rbind(ID=colnames(out), out)
fwrite(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F,row.names = T)

#运行主代码，得到免疫细胞浸润的结果
source("xxdimmune.R",encoding ="utf-8" )
results=CIBERSORT("scimmune.txt", "uniq.symbol.txt", perm=2, QN=F)
immune=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<0.05,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])
data=t(immune)
col=rainbow(nrow(data),s=0.7,v=0.7)
pdf("1.反卷积样本总览.pdf",height=10,width=20)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=60,xpd=T);text(a1,-0.03,colnames(data),adj=1,cex=0.18);par(srt=0)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()

gene="CLEC4M"
sameSample=intersect(rownames(immune), colnames(Tumordata))
immune=immune[sameSample,]
Tumordata=Tumordata[,sameSample]
rt=data.frame(immune,gene=as.numeric(Tumordata[gene,]))
rt$gene=ifelse(rt$gene>=median(rt$gene),"High","Low")
data=rt[,-(ncol(rt)-1)]
data=reshape2::melt(data,id.vars=c("gene"))
colnames(data)=c("gene", "Immune", "Expression")
#绘制箱线图
group=levels(factor(data$gene))
data$gene=factor(data$gene, levels=c("Low","High"))
bioCol=c("#223D6C","#FFD121")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="gene",
                  xlab="",
                  ylab="Fraction",
                  legend.title=gene,
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(50)+
  stat_compare_means(aes(group=gene),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#输出图片
pdf(file="2.反卷积箱线图.pdf", width=7, height=6)
print(boxplot)
dev.off()
####################################################
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
#######欢迎关注《叉叉滴同学的生信笔记》#########
####################################################

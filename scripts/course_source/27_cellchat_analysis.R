# Course-derived reference script
# English filename: 27_cellchat_analysis.R
# Original course path: see source_manifest.csv
# Role: cell communication reference
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
library(NMF)
library(presto)
library(ggalluvial)
scRNA=qread("deepseek注释.qs",nthreads = detectCores())#读取数据
###拟时序分析
##随机抽取5%个细胞，仅作为演示，实际操作或高性能工作站直接运行scRNA_tpm<-scRNA
N=length(colnames(scRNA))/20
N=round(N)
scRNA<-scRNA[,sample(x=colnames(scRNA),size = N,replace=F)]
#获取表达矩阵
cellchat <- createCellChat(scRNA@assays$RNA@data)
#获取注释细胞标签
meta <- data.frame(cellType = scRNA$cellType, row.names =  Cells(scRNA))
#把metadata信息加到CellChat对象中，添加细胞标签
cellchat <- addMeta(cellchat, meta = meta, meta.name = "cellType")
#把细胞标签设置成默认的ID
cellchat <- setIdent(cellchat, ident.use = "cellType") 
#统计每个细胞亚群中的细胞数目
groupSize <- as.numeric(table(cellchat@idents)) 

#CellChat提供了人和小鼠的配受体数据库，分别可以用CellChatDB.human,CellChatDB.mouse来导入
CellChatDB <- CellChatDB.human


#择特定的信息描述细胞间的相互作用，比用一个大的配体库更精细
#用Secreted Signaling来分析细胞通信
#查看有哪些可以选择的子数据库
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")


#将使用的数据库信息写入到cellchat对象中
cellchat@DB <- CellChatDB.use

#抽取信号通路相关基因的表达矩阵
cellchat <- subsetData(cellchat) 


#对表达数据进行预处理，用于细胞间的通信分析。
#首先在一个细胞组中识别过表达的配体或受体，然后将基因表达数据投射到蛋白-蛋白相互作用(PPI)网络上。
#如果配体或受体过表达，则识别过表达配体和受体之间的相互作用。
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)

#相互作用推断
#为每个相互作用分配一个概率值并进行置换检验来推断生物意义上的细胞-细胞通信
cellchat <- computeCommunProb(cellchat)

#通过计算与每个信号通路相关的所有配体-受体相互作用的通信概率来推断信号通路水平上的通信概率
#通过计算链路的数量或汇总通信概率来计算细胞间的聚合通信网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

#可视化聚合的通讯网络
groupSize <- as.numeric(table(cellchat@idents))
pdf("1.细胞通讯网络_细胞交互次数.pdf")
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf("2.细胞通讯网络_细胞交互权重.pdf")
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

#比较不同细胞通讯网络网络之间的权重
mat <- cellchat@net$weight
pdf("3.不同细胞通讯网络权重.pdf")
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
#显示重要通信的信号路径
cellchat@netP$pathways
levels(cellchat@idents) 
#在层次绘图的时候，第一列显示的细胞类型的数目，一般有几个细胞就选几个
vertex.receiver = seq(1,5) 
#显示的信号通路
pathways.show <- "CypA"
pdf("4.重要通路通讯_贝壳图.pdf")
par(mfrow=c(1,1),xpd=TRUE)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
pdf("5.重要通路通讯_弦图.pdf")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()
#绘制通讯热图
pdf("6.通讯热图.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()

#绘制配体受体对整个信号通路的贡献度
pdf("7.配体受体贡献度.pdf")
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

#可视化单个配体受体通讯网络
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] 
vertex.receiver = seq(1,5) 
pdf("8.单个配体受体通讯网络_贝壳图.pdf")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
dev.off()

pdf("9.单个配体受体通讯网络_弦图.pdf")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()

#可视化由多种配体-受体或信号通路介导的细胞间通讯
levels(cellchat@idents)
pdf("10.多配体受体介导气泡图.pdf")
netVisual_bubble(cellchat, sources.use = 5, remove.isolate = FALSE)
dev.off()
pdf("11.多配体受体介导弦图.pdf")
netVisual_chord_gene(cellchat, sources.use = 5, lab.cex = 0.5,legend.pos.y = 30)
dev.off()

#识别细胞群体中的信号作用（例如，显性发送者、接收者）以及主要贡献的信号
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf("12.中心性得分.pdf")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
dev.off()

pdf("13.发送者和接收者散点图.pdf")
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

#识别对某些细胞群体出向或入向信号传递贡献最大的信号
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
pdf("14.信号通路在细胞通讯中的作用.pdf",width = 12,height = 6)
ht1 + ht2
dev.off()

#识别并可视化分泌细胞的输出通信模式
#推断模式数量
pdf("15.selectK推断_输出.pdf")
selectK(cellchat, pattern = "outgoing")
dev.off()
#输出模式数量为4时曲线骤降
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
pdf("16.通讯模式热图_输出.pdf")
identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.off()
pdf("17.通讯模式桑基图_输出.pdf")
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
pdf("18.通讯模式点图_输出.pdf")
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

#识别并可视化分泌细胞的输入通信模式
pdf("19.selectK推断_输入.pdf")
selectK(cellchat, pattern = "incoming")
dev.off()
#输出模式数量为4时曲线骤降
nPatterns = 7
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
pdf("20.通讯模式热图_输入.pdf")
identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
dev.off()
pdf("21.通讯模式桑基图_输入.pdf")
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()
pdf("22.通讯模式点图_输入.pdf")
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()

#信号网络的多维和分类学习分析
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional",umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "functional")
pdf("23.二维分类图.pdf")
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
dev.off()

#根据相似性识别信号组
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural",umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "structural")
pdf("24.二维信号组图.pdf")
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
dev.off()
pdf("25.分组二维信号图.pdf")
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()
qsave(cellchat, file = "cellchat.qs",nthreads = detectCores() )
####################################################
####################################################

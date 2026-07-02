# Course-derived reference script
# English filename: 24_go_kegg_enrichment.R
# Original course path: see source_manifest.csv
# Role: enrichment reference
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
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(DOSE)
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

#挑选每个细胞亚群中高表达的前200个基因
sample.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top200 <- sample.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC) 
top200_gene<-unstack(top200,gene~cluster)
#将gene symbol转换成entriz gene id
#做KEGG富集分析需要
top200_entrez <- lapply(X = top200_gene, FUN = function(x) {
  x <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  x <- x[,-1]
})

######################################
#对所有亚群做GO富集分析 Gene ontoloty
#####################################
######################################
# CC cellular component
######################################
all_cc <- compareCluster(top200_entrez, 
                         fun='enrichGO',
                         ont= 'CC',
                         OrgDb='org.Hs.eg.db')
write.csv(as.data.frame(all_cc),"1.所有亚群GO分析_cc.csv")

#绘制CC富集图
pdf('2.所有亚群GO分析_CC.pdf', width = 15, height = 6)
dotplot(all_cc, includeAll=FALSE,label_format=100)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

##################################################
# MF molecular function
#################################################
all_mf <- compareCluster(top200_entrez, 
                         fun='enrichGO',
                         ont= 'MF',
                         OrgDb='org.Hs.eg.db')
write.csv(as.data.frame(all_mf),"3.所有亚群GO分析_mf.csv")

#绘制MF富集图
pdf('4.所有亚群GO分析_MF.pdf', width = 15, height = 6)
dotplot(all_mf, includeAll=FALSE,label_format=100)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()
#############################################
# BP biological process
#############################################
all_bp <- compareCluster(top200_entrez, 
                         fun='enrichGO',
                         ont= 'BP',
                         OrgDb='org.Hs.eg.db')
write.csv(as.data.frame(all_bp),"5.所有亚群GO分析_bp.csv")

#绘制BP富集图
pdf('6.所有亚群GO分析_BP.pdf', width = 15, height = 6)
dotplot(all_bp, includeAll=FALSE,label_format=100)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()


#批量对每个细胞进行GO富集分析
lapply(1:length(top200_entrez),function(i){
all_go_cell  <- enrichGO(gene = top200_entrez[[i]],  
                      OrgDb = 'org.Hs.eg.db',
                      ont = 'ALL',
                      pvalueCutoff = 0.05,
                      qvalueCutoff =0.2)

all_go_cell<-setReadable(all_go_cell, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(as.data.frame(all_go_cell),paste0("7.",names(top200_entrez[i]),"GO分析.csv"))
plotGO<-dotplot(all_go_cell, split="ONTOLOGY",label_format=100) + facet_grid(ONTOLOGY~., scale="free")
pdf(paste0('8.',names(top200_entrez[i]),'GO分析.pdf'), width = 16, height = 8)
print(plotGO)
dev.off()
})
###################################
#对所有亚群做KEGG富集分析
###################################
all_kegg <- compareCluster(top200_entrez,
                           fun='enrichKEGG',
                           pvalueCutoff=0.05, 
                           organism="hsa"
)

write.csv(as.data.frame(all_kegg),"9.所有亚群_KEGG.csv")
pdf('10.所有亚群_KEGG.pdf', width = 15, height = 8)
dotplot(all_kegg, showCategory=5, includeAll=FALSE,label_format=100)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

#批量对每个细胞进行KEGG富集分析
lapply(1:length(top200_entrez),function(i){
kegg_cell <- enrichKEGG(gene = top200_entrez[[i]],
                     organism  = 'hsa', 
                     pvalueCutoff = 1,
                     qvalueCutoff =1)
kegg_cell<-setReadable(kegg_cell, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(kegg_cell,paste0("11.",names(top200_entrez[i]),"特定亚群KEGG.csv"))
plotKEGG<-dotplot(kegg_cell,label_format=100)
pdf(paste0('12.',names(top200_entrez[i]),'特定亚群KEGG.pdf'), width = 12, height = 6)
print(plotKEGG)
dev.off()
})
###################################################
####################################################
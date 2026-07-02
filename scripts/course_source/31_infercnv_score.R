# Course-derived reference script
# English filename: 31_infercnv_score.R
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
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(circlize)
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

#读入数据
scRNA=qread("手动注释.qs",nthreads = detectCores())#读取数据
#读取cnv结果
infercnv_obj = readRDS("run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices%>%unlist()
tumor_loc <- infercnv_obj@observation_grouped_cell_indices%>%unlist()
#定义细胞结果
anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc],colnames(expr)[tumor_loc]),
  class=c(rep("normal",length(normal_loc)),rep("tumor",length(tumor_loc)))
)
gn <- rownames(expr)
#读取基因文件
geneFile <- read.table("geneFile.txt",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile) <- geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
#提取cnv结果基因表达量
expr <- expr[intersect(gn,geneFile$V1),]
#提取结果细胞
samecells<-match(colnames(scRNA),anno.df$CB)
anno.df<-anno.df[samecells,]
scRNA<-scRNA[,anno.df$CB]
scRNA$tissue<-anno.df$class
#infercnv结果可视化
pdf(file = "1.infercnv结果UMAP.pdf",width =7.5,height = 5.5)
DimPlot(scRNA,group.by = "tissue", reduction = "umap", label = T,label.box = T, label.size = 3.5,pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "2.infercnv结果TSEN.pdf",width =7.5,height = 5.5)
DimPlot(scRNA,group.by = "tissue", reduction = "tsne", label = T,label.box = T, label.size = 3.5,pt.size = 0.01)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()

#进行聚类
kmeans.result <- kmeans(t(expr), 5)
kmeans_df <- data.frame(Cluster=paste0('C',kmeans.result$cluster))
kmeans_df$CB <- colnames(expr)
kmeans_df <- kmeans_df%>%inner_join(anno.df,by="CB") 
kmeans_df_s <- arrange(kmeans_df,Cluster)
rownames(kmeans_df_s) <- kmeans_df_s$CB
kmeans_df_s$CB <- NULL
kmeans_df_s$Cluster <- as.factor(kmeans_df_s$Cluster) 
kmeans_df_s$cellType<-as.factor(scRNA$cellType)
#绘制热图
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:5]
names(color_v)=paste0('C',1:5)
left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("tumor"="red","normal" = "blue"),Cluster=color_v))
pdf("3.CNV聚类热图.pdf",width = 15,height = 10)
ht = Heatmap(t(expr)[rownames(kmeans_df_s),], #绘图数据的CB顺序和注释CB顺序保持一致
             col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")), #如果是10x的数据，这里的刻度会有所变化
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), #这一步可以控制染色体顺序，即使你的基因排序文件顺序是错的
             column_gap = unit(2, "mm"),
             heatmap_legend_param = list(title = "Scores",
                                         at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
             top_annotation = top_anno,left_annotation = left_anno, #添加注释
             row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()

#计算CNV打分
cnvScore <- function(data){
  require(tidyverse)
  require(scales)
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    scales::rescale(to=c(-1, 1)) %>% 
    t()
  
  cnv_score <- data.frame(ID=colnames(data),Score=colSums(data * data),row.names = NULL)
  return(cnv_score)
}
cnv_score <- cnvScore(expr)
data <- merge(kmeans_df_s,cnv_score,by.x=0,by.y=1)

#绘制CNV打分箱线图
#聚类肿瘤分组CNV箱线图
ggplot(data,aes(Cluster,Score))+
  geom_boxplot(aes(fill=class),outlier.colour = 'grey30',outlier.size = 0.3)+
  labs(x=NULL,y='CNV score',fill=NULL)+
  ggsci::scale_fill_jco()+theme_bw()
ggsave('4.聚类肿瘤分组CNV箱线图.pdf',width = 5,height = 3)

#聚类CNV箱线图
ggplot(data,aes(Cluster,Score))+
  geom_boxplot(aes(fill=Cluster),outlier.colour = 'grey30',outlier.size = 0.3)+
  labs(x=NULL,y='CNV score',fill=NULL)+
  ggsci::scale_fill_jco()+theme_bw()
ggsave('5.聚类CNV箱线图.pdf',width = 5,height = 3)

#细胞CNV箱线图
ggplot(data,aes(cellType,Score))+
  geom_boxplot(aes(fill=cellType),outlier.colour = 'grey30',outlier.size = 0.3)+
  labs(x=NULL,y='CNV score',fill=NULL)+
  ggsci::scale_fill_jco()+theme_bw()
ggsave('6.细胞CNV箱线图.pdf',width = 5,height = 3)

#细胞分组CNV箱线图
ggplot(data,aes(cellType,Score))+
  geom_boxplot(aes(fill=class),outlier.colour = 'grey30',outlier.size = 0.3)+
  labs(x=NULL,y='CNV score',fill=NULL)+
  ggsci::scale_fill_jco()+theme_bw()
ggsave('7.细胞分组CNV箱线图.pdf',width = 5,height = 3)
####################################################
####################################################
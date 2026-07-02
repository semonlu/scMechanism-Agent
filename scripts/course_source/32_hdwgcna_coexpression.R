# Course-derived reference script
# English filename: 32_hdwgcna_coexpression.R
# Original course path: see source_manifest.csv
# Role: coexpression reference
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
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(qs)
library(igraph)
library(parallel)
library(scCustomize)
library(doParallel)
#使用ggplot的cowplot主题画图用
theme_set(theme_cowplot())
options(future.globals.maxSize = 1024^10000)
#设置种子保证结果可以复现
set.seed(12345)

#启用多线程
enableWGCNAThreads(nThreads = 6)

#加载数据
seurat_obj <- qread("手动注释.qs",nthreads = detectCores())#读取数据
##随机抽取10%个细胞，仅作为演示，实际操作或高性能工作站直接运行
N=length(colnames(seurat_obj))/10
N=round(N)
seurat_obj<-seurat_obj[,sample(x=colnames(seurat_obj),size = N,replace=F)]
#fraction ：选择至少在一部分细胞中表达的基因；此时需要追加fraction 参数，如fraction =0.01 表示选择至少在1%的细胞中表达的基因
seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "fraction",
    fraction = 0,
    wgcna_name = "hdWGCNA_project"
)

# 构造metacell用于WGCNA
seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("cellType", "Type"),
    reduction = "harmony", 
    k = 25, # KNN参数
    max_shared = 10, 
    ident.group = "cellType"  
)

# 标准化metacell表达矩阵
seurat_obj <- NormalizeMetacells(seurat_obj)
#创建表达矩阵用于WGCNA
#该部分展示如何在细胞上执行共表达网络分析
seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = c("T cells"), #输入想要分析的细胞
    assay = "RNA",
    layer = "data"
)

#此部分用于计算合适的软阈值，以创建基因共表达网络。
#在WGCNA中，networkType参数用于指定构建基因共表达网络的类型，主要有以下几个可选参数：
#1. unsigned：无符号网络，即不考虑基因表达的正负号，只考虑它们之间的关联性。
#2. signed：带符号网络，即同时考虑基因表达的正负号和大小，可以反映基因之间的正负调控关系。
#3. hybrid：混合网络，即将无符号网络和带符号网络结合起来，利用它们之间的优势来提高网络分析的准确性和可靠性。
#一般选择"signed"参数。
seurat_obj <- TestSoftPowers(
    seurat_obj,
    networkType = "signed" # 此外，还可选择“unsigned”或“signed hybrid”参数
)
plot_list <- PlotSoftPowers(seurat_obj)
pdf("1.最佳软阈值.pdf")
wrap_plots(plot_list, ncol=2)
dev.off()
#hdWGCNA将共表达网络表示为拓扑重叠矩阵（Topoligcal Overlap Matrix, TOM），值代表基因之间的拓扑重叠，TOM矩阵可以简单认为是基因相关性矩阵的升级版。
seurat_obj <- ConstructNetwork(
    seurat_obj, soft_power = 7, #最佳软阈值
    setDatExpr = FALSE,
    tom_name = 'T cells' # TOM矩阵名称，运行后会在当前目录生成一个TOM文件夹保存TOM矩阵
)

pdf("2.模块聚类.pdf")
PlotDendrogram(seurat_obj, main='T cells hdWGCNA Dendrogram')
dev.off()
#计算模块特征值ME
#ME: Module Eigengenes，模块特征基因（值），一般指每个模块的基因表达矩阵的PCA降维后的第一个主成分（PC）。
# harmony前需要进行正态化（Scale）
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

# 计算所有的MEs，比较耗时
seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars = "cellType"  
)

#计算模块连接性
#hub genes：那些与每个模块高度连接的基因
#kME：eigengene-based connectivity，代表着每个基因基于特征基因的连通性，用于确定hubgenes。
#通过ModuleConnectivity 函数在整个单细胞数据集上计算kME values，kME values值越高，该基因为hub genes的可能性越大。
seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = 'cellType', group_name = 'T cells'  # 这里继续计算T cells（感兴趣的细胞类型）的kME
)

# 重命名module，以表明这些module由T cells计算所得
seurat_obj <- ResetModuleNames(
    seurat_obj,
    new_name = "T cells-M"
)
#获取模块分配表
modules <- GetModules(seurat_obj)
# 可视化每个模块中的基因的kME
pdf("3.可视化每个模块中基因的kME.pdf",width = 15,height = 15)
PlotKMEs(seurat_obj)
dev.off()
# 通过GetHubGenes函数抽取每个模块的topN hub节点，即每个模块中kME高的topN个基因
hub_df <- GetHubGenes(seurat_obj = seurat_obj, n_hubs = 10)

#至此，基本完成了hdWGCNA分析
qsave(seurat_obj,"hdWGCNA.qs",nthreads = detectCores())
#利用hub genes进行打分
# 根据每个模块的topN个基因对细胞进行打分
seurat_obj <- ModuleExprScore(
    seurat_obj,
    n_genes = 25, # topN hub genes
    method = "Seurat" # 打分方法，分为Seurat、UCell
)

#模块分数可视化
#MEs：模块特征值
#hMEs：harmony校正过的模块特征值
#scores：模块分数，ModuleExprScore函数添加
#average： GetAvgModuleExpr(seurat_obj) 计算ModuleScore均值
# 根据hMEs进行FeaturePlot
plot_list <- ModuleFeaturePlot(
    seurat_obj,
    features = 'hMEs', # 可选择MEs、hMEs、scores、average
    order = TRUE # 排序，使hME最高的点位于顶部
)

pdf("4.模块分数可视化.pdf",width = 15,height = 15)
wrap_plots(plot_list, ncol = 6)
dev.off()
# 模块相关性
ModuleCorrelogram(seurat_obj)

# 利用Seurat自带的可视化方法
#从seurat对象获取hME
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

#将hME添加到Seurat元数据
#细胞模块相关性
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
p <- DotPlot(seurat_obj, features = mods, group.by = "cellType")
p <- p + 
    RotatedAxis() + # 旋转坐标轴标签
    scale_color_gradient2(high = "red", mid = "grey95", low = "blue")  # 改颜色
pdf("5.细胞模块相关性点图.pdf",,width = 15,height = 15)
p
dev.off()
#模块细胞相关性小提琴图
p <- VlnPlot(
    seurat_obj,
    features = "T cells-M13",
    group.by = "cellType",
    pt.size = 0
)
p <- p + geom_boxplot(width = .25, fill = "white")
p <- p + xlab("") + ylab("hME") + NoLegend()
pdf("6.模块细胞相关性小提琴图.pdf")
p
dev.off()
theme_set(theme_cowplot())
set.seed(12345)
#绘制特定模块基因网络
ModuleNetworkPlot(seurat_obj = seurat_obj, mods = "T cells-M14")
options(repr.plot.height = 10, repr.plot.width = 10)
#计算hub基因网络
pdf("7.hub基因网络.pdf",width = 15,height = 15)
HubGeneNetworkPlot(
    seurat_obj,
    n_hubs = 3, # 用于可视化的hub gene
    n_other = 5, # 随机选取的gene
    edge_prop = 0.75, # 采样的边数
    mods = "all"
)
dev.off()
#计算hub基因在UMAP中的空间关系

seurat_obj <- RunModuleUMAP(
    seurat_obj, 
    n_hubs = 10, # 用于UMAP嵌入的hub genes
    n_neighbors = 15, # UMAP参数
    min_dist = 0.1 # 两个点在UMAP空间中的最短距离
)

umap_df <- GetModuleUMAP(seurat_obj)
pdf("8.hub基因在UMAP中的空间关.pdf",width = 15,height = 15)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) + # 设置坐标轴
    geom_point(
        color = umap_df$color, # 按module对点（gene）着色
        size = umap_df$kME * 2 # 点的大小与其kME值相关
    ) + 
    umap_theme()
dev.off()

# 需要与模块进行关联的特征，该特征从Seurat_obj@meta.data抽取
cur_traits <- c('nCount_RNA', 'nFeature_RNA')#可以定义临床信息

# 将基因模块与特征进行关联
seurat_obj <- ModuleTraitCorrelation(
	seurat_obj,
	traits = cur_traits,  # 需要与模块进行关联的特征
	group.by = 'cellType'  # 该参数表示以cell_type进行分组，以cell_type为单位计算每个cell_type的模块-特征相关性矩阵
)

# 查看相关性矩阵
mt_cor <- GetModuleTraitCorrelation(seurat_obj)
# 模块-特征相关性热图
pdf("9.模块相关性热图.pdf",width = 15,height = 15)
PlotModuleTraitCorrelation(
	seurat_obj,
	label = 'fdr',  # 可选pval、fdr作为显著性分数
	label_symbol = 'stars',  # 以*号作为显著性标记，numeric则显示数值
	text_size = 2,
	text_digits = 2,
	text_color = 'white',
	high_color = 'yellow',
	mid_color = 'black',
	low_color = 'purple',
	plot_max = 0.2,
	combine=TRUE
)
dev.off()
####################################################
####################################################

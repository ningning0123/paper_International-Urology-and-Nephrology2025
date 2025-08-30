# ------------------- 预加载R包及依赖 --------------------
library(Seurat)       # 主要的单细胞分析包
library(dplyr)        # 数据整理、管道操作
library(ggplot2)      # 画图
library(magrittr)     # 管道符 %>%
library(RColorBrewer) # 颜色方案
library(limma)        # 差异分析
library(tidyr)        # 数据整理（长宽表互转等）
library(NMF)          # 非负矩阵分解相关（部分包依赖）
library(CellChat)     # 细胞通讯分析
library(ggalluvial)   # Alluvial 图
library(svglite)      # 导出svg格式图
library(celldex)      # SingleR注释参考数据库
library(SingleR)      # 自动细胞类型注释
library(monocle)      # 单细胞轨迹分析

# ----------------- 设置工作目录 -----------------
workDir <- "F:\\MPM\\03.单细胞分析"  # 工作目录
setwd(workDir)  # 切换到工作目录

# ----------- 1. 参数设置区 ------------
logFC_filter       <- 1        # 差异分析logFC阈值
p_adj_filter       <- 0.05     # 差异分析校正p值阈值
min_cells_gene     <- 5        # 至少在5个细胞出现的基因参与分析
min_genes_per_cell <- 200      # 每个细胞至少检测到的基因数
post_filter_cells  <- 300      # 二次过滤时，细胞内基因数大于此值
post_filter_mito   <- 20       # 线粒体比例不得高于此值（%）
n_top_var_features <- 2500     # 选取前2500个高度变异基因
n_pcs              <- 22        # PCA降维时使用的主成分数
n_topmarker_heat   <- 10       # 每簇热图展示top10 marker
cluster_resolution <- 0.6      # 聚类分辨率
neighbor_dims      <- 15       # 聚类与UMAP时用多少PC
qc_vln_width       <- 15       # 质控小提琴图宽度
qc_vln_height      <- 7        # 质控小提琴图高度
input_expr_file    <- "single_cell_data.csv"  # 输入表达数据文件名
output_dir         <- "analysis_results"      # 结果输出文件夹

# ----------- 2. 文件夹确认 -------
if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)  # 如果输出目录不存在则创建
cat("输出目录设置为：", output_dir, "\n")  # 输出目录提示

# ----------- 3. 数据读取与初处理 -------
cat("Step 1: 读取表达数据...\n")  # 流程进度提示
if(!file.exists(input_expr_file)) stop("表达数据文件不存在！")  # 检查输入文件是否存在
expr_table <- read.csv(input_expr_file, header=TRUE, sep=",", stringsAsFactors=FALSE, check.names=FALSE) # 读取表达矩阵
if(nrow(expr_table)==0) stop("表达矩阵无内容！")  # 检查表格非空
rownames(expr_table) <- expr_table[,1]   # 首列为基因名，设为行名
gene_names <- expr_table[,1]             # 保存基因名向量
expr_values <- as.matrix(expr_table[,-1]) # 去掉第一列后的原始表达值
expr_matrix <- expr_values                # 表达矩阵后续用于Seurat

# ----------- 4. 创建Seurat对象及QC小提琴图 ---------
cat("Step 2: 创建Seurat对象...\n") # 流程提示
scObject <- CreateSeuratObject(
  counts = as.matrix(expr_matrix),         # 原始表达量
  project = "scRNAseqProject",             # 项目名
  min.cells = min_cells_gene,              # 参与分析的基因必须在多少细胞中存在
  min.features = min_genes_per_cell,       # 细胞需检测多少基因
  names.delim = "_"                        # 分隔符
)
scObject[["percent.mito"]] <- PercentageFeatureSet(scObject, pattern = "^MT-") # 线粒体基因比例
pdf(file.path(output_dir, "QC_violin_basicMetrics_pubstyle.pdf"), width=12, height=6)
p <- VlnPlot(
  scObject,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
  pt.size = 0,
  group.by = "orig.ident"
) +
  theme_classic(base_size = 16) +  # 高级美化
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x  = element_text(size = 12, face = "bold", angle = 45, vjust=1, hjust=1),
    axis.text.y  = element_text(size = 12, face = "bold"),
    strip.text   = element_text(size = 15, face = "bold", color = "black")
  )
print(p)  # 必须print
dev.off()

# ----------- 5. 二次过滤与QC散点 -------
cat("Step 3: 二次细胞质控...\n") # 流程提示
cell_num_before <- ncol(scObject) # 过滤前细胞数
scObject <- subset(scObject, subset = nFeature_RNA > post_filter_cells & percent.mito < post_filter_mito) # 基于基因数和线粒体比过滤
cell_num_after <- ncol(scObject) # 过滤后细胞数
cat(sprintf("细胞过滤: 原%d，剩%d。\n", cell_num_before, cell_num_after)) # 输出细胞过滤信息
if(cell_num_after < 10) stop("过滤后过少细胞，检查阈值！") # 过少细胞则终止

pdf(file.path(output_dir, "QC_scatter_metrics.pdf"), width=13, height=7) # PDF输出
plot1 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size = 1.5)   # UMI-线粒体比例
plot2 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 1.5)   # UMI-基因数
CombinePlots(plots = list(plot1, plot2)) # 合并图
dev.off()

# ----------- 6. 归一化及变异基因 -------
cat("Step 4: 归一化+高变基因...\n") # 流程提示
scObject <- NormalizeData(scObject, normalization.method="LogNormalize", scale.factor=10000) # 归一化
scObject <- FindVariableFeatures(scObject, selection.method="vst", nfeatures=n_top_var_features) # 高变基因
pdf(file.path(output_dir, "VarGenes_overview.pdf"), width=10, height=6) # 输出高变基因分布
VariableFeaturePlot(scObject)
dev.off()

# ----------- 7. PCA降维与可视化 ----------
cat("Step 5: PCA降维...\n") # 流程提示
scObject <- ScaleData(scObject) # 标准化数据（均值0方差1）
scObject <- RunPCA(scObject, npcs=n_pcs, features=VariableFeatures(scObject)) # PCA降维，n_pcs主成分
pdf(file.path(output_dir, "PCA_geneLoadings.pdf")) # 主成分主要基因可视化
VizDimLoadings(scObject, dims=1:3, reduction="pca", nfeatures=25)
dev.off()
pdf(file.path(output_dir, "PCA_heatmap_topGenes.pdf"), width=8, height=7) # 前三主成分主要基因热图
DimHeatmap(scObject, dims=1:3, cells=400, balanced=TRUE)
dev.off()

# ----------- 8. 聚类与UMAP降维 ----------
cat("Step 6: 聚类与UMAP降维...\n") # 流程提示
scObject <- FindNeighbors(scObject, dims=1:neighbor_dims)      # 构建K近邻图
scObject <- FindClusters(scObject, resolution=cluster_resolution) # 聚类
scObject <- RunUMAP(scObject, dims=1:neighbor_dims)            # UMAP降维
nColors <- length(unique(scObject$seurat_clusters))             # 统计聚类类别数量
myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nColors)  # 自定义颜色
pdf(file.path(output_dir, "UMAP_clustered_samples.pdf"), width=7, height=5)
DimPlot(scObject, reduction="umap", label=TRUE, cols=myColors) # 绘制UMAP并标注聚类
dev.off()
write.csv(data.frame(Cell=colnames(scObject), Cluster=scObject$seurat_clusters), 
          file=file.path(output_dir, "CellCluster_UMAP_assignments.csv"), row.names=FALSE) # 保存每个细胞的聚类信息

# ----------- 9. 差异基因查找/热图 -----------
cat("Step 7: 差异表达分析...\n") # 流程提示
marker_df <- tryCatch({
  FindAllMarkers(scObject, only.pos=TRUE, min.pct=0.2, logfc.threshold=logFC_filter) # 所有聚类差异表达基因
}, error=function(e){cat("Marker分析异常\n"); NULL})
if(!is.null(marker_df) && nrow(marker_df)>0){
  sigMarkers <- marker_df[abs(marker_df$avg_log2FC)>logFC_filter & marker_df$p_val_adj<p_adj_filter,] # 显著marker筛选
  write.csv(sigMarkers, file=file.path(output_dir, "Cluster_Markers_DEGs.csv"), row.names=FALSE)      # 保存marker
  topmarker <- marker_df %>% group_by(cluster) %>% top_n(n_topmarker_heat, avg_log2FC)                # 每类top marker
  pdf(file.path(output_dir, "Markers_DoHeatmap.pdf"))
  DoHeatmap(scObject, features=topmarker$gene, size=4) + NoLegend()  # marker热图
  dev.off()
} else {
  cat('未检出显著marker。\n') # 若未检测到显著基因
}

# ----------- 10. SingleR细胞自动注释+合并marker表(紧跟marker后，防变量丢失) -----------
cat("Step 8: SingleR细胞类型自动注释...\n") # 流程提示
expr4singler <- GetAssayData(scObject, slot = "data")                # 获取表达数据用于SingleR
clusters4singler <- scObject$seurat_clusters                         # 获取聚类标签
ref_hpa <- celldex::HumanPrimaryCellAtlasData()                      # 加载人类细胞参考库
singler_result <- SingleR(
  test = expr4singler,
  ref = ref_hpa,
  labels = ref_hpa$label.main,
  clusters = clusters4singler
)   # 运行SingleR注释
scObject$SingleR_celltype <- singler_result$labels[as.numeric(scObject$seurat_clusters)+1]  # 将注释写回meta
write.csv(
  data.frame(Cluster=rownames(singler_result), CellType=singler_result$labels),
  file = file.path(output_dir, "CellCluster_AutoAnno_SingleR.csv"), row.names=FALSE
)   # 输出注释表
if(!is.null(marker_df) && !is.null(singler_result)) {
  cluster_anno <- data.frame(
    cluster = as.character(rownames(singler_result)),
    CellType = singler_result$labels,
    stringsAsFactors = FALSE
  )  # 建立聚类注释表
  marker_df$cluster <- as.character(marker_df$cluster)
  marker_df_anno <- marker_df %>% left_join(cluster_anno, by = "cluster") # 合并marker和注释
  write.csv(marker_df_anno, file=file.path(output_dir,"Cluster_Markers_DEGs_with_CellType.csv"), row.names=FALSE)
}

nCellTypes <- length(unique(scObject$SingleR_celltype))  # 统计细胞类型数
celltypeColors <- colorRampPalette(brewer.pal(12, "Set3"))(nCellTypes) # 颜色
pdf(file.path(output_dir, "UMAP_celltype_autoAnnot_Set3.pdf"), width = 8, height = 6)
DimPlot(scObject, group.by = "SingleR_celltype", reduction = "umap", label = TRUE, cols = celltypeColors) # 画UMAP并标注细胞类型
dev.off()

# ----------- 11. 细胞类型统计可视化及marker数 Alluvial等 -----------
cat("Step 9: 各种统计可视化输出...\n") # 流程提示
meta_df <- scObject@meta.data
celltype_stats <- as.data.frame(table(
  Sample = meta_df$orig.ident, 
  CellType = meta_df$SingleR_celltype
))         # 按样本和细胞类型做数量统计
celltype_stats <- celltype_stats %>%
  group_by(Sample) %>%
  mutate(Ratio = Freq / sum(Freq)) %>%
  ungroup()  # 计算各样本细胞类型的比例
celltype_stats$Sample <- factor(celltype_stats$Sample, levels = unique(celltype_stats$Sample)) # 因子化
celltype_stats$CellType <- factor(celltype_stats$CellType, levels = unique(celltype_stats$CellType)) # 因子化
cell_types <- levels(celltype_stats$CellType) # 获得全部细胞类型
nColors <- length(cell_types) # 类型数
myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nColors) # 颜色
names(myColors) <- cell_types # 命名

pdf(file.path(output_dir, "Sample_CellType_Composition_barplot.pdf"), width=7, height=6) # 打开pdf
ggplot(celltype_stats, aes(x = Sample, y = Ratio, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7, color = "grey70") +
  scale_y_continuous(expand = c(0,0), labels = scales::percent, name = "Proportion") +
  scale_fill_manual(values = myColors, name="Cell Type") +
  labs(x = "Sample") +
  theme_classic(base_size=16) +
  coord_flip() +  
  theme(
    legend.title=element_text(size=14),
    legend.text=element_text(size=13),
    axis.text = element_text(size=13),
    axis.title = element_text(size=14)
  ) # 绘制每个样本各细胞类型比例的条形图
dev.off()

cell_counts <- as.data.frame(table(Sample = meta_df$orig.ident, CellType = meta_df$SingleR_celltype)) # 统计数量
write.csv(cell_counts, file.path(output_dir, "CellNumber_perSample_CellType.csv"), row.names = FALSE) # 输出分组计数表
celltype_matrix <- celltype_stats %>% select(Sample, CellType, Ratio) %>%
  tidyr::spread(CellType, Ratio, fill = 0)
write.csv(celltype_matrix, file.path(output_dir, "CellType_Proportion_Matrix.csv"), row.names=FALSE) # 输出矩阵
cell_annot <- data.frame(CellBarcode = rownames(meta_df), Sample = meta_df$orig.ident,
                         Cluster = meta_df$seurat_clusters, CellType = meta_df$SingleR_celltype)
write.csv(cell_annot, file.path(output_dir, "Cell_Full_Annotation.csv"), row.names=FALSE) # 输出所有细胞信息注释表

# marker计数及绘图
if (exists("marker_df_anno")) {
  gene_count_by_celltype <- marker_df_anno %>%
    group_by(CellType) %>%
    summarise(gene_number = n_distinct(gene)) %>%
    arrange(desc(gene_number))   # 每类marker基因数
  write.csv(gene_count_by_celltype, file = file.path(output_dir, "GeneNumber_per_CellType.csv"), row.names = FALSE) # 保存
  pdf(file.path(output_dir, "MarkerGeneCount_perCellType_bar.pdf"), width=8, height=5)
  ggplot(gene_count_by_celltype, aes(x = reorder(CellType, gene_number), y = gene_number, fill=CellType)) +
    geom_bar(stat="identity", width=0.7, color="grey40") +
    scale_fill_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(nrow(gene_count_by_celltype))) +
    coord_flip() +
    labs(x="Cell Type", y="Number of Marker Genes", title="Number of Marker Genes per Cell Type") +
    theme_classic(base_size=16) +
    theme(legend.position = "none")
  dev.off()
  # Alluvial流图
  top_genes_per_CT <- marker_df_anno %>%
    group_by(CellType) %>% top_n(5, abs(avg_log2FC)) %>% ungroup() # 每类top5 marker
  pdf(file.path(output_dir, "MarkerGene_CellType_Alluvial.pdf"), width=9, height=6)
  ggplot(top_genes_per_CT, aes(axis1 = CellType, axis2 = gene, y = abs(avg_log2FC))) +
    geom_alluvium(aes(fill=CellType), width=1/12) +
    geom_stratum(width=1/6, fill="grey90", color="black") +
    geom_text(stat="stratum", aes(label=after_stat(stratum)), size=3) +
    scale_x_discrete(limits = c("CellType","MarkerGene"), expand = c(.05, .05)) +
    theme_minimal(base_size=14) +
    ggtitle("Top Marker Genes Alluvial Plot")
  dev.off()
}

# ----------- PCA方差解释与推荐主成分数输出 -----------
pca_stdev <- scObject[["pca"]]@stdev
pca_var_explained <- (pca_stdev^2) / sum(pca_stdev^2)
cumulative_var <- cumsum(pca_var_explained)
pca_table <- data.frame(
  PC = seq_along(pca_var_explained),
  Variance_Explained = pca_var_explained,
  Cumulative_Variance = cumulative_var
)
write.csv(pca_table, file = file.path(output_dir, "PCA_Variance_Explained.csv"), row.names=FALSE)

opt_cutoff <- 0.8  # 推荐主成分累计解释率阈值80%
best_n_pcs <- which(cumulative_var >= opt_cutoff)[1]
write.csv(
  data.frame(Recommended_n_PCs = best_n_pcs, Cumulative_Variance = cumulative_var[best_n_pcs]),
  file = file.path(output_dir, "PCA_Recommended_nPCs.csv"),
  row.names=FALSE
)


cat("所有结果/表格/图形均已输出到文件夹：", output_dir, "\n") # 流程结束提示

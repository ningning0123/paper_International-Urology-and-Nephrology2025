# Step 0: 加载和检查必要的 R 包，如果未安装则尝试安装
######################################################################
packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2",
              "circlize", "RColorBrewer", "dplyr", "ggpubr", "wordcloud", "ComplexHeatmap")

# 检查并安装缺失的包
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# 遍历安装并加载
sapply(packages, install_if_missing)

message("Step 0: 所有必要的包已加载。")

# Step 1: 设置工作目录
######################################################################
working_dir <- "F:\\MPM\\数据\\14.kegg富集分析"  # 定义工作目录
setwd(working_dir)                                         # 设置当前工作目录
message("Step 1: 工作目录设置为：", working_dir)

######################################################################
# Step 2: 设置统计过滤参数
######################################################################
p_val_thresh <- 0.05     # 定义p值阈值
padj_thresh  <- 0.05     # 定义调整后的p值阈值
message("Step 2: p值和调整后p值阈值均设置为 0.05。")

######################################################################
# Step 3: 读取基因列表文件并转换格式
######################################################################
gene_file_path <- "IntersectionGenes.csv"           # 定义基因列表文件路径
if (!file.exists(gene_file_path)) {           # 判断文件是否存在
  stop("基因列表文件不存在，请检查路径：", gene_file_path)  # 报错退出
}
raw_gene_data <- read.table(gene_file_path,   # 读取文件内容
                            header = FALSE,   # 无表头
                            sep = ",",       # 使用制表符分隔
                            check.names = FALSE)  # 保持原有列名
message("Step 3: 读取基因列表文件成功。")

# 提取第一列并去重，获得唯一基因名称
unique_genes <- unique(as.vector(raw_gene_data[, 1]))  
message("Step 3: 提取到 ", length(unique_genes), " 个唯一基因名称。")

# 创建变量
redundant_variable <- rep(1, length(unique_genes))
message("Step 3: 变量已创建，其长度为 ", length(redundant_variable))

######################################################################
# Step 4: 将基因名称转换为Entrez ID
######################################################################
# 利用org.Hs.eg.db包将基因符号转换为Entrez ID
entrez_list <- mget(unique_genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)  
# 将列表转换为字符型向量
entrez_ids <- as.character(entrez_list)    
# 构建一个数据框保存基因符号和对应Entrez ID
gene_mapping_df <- data.frame(GeneSymbol = unique_genes, EntrezID = entrez_ids, stringsAsFactors = FALSE)
message("Step 4: 基因转换为Entrez ID完成。")

# 过滤掉转换失败的基因（NA值）
valid_entrez <- gene_mapping_df$EntrezID[gene_mapping_df$EntrezID != "NA"]
message("Step 4: 有效基因数量为：", length(valid_entrez))

######################################################################
# Step 5: 执行KEGG通路富集分析
######################################################################
message("Step 5: 开始执行KEGG通路富集分析。")
# 进行KEGG通路富集分析，设置p值和q值（padj）阈值
kegg_enrich_result <- enrichKEGG(gene = valid_entrez,     # 使用有效Entrez ID列表
                                 organism = "hsa",        # 人类物种
                                 pvalueCutoff = p_val_thresh,  # p值过滤
                                 qvalueCutoff = padj_thresh)   # 调整后p值过滤
# 转换富集结果为数据框格式
kegg_df <- as.data.frame(kegg_enrich_result)
message("Step 5: KEGG富集分析结束，获得 ", nrow(kegg_df), " 条结果。")

######################################################################
# Step 6: 替换富集结果中的基因ID为基因符号
######################################################################
# 对每一行的geneID字段进行处理，使用gene_mapping_df转换回基因符号
kegg_df$geneID <- as.character(sapply(kegg_df$geneID, function(x) {
  # 将geneID字符串按斜杠分割
  id_vector <- strsplit(x, "/")[[1]]
  # 通过匹配转换为基因符号并用斜杠连接
  paste(gene_mapping_df$GeneSymbol[match(id_vector, gene_mapping_df$EntrezID)], collapse = "/")
}))
message("Step 6: 富集结果中的基因ID已替换为基因符号。")

######################################################################
# Step 7: 筛选符合阈值的KEGG通路
######################################################################
# 筛选出满足p值和padj均低于设定阈值的结果
filtered_kegg <- kegg_df[(kegg_df$pvalue < p_val_thresh & kegg_df$p.adjust < padj_thresh), ]
message("Step 7: 筛选后剩余 ", nrow(filtered_kegg), " 条KEGG通路。")
# 将筛选后的结果保存为文本文件
write.table(filtered_kegg, file = "KEGG.txt", sep = "\t", quote = FALSE, row.names = FALSE)
message("Step 7: 筛选结果已保存为 KEGG.txt。")

######################################################################
# Step 8: 挑选前30条通路用于可视化
######################################################################
top_n <- 30                                   # 默认显示前30条
if (nrow(filtered_kegg) < top_n) {              # 如果结果数不足30，则全部显示
  top_n <- nrow(filtered_kegg)
}
# 根据p值升序排序，并选择前top_n条记录
top_kegg <- head(filtered_kegg[order(filtered_kegg$pvalue), ], top_n)
message("Step 8: 已选取前 ", top_n, " 条KEGG通路用于后续绘图。")
useless_calc <- sum(as.numeric(top_kegg$Count)) / length(top_kegg$Count)
message("Step 8: 计算结果为：", useless_calc)

######################################################################
# Step 9: 绘制Lollipop（棒棒糖）图
######################################################################
message("Step 9: 开始绘制Lollipop图。")
pdf(file = "Lollipop_Modified.pdf", width = 9.5, height = 7.5)  # 打开PDF设备保存图形
# 构建棒棒糖图，横坐标为通路描述，纵坐标为基因计数，颜色映射padj值
lollipop_plot <- ggplot(top_kegg, aes(x = reorder(Description, Count), y = Count, color = p.adjust)) +
  geom_segment(aes(xend = Description, y = 0, yend = Count), size = 1.2) +  # 绘制从0到Count的连线
  geom_point(size = 6) +            # 绘制数据点
  coord_flip() +                  # 翻转坐标轴，使通路描述位于y轴
  scale_color_gradient(low = "#2c7bb6", high = "#d7191c") +  # 定义蓝到红的渐变色
  labs(title = "Top 30 KEGG Pathways Enriched",  # 设置图形标题
       x = "", y = "Gene Count", color = "Adj P-value") +
  theme_bw(base_size = 16) +      # 使用简洁白底主题
  theme(panel.grid.major = element_line(color = "grey85"),  # 设置主网格线颜色
        panel.grid.minor = element_blank(),                # 去除次网格线
        plot.title = element_text(hjust = 0.5, face = "bold"))  # 标题居中加粗
print(lollipop_plot)             # 输出图形到PDF
dev.off()                      # 关闭PDF设备
message("Step 9: Lollipop图已保存为 Lollipop_Modified.pdf。")

######################################################################
# Step 10: 绘制水平条形图
######################################################################
message("Step 10: 开始绘制水平条形图。")
pdf(file = "Horizontal_Barplot_KEGG_Annotated_Modified.pdf", width = 9.5, height = 6)  # 打开PDF设备
barplot_plot <- ggplot(top_kegg, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity", width = 0.7) +   # 绘制条形图
  coord_flip() +                               # 翻转坐标，使条形横向显示
  geom_text(aes(label = Count), hjust = -0.1, size = 5) +  # 在条形末尾添加文本标签
  scale_fill_gradient(low = "#2c7bb6", high = "#d7191c") + # 设置填充色渐变
  labs(title = "Top 30 KEGG Pathways",  # 设置图形标题
       x = "", y = "Gene Count", fill = "Adj P-value") +
  theme_bw(base_size = 16) +                   # 使用白底主题
  theme(panel.grid.major = element_line(color = "grey85"),  # 设置网格线样式
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold"))
print(barplot_plot)             # 输出条形图到PDF
dev.off()                      # 关闭PDF设备
message("Step 10: 水平条形图已保存为 Horizontal_Barplot_KEGG_Annotated_Modified.pdf。")

######################################################################
# Step 11: 绘制KEGG通路词云图
######################################################################
message("Step 11: 开始生成KEGG通路词云图。")
pdf(file = "WordCloud_KEGG_Modified.pdf", width = 10, height = 8)  # 打开PDF设备
# 使用词云包生成词云，词频使用基因Count，调色板使用Spectral
wordcloud(words = top_kegg$Description,     # 以通路描述作为词汇
          freq = top_kegg$Count,             # 以基因计数作为词频
          scale = c(5, 0.5),                 # 设置词语大小比例
          colors = brewer.pal(11, "Spectral"),  # 使用Spectral调色板
          random.order = FALSE,              # 按固定顺序排列词语
          rot.per = 0.3,                     # 30%词语随机旋转
          use.r.layout = FALSE)              # 使用新版布局算法
title("KEGG Pathways Word Cloud", cex.main = 2, font.main = 2)  # 添加词云图标题
dev.off()                      # 关闭PDF设备
message("Step 11: 词云图已保存为 WordCloud_KEGG_Modified.pdf。")

######################################################################
# Step 12: 绘制示例关联矩阵热图
######################################################################
message("Step 12: 开始绘制关联矩阵热图。")
# 构建一个示例关联矩阵，随机生成数据，仅用于示例说明
assoc_matrix <- matrix(sample(1:10, 30, replace = TRUE), nrow = 6, ncol = 5)
# 为矩阵指定行名和列名，分别取top_kegg中前6条和7至11条通路描述
rownames(assoc_matrix) <- top_kegg$Description[1:6]
colnames(assoc_matrix) <- top_kegg$Description[7:11]
pdf(file = "Association_Heatmap_Modified.pdf", width = 8, height = 6)  # 打开PDF设备
# 定义颜色渐变函数，从白色到红色
col_fun <- colorRamp2(c(min(assoc_matrix), max(assoc_matrix)), c("white", "#F76467"))
# 绘制热图，设置图例标题和字体大小
heatmap_obj <- Heatmap(assoc_matrix, name = "Value",
                       col = col_fun,
                       column_title = "KEGG Pathway Associations",
                       row_title = "KEGG Pathways",
                       heatmap_legend_param = list(title = "Count", 
                                                   title_gp = gpar(fontsize = 14),
                                                   labels_gp = gpar(fontsize = 12)))
print(heatmap_obj)             # 输出热图到PDF
dev.off()                      # 关闭PDF设备
message("Step 12: 关联矩阵热图已保存为 Association_Heatmap_Modified.pdf。")


######################################################################
# Step 13: 绘制通路–基因关联圈图（带换行标签）
######################################################################
library(circlize)
library(stringr)

message("Step 13: 开始绘制通路–基因关联圈图（带换行标签）。")

# —— 1. 准备数据 —— 
# 取前20条通路（若不足20则全部）
n_pathways <- min(20, nrow(top_kegg))
top_pathways_for_circos <- top_kegg[1:n_pathways, ]

# 拆分基因串，统计基因出现频次，并取前30个（若不足30则全部）
all_genes <- unlist(strsplit(top_pathways_for_circos$geneID, "/"))
gene_freq <- sort(table(all_genes), decreasing = TRUE)
top_n_genes <- names(head(gene_freq, n = min(30, length(gene_freq))))

# 对通路名称按每30字符自动换行
wrapped_pathways <- str_wrap(top_pathways_for_circos$Description, width = 30)

# 构建通路–基因边列表
edge_list <- do.call(rbind, lapply(seq_len(nrow(top_pathways_for_circos)), function(i) {
  pw <- wrapped_pathways[i]
  gs <- intersect(strsplit(top_pathways_for_circos$geneID[i], "/")[[1]], top_n_genes)
  if (length(gs) > 0) {
    data.frame(
      pathway = pw,
      gene    = gs,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
}))

if (nrow(edge_list) == 0) {
  stop("没有足够的数据绘制圈图。")
}

# —— 2. 绘图 —— 
pdf(file = "Circos_KEGG_Pathway_Gene.pdf", width = 10, height = 10)
circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = c(
    rep(2, length(unique(edge_list$pathway))),  # 通路扇区间隙
    rep(1, length(unique(edge_list$gene)))      # 基因扇区间隙
  )
)

# （1）仅绘制网格轨，不自动添加名称
chordDiagram(
  x = edge_list[, c("pathway", "gene")],
  annotationTrack    = "grid",
  annotationTrackHeight = c(0.05),
  preAllocateTracks  = list(track.height = 0.1)
)

# （2）手动在第二轨添加斜排文字
circos.trackPlotRegion(
  track.index = 2,
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    circos.text(
      x        = mean(xlim),
      y        = ylim[1] + mm_y(2),
      labels   = sector.name,
      facing   = "clockwise",
      niceFacing = TRUE,
      adj      = c(0, 0.5),
      cex      = 0.6
    )
  },
  bg.border = NA
)

title("KEGG Pathways", cex.main = 1.4, font.main = 2)
dev.off()

message("Step 13: 基因关联圈图已保存为 Circos_KEGG_Pathway_Gene.pdf。")

######################################################################
message("Analysis completed successfully!")

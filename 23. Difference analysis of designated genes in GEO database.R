# 安装和加载必要的R包
# if (!requireNamespace("BiocManager", quietly = TRUE)) 
#     install.packages("BiocManager")
# BiocManager::install("limma")
# install.packages("ggpubr")

# 加载必要的R包
library(limma)   # 用于微阵列数据分析
library(ggpubr)  # 用于绘制美观的统计图形
library(pROC)    # 用于绘制和分析ROC曲线
library(ggplot2) # 用于图形的增强

# 设置工作目录，替换为您本地的路径
setwd("F:\\MPM\\数据\\22.外部验证-差异和ROC曲线")

# 读取标准化后的表达矩阵文件
expFile = "Sample Type Matrix.csv"  # 文件路径
rt = read.csv(expFile, header = TRUE, sep = ",", check.names = FALSE, row.names = 1)

# 提取样本的分组标签，并转换为0和1
y = gsub("(.*)\\_(.*)", "\\2", colnames(rt))  # 从列名提取分组标签
y = ifelse(y == "con", 0, 1)  # 对照组为0，实验组为1

# 读取基因列表文件
geneFile = "IntersectionGenes.csv"  # 文件路径
geneRT = read.csv(geneFile, header = FALSE, sep = ",", check.names = FALSE)

# 从基因列表中提取所需的基因
selectedGenes = as.vector(geneRT[, 1])  # 提取基因名
rt_filtered = rt[selectedGenes, , drop = FALSE]  # 仅保留基因列表中的基因

# 设置比较组（对照组和实验组）样本数量
conNum = sum(y == 0)
treatNum = sum(y == 1)
Type = c(rep("con", conNum), rep("treat", treatNum))  # 对照组与实验组标签

# 设置调色板
colors <- c("con" = "#FF6347", "treat" = "#4682B4")  # 红色和蓝色

# 创建一个数据框来存储基因的差异分析结果
result_df <- data.frame(
  Gene = character(),
  P_Value = numeric(),
  AUC = numeric(),
  AUC_Lower_CI = numeric(),
  AUC_Upper_CI = numeric(),
  SE = numeric(),
  Mean_Con = numeric(),
  Mean_Treat = numeric(),
  stringsAsFactors = FALSE
)

# 设置进度条来跟踪分析进度
cat("开始基因差异分析...\n")
total_genes = nrow(rt_filtered)
pb <- txtProgressBar(min = 0, max = total_genes, style = 3)  # 初始化进度条

# 遍历基因列表中的每个基因，绘制箱线图并保存
for (i in rownames(rt_filtered)) {
  
  # 获取该基因的表达数据，并创建数据框
  rt1 = data.frame(expression = as.numeric(rt_filtered[i, ]), Type = Type)
  
  # 执行t检验
  t_test_result = t.test(as.numeric(rt_filtered[i, ]) ~ Type)
  
  # 提取t检验的p值和标准误（SE）
  p_value = t_test_result$p.value
  SE = t_test_result$stderr
  
  # 计算每组的均值
  mean_con = mean(rt1[rt1$Type == "con", "expression"])
  mean_treat = mean(rt1[rt1$Type == "treat", "expression"])
  
  # 计算并绘制ROC曲线
  roc1 = roc(y, as.numeric(rt_filtered[i, ]))
  ci1 = ci.auc(roc1, method = "bootstrap")  # 计算AUC的95%置信区间
  ciVec = as.numeric(ci1)
  
  # 将每个基因的统计结果添加到result_df中
  result_df = rbind(result_df, data.frame(
    Gene = i,
    P_Value = p_value,
    AUC = ciVec[2],
    AUC_Lower_CI = ciVec[1],
    AUC_Upper_CI = ciVec[3],
    SE = SE,
    Mean_Con = mean_con,
    Mean_Treat = mean_treat
  ))
  
  # 绘制箱线图并添加显著性比较
  boxplot = ggplot(rt1, aes(x = Type, y = expression, fill = Type)) +
    geom_boxplot(outlier.colour = "red", outlier.size = 3, width = 0.5, alpha = 0.7) +
    scale_fill_manual(values = colors) +
    geom_jitter(color = "black", size = 2, width = 0.2) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      legend.title = element_blank(),
      legend.position = "none"
    ) +
    labs(y = paste(i, "expression")) +
    stat_compare_means(method = "t.test")  # 显示统计比较
  
  # 保存箱线图到PDF文件
  pdf(file = paste0("boxplot.", i, ".pdf"), width = 3.4, height = 4.5)
  print(boxplot)
  dev.off()
  
  # 绘制并保存ROC曲线
  pdf(file = paste0("ROC.", i, ".pdf"), width = 5, height = 5)
  plot(roc1, print.auc = TRUE, col = "orange", legacy.axes = TRUE, main = i)  # 绘制ROC曲线
  
  # 添加AUC的置信区间阴影
  polygon(c(roc1$specificities, rev(roc1$specificities)), 
          c(roc1$sensitivities, rep(0, length(roc1$specificities))),
          col = rgb(1, 0.647, 0, 0.3), border = NA)  # 设置阴影的颜色和透明度
  
  # 在图中添加置信区间文本
  text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), col = "orange")
  dev.off()
  
  # 绘制并保存密度图
  density_plot = ggplot(rt1, aes(x = expression, fill = Type)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    labs(title = paste(i, "Density Plot"), x = "Expression", y = "Density")
  
  # 保存密度图
  pdf(file = paste0("density_plot.", i, ".pdf"), width = 6, height = 6)
  print(density_plot)
  dev.off()
  
  # 更新进度条
  setTxtProgressBar(pb, which(rownames(rt_filtered) == i))  # 更新进度条
}

# 关闭进度条
close(pb)

# 将结果保存为CSV文件
write.csv(result_df, file = "gene_analysis_results.csv", row.names = FALSE)

# 输出结果表格
cat("差异分析结果已保存\n")
print(result_df)  # 打印最终结果

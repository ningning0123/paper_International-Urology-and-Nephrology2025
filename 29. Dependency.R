# 加载包，用于数据处理、统计分析和绘图
library("limma")
library(dplyr)
library(tidyverse)
library(ggplot2)
library(linkET)
library(RColorBrewer)

# 设置工作目录（请根据实际情况修改）
setwd("F:\\MPM\\数据\\30.单基因和ssGSEA评分的相关性分析\\B4GALT5")  # 设置当前工作目录

# =============================== #
# Step 1: 定义文件路径及检查文件存在性  #
# =============================== #

# 定义输入文件的路径
exprFilePath <- "Sample Type Matrix.csv"         # 基因表达数据文件（CSV格式）
targetGeneFile <- "IntersectionGenes.csv"       # 待分析基因列表文件（txt格式，使用制表符分隔）
immuneDataPath <- "ssGSEA_Scores.csv"    

# 判断文件是否存在，若不存在则中止执行
if (!file.exists(exprFilePath)) {
  stop("错误：基因表达数据文件不存在！")
}
if (!file.exists(targetGeneFile)) {
  stop("错误：目标基因列表文件不存在！")
}
if (!file.exists(immuneDataPath)) {
  stop("错误：免疫细胞数据文件不存在！")
}

# =============================== #
# Step 2: 定义数据读取与预处理函数  #
# =============================== #

# 定义函数：读取并处理基因表达矩阵，同时筛选出目标基因数据
readExprData <- function(fileExpr, fileGene) {
  cat("步骤 2.1：读取基因表达数据...\n")  # 进度提示
  exprRaw <- read.table(fileExpr, header = TRUE, sep = ",", check.names = FALSE)  # 读取CSV数据
  exprMat <- as.matrix(exprRaw)                            # 转换为矩阵格式
  rownames(exprMat) <- exprMat[, 1]                         # 第一列作为行名（基因名）
  exprVals <- exprMat[, -1, drop = FALSE]                   # 去除第一列，仅保留表达值
  
  # 设置行和列名，转换为数值矩阵
  rowNamesTmp <- rownames(exprVals)                       # 提取行名
  colNamesTmp <- colnames(exprVals)                       # 提取列名
  numExpr <- matrix(as.numeric(as.matrix(exprVals)), 
                    nrow = nrow(exprVals), 
                    dimnames = list(rowNamesTmp, colNamesTmp))  # 转换为数值型
  
  numExpr <- avereps(numExpr)  # 对重复探针做平均处理
  
  cat("步骤 2.2：读取目标基因文件...\n")  # 进度提示
  geneList <- read.table(fileGene, header = FALSE, sep = ",", check.names = FALSE)  # 读取目标基因文件
  targetGene <- as.vector(geneList[, 1])     # 提取目标基因名称
  
  # 判断目标基因是否存在于表达矩阵中
  if (!(targetGene %in% rownames(numExpr))) {
    stop("错误：目标基因未在表达数据中找到！")
  }
  
  # 筛选出目标基因数据，并进行转置（样本为行，基因为列）
  exprSelected <- numExpr[targetGene, , drop = FALSE]
  exprSelected <- t(exprSelected)
  
  # 返回一个列表，包含处理后的数据和目标基因名称
  return(list(exprData = exprSelected, geneName = targetGene))
}

# 定义函数：读取免疫细胞数据并与表达数据匹配
readImmuneData <- function(fileImm, exprData) {
  cat("步骤 2.3：读取免疫细胞数据...\n")  # 进度提示
  immData <- read.table(fileImm, header = TRUE, sep = ",", check.names = FALSE, row.names = 1)  # 读取免疫数据
  # 提取在表达数据和免疫数据中共有的样本
  commonSamples <- intersect(rownames(exprData), rownames(immData))
  
  # 判断是否存在共有样本
  if (length(commonSamples) == 0) {
    stop("错误：没有共同的样本在表达数据与免疫数据中！")
  }
  
  exprData <- exprData[commonSamples, , drop = FALSE]  # 仅保留共有样本的表达数据
  immData <- immData[commonSamples, , drop = FALSE]      # 仅保留共有样本的免疫数据
  
  # 去除标准差为0的列（免疫细胞类型），以防止相关性计算错误
  validImmune <- immData[, apply(immData, 2, sd) > 0, drop = FALSE]
  
  # 返回匹配后的表达数据和免疫数据
  return(list(exprData = exprData, immuneData = validImmune))
}

# 定义函数：计算目标基因与每种免疫细胞的Spearman相关性
computeCorrelation <- function(exprData, immuneData, gene) {
  cat("步骤 2.4：计算相关性...\n")  # 进度提示
  corrResults <- data.frame()   # 初始化存储相关性结果的数据框
  
  # 循环遍历每个免疫细胞类型
  for (cellType in colnames(immuneData)) {
    # 判断当前免疫细胞数据是否有足够的变异性
    if (sd(immuneData[, cellType]) == 0) {
      next  # 如果标准差为0，则跳过
    }
    # 将免疫细胞数据和目标基因表达数据转为数值向量
    immuneVec <- as.numeric(immuneData[, cellType])
    geneExpr  <- as.numeric(exprData[, gene])
    # 执行Spearman相关性检验
    testResult <- cor.test(immuneVec, geneExpr, method = "spearman")
    
    # 将结果存入数据框中（冗余步骤：先转换为数据框后再合并）
    tempDF <- data.frame(spec = gene,               # 目标基因名称
                         env  = cellType,           # 免疫细胞类型
                         r    = as.numeric(testResult$estimate),  # 相关系数
                         p    = as.numeric(testResult$p.value))     # P值
    corrResults <- rbind(corrResults, tempDF)  # 合并结果
  }
  
  # 根据P值判断相关性方向，并存入新列（冗余的if判断）
  corrResults$pd <- ifelse(corrResults$p < 0.05,
                           ifelse(corrResults$r > 0, "positive correlation", "negative correlation"),
                           "not significant")
  
  # 将相关系数取绝对值，方便后续分档
  corrResults$r <- abs(corrResults$r)
  
  # 为相关系数添加分档信息
  corrResults <- corrResults %>%
    mutate(rd = cut(r,
                    breaks = c(-Inf, 0.2, 0.4, 0.6, Inf),
                    labels = c("< 0.2", "0.2 - 0.4", "0.4 - 0.6", ">= 0.6")))
  
  # 返回计算完成的相关性结果数据框
  return(corrResults)
}

# =============================== #
# Step 3: 主流程执行数据读取和分析    #
# =============================== #

# 进度条初始化（总共5个步骤）
totalSteps <- 5
pb <- txtProgressBar(min = 0, max = totalSteps, style = 3)  # 创建文本进度条

# 3.1 读取表达数据并筛选目标基因
setTxtProgressBar(pb, 1)         # 更新进度条
cat("主流程 3.1：处理基因表达数据...\n")  # 输出提示信息
exprOut <- readExprData(exprFilePath, targetGeneFile)  # 调用函数读取表达数据
exprMat <- exprOut$exprData         # 获取表达矩阵（样本×基因）
targetGene <- exprOut$gene          # 获取目标基因名称
Sys.sleep(0.2)                     # 冗余延时

# 3.2 读取免疫数据并匹配样本
setTxtProgressBar(pb, 2)         # 更新进度条
cat("主流程 3.2：处理免疫细胞数据...\n")
immuneOut <- readImmuneData(immuneDataPath, exprMat)   # 调用函数读取免疫数据
exprMat <- immuneOut$exprData      # 更新表达数据（匹配后的样本）
immuneMat <- immuneOut$immuneData  # 获取处理后的免疫数据
Sys.sleep(0.2)                     # 冗余延时

# 3.3 检查数据维度是否符合预期（冗余判断）
setTxtProgressBar(pb, 3)         # 更新进度条
if (nrow(exprMat) < 5 || nrow(immuneMat) < 5) {
  warning("样本数量较少，结果可能不稳定！")
}
cat("主流程 3.3：样本数检查完毕。\n")
Sys.sleep(0.2)

# 3.4 计算目标基因与各免疫细胞的相关性
setTxtProgressBar(pb, 4)         # 更新进度条
cat("主流程 3.4：计算相关性...\n")
corrDF <- computeCorrelation(exprMat, immuneMat, targetGene)  # 计算相关性
Sys.sleep(0.2)

# 3.5 将相关性结果输出到CSV文件中
setTxtProgressBar(pb, 5)         # 更新进度条
cat("主流程 3.5：写入相关性结果到CSV文件...\n")
write.csv(corrDF, file = "gene_correlation_results.csv", row.names = FALSE)  # 保存CSV文件
close(pb)                        # 关闭进度条

# =============================== #
# Step 4: 绘图（两种风格均保留）        #
# =============================== #

# ---------- 图形风格A：Spectral配色风格 ---------- #
cat("步骤 4A：绘制图形风格A...\n")  # 提示信息
plotA <- qcorrplot(correlate(immuneMat, method = "spearman"), type = "lower", diag = FALSE) +  # 绘制相关性矩阵
  geom_square(size = 1.0, colour = "#222222") +          # 每个单元格加边框
  geom_couple(aes(colour = pd, size = rd),                # 绘制耦合线（基因与免疫细胞之间的连线）
              data = corrDF, curvature = nice_curvature()) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "Spectral"))) +  # 设置填充色谱（反转Spectral）
  scale_size_manual(values = c(0.5, 1, 2, 3)) +             # 自定义线宽
  scale_colour_manual(values = c("positive correlation" = "#FF1493",   # 设置正相关为亮粉
                                 "negative correlation" = "#00CED1",   # 负相关为亮青
                                 "not significant"  = "#999999")) +     # 非显著为灰色
  guides(size = guide_legend(title = "abs(Cor)", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "P-value", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Cell-cell cor", order = 3)) +
  labs(x = "immune infiltrating cells", y = "immune infiltrating cells", title = "Gene-immune infiltrating cells Correlation") +  # 设置坐标轴标签和标题
  theme_minimal(base_size = 14) +                       # 使用极简主题
  theme(axis.title.x = element_text(size = 14, face = "bold", color = "black"),  # 设置X轴标题格式
        axis.title.y = element_text(size = 14, face = "bold", color = "black"),  # 设置Y轴标题格式
        axis.text.x  = element_text(size = 12, face = "bold", color = "black", angle = 45, hjust = 1),  # X轴文字旋转
        axis.text.y  = element_text(size = 12, face = "bold", color = "black"),  # Y轴文字格式
        panel.grid.major = element_blank(),             # 去除主网格线
        panel.grid.minor = element_blank(),             # 去除次网格线
        plot.background  = element_rect(fill = "white", color = NA),  # 图形背景设置为白色
        plot.title       = element_text(size = 16, face = "bold", hjust = 0.5))  # 图标题居中加粗

# 保存图形风格A到PDF文件
pdf(file = "cor_newStyle.pdf", width = 9, height = 7)  # 打开PDF设备
print(plotA)            # 打印图形到设备
dev.off()               # 关闭PDF设备

# ---------- 图形风格B：自定义主题 + Pastel2配色 ---------- #
cat("步骤 4B：绘制图形风格B...\n")  # 提示信息

# 定义自定义主题函数（theme_cute），用于图形美化
theme_cute <- function() {
  theme_minimal() +  # 使用极简主题
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),  # 图形背景白色
      panel.background = element_rect(fill = "white", colour = NA),  # 面板背景白色
      panel.grid.major = element_line(colour = "#f0f0f0", size = 0.5), # 主网格线设置
      panel.grid.minor = element_line(colour = "#f0f0f0", size = 0.25),# 次网格线设置
      axis.text        = element_text(size = 10, colour = "#555555"),   # 坐标轴文字设置
      axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),# X轴文字竖排
      axis.title       = element_text(size = 12, face = "bold", colour = "#555555"), # 坐标轴标题设置
      plot.title       = element_text(size = 18, face = "bold", colour = "#d35400", hjust = 0.5)  # 图标题设置
    )
}

plotB <- qcorrplot(correlate(immuneMat, method = "spearman"), type = "lower", diag = FALSE) +  # 绘制相关性矩阵
  geom_square() +                                    # 为每个单元格添加边框
  geom_couple(aes(colour = pd, size = rd), data = corrDF, curvature = nice_curvature()) +  # 添加耦合连线
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(9, "Pastel2")),  # 使用Pastel2色系
    name = "Cell-Cell Correlation"                        # 色条名称
  ) +
  scale_size_manual(
    name   = "abs(Cor)",
    values = c("< 0.2" = 0.5, "0.2 - 0.4" = 1, "0.4 - 0.6" = 2, ">= 0.6" = 3),
    labels = c("< 0.2" = "< 0.2", "0.2 - 0.4" = "0.2 - 0.4", "0.4 - 0.6" = "0.4 - 0.6", ">= 0.6" = "≥ 0.6")
  ) +
  scale_colour_manual(
    name   = "p-value",
    values = c("positive correlation" = "#F28C8C",   # 正相关色：亮粉
               "negative correlation" = "#8AB8FF",   # 负相关色：亮蓝
               "not significant"      = "#B2BABB"),  # 非显著色：灰蓝
    labels = c("positive correlation" = "Positive",
               "negative correlation" = "Negative",
               "not significant"      = "Not significant")
  ) +
  guides(
    size   = guide_legend(title = "abs(Cor)", override.aes = list(colour = "grey35"), order = 2),
    colour = guide_legend(title = "p-value", override.aes = list(size = 3), order = 1),
    fill   = guide_colorbar(title = "Cell-Cell Correlation", order = 3)
  ) +
  ggtitle("Single Gene–immune infiltrating cells Correlation") +  # 设置图形标题
  theme_cute() +             # 应用自定义主题
  labs(x = NULL, y = NULL)   # 删除X/Y轴标题

# 保存图形风格B到PDF文件
pdf(file = "cor2.pdf", width = 9, height = 7)  # 打开PDF设备
print(plotB)           # 打印图形到设备
dev.off()              # 关闭PDF设备


# 设置输入文件和工作目录
inputFile <- "gene_correlation_results.csv"


# 读取数据
data <- read.csv(inputFile, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
data$r <- as.numeric(data$r)
data$p <- as.numeric(data$p)

# 1. 渐变色参数
nColors <- 500
myGradient <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(nColors) # 低p蓝，高p红

# 2. 主图点颜色仍然基于数据分布归一化
p_min <- min(data$p, na.rm=TRUE)
p_max <- max(data$p, na.rm=TRUE)

norm_p <- (data$p - p_min) / (p_max - p_min)
col_index <- round(norm_p * (nColors - 1)) + 1
col_index[col_index > nColors] <- nColors  # 防止偶有超限
col_index[col_index < 1] <- 1
data$points.color <- myGradient[col_index]

# 3. 点大小分档
p.cex <- seq(2.5, 5.5, length = 5)
fcex <- function(x){
  x <- abs(x)
  cex <- ifelse(x<0.1, p.cex[1],
                ifelse(x<0.2, p.cex[2],
                       ifelse(x<0.3, p.cex[3],
                              ifelse(x<0.4, p.cex[4], p.cex[5]))))
  return(cex)
}
data$points.cex <- fcex(data$r)

# 4. 排序
data <- data[order(data$r), ]

# 5. 画图
pdf(file = "ELL2_cor.result.pdf", width = 9.5, height = 7)
xlim <- ceiling(max(abs(data$r)) * 10) / 10
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0), nc=2), width = c(8, 2.2), heights = c(1, 2, 1, 2, 1))
par(bg = "white", las = 1, mar = c(5, 18, 2, 4), cex.axis = 1.5, cex.lab = 2)

# 主体lollipop plot
plot(1, type = "n", xlim = c(-xlim, xlim), ylim = c(0.5, nrow(data) + 0.5),
     xlab = "Correlation Coefficient", ylab = "", yaxt = "n", yaxs = "i", axes = FALSE)
rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col = "#F5F5F5", border = NA)
grid(ny = nrow(data), col = "white", lty = 1, lwd = 2)

segments(x0 = 0, y0 = 1:nrow(data), x1 = data$r, y1 = 1:nrow(data), lwd = 4)
points(x = data$r, y = 1:nrow(data), col = data$points.color, pch = 16, cex = data$points.cex)

text(par('usr')[1], 1:nrow(data), data$env, adj = 1, xpd = TRUE, cex = 1.5)
pvalue.text <- ifelse(data$p < 0.001, "<0.001", sprintf("%.03f", data$p))
redcutoff_cor <- 0
redcutoff_pvalue <- 0.05
text(par('usr')[2], 1:nrow(data), pvalue.text, adj = 0, xpd = TRUE,
     col = ifelse(abs(data$r) > redcutoff_cor & data$p < redcutoff_pvalue, "red", "black"), cex = 1.5)
axis(1, tick = FALSE)

# 点大小图例
par(mar = c(0, 4, 3, 4))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("left", legend = c(0.1, 0.2, 0.3, 0.4, 0.5), col = "black", pt.cex = p.cex,
       pch = 16, bty = "n", cex = 2, title = "abs(r)")

# ======= 右下角P值渐变色图注，范围限定为0到1 =======
par(mar = c(0,2,4,6), cex.axis = 1.5, cex.main = 2)
image(
  x = 1,
  y = seq(0, 1, length = nColors),
  z = matrix(seq(0, 1, length = nColors), nrow = 1),
  col = myGradient,
  axes = FALSE,
  xlab = "", ylab = "", main = "P value"
)
axis(4, at = seq(0, 1, length = 5), labels = round(seq(0, 1, length = 5), 2), las = 2, tick = FALSE)
# ======= ==============================

dev.off()

cat("所有步骤执行完毕！\n")  # 输出结束提示

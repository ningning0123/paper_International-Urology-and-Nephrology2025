# 加载limma包，用于生物信息学中的差异表达分析
library(limma)  # 加载limma包

# 定义总步骤数（用于进度条显示）
total_steps <- 8  # 定义总步骤数为8
# 初始化进度条，min为0，max为total_steps，style=3为进度条样式
pb <- txtProgressBar(min = 0, max = total_steps, style = 3)  # 初始化进度条

# ---------------------- 步骤1：设置工作目录 ---------------------- 
setwd("F:\\MPM\\数据\\15lasso回归\\01.提取交集基因的矩阵文件")  # 设置工作目录
Sys.sleep(0.5)  # 暂停0.5秒，模拟处理延迟
setTxtProgressBar(pb, 1)  # 更新进度条到步骤1

# ---------------------- 步骤2：检查输入文件 ----------------------
# 判断表达矩阵数据文件
if (!file.exists("Sample Type Matrix.csv")) {  # 如果文件不存在
  stop("错误:Sample Type Matrix.csv 文件不存在！")  # 则终止程序并报错
} else {  
  cat("Sample Type Matrix.csv文件已找到。\n")  # 否则输出提示信息
}
Sys.sleep(0.5)  # 暂停0.5秒
setTxtProgressBar(pb, 2)  # 更新进度条到步骤2

# ---------------------- 步骤3：读取表达矩阵数据文件 ----------------------
rawData <- read.table("Sample Type Matrix.csv", sep = ",", header = TRUE, check.names = FALSE)  # 读取基因表达矩阵数据
Sys.sleep(0.5)  # 暂停0.5秒
setTxtProgressBar(pb, 3)  # 更新进度条到步骤3

# ---------------------- 步骤4：数据预处理与转换 ----------------------
rawDataMatrix <- as.matrix(rawData)  # 将数据框转换为矩阵格式
rownames(rawDataMatrix) <- rawDataMatrix[, 1]  # 将矩阵第一列设为行名（假定第一列为基因名称）
expData <- rawDataMatrix[, 2:ncol(rawDataMatrix)]  # 提取表达数据，从第二列到最后一列
backupExpData <- expData  # 冗余代码：备份一份原始表达数据（用于后续可能的调试）
dimnamesList <- list(rownames(expData), colnames(expData))  # 定义表达矩阵的行名和列名
numericData <- matrix(as.numeric(as.matrix(expData)), nrow = nrow(expData), dimnames = dimnamesList)  # 将表达数据转换为数值型矩阵
Sys.sleep(0.5)  # 暂停0.5秒
setTxtProgressBar(pb, 4)  # 更新进度条到步骤4

# ---------------------- 步骤5：处理重复值及低表达基因过滤 ----------------------
averagedData <- avereps(numericData)  # 对重复值进行平均处理
filteredData <- averagedData[rowMeans(averagedData) > 0, ]  # 过滤掉行均值为0或以下的基因
# 判断过滤后的数据是否为空
if (nrow(filteredData) == 0) {  # 如果过滤后的数据行数为0
  stop("错误: 过滤后没有剩余数据！")  # 则终止程序并报错
}
Sys.sleep(0.5)  # 暂停0.5秒
setTxtProgressBar(pb, 5)  # 更新进度条到步骤5

# ---------------------- 步骤6：检查并读取基因列表文件 ----------------------
# 判断基因列表文件IntersectionGenes.csv是否存在
if (!file.exists("IntersectionGenes.csv")) {  # 如文件不存在
  stop("错误: IntersectionGenes.csv文件不存在！")  # 则终止程序并报错
} else {  
  cat("IntersectionGenes.csv 文件已找到。\n")  # 否则输出提示信息
}
geneList <- read.table("IntersectionGenes.csv", header = FALSE, check.names = FALSE, sep = "\t")  # 读取基因列表数据
Sys.sleep(0.5)  # 暂停0.5秒
setTxtProgressBar(pb, 6)  # 更新进度条到步骤6

# ---------------------- 步骤7：寻找共同基因并提取表达数据 ----------------------
geneNames <- as.vector(geneList[, 1])  # 提取基因列表中第一列为基因名称
commonGenes <- intersect(geneNames, rownames(filteredData))  # 获取表达矩阵与基因列表中的共同基因
# 判断是否存在共同基因
if (length(commonGenes) == 0) {  # 如果没有共同基因
  stop("错误: 没有找到共同基因！")  # 则终止程序并报错
} else {
  cat("找到", length(commonGenes), "个共同基因。\n")  # 否则输出找到的共同基因数量
}
geneExpression <- filteredData[commonGenes, ]  # 提取共同基因的表达数据
Sys.sleep(0.5)  # 暂停0.5秒
setTxtProgressBar(pb, 7)  # 更新进度条到步骤7

# ---------------------- 步骤8：将结果输出到CSV文件 ----------------------
outputMatrix <- rbind(ID = colnames(geneExpression), geneExpression)  # 在结果矩阵中添加列名作为第一行
write.table(outputMatrix, file = "geneexp.csv", sep = ",", quote = FALSE, col.names = FALSE)  # 将结果写入CSV文件
Sys.sleep(0.5)  # 暂停0.5秒
setTxtProgressBar(pb, 8)  # 更新进度条到步骤8

# 关闭进度条
close(pb)  # 关闭进度条
cat("\n所有步骤完成，结果已保存至 geneexp.csv 文件。\n")  # 输出完成提示信息

# 加载必要的库
library(dplyr)

# 设置工作目录
setwd("F:\\MPM\\数据\\24miRNA预测lncRNA\\B4GALT5")  # 确保路径中没有转义字符问题

# 读取TSV文件
filtered_results <- read.table("filtered_results.tsv", header = TRUE, sep = "\t")
sponge_scan <- read.table("spongeScan.tsv", header = TRUE, sep = "\t")

# 提取filtered_results的第二列 (miRNA)
miRNA_list <- unique(filtered_results$miRNA)

# 过滤sponge_scan以找到每个miRNA对应的lncRNA
result <- sponge_scan %>%
  filter(miRNAs %in% miRNA_list) %>%
  select(miRNAs, lncRNA)

# 将结果写入一个文本文件
output_file <- "miRNA_lncRNA_output.txt"
write.table(result, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

# 打印确认信息
cat("输出已写入", output_file, "\n")

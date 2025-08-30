setwd("F:\\MPM\\数据\\23mRNA预测miRNA\\B4GALT5")     # 设置工作目录

# 加载必要的包
library(dplyr)
library(readr)

# 读取数据
miRDB <- read_tsv("miRDB.tsv", col_names = TRUE)
TargetScan <- read_tsv("TargetScan.tsv", col_names = TRUE)
miRanda <- read_tsv("miRanda.tsv", col_names = TRUE)

# 打印列名以检查它们是否一致
print(colnames(miRDB))
print(colnames(TargetScan))
print(colnames(miRanda))
 
# 确定每个数据集的列名
# 假设每个数据库文件的格式如下：
# miRDB: hsa-miR-?, Gene
# TargetScan: hsa-miR-?, Gene
# miRanda: hsa-miR-?, Gene

# 如果列名不一致，可以使用以下代码进行重命名：
colnames(miRDB) <- c("miRNA", "Gene")
colnames(TargetScan) <- c("miRNA", "Gene")
colnames(miRanda) <- c("miRNA", "Gene")

# 指定基因
gene <- "ELL2"

# 添加标志列
miRDB <- miRDB %>% mutate(miRDB = 1)
TargetScan <- TargetScan %>% mutate(TargetScan = 1)
miRanda <- miRanda %>% mutate(miRanda = 1)

# 合并数据
merged_data <- full_join(miRDB, TargetScan, by = c("Gene", "miRNA")) %>%
  full_join(., miRanda, by = c("Gene", "miRNA"))

# 用0填充NA值
merged_data[is.na(merged_data)] <- 0

# 计算Sum列
merged_data <- merged_data %>% mutate(Sum = miRDB + TargetScan + miRanda)

# 过滤指定基因的数据并只保留在所有数据库中都存在的条目
filtered_data <- merged_data %>% filter(Gene == gene & miRDB == 1 & TargetScan == 1 & miRanda == 1)

# 重新排列列顺序，第一列是Gene，第二列是miRNA
final_data <- filtered_data %>% select(Gene, miRNA, miRanda, miRDB, TargetScan, Sum)

# 输出结果
write_tsv(final_data, "filtered_results.tsv")

# 显示前几行结果
head(final_data)

# 如果尚未安装，请取消注释以下安装包的代码
# install.packages("reshape2")
# install.packages("ggplot2")
# install.packages("RColorBrewer")
# install.packages("ggpubr")
# install.packages("dplyr")
# install.packages("broom")

# ------------------ 加载所需包 ------------------
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(broom)

# ------------------ 参数设置 ------------------
inputFile <- "ssGSEA_Scores.csv"
outdir <- "F:\\MPM\\数据\\28.ssGSEA评分可视化1.热图"
setwd(outdir)
gap_label <- "gap" # gap插入名
pdf_outfile <- "barplot_two_lines_set3_gap.pdf"
csv_outfile <- "barplot_data_long.csv"

# ------------------ 数据准备 ------------------
# 假设Score.csv行为样本、列为免疫细胞，行为SampleID
rt <- read.table(inputFile, header = TRUE, sep = ",", 
                 check.names = FALSE, row.names = 1)
rt <- t(rt) # 需转置：行为细胞，列为样本

# 组别判断：行名含 _con 为"Control"，_tre 为"Treat"
sample_names <- rownames(rt)
Group <- ifelse(grepl("_con$", sample_names, ignore.case = TRUE), "Control",
                ifelse(grepl("_tre$", sample_names, ignore.case = TRUE), "Treat", NA))
rt <- as.data.frame(rt)
rt$Group <- Group

# 获取各组样本顺序
control_samples <- sample_names[Group == "Control"]
treat_samples   <- sample_names[Group == "Treat"]
# 可扩展至N组——此处手写两组
all_samples_ordered <- c(control_samples, gap_label, treat_samples)

# ------------------ 数据转换（长格式） ------------------
rt$Sample <- rownames(rt)
data_long <- melt(rt, id.vars = c("Sample", "Group"), 
                  variable.name = "Immune", 
                  value.name = "Fraction")
data_long$Sample <- factor(data_long$Sample, levels = all_samples_ordered)

# 导出长格式数据表
write.csv(data_long, file = csv_outfile, row.names = FALSE)
cat("Barplot 长格式数据已保存为：", csv_outfile, "\n")

# ------------------ 自适应配色 ------------------
immune_types <- unique(data_long$Immune)
nColors <- length(immune_types)
myColors <- colorRampPalette(brewer.pal(12, "Set3"))(nColors)

# ------------------ 条形图、gap分割及分组注释 ------------------
control_count <- length(control_samples)
treat_count   <- length(treat_samples)

p <- ggplot(data_long, aes(x = Sample, y = Fraction, fill = Immune)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = myColors) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal(base_size = 18) +
  theme(
    text = element_text(face = "bold"),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = NULL, 
    y = "Relative Percent", 
    fill = "Immune\nCell Type",
    title = "Immune Cell Distribution",
    subtitle = "Control vs. Treat"
  ) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.05)))

# 分组线与标签
p_annot <- p +
  # 控制组
  annotate("segment", x = 0.5, xend = control_count + 0.5, 
           y = -0.04, yend = -0.04, color = "#D65DB1", size = 5) +
  annotate("text", x = (control_count)/2 + 0.5, y = -0.08, 
           label = "Control", color = "#D65DB1", size = 7, fontface = "bold") +
  # 实验组
  annotate("segment", x = control_count + 1.5, 
           xend = control_count + treat_count + 1.5, 
           y = -0.04, yend = -0.04, color = "#0089BA", size = 5) +
  annotate("text", x = control_count + (treat_count)/2 + 1.5, y = -0.08, 
           label = "Treat", color = "#0089BA", size = 7, fontface = "bold")

# 导出PDF图片
ggsave(pdf_outfile, p_annot, width = 12, height = 7.5)
cat("堆叠条形图已保存为：", pdf_outfile, "\n")

# ------------------ （可选）查看图
# print(p_annot)

cat("全部完成！PDF图片和CSV长表都已保存。\n")

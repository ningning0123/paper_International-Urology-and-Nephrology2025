# 如果尚未安装，请取消注释以下安装包的代码
# install.packages("reshape2")
# install.packages("ggplot2")
# install.packages("RColorBrewer")
# install.packages("ggpubr")
# install.packages("dplyr")
# install.packages("broom")
# install.packages("ggsci")
# install.packages("corrplot")
# install.packages("ggridges")

# ------------------ 加载所需包 ------------------
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(broom)
library(ggsci)
library(corrplot)
library(ggridges)

# ------------------ 设置路径和参数 ------------------
setwd("F:\\MPM\\数据\\29.ssGSEA评分可视化2.差异")
inputFile <- "ssGSEA_Scores.csv"

# ------------------ 数据读取与转置 + 分组赋值 ------------------
rt <- read.table(inputFile, header = TRUE, sep = ",", check.names = FALSE, row.names = 1)
rt <- t(rt)
sample_names <- rownames(rt)
rt <- as.data.frame(rt)
rt$Group <- ifelse(grepl("_con$", sample_names, ignore.case = TRUE), "Control",
                   ifelse(grepl("_tre$", sample_names, ignore.case = TRUE), "Treat", NA))
rt$Sample <- sample_names

# 控制组与实验组样本名排序
control_samples <- sample_names[rt$Group == "Control"]
treat_samples   <- sample_names[rt$Group == "Treat"]
all_samples_ordered <- c(control_samples, "gap", treat_samples)

# ------------------ 数据转换（长表） ------------------
data_long <- melt(rt, id.vars = c("Sample", "Group"), 
                  variable.name = "Immune", value.name = "Fraction")
data_long$Sample <- factor(data_long$Sample, levels = all_samples_ordered)

# 导出长格式表
write.csv(data_long, "barplot_data_long.csv", row.names = FALSE)
cat("1. Barplot 长格式数据表: barplot_data_long.csv\n")

# ------------------ 全组混合箱线图（含样本点） ------------------
countControl <- sum(data_long$Group == "Control", na.rm=TRUE) / length(unique(data_long$Immune))
countTreat   <- sum(data_long$Group == "Treat", na.rm=TRUE) / length(unique(data_long$Immune))
myCuteColors <- c("Control" = "#FFC0CB", "Treat" = "#87CEFA")

boxplot_cute <- ggboxplot(
  data_long,
  x = "Immune", 
  y = "Fraction",
  fill = "Group",           
  palette = myCuteColors,  
  xlab = "",
  ylab = "Fraction",
  legend.title = "Group",
  notch = FALSE,
  width = 0.8
) +
  stat_compare_means(
    aes(group = Group),
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns")
    )
  ) +
  geom_jitter(shape = 21, color = "black", alpha = 0.7, width = 0.15, size = 1.8) +  # 添加样本点
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 1)
  ) +
  labs(
    title = "Immune Cell Comparison",
    subtitle = paste0("(Control n=", countControl, ", Treat n=", countTreat, ")")
  )

ggsave("immune_diff-points_n.pdf", boxplot_cute, width = 8, height = 6)
cat("2. 全组免疫细胞箱线图: immune_diff-points_n.pdf\n")

# ------------------ 汇总统计表输出 ------------------
summary_table <- data_long %>%
  filter(!is.na(Group)) %>%
  group_by(Immune, Group) %>%
  summarise(
    MeanFraction = mean(Fraction, na.rm = TRUE),
    MedianFraction = median(Fraction, na.rm = TRUE),
    SD = sd(Fraction, na.rm = TRUE),
    Count = n(),
    .groups = "drop"
  )
write.csv(summary_table, "boxplot_summary_table.csv", row.names = FALSE)
cat("3. 免疫细胞箱线图统计表: boxplot_summary_table.csv\n")

# ------------------ Wilcox检验p值表兑现 ------------------
pvalue_table <- data_long %>%
  filter(!is.na(Group)) %>%
  group_by(Immune) %>%
  do(tidy(wilcox.test(Fraction ~ Group, data = .))) %>%
  select(Immune, p.value)
write.csv(pvalue_table, "immune_pvalues.csv", row.names = FALSE)
cat("4. Wilcox检验p值表: immune_pvalues.csv\n")

# ------------------ 按细胞类型单独输出SCI风箱线图+统计 ------------------
output_folder <- "immune_boxplots_SCI"
if(!dir.exists(output_folder)) dir.create(output_folder)
immune_types <- unique(data_long$Immune)
for (immune_cell in immune_types) {
  data_subset <- subset(data_long, Immune == immune_cell & !is.na(Group))
  y_pos <- max(data_subset$Fraction, na.rm = TRUE) * 1.1
  if (y_pos == 0) y_pos <- 0.001
  p <- ggplot(data_subset, aes(x = Group, y = Fraction, fill = Group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, color = "black", alpha = 0.8) +
    geom_jitter(shape = 21, color = "black", alpha = 0.7, width = 0.15, size = 2) +
    scale_fill_npg() +
    stat_compare_means(aes(group = Group), label = "p.format", method = "wilcox.test",
                       label.x = 1.5, label.y = y_pos) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
      plot.subtitle = element_text(face = "italic", size = 14, hjust = 0.5, color = "gray30"),
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 14, color = "black"),
      axis.line = element_line(color = "black"),
      legend.position = "none"
    ) +
    labs(
      title = paste("Immune Cell:", immune_cell),
      subtitle = paste0("Control n=", sum(data_subset$Group == "Control"),
                        ", Treat n=", sum(data_subset$Group == "Treat")),
      x = NULL,
      y = "Fraction"
    )
  ggsave(file.path(output_folder, paste0(immune_cell, "_boxplot_SCI.pdf")), p, width = 6, height = 8)
  cat("5. 单细胞PDF箱线图输出：", file.path(output_folder, paste0(immune_cell, "_boxplot_SCI.pdf")), "\n")
}

# ------------------ 相关性热图输出（PDF） ------------------
immuno_mat <- rt[, !(colnames(rt) %in% c("Group", "Sample")), drop=FALSE]
cor_mat <- cor(immuno_mat, method = "spearman", use = "pairwise.complete.obs")
pdf("immune_corrplot.pdf", width = 11, height = 11)
corrplot(cor_mat, method = "circle", type = "upper", tl.col = "black",
         tl.srt = 45, addCoef.col = "black", number.cex = 0.7)
dev.off()
cat("6. 相关性热图 PDF 已保存：immune_corrplot.pdf\n")

# ------------------ 叠加山脊图输出（PDF） ------------------
p_ridges <- ggplot(data_long, aes(x = Fraction, y = Immune, fill = Group, color = Group)) +
  geom_density_ridges(alpha = 0.5, position = "identity", scale = 0.8, size = 0.8) +
  scale_fill_npg() + scale_color_npg() +
  theme_minimal(base_size = 15) +
  labs(title = "Ridgeline Plot of Immune Cells (Overlayed)",
       x = "Fraction", y = "Immune Cell")
ggsave("immune_ridges_overlay.pdf", p_ridges, width = 10, height = 8)
cat("7. Overlayed Ridgeline Plot PDF 输出：immune_ridges_overlay.pdf\n")

# ------------------ Facet分面箱线+样本点 ------------------
p_box_with_points <- ggplot(data_long, aes(x = Group, y = Fraction, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(aes(color = Group), 
              position = position_jitter(width = 0.2), 
              size = 2, alpha = 0.7) +
  stat_compare_means(label = "p.format", method = "wilcox.test", 
                     label.y.npc = 0.8) +
  scale_fill_npg() + scale_color_npg() +
  theme_classic(base_size = 14) +
  facet_wrap(~ Immune, scales = "free") +
  labs(title = "Boxplot with Sample Distribution Points",
       x = "Group", y = "Fraction")
ggsave("boxplot_with_points_facet.pdf", p_box_with_points, width = 13, height = 15)
cat("8. Facet分面箱线图 PDF 已保存：boxplot_with_points_facet.pdf\n")

cat("全部可视化与统计表格已输出（PDF/CSV），请查收结果文件。\n")

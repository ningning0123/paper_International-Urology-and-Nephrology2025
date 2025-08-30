#############################
### Step 0: 导入R包与工具 ###
#############################

cat("Step 0/6: 检查并加载R包......\n")  # 步骤进度提示

# 加载基础包
suppressPackageStartupMessages({
  if(!require(ggplot2)) install.packages("ggplot2")
  if(!require(reshape2)) install.packages("reshape2")
  if(!require(dplyr)) install.packages("dplyr")
  if(!require(progress)) install.packages("progress")
})

library(ggplot2)
library(reshape2)
library(dplyr)
library(progress)

################################################
### Step 1: 读取与检查原始数据 #################
################################################

cat("Step 1/6: 加载原始表格，执行文件检查......\n")

# 指定文件路径
proj_dir <- "F:\\MPM\\数据\\32.分型后ssGSEA差异分析"
setwd(proj_dir)
clust_name <- "geneCluster.txt"
score_name <- "ssGSEA_Scores.csv"

# 检查文件是否存在
if (!file.exists(clust_name)) stop(paste("未找到分型文件:", clust_name))
if (!file.exists(score_name)) stop(paste("未找到打分文件:", score_name))

# 读取数据
labels_df <- read.table(clust_name, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
colnames(labels_df) <- "TypeGroup"   # 防止与原代码冲突
score_df  <- read.table(score_name,  sep=",", header=TRUE, row.names=1, check.names=FALSE)

# 判断数据体积
if(nrow(labels_df)<2 | nrow(score_df)<2) stop("分型文件或打分文件数据量过小。")

############################################
### Step 2: 样本对齐与合并 #################
############################################

cat("Step 2/6: 样本名取交集并合并......\n")

common_samples <- intersect(rownames(labels_df), rownames(score_df))
if(length(common_samples)<2) stop("两个文件交集样本太少，后续分析毫无意义。")
labels_df <- labels_df[common_samples,,drop=FALSE]
score_df  <- score_df[common_samples,,drop=FALSE]

combo_data <- cbind(score_df, labels_df)

########################################################
### Step 3: 数据整理、转换成长表，并基本检查 ###########
########################################################

cat("Step 3/6: 转成长表格式，检查异常值和缺失......\n")

longdat <- melt(combo_data, id.vars="TypeGroup")  # "value"为分数
colnames(longdat) <- c("TypeGroup", "ImmuneType", "ScoreVal")
longdat <- na.omit(longdat)                       # 丢弃缺失

# 强制因子化
longdat$TypeGroup <- as.factor(longdat$TypeGroup)
longdat$ImmuneType <- as.factor(longdat$ImmuneType)

if(nlevels(longdat$TypeGroup)<2) stop("你的分型类别数太少（应为2类以上）！")
if(any(table(longdat$TypeGroup, longdat$ImmuneType)==0)) warning("有的分组与免疫类型组别数量极少或无数据。")

###########################################################
### Step 4: 逐类做差异检验，自动生成星号/输出差异表#########
###########################################################

cat("Step 4/6: 各类免疫类型进行双组检验，生成差异星号表...\n")

# 准备一个差异统计结果表
star_levels <- function(p) {  # P值转星号
  if (is.na(p)) return("")
  if (p <= 0.001) return("***")
  if (p <= 0.01 ) return("**")
  if (p <= 0.05 ) return("*")
  return("")
}
pb <- progress_bar$new(total=length(levels(longdat$ImmuneType)), format="  [:bar] :percent  :message")
res_list <- list()
for(i in levels(longdat$ImmuneType)) {
  curdat <- subset(longdat, ImmuneType==i)
  pval <- tryCatch(wilcox.test(ScoreVal ~ TypeGroup, data=curdat)$p.value, error=function(e) NA)
  mean1 <- mean(curdat$ScoreVal[curdat$TypeGroup==levels(curdat$TypeGroup)[1]])
  mean2 <- mean(curdat$ScoreVal[curdat$TypeGroup==levels(curdat$TypeGroup)[2]])
  res_list[[i]] <- data.frame(
    ImmuneType = i,
    Group1 = levels(curdat$TypeGroup)[1],
    Group2 = levels(curdat$TypeGroup)[2],
    Mean1 = mean1,
    Mean2 = mean2,
    P_value = ifelse(is.na(pval), NA, signif(pval,4)),
    Stars = star_levels(pval),
    stringsAsFactors = FALSE
  )
  pb$tick(tokens = list(message = paste0(i," 完成")))
}
diffstat_tab <- do.call(rbind, res_list)

# 保存比较结果表到本地
write.csv(diffstat_tab, "ImmuneGroup_comparison.csv", row.names=FALSE)
cat("差异分析表已保存为 'ImmuneGroup_comparison.csv'\n")

##########################################################
### Step 5: 配置箱线图颜色，绘图，并加P值星号 ###########
##########################################################

cat("Step 5/6: 绘图并加星号......\n")

# 只用你要求的颜色
mycol <- c("#99CC99", "#FF9966")[1:length(levels(longdat$TypeGroup))]

# 画图
g <- ggplot(longdat, aes(x=ImmuneType, y=ScoreVal, fill=TypeGroup)) +
  geom_boxplot(position=position_dodge(0.8), width=0.72, outlier.size=0.72, outlier.alpha=0.6, lwd=0.5) +  # 箱线设置
  scale_fill_manual(values=mycol, name="Subtype") +       # 手动色板
  # 标星号层：只保留有显著性的星号
  geom_text(
    data = diffstat_tab,
    aes(x = ImmuneType, y = max(longdat$ScoreVal[longdat$ImmuneType == ImmuneType])*1.07, label = Stars, group = NULL),
    inherit.aes=FALSE,
    size=5,
    color="black",
    vjust=0
  ) +
  labs(x="Immune Cell Type", y="Score&Fraction") +           # 轴标题
  theme_bw(base_size=16) +
  theme(
    axis.text.x = element_text(angle=48, hjust=1, vjust=1, size=12),
    axis.text.y = element_text(size=13),
    axis.title = element_text(face="bold", size=18),
    legend.position = "right",
    legend.title = element_text(face="bold", size=14),
    legend.text = element_text(size=13)
  )

#####################################
### Step 6: 输出最终图片文件 #########
#####################################

cat("Step 6/6: 导出PDF文件：boxplot_stars.pdf......\n")
pdf("boxplot_stars.pdf", width=14, height=6)
print(g)
dev.off()
cat("全部流程已完成！\n")

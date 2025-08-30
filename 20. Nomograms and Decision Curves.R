# ==== 基因分析pipeline 一体化完整R脚本 ====

# Step 0: 设置分析主目录
working_directory <- "F:\\MPM\\数据\\19.列线图及校准曲线"
if (!dir.exists(working_directory)) {
  stop(paste("指定的工作目录不存在：", working_directory))
}
setwd(working_directory)

# 包检查与加载
required_pkgs <- c("rms", "regplot", "rmda", "pROC", "tools")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
suppressMessages({
  library(rms)
  library(regplot)
  library(rmda)
  library(pROC)
  library(tools)
})

# Step 1: 进度条与分析步骤
cat("==== 基因分析pipeline启动 ====\n")
steps <- c(
  "创建输出目录", "读取表达数据", "读取重要基因表", "筛选基因",
  "样本结构化", "特征准备", "建模/列线图", "校准曲线", "决策曲线"
)
pb <- txtProgressBar(min=0, max=length(steps), width=40, style=3)

# Step 2: 创建输出文件夹
result_dir <- "Analysis_Result"
if (!dir.exists(result_dir)) dir.create(result_dir)
setwd(result_dir)
setTxtProgressBar(pb, 1)

# Step 3: 读取数据
input_expr <- file.path("..", "geneexp.csv")
gene_file  <- file.path("..", "IntersectionGenes.csv")
if (!file.exists(input_expr)) stop(sprintf("未找到文件: %s", input_expr))
if (!file.exists(gene_file)) stop(sprintf("未找到文件: %s", gene_file))
dat_expr <- tryCatch(
  read.table(input_expr, sep = ",", header = TRUE, row.names = 1, check.names = FALSE), 
  error = function(e) stop("表达文件读取失败")
)
setTxtProgressBar(pb, 2)
gene_tbl <- tryCatch(
  read.table(gene_file, sep = ",", header = TRUE, check.names = FALSE),
  error = function(e) stop("基因列表读取失败")
)
if (ncol(gene_tbl)<1) stop("基因文件无有效列")
feature_genes <- as.character(gene_tbl[,1])
if (length(feature_genes)<2) stop("特征基因数过少")
setTxtProgressBar(pb, 3)

# Step 4: 筛选表达矩阵
dat_expr <- dat_expr[gsub("-", "_", rownames(dat_expr)), , drop=FALSE] # 格式统一
found_genes <- feature_genes %in% rownames(dat_expr)
if (sum(found_genes)==0) stop("IntersectionGenes.csv中基因在表达数据中均找不到")
if (any(!found_genes)) warning(paste0("以下基因不在表达矩阵中：", paste(feature_genes[!found_genes],collapse=",")))
dat_expr_filt <- dat_expr[feature_genes[found_genes],,drop=FALSE]
setTxtProgressBar(pb, 4)

# Step 5: 转置+分组变量赋值
df_expr <- as.data.frame(t(dat_expr_filt))
sample_names <- rownames(df_expr)
groups <- gsub("(.*)_([A-Za-z0-9]+)$", "\\2", sample_names)
df_expr$GroupType <- groups
if (length(unique(groups)) < 2) warning("分组数不足二，分析或有问题")
setTxtProgressBar(pb, 5)

# Step 6: 建模环境准备
ddinfo <- datadist(df_expr)
options(datadist="ddinfo")
setTxtProgressBar(pb, 6)

# Step 7: 组装回归公式
model_vars <- setdiff(colnames(df_expr), "GroupType")
reg_formula <- as.formula(paste("GroupType ~", paste(model_vars, collapse=" + ")))
setTxtProgressBar(pb, 7)

# Step 8: 模型拟合、列线图绘制
lrm_fit <- lrm(reg_formula, data=df_expr, x=TRUE, y=TRUE)
nomo_obj <- nomogram(
  lrm_fit, fun = plogis,
  fun.at = c(0.001,0.1,0.3,0.5,0.7,0.9,0.99), lp=FALSE, funlabel="Disease Risk"
)
pdf("Nomogram_Plot.pdf", width=11, height=6)
plot(nomo_obj)
dev.off()
cat(">> 列线图输出至: Nomogram_Plot.pdf\n")
setTxtProgressBar(pb, 8)

# Step 9: 校准曲线
calibrate_obj <- calibrate(lrm_fit, method="boot", B=800)
pdf("Calibration_Curve.pdf", width=5.5, height=5.5)
plot(calibrate_obj, xlab="Predicted", ylab="Observed", sub=FALSE)
dev.off()
cat(">> 校准曲线输出至: Calibration_Curve.pdf\n")

# ==== Step 10: 决策曲线分析 (DCA) ====
df_expr$GroupType <- as.numeric(factor(df_expr$GroupType)) - 1
if(!all(df_expr$GroupType %in% c(0, 1))) stop("GroupType不是严格的0/1!")
set.seed(123)
dca_obj <- decision_curve(
  formula = reg_formula,
  data = df_expr,
  thresholds = seq(0, 1, by = 0.01),
  family = binomial(link = "logit"),
  bootstraps = 100
)
if(is.null(dca_obj$derived.data)) stop("dca_obj生成失败，请检查数据！")
pdf(file = "DCA.pdf", width = 5.5, height = 5.5)
plot_decision_curve(
  dca_obj,
  xlab = "Threshold Probability",
  col = "orange",
  confidence.intervals = TRUE,    # TRUE显示置信区间
  standardize = TRUE,
  cost.benefit.axis = TRUE
)
dev.off()

cat(">> DCA决策曲线输出至 DCA.pdf\n")




# Step 11: 输出主要分析表格
write.csv(df_expr, file="Filtered_Expression_Matrix.csv", quote=F)
write.csv(as.data.frame(dca_obj$derived.data), file="Decision_Curve_Data.csv")
pred_probs <- predict(lrm_fit, newdata=df_expr, type="fitted")
df_out <- data.frame(Sample=rownames(df_expr), GroupType=df_expr$GroupType, Predicted_Prob=pred_probs)
write.csv(df_out, file="Sample_Predicted_Probabilities.csv")

# Step 12: 模型系数/OR/CI结果表
coef_df <- as.data.frame(summary(lrm_fit))
if("Effect" %in% names(coef_df) & "S.E." %in% names(coef_df)) {
  coef_df$Effect <- as.numeric(as.character(coef_df$Effect))
  coef_df$`S.E.` <- as.numeric(as.character(coef_df$`S.E.`))
  coef_df$OR <- exp(coef_df$Effect)
  coef_df$OR_low <- exp(coef_df$Effect - 1.96 * coef_df$`S.E.`)
  coef_df$OR_high <- exp(coef_df$Effect + 1.96 * coef_df$`S.E.`)
  write.csv(coef_df, file="Model_Coefficients.csv", row.names=FALSE)
  if ("Lower 0.95" %in% names(coef_df) & "Upper 0.95" %in% names(coef_df)) {
    write.csv(coef_df, file="Model_Coefficients_With_CI.csv", row.names=FALSE)
  }
  or_selected <- coef_df[!is.na(coef_df$OR) & coef_df$OR > 1, ]
  write.csv(or_selected, file="Model_OR_GT1.csv", row.names=FALSE)
  if("P" %in% names(coef_df)) {
    sig_coef <- coef_df[!is.na(coef_df$P) & coef_df$P < 0.05, ]
    write.csv(sig_coef, file="Model_Significant_Coefficients.csv", row.names=FALSE)
  }
}

# Step 13：增强标注列线图（高亮中位数虚拟患者）
median_point <- sapply(df_expr[,model_vars,drop=FALSE], median)
median_df <- as.data.frame(t(median_point)); median_df$GroupType <- 1
rownames(median_df) <- "AllMedian"
regplot(
  lrm_fit, showP=TRUE, rank="sd", distribution=TRUE,
  observation=median_df, title="Prediction Nomogram"
)
# ==== 联合模型 ROC（主模型红色，风格与单基因统一）====


roc_y <- as.numeric(df_expr$GroupType)
roc_pred <- pred_probs
roc_obj <- roc(roc_y, roc_pred, levels = c(0,1), direction = "<")

# 主模型/单基因统一色板（前几个色都与后面单基因一致）
my_cols <- c("red", "deepskyblue", "forestgreen", "orange", 
             "purple", "gray40", "black", "magenta", "gold", "brown")
while(length(my_cols) < length(model_vars) + 1) { 
  my_cols <- c(my_cols, rainbow(length(model_vars) + 1 - length(my_cols)))
}

## ====== 联合主模型ROC ======
pdf("CombinedGenes_ROC_SCI.pdf", width = 6, height = 6)
par(mar = c(5,6,4,2)+0.1, cex = 1.3)
plot(1-roc_obj$specificities, roc_obj$sensitivities, type = "l",
     col = my_cols[1], lwd = 4, lty = 1,
     xlab = expression("1 - Specificity"), ylab = "Sensitivity",
     main = " Model ROC Curve", 
     cex.lab = 1.4, cex.axis = 1.15, cex.main = 1.45,
     xlim = c(0,1), ylim = c(0,1))
abline(0, 1, lty = 2, col = "gray70", lwd = 2)
auc_val <- as.numeric(auc(roc_obj))
legend("bottomright", legend = sprintf("AUC = %.3f", auc_val), col = my_cols[1], 
       lwd = 4, lty=1, bty = "n", cex = 1.2)
dev.off()
cat(">> 合并主模型 ROC 曲线 (红色) 已输出至 CombinedGenes_ROC_SCI.pdf\n")

## ====== 单基因多色ROC ======
roc_list <- list(); auc_list <- c(); leglab <- c()
pdf("IndividualGenes_ROC_SCI.pdf", width = 6, height = 6)
par(mar = c(5, 6, 4, 2)+0.1, cex = 1.3)
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
     xlab = expression("1 - Specificity"), ylab = "Sensitivity",
     main = " Gene ROC Curves", 
     cex.lab = 1.4, cex.axis = 1.15, cex.main = 1.45
)
abline(0, 1, lty = 2, col = "gray70", lwd = 2)

for(i in seq_along(model_vars)) {
  g <- model_vars[i]
  cur_roc <- roc(roc_y, df_expr[[g]], levels = c(0,1), direction = "<")
  # 横轴: 1-specificity，纵轴: sensitivity
  lines(1-cur_roc$specificities, cur_roc$sensitivities, col = my_cols[i+1], lwd = 3)
  roc_list[[g]] <- cur_roc
  cur_auc <- as.numeric(auc(cur_roc))
  auc_list <- c(auc_list, cur_auc)
  leglab <- c(leglab, sprintf("%s   AUC=%.3f", g, cur_auc))
}
legend("bottomright", legend = leglab, col = my_cols[2:(length(model_vars)+1)], 
       lwd = 3, lty = 1, bty = "n", cex = 1.15, y.intersp = 1.18)
dev.off()

# AUC同步输出
gene_auc_df <- data.frame(Gene = model_vars, AUC = auc_list)
write.csv(gene_auc_df, file = "IndividualGenes_AUC.csv", row.names = FALSE)
cat(">> 单基因 ROC 曲线和 AUC 已输出至 IndividualGenes_ROC_SCI.pdf / IndividualGenes_AUC.csv\n")


# 收尾
cat("==== 全流程执行完毕！====\n")

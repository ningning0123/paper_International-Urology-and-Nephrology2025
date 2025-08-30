#---------------------#
# 1. 包加载函数 #
# 分析包（Bioconductor）
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
#------------------------------#
# 2. 初始化参数/环境与进度条    #
#------------------------------#
wdPath <- "F:\\MPM\\数据\\13.GO富集分析"
if (!dir.exists(wdPath)) dir.create(wdPath, recursive = TRUE)
setwd(wdPath)

updateProgress <- function(pb, value, message) {
  setTxtProgressBar(pb, value)
  cat(message, "\n")
}
totalSteps <- 9
pb <- txtProgressBar(min = 0, max = totalSteps, style = 3)

pvalThreshold <- 0.05
padjThreshold <- 0.05
colorParameter <- if (padjThreshold > 0.05) "pvalue" else "p.adjust"
updateProgress(pb, 1, "第1步：设置参数和工作环境完成")

#------------------------------#
# 3. 读取基因数据               #
#------------------------------#
geneFile <- "IntersectionGenes.csv"
if (!file.exists(geneFile)) stop("错误：基因文件 IntersectionGenes.csv 不存在")
geneData <- read.table(geneFile, header = FALSE, sep = ",", check.names = FALSE)
updateProgress(pb, 2, "第2步：读取基因数据完成")

geneSymbols <- unique(as.vector(geneData[, 1]))
if (length(geneSymbols) == 0) stop("错误：未找到有效的基因符号")

#------------------------------#
# 4. 基因符号转Entrez ID       #
#------------------------------#
entrezMapping <- mget(geneSymbols, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezMapping)
validGenes <- entrezIDs[!is.na(entrezIDs) & entrezIDs!="NA"]
if (length(validGenes) == 0) stop("错误：无有效的 Entrez 基因 ID")
updateProgress(pb, 3, "第3步：基因符号转换为 EntrezID完成")

#------------------------------#
# 5. GO富集与过滤              #
#------------------------------#
goAnalysis <- enrichGO(
  gene = validGenes,
  OrgDb = org.Hs.eg.db,
  pvalueCutoff = 1, qvalueCutoff = 1,
  ont = "all",
  readable = TRUE
)
goResult <- as.data.frame(goAnalysis)
if (nrow(goResult)==0) stop("警告：未检测到任何富集结果")
filteredGO <- goResult[goResult$pvalue < pvalThreshold & goResult$p.adjust < padjThreshold, ]
updateProgress(pb, 4, "第4步：GO 富集分析及结果过滤完成")

#------------------------------#
# 6. 结果输出                  #
#------------------------------#
outputFile <- "GO_results.txt"
write.table(filteredGO, file = outputFile, sep = "\t", quote = FALSE, row.names = FALSE)
updateProgress(pb, 5, "第5步：富集结果写入文件完成")

#------------------------------#
# 7. 柱状图/气泡图 PDF绘制     #
#------------------------------#
# 柱状图
pdf("GO_barplot.pdf", width=8, height=10)
barPlot <- barplot(
  goAnalysis,
  drop=TRUE, showCategory=10, label_format=50,
  split="ONTOLOGY", color=colorParameter
) +
  facet_grid(ONTOLOGY ~ ., scale = 'free') +
  scale_fill_gradientn(colors = c("#FF6666", "#FFB266", "#FFFF99", "#99FF99", "#6666FF", "#7F52A0", "#B266FF"))
print(barPlot)
dev.off()

# 气泡图
pdf("GO_bubble.pdf", width=8, height=10)
bubblePlot <- dotplot(
  goAnalysis,
  showCategory = 10,
  orderBy = "GeneRatio",
  label_format = 50,
  split = "ONTOLOGY",
  color = colorParameter
) +
  facet_grid(ONTOLOGY ~ ., scale = 'free') +
  scale_color_gradientn(colors = c("#FFB266", "#FFFF99", "#99FF99", "#6666FF", "#7F52A0", "#B266FF"))
print(bubblePlot)
dev.off()
updateProgress(pb, 6, "第6步：柱状图和气泡图生成完成")

#------------------------------#
# 8. 分组条形图 ggbarplot       #
#------------------------------#
topGO <- filteredGO %>% group_by(ONTOLOGY) %>% slice_head(n = 10)
pdf("GO_grouped_barplot.pdf", width=11, height=8)
groupBarPlot <- ggbarplot(
  topGO,
  x="Description", y="Count", fill="ONTOLOGY", color="white",
  xlab="", palette="aaas",
  legend="right", sort.val="desc", sort.by.groups=TRUE,
  position=position_dodge(0.9)
) +
  rotate_x_text(75) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size=10, color="black")) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  geom_text(
    aes(label=Count),
    position=position_dodge(0.9), vjust=-0.3, size=3
  )
print(groupBarPlot)
dev.off()
updateProgress(pb, 7, "第7步：分组条形图生成完成")

#------------------------------#
# 9. chord弦图 PDF绘制         #
#------------------------------#
pdf("GO_chord_diagram.pdf", width=12, height=12)
go <- read.delim("GO_results.txt", header=TRUE, stringsAsFactors=FALSE)
top_terms <- go %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup()
insert_linebreak <- function(text, line_length=35) {
  if(nchar(text) <= line_length) return(text)
  paste(strwrap(text, width=line_length), collapse = "\n")
}
top_terms$Description_new <- sapply(top_terms$Description, insert_linebreak)
mat <- table(
  factor(top_terms$ONTOLOGY, levels=c("BP","CC","MF")),
  factor(top_terms$Description_new, levels=unique(top_terms$Description_new))
)
n_ont <- 3
n_term <- ncol(mat)
gap.deg <- c(rep(1, n_ont-1), 10, rep(1, n_term-1), 10)
grid.col <- c(BP="#E69F00", CC="#56B4E9", MF="#009E73",
              setNames(rep("#BBBBBB", n_term), colnames(mat)))
circos.clear()
circos.par(gap.degree = gap.deg, start.degree = 90)
chordDiagram(
  mat,
  grid.col = grid.col,
  transparency = 0.4,
  annotationTrack = c("", "grid"),
  preAllocateTracks = list(track.height = 0.2)
)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1] + 0.1, sector.name,
              facing="clockwise", niceFacing=TRUE,
              adj=c(0,0.5), cex=0.60)
}, bg.border=NA)
title("GO Ontology")
circos.clear()
dev.off()
updateProgress(pb, 8, "第8步：chord弦图生成完成")

#------------------------------#
# 10. 关闭进度条               #
#------------------------------#
close(pb)
cat('全部分析及绘图完成！\n')

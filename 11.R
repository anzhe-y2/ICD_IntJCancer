library(limma)
library(pheatmap)

expFile <- "NRGLncExp.txt"
logFCfilter <- 1.5
fdrFilter <- 0.05
setwd("D:\\桌面文件夹\\3\\11 diff lncRNA expression")

rt <- read.table(expFile, header = T, sep = "\t", check.names = F)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)
data <- data[rowMeans(data) > 0, ]

group <- sapply(strsplit(colnames(data), "\\-"), "[", 4)
group <- sapply(strsplit(group, ""), "[", 1)
group <- gsub("2", "1", group)
conNum <- length(group[group == 1])
treatNum <- length(group[group == 0])
grade <- c(rep(1, conNum), rep(2, treatNum))

outTab <- data.frame()
for (i in row.names(data)) {
    geneName <- unlist(strsplit(i, "\\|", ))[1]
    geneName <- gsub("\\/", "_", geneName)
    rt <- rbind(expression = data[i, ], grade = grade)
    rt <- as.matrix(t(rt))
    wilcoxTest <- wilcox.test(expression ~ grade, data = rt)
    pvalue <- wilcoxTest$p.value
    conGeneMeans <- mean(data[i, 1:conNum])
    treatGeneMeans <- mean(data[i, (conNum + 1):ncol(data)])
    logFC <- log2(treatGeneMeans) - log2(conGeneMeans)
    conMed <- median(data[i, 1:conNum])
    treatMed <- median(data[i, (conNum + 1):ncol(data)])
    diffMed <- treatMed - conMed
    if (((logFC > 0) & (diffMed > 0)) | ((logFC < 0) & (diffMed < 0))) {
        outTab <- rbind(outTab, cbind(gene = i, conMean = conGeneMeans, treatMean = treatGeneMeans, logFC = logFC, pValue = pvalue))
    }
}
pValue <- outTab[, "pValue"]
fdr <- p.adjust(as.numeric(as.vector(pValue)), method = "fdr")
outTab <- cbind(outTab, fdr = fdr)

write.table(outTab, file = "all.xls", sep = "\t", row.names = F, quote = F)

outDiff <- outTab[(abs(as.numeric(as.vector(outTab$logFC))) > logFCfilter & as.numeric(as.vector(outTab$fdr)) < fdrFilter), ]
write.table(outDiff, file = "diff.xls", sep = "\t", row.names = F, quote = F)

heatmap <- rbind(ID = colnames(data[as.vector(outDiff[, 1]), ]), data[as.vector(outDiff[, 1]), ])
write.table(heatmap, file = "diffLncExp.txt", sep = "\t", col.names = F, quote = F)

geneNum <- 50
outDiff <- outDiff[order(as.numeric(as.vector(outDiff$logFC))), ]
diffGeneName <- as.vector(outDiff[, 1])
diffLength <- length(diffGeneName)
hmGene <- c()
if (diffLength > (2 * geneNum)) {
    hmGene <- diffGeneName[c(1:geneNum, (diffLength - geneNum + 1):diffLength)]
} else {
    hmGene <- diffGeneName
}
hmExp <- log2(data[hmGene, ] + 0.01)
Type <- c(rep("Normal", conNum), rep("LUAD", treatNum))
names(Type) <- colnames(data)
Type <- as.data.frame(Type)
pdf(file = "heatmap.pdf", width = 10, height = 7)
pheatmap(hmExp,
    annotation = Type,
    color = colorRampPalette(c(rep("orange", 3), "white", rep("blue", 3)))(50),
    cluster_cols = F,
    show_colnames = F,
    scale = "row",
    fontsize = 8,
    fontsize_row = 5,
    fontsize_col = 8
)
dev.off()

pdf(file = "vol.pdf", width = 5, height = 5)
xMax <- 6
yMax <- max(-log10(outTab$fdr)) + 1
plot(as.numeric(as.vector(outTab$logFC)), -log10(outTab$fdr),
    xlab = "logFC", ylab = "-log10(fdr)",
    main = "Volcano", ylim = c(0, yMax), xlim = c(-xMax, xMax), yaxs = "i", pch = 20, cex = 0.8
)
diffSub <- subset(outTab, fdr < fdrFilter & as.numeric(as.vector(logFC)) > logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch = 20, col = "orange", cex = 1)
diffSub <- subset(outTab, fdr < fdrFilter & as.numeric(as.vector(logFC)) < (-logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch = 20, col = "blue", cex = 1)
abline(v = 0, lty = 2, lwd = 3)
dev.off()

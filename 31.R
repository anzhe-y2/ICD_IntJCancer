library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
set.seed(12345)

pFilter <- 0.001
expFile <- "symbol.txt"
riskFile <- "risk.all.txt"
setwd("C:\\Users\\Zhang TONGTONG\\Desktop\\3\\31")
allDrugs <- c("Gemcitabine")

rt <- read.table(expFile, header = T, sep = "\t", check.names = F)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)
data <- data[rowMeans(data) > 0.5, ]

group <- sapply(strsplit(colnames(data), "\\-"), "[", 4)
group <- sapply(strsplit(group, ""), "[", 1)
group <- gsub("2", "1", group)
data <- data[, group == 0]
data <- t(data)
rownames(data) <- gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
data <- avereps(data)
data <- t(data)

riskRT <- read.table(riskFile, header = T, sep = "\t", check.names = F, row.names = 1)

for (drug in allDrugs) {
    senstivity <- pRRopheticPredict(data, drug, selection = 1)
    senstivity <- senstivity[senstivity != "NaN"]
    # senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)

    sameSample <- intersect(row.names(riskRT), names(senstivity))
    risk <- riskRT[sameSample, "risk", drop = F]
    senstivity <- senstivity[sameSample]
    rt <- cbind(risk, senstivity)

    rt$risk <- factor(rt$risk, levels = c("low", "high"))
    type <- levels(factor(rt[, "risk"]))
    comp <- combn(type, 2)
    my_comparisons <- list()
    for (i in 1:ncol(comp)) {
        my_comparisons[[i]] <- comp[, i]
    }

    test <- wilcox.test(senstivity ~ risk, data = rt)

    if (test$p.value < pFilter) {
        boxplot <- ggboxplot(rt,
            x = "risk", y = "senstivity", fill = "risk",
            xlab = "Risk",
            ylab = paste0(drug, " senstivity (IC50)"),
            legend.title = "Risk",
            palette = c("#0066FF", "#FF0000")
        ) +
            stat_compare_means(comparisons = my_comparisons)
        pdf(file = paste0("durgSenstivity.", drug, ".pdf"), width = 5, height = 4.5)
        print(boxplot)
        dev.off()
    }
}

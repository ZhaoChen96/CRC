rm(list=ls())
setwd('/data3/zhaochen/project/colon_cancer/colon_chip/htseqCount')
library(DESeq2)
require(BiocParallel)
library(ggplot2)

inputf <- read.table('/data3/zhaochen/project/colon_cancer/colon_chip/htseqCount/chip_readcounts.txt',
                     header = T, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
Times <- c("control","2weeks","4weeks","7weeks","10weeks")
Markers <- c('H3K27ac','H3K4me1','H3K4me3','H3K27me3','H3K9me3')
sample_vector <- c()
for (marker in Markers){
  for (time in Times){
    for (rep in 1:3){
      sample <- paste(time,rep,marker,sep = "-")
      sample_vector <- c(sample_vector,sample)
    }
  }
}
names(data) <- sample_vector
condition <- factor(rep(rep(Times,each=3),levels = Times,5))
histone <- rep(Markers,each=15)
cData <- data.frame(row.names = colnames(data),condition,histone)
cData
dds <- DESeqDataSetFromMatrix(countData = data, colData = cData, design = ~ condition)
dds <- DESeq(dds)
class(dds)

a <- show(dds)
save(a,file = 'chip-seq read count DESeq2 dds.Rdata')
res <- result(dds)
write.table(res,file = 'chip-seq read count DESeq2 res.txt',sep = "\t")

vsd <- getVarianceStabilizedData(dds)
pdf('Correlation heatmap of ChIP-seq Data')
heatmap(cor(vsd),cexCol = 0.75,cexCol = 0.75,margins = c(5,5),)
dev.off()

library(RColorBrewer)
pr <- prcomp(t(vsd))
class(pr)
write.table(pr,file = '')
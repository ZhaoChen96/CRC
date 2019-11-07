library("maSigPro")

# enhancer ---------------------------------------------------------
#construct matrix 
Time = c(rep(0, 9), rep(2, 9), rep(4, 9), rep(7, 9), rep(10, 9)) 
Replicates = rep(1:15, each=3)
Control = rep(c(rep(1,3),rep(0,6)),5)
H3K27ac = rep(c(rep(0,3), rep(1,3),rep(0,3)), 5)
H3K4me1 = rep(c(rep(0,6), rep(1,3)), 5)
CRC.design = cbind(Time,Replicates,Control,H3K27ac,H3K4me1)
sample_vector = c()
for (we in c("ctrl", "2weeks", "4weeks", "7weeks", "10weeks")){
  for (mark in c("Input","H3K27ac", "H3K4me1")){
    for ( rep in c("1", "2" ,"3")){
      sample = paste(paste(we, rep, sep="_"), mark, sep="_")
      sample_vector = c(sample_vector, sample)
    }
  }
}
rownames(CRC.design) = sample_vector
d <- make.design.matrix(CRC.design, degree = 4) #degree 5个时间点,自由度是4
d

df <- read.csv('/home/chenzhao/CRC/Chip_analysis/MaSigPro/chip/enhancer/enhancer_total_readCount.txt',sep='\t',
               row.names = 1,check.names = FALSE)
df <- data.matrix(df)
df <- scale(df,center = TRUE,scale = TRUE) #scale()函数,先中心化:df-均值,再标准化:除以标准差
fit <- p.vector(df, d, Q = 0.05, MT.adjust = "BH", min.obs = 20)
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
save(tstep,  file = "/home/chenzhao/test/enhancer/no_control_enhancer.RData")

get <- get.siggenes(tstep, rsq=0.75, vars="all")
pdf("/home/chenzhao/test/enhancer/new_chip_rsq0.75.pdf")
cluster_result = see.genes(get$sig.genes, k = 9, newX11 = FALSE)
dev.off()

i = 1
clu1 = cluster_result$cut[cluster_result$cut == i]
#clu9 = cluster_result$cut[cluster_result$cut == 9]
#clu15 = cluster_result$cut[cluster_result$cut == 15]
#deg_list = c(names(clu8), names(clu9), names(clu15))

library(org.Mm.eg.db)
library("clusterProfiler")
ego <- enrichGO(gene          = names(clu1),
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

table_name = paste("/home/chenzhao/test/enhancer/chip_", i , "KEGG.txt", sep="")
write.table(ego@result, table_name, sep="\t", quote=F, row.names=F, col.names = T, na="NA", eol="\n")
pdf_name = paste("/home/chenzhao/test/enhancer/pdf_", i , ".pdf", sep="")
pdf(pdf_name, width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()

ego <- enrichKEGG(gene          = names(clu1),
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

table_name = paste("/home/chenzhao/test/enhancer/chip_", i , "KEGG.txt", sep="")
write.table(ego@result, table_name, sep="\t", quote=F, row.names=F, col.names = T, na="NA", eol="\n")
pdf_name = paste("/home/chenzhao/test/enhancer/pdf_", i , "KEGG.pdf", sep="")
pdf(pdf_name, width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()

i = 2
clu1 = cluster_result$cut[cluster_result$cut == i]
#clu9 = cluster_result$cut[cluster_result$cut == 9]
#clu15 = cluster_result$cut[cluster_result$cut == 15]
#deg_list = c(names(clu8), names(clu9), names(clu15))

library(org.Mm.eg.db)
library("clusterProfiler")
ego <- enrichGO(gene          = names(clu1),
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

table_name = paste("/home/chenzhao/test/enhancer/chip_", i , ".txt", sep="")
write.table(ego@result, table_name, sep="\t", quote=F, row.names=F, col.names = T, na="NA", eol="\n")
pdf_name = paste("/home/chenzhao/test/enhancer/pdf_", i , ".pdf", sep="")
pdf(pdf_name, width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()

ego <- enrichKEGG(gene          = names(clu1),
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

table_name = paste("/home/chenzhao/test/enhancer/chip_", i , "KEGG.txt", sep="")
write.table(ego@result, table_name, sep="\t", quote=F, row.names=F, col.names = T, na="NA", eol="\n")
pdf_name = paste("/home/chenzhao/test/enhancer/pdf_", i , ".pdf", sep="")
pdf(pdf_name, width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()

i = 3
clu1 = cluster_result$cut[cluster_result$cut == i]
#clu9 = cluster_result$cut[cluster_result$cut == 9]
#clu15 = cluster_result$cut[cluster_result$cut == 15]
#deg_list = c(names(clu8), names(clu9), names(clu15))

library(org.Mm.eg.db)
library("clusterProfiler")
ego <- enrichGO(gene          = names(clu1),
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

table_name = paste("/home/chenzhao/test/enhancer/chip_", i , ".txt", sep="")
write.table(ego@result, table_name, sep="\t", quote=F, row.names=F, col.names = T, na="NA", eol="\n")
pdf_name = paste("/home/chenzhao/test/enhancer/pdf_", i , ".pdf", sep="")
pdf(pdf_name, width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()


i = 4
clu1 = cluster_result$cut[cluster_result$cut == i]
#clu9 = cluster_result$cut[cluster_result$cut == 9]
#clu15 = cluster_result$cut[cluster_result$cut == 15]
#deg_list = c(names(clu8), names(clu9), names(clu15))

library(org.Mm.eg.db)
library("clusterProfiler")
ego <- enrichGO(gene          = names(clu1),
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

table_name = paste("/home/chenzhao/test/enhancer/chip_", i , ".txt", sep="")
write.table(ego@result, table_name, sep="\t", quote=F, row.names=F, col.names = T, na="NA", eol="\n")
pdf_name = paste("/home/chenzhao/test/enhancer/pdf_", i , ".pdf", sep="")
pdf(pdf_name, width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()


i = 5
clu1 = cluster_result$cut[cluster_result$cut == i]
#clu9 = cluster_result$cut[cluster_result$cut == 9]
#clu15 = cluster_result$cut[cluster_result$cut == 15]
#deg_list = c(names(clu8), names(clu9), names(clu15))

library(org.Mm.eg.db)
library("clusterProfiler")
ego <- enrichGO(gene          = names(clu1),
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

table_name = paste("/home/chenzhao/test/enhancer/chip_", i , ".txt", sep="")
write.table(ego@result, table_name, sep="\t", quote=F, row.names=F, col.names = T, na="NA", eol="\n")
pdf_name = paste("/home/chenzhao/test/enhancer/pdf_", i , ".pdf", sep="")
pdf(pdf_name, width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()


i = 6
clu1 = cluster_result$cut[cluster_result$cut == i]
#clu9 = cluster_result$cut[cluster_result$cut == 9]
#clu15 = cluster_result$cut[cluster_result$cut == 15]
#deg_list = c(names(clu8), names(clu9), names(clu15))

library(org.Mm.eg.db)
library("clusterProfiler")
ego <- enrichGO(gene          = names(clu1),
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

table_name = paste("/home/chenzhao/test/enhancer/chip_", i , ".txt", sep="")
write.table(ego@result, table_name, sep="\t", quote=F, row.names=F, col.names = T, na="NA", eol="\n")
pdf_name = paste("/home/chenzhao/test/enhancer/pdf_", i , ".pdf", sep="")
pdf(pdf_name, width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()


i = 7
clu1 = cluster_result$cut[cluster_result$cut == i]
#clu9 = cluster_result$cut[cluster_result$cut == 9]
#clu15 = cluster_result$cut[cluster_result$cut == 15]
#deg_list = c(names(clu8), names(clu9), names(clu15))

library(org.Mm.eg.db)
library("clusterProfiler")
ego <- enrichGO(gene          = names(clu1),
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

table_name = paste("/home/chenzhao/test/enhancer/chip_", i , ".txt", sep="")
write.table(ego@result, table_name, sep="\t", quote=F, row.names=F, col.names = T, na="NA", eol="\n")
pdf_name = paste("/home/chenzhao/test/enhancer/pdf_", i , ".pdf", sep="")
pdf(pdf_name, width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()


i = 8
clu1 = cluster_result$cut[cluster_result$cut == i]
#clu9 = cluster_result$cut[cluster_result$cut == 9]
#clu15 = cluster_result$cut[cluster_result$cut == 15]
#deg_list = c(names(clu8), names(clu9), names(clu15))

library(org.Mm.eg.db)
library("clusterProfiler")
ego <- enrichGO(gene          = names(clu1),
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

table_name = paste("/home/chenzhao/test/enhancer/chip_", i , ".txt", sep="")
write.table(ego@result, table_name, sep="\t", quote=F, row.names=F, col.names = T, na="NA", eol="\n")
pdf_name = paste("/home/chenzhao/test/enhancer/pdf_", i , ".pdf", sep="")
pdf(pdf_name, width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()


i = 9
clu1 = cluster_result$cut[cluster_result$cut == i]
#clu9 = cluster_result$cut[cluster_result$cut == 9]
#clu15 = cluster_result$cut[cluster_result$cut == 15]
#deg_list = c(names(clu8), names(clu9), names(clu15))

library(org.Mm.eg.db)
library("clusterProfiler")
ego <- enrichGO(gene          = names(clu1),
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

table_name = paste("/home/chenzhao/test/enhancer/chip_", i , ".txt", sep="")
write.table(ego@result, table_name, sep="\t", quote=F, row.names=F, col.names = T, na="NA", eol="\n")
pdf_name = paste("/home/chenzhao/test/enhancer/pdf_", i , ".pdf", sep="")
pdf(pdf_name, width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()
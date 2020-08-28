# chip-seq histone modification maSigPro and RNA-seq maSigPro clusters find overlap gene
# RNA-seq nonfilter cluster1: cancer 1049 gene   cluster8: immunity 414 gene
# RNA-seq filter cluster1: cancer 670 gene

rnaDir <- "/data3/zhaochen/project/colon_cancer/RNA-seq/maSigPro/nonfilter/RNA-seq_masigpro_genelist_"
rna_cluster1 <- read.table(paste(rnaDir,"cluster1.txt",sep = ""),check.names = FALSE,header = TRUE)
rna_cluster8 <- read.table(paste(rnaDir,"cluster8.txt",sep = ""),check.names = FALSE,header = TRUE)

chipDir <- "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro"
H3K27ac_cluster3 <- read.table(paste(chipDir,H3K27ac,H3K27ac_cluster3.txt,sep = "/"),check.names = FALSE,header = TRUE)
H3K27ac_cluster4 <- read.table(paste(chipDir,H3K27ac,H3K27ac_cluster4.txt,sep = "/"),check.names = FALSE,header = TRUE)
H3K27ac_cluster8 <- read.table(paste(chipDir,H3K27ac,H3K27ac_cluster8.txt,sep = "/"),check.names = FALSE,header = TRUE)

a <- intersect(rna_cluster1,H3K27ac_cluster3)
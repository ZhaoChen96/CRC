#this script is used to get specific region related to a gene


rm(list=ls())
library(GenomicFeatures)
library(GenomicRanges)
setwd("/home/zhihl/Project/CRC_2021/step51_annotation")
#make txdb
txdb = makeTxDbFromGFF("/home/zhihl/Project/CRC_2021/reference/gencode.vM20.annotation.gtf", format="gtf")

#1. get TSS upstream and downstream
#get TSS region
get_TSS = function (txdb, upstream, downstream, output_file){ 

PR <- promoters(genes(txdb), upstream=upstream, downstream=downstream)
#get bed dataframe
df <- data.frame(chr=seqnames(PR),
                 start=start(PR),
                 end=end(PR),
                 names=names(PR),
                 scores=c(rep(".", length(PR))),
                 strands=strand(PR))

#remove gene version
df$names = gsub("\\..*","",df$names)
#write dataframe
write.table(df, file=output_file, quote=F, sep="\t", row.names=F, col.names=F)
}

get_TSS(txdb = txdb, upstream = 3000, downstream = 3000, output_file = paste("TSS_","3000", ".bed",sep=""))
get_TSS(txdb = txdb, upstream = 4000, downstream = 4000, output_file = paste("TSS_","4000", ".bed",sep=""))
get_TSS(txdb = txdb, upstream = 10000, downstream = 10000, output_file = paste("TSS_","10000", ".bed",sep=""))


#2. get genebody and upstream and down stream
get_genes_resion = function(txdb, upstream, downstream, output_file){
gene <- genes(txdb)
df <- data.frame(chr=seqnames(gene),
                 start=start(gene) - upstream,
                 end=end(gene) + downstream,
                 names=names(gene),
                 scores=c(rep(".", length(gene))),
                 strands=strand(gene))
df$names = gsub("\\..*","",df$names)
write.table(df, file=output_file, quote=F, sep="\t", row.names=F, col.names=F)
}
get_genes_resion(txdb = txdb, upstream = 2000, downstream = 2000, output_file = paste("gene_","2000", ".bed",sep=""))
get_genes_resion(txdb = txdb, upstream = 4000, downstream = 4000, output_file = paste("gene_","4000", ".bed",sep=""))
get_genes_resion(txdb = txdb, upstream = 10000, downstream = 10000, output_file = paste("gene_","10000", ".bed",sep=""))
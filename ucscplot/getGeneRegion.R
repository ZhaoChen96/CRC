rm(list=ls())
library(GenomicFeatures)
library(GenomicRanges)
library(biomaRt)

# make txdb
txdb <- makeTxDbFromGFF("/data3/zhaochen/reference/mm10/mm10_annotation.gtf",format = "gtf")
genes(txdb)
# get gene body and upstream and down stream
get_gene_region <- function (txdb,upstream,downstream,output_file) {
  gene <- genes(txdb)
  df <- data.frame(chr = seqnames(gene),
                   start = start(gene) - upstream,
                   end = end(gene) + downstream,
                   names = names(gene),
                   score = c(rep(".",length(gene))),
                   strands = strand(gene)
  )
  df$names <- gsub("\\..*","",df$names)
  write.table(df,file = output_file,quote = F,sep = "\t",row.names = F,col.names = F)
}

get_gene_region(txdb = txdb,upstream = 10000,downstream = 10000,output_file = "mm10_genebody_10000.bed")


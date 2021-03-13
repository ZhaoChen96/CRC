rm(list=ls())
library(ggplot2)
library(reshape2)

genecountDir <- "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/genecount/"

plotSaccter <- function (filename,rnafile,marker,region,time,cols) {
  # chip-seq rpm
  file <- read.delim(filename,sep = "\t",check.names=FALSE)

  # rna-seq fpkm
  rna <- read.delim(rnafile,sep = ",",check.names=FALSE)

  rna <- rna[,c(22,2:16)]

  #a <- basename(filename)
  #region <- paste(unlist(strsplit(a,split = "_"))[2],unlist(strsplit(a,split = "_"))[3],sep = "_")
  #marker <- unlist(strsplit(a,split = "_"))[4]

  chip_df <- melt(file[,cols],id.vars = "gene_name", variable.name ="rep",value.name = "rpm")
  rna_df <- melt(rna[,cols],id.vars = "gene_name",value.name = "fpkm")

  data <- merge(x = rna_df,y = chip_df)
  #write.table(data,file = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/genecount/chip/H3K27ac/mm10_genebody_0bp_chip_rna.txt",sep = "\t")
  data$fpkm <- log2(data$fpkm + 1)
  data$rpm <- log2(data$rpm)

  correlation <- cor(data$fpkm,data$rpm,method = "spearman")
  correlation

  ggplot(data,aes(y = fpkm,x = rpm)) +
    geom_point(position = position_jitter(width = 0.3,height = 0.06),alpha=0.4,shape=19,size=1.5,colour="grey60")  +
  #  stat_density2d(aes(alpha=..density..),geom = "tile",contour = FALSE) +
    stat_smooth(formula = data$rpm ~ data$fpkm,method = lm,se = FALSE,colour="black") +
    theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) +
    #annotate("text",label=paste("cor:",round(correlation,4),sep = " "),x = 500,y = 1500,parse=TRUE,hjust=0.5,vjust=0.5,size=8) +
    scale_x_continuous(expand = c(0.03,0)) +
    scale_y_continuous(expand = c(0.03,0)) +
    labs(title = paste(time,region,"\ncor:",round(correlation,4),sep = " "),y="RNA-seq log2(FPKM+1)",x=paste(marker ,"log2(RPM)",sep = " "))  +
    theme(plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
          axis.title = element_text(size = 18,colour = "black"),
          axis.text = element_text(size = 18,colour = "black"))

  ggsave(paste("/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/genecount/chip/",marker,"/",marker,"_",time,"_",region,"_log2.png",sep = ""))
}

Markers <- c("H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me3")
rnafiles <- list.files(path = genecountDir,pattern = "merge_rna.csv")

# need change time and cols
for (marker in Markers) {
  for (rnafile in rnafiles) {
    a <- basename(rnafile)
    region <- paste(unlist(strsplit(a,split = "_"))[2],unlist(strsplit(a,split = "_"))[3],sep = "_")
    rnafile <- paste(genecountDir,rnafile,sep = "")
    filename <- paste(genecountDir,"chip/",marker,"/","mm10_",region,"_",marker,"_rpm.txt",sep = "")
    plotSaccter(filename = filename,rnafile = rnafile,marker = marker,region = region,time = "0week",cols=1:4)
  }
}

# for (marker in Markers) {
#   for (rnafile in rnafiles) {
#     a <- basename(rnafile)
#     region <- paste(unlist(strsplit(a,split = "_"))[2],unlist(strsplit(a,split = "_"))[3],sep = "_")
#     rnafile <- paste(genecountDir,rnafile,sep = "")
#     filename <- paste(genecountDir,"chip/",marker,"/","mm10_",region,"_",marker,"_rpm.txt",sep = "")
#     plotSaccter(filename = filename,rnafile = rnafile,marker = marker,region = region,time = "2weeks",cols=c(1,5:7))
#   }
# }
#
# for (marker in Markers) {
#   for (rnafile in rnafiles) {
#     a <- basename(rnafile)
#     region <- paste(unlist(strsplit(a,split = "_"))[2],unlist(strsplit(a,split = "_"))[3],sep = "_")
#     rnafile <- paste(genecountDir,rnafile,sep = "")
#     filename <- paste(genecountDir,"chip/",marker,"/","mm10_",region,"_",marker,"_rpm.txt",sep = "")
#     plotSaccter(filename = filename,rnafile = rnafile,marker = marker,region = region,time = "4weeks",cols=c(1,8:10))
#   }
# }
#
# for (marker in Markers) {
#   for (rnafile in rnafiles) {
#     a <- basename(rnafile)
#     region <- paste(unlist(strsplit(a,split = "_"))[2],unlist(strsplit(a,split = "_"))[3],sep = "_")
#     rnafile <- paste(genecountDir,rnafile,sep = "")
#     filename <- paste(genecountDir,"chip/",marker,"/","mm10_",region,"_",marker,"_rpm.txt",sep = "")
#     plotSaccter(filename = filename,rnafile = rnafile,marker = marker,region = region,time = "7weeks",cols=c(1,11:13))
#   }
# }
#
# for (marker in Markers) {
#   for (rnafile in rnafiles) {
#     a <- basename(rnafile)
#     region <- paste(unlist(strsplit(a,split = "_"))[2],unlist(strsplit(a,split = "_"))[3],sep = "_")
#     rnafile <- paste(genecountDir,rnafile,sep = "")
#     filename <- paste(genecountDir,"chip/",marker,"/","mm10_",region,"_",marker,"_rpm.txt",sep = "")
#     plotSaccter(filename = filename,rnafile = rnafile,marker = marker,region = region,time = "10weeks",cols=c(1,14:16))
#   }
# }
# png("/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/genecount/chip/H3K27ac/H3K27ac_0week_genebody_log2_scatter.png")
# smoothScatter(data$rpm, data$fpkm,pch = 19,
#               transformation = function(x) x ^ 0.5 # Scale
#               )
# dev.off()

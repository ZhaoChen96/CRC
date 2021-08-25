rm(list = ls())
library(Hmisc)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(cowplot)

fpkmfile <- read.table("~/project/colon cancer/RNA-seq/CRC_mouse_rna_genes_fpkm.txt",sep = "\t",header = TRUE,check.names = FALSE)
tf_fpkm <- function(gene,ymax) {
  genedf <- subset(fpkmfile,external_gene_name == gene)
  df <- melt(genedf,id.vars="external_gene_name",variable.name = "time",value.name = "fpkm")
  times <- factor(c("control","2weeks","4weeks","7weeks","10weeks"),levels =c("control","2weeks","4weeks","7weeks","10weeks"))
  timepoint <- rep(times,each=3)
  data <- cbind(df,timepoint)
  p <- ggplot(data,aes(x = factor(timepoint),y = fpkm,size=fpkm)) +
    geom_point(shape=21,colour="white",fill="#1B9E77",alpha=0.6) + #,position = position_jitter(width = 0.25,height = 0)
    theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) +
    labs(title = gene,x="week(s)",y="RNA-seq gene expression level (FPKM)",size=NULL) +
    scale_y_continuous(expand = c(0,0.01),limits = c(0,ymax)) +
    scale_x_discrete(breaks=times,labels=c(0,2,4,7,10)) +
    theme(aspect.ratio = 1/0.618,
          plot.margin = margin(t=2,b=2,l=0.5,r=0.5,unit = "cm"),
          plot.title = element_text(size = 18,face = "bold",margin = margin(b=0,unit = "cm"),hjust = 0.5,vjust = 0.5),
          axis.title.y = element_text(margin = margin(t = 0,b = 0,l = 0,r = 0,unit = "cm")),
          axis.text = element_text(colour = "black",size = 18,margin = margin(t = 0.1,b = 0.1,l = 0.1,r = 0.1,unit = "cm")),
          axis.text.x = element_text(angle = 0,hjust = 0.5)) +
    guides(size=FALSE)
  print(p)
  #ggsave(paste("~/project/colon cancer/RNA-seq/writer_FPKM/",gene,"_FPKM.pdf",sep = ""),width = 6,height = 6) 
}
tf_fpkm(gene = "Otx2",ymax = 0.3)
tf_fpkm(gene = "Axin2",ymax = 120)
tf_fpkm(gene = "Wnt3",ymax = 120)

# AOM/DSS mouse gene FPKM point plot --------------------------------------
rm(list = ls())
library(ggpubr)
library(ggplot2)
library(reshape2)
fpkm <- read.table("~/project/colon cancer/chip-seq/luo_figure/genes.fpkm_table",sep = "\t",header = TRUE,check.names = FALSE)
names(fpkm)[1] <- "ensembl"
gene_detail = read.delim("~/project/colon cancer/RNA-seq/mart convert gene/mart_export_mouse_geneid_symbol.txt", sep="\t", header=T)
symbol_to_ensembl = function(symbol_list){
  ensembl = gene_detail[match(symbol_list, gene_detail$ensembl_gene_id),c("ensembl_gene_id", "external_gene_name")]
  colnames(ensembl) = c("ensembl", "symbol")
  ensembl = unique(na.omit(ensembl))
  rownames(ensembl) <- ensembl[,1]
  #ensembl = ensembl[,2]
  return(ensembl) 
}
library(dplyr)
fpkmfile = symbol_to_ensembl(symbol_list = fpkm$ensembl)
fpkmfile <- inner_join(fpkmfile,fpkm,by="ensembl")

fpkmfile <- read.csv("~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody2000bp/mm10_genebody_2000bp_merge_rna.csv",check.names = FALSE)
fpkmfile <- fpkmfile[,c(2:16,22)]
tf_fpkm <- function(gene,ymax) {
  genedf <- subset(fpkmfile,gene_name == gene)
  #genedf <- genedf[,-1]
  df <- melt(genedf,id.vars="gene_name",variable.name = "time",value.name = "fpkm")
  times <- factor(c("control","2weeks","4weeks","7weeks","10weeks"),levels =c("control","2weeks","4weeks","7weeks","10weeks"))
  timepoint <- rep(times,each=3)
  data <- cbind(df,timepoint)
  p <- ggplot(data,aes(x = factor(timepoint),y = fpkm,size=fpkm)) +
    geom_point(shape=21,colour="white",fill="#1B9E77",alpha=0.6) + #,position = position_jitter(width = 0.25,height = 0)
    theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) +
    labs(title = gene,x="week(s)",y="Expression (FPKM)",size=NULL) +
    scale_y_continuous(expand = c(0,0.01),limits = c(0,ymax)) +
    scale_x_discrete(breaks=times,labels=c(0,2,4,7,10)) +
    theme(aspect.ratio = 1/0.618,
          plot.margin = margin(t=0,b=0,l=0,r=0,unit = "cm"),
          plot.title = element_text(size = 18,face = "bold",margin = margin(b=0.3,unit = "cm"),hjust = 0.5,vjust = 0.5),
          axis.title.y = element_text(margin = margin(t = 0,b = 0,l = 0.1,r = 0.1,unit = "cm"),size = 18),
          axis.text = element_text(colour = "black",size = 18,margin = margin(t = 0.1,b = 0.1,l = 0.2,r = 0.1,unit = "cm")),
          axis.text.x = element_text(angle = 0,hjust = 0.5,size=18)) +
    guides(size=FALSE)
  ggsave(paste("~/project/colon cancer/RNA-seq/writer_FPKM/",gene,"_FPKM.png",sep = ""),p,width = 3,height = 3.7) 
  p
}

tf_fpkm(gene = "Otx2",ymax = 0.3)
# wnt 
tf_fpkm(gene = "Axin2",ymax = 120)
tf_fpkm(gene = "Wnt3",ymax = 3)
tf_fpkm(gene = "Nlk",ymax = 25)
tf_fpkm(gene = "Wif1",ymax = 400)
tf_fpkm(gene = "Lef1",ymax = 20)
tf_fpkm(gene = "Wnt7a",ymax = 1.5)
tf_fpkm(gene = "Lgr6",ymax = 6)
tf_fpkm(gene = "Sox9",ymax = 40)
tf_fpkm(gene = "Sox4",ymax = 60)

# H3K27ac cluster3
tf_fpkm(gene = "Pabpc1",ymax = 800)
tf_fpkm(gene = "Thap6",ymax = 6)
tf_fpkm(gene = "Mpzl3",ymax = 40)
tf_fpkm(gene = "Ywhaq",ymax = 120)
# H3K27ac cluster4
tf_fpkm(gene = "Fzd2",ymax = 15)
tf_fpkm(gene = "Mcc",ymax = 3)
tf_fpkm(gene = "Sqstm1",ymax = 300)
tf_fpkm(gene = "Ptx3",ymax = 1.5)
tf_fpkm(gene = "Nipal1",ymax = 12)
tf_fpkm(gene = "Il1r1",ymax = 12)
tf_fpkm(gene = "Cxcl5",ymax = 800)
tf_fpkm(gene = "Cxcl1",ymax = 80)
tf_fpkm(gene = "Fst",ymax = 15)
tf_fpkm(gene = "Lhfpl2",ymax = 40)
tf_fpkm(gene = "Ptpn22",ymax = 6)
tf_fpkm(gene = "Trim12c",ymax = 15)
tf_fpkm(gene = "Ets2",ymax = 200)
# OTX2 target gene
tf_fpkm(gene = "Ifitm1",ymax = 600)
tf_fpkm(gene = "Ifitm2",ymax = 800)
tf_fpkm(gene = "Ifitm3",ymax = 1200)
tf_fpkm(gene = "Ifitm5",ymax = 2)
tf_fpkm(gene = "Ifitm6",ymax = 20)
tf_fpkm(gene = "Mta3",ymax = 30)
tf_fpkm(gene = "Clcn3",ymax = 50)

# IFITM family FPKM
ifitm <- c("Ifitm1","Ifitm2","Ifitm3","Ifitm5","Ifitm6")
fpkmfile <- fpkmfile[,-1]
data <- fpkmfile[fpkmfile$symbol %in% ifitm,]
df <- melt(data = data,id.vars = "symbol",variable.name = "sample",value.name = "fpkm")
times <- factor(c("control","2weeks","4weeks","7weeks","10weeks"),levels =c("control","2weeks","4weeks","7weeks","10weeks"))
df$time <- rep(times,each=15)
df <- na.omit(df)

f <- log2(df$fpkm + 1)
summary(f[61:75])

library(ggpubr)
mycomparsion <- list(c("7weeks","10weeks"),c("4weeks","10weeks"),c("2weeks","10weeks"),c("control","10weeks"))
ggplot(data = df,aes(x = factor(time),y = log2(fpkm + 1),fill=log2(fpkm+1))) +
  geom_point(shape=21,colour="white",fill="#1B9E77",alpha=0.6,size=4) + #,position = position_jitter(width = 0.2)
  theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) +
  stat_compare_means(method = "wilcox",paired = FALSE,comparisons = mycomparsion,label = "p.signif",label.y = c(11,12,13,14)) +
  geom_segment(aes(x = 0.75,y = 4.82,xend = 1.25,yend = 4.82),lineend = "square",linejoin = "round",size=1.2) +
  geom_segment(aes(x = 1.75,y = 6.67,xend = 2.25,yend = 6.67),lineend = "square",linejoin = "round",size=1.2) +
  geom_segment(aes(x = 2.75,y = 6.4,xend = 3.25,yend = 6.4),lineend = "square",linejoin = "round",size=1.2) +
  geom_segment(aes(x = 3.75,y = 8.39,xend = 4.25,yend = 8.39),lineend = "square",linejoin = "round",size=1.2) +
  geom_segment(aes(x = 4.75,y = 8.62,xend = 5.25,yend = 8.62),lineend = "square",linejoin = "round",size=1.2) +
  labs(title = "IFITM family",x="week(s)",y="RNA expression log2(FPKM+1)",size=NULL) +
  scale_y_continuous(expand = c(0,0.01),limits = c(0,15)) +
  scale_x_discrete(breaks=times,labels=c(0,2,4,7,10)) +
  theme(aspect.ratio = 1/0.618,
        #plot.margin = margin(t=2,b=2,l=0.5,r=0.5,unit = "cm"),
        plot.title = element_text(size = 18,face = "bold",margin = margin(b=0,unit = "cm"),hjust = 0.5,vjust = 0.5),
        axis.title.y = element_text(margin = margin(t = 0,b = 0,l = 0.1,r = 0.1,unit = "cm")),
        axis.text = element_text(colour = "black",size = 18,margin = margin(t = 0.1,b = 0.1,l = 0.1,r = 0.1,unit = "cm")),
        axis.text.x = element_text(angle = 0,hjust = 0.5)) +
  guides(size=FALSE)
ggsave("~/project/colon cancer/RNA-seq/writer_FPKM/IFITM_family_FPKM.png") #,width = 3,height = 4.2) 


# TCGA boxplot ------------------------------------------------------------
tcga <- read.csv("~/project/colon cancer/TCGA/total_COAD_FPKM.txt",sep = "\t",check.names = FALSE)
library(stringr)
tcga$symbol <- unlist(str_split(tcga$symbol,"[.]",simplify=T))[,1]
tcga$ensembl <- unlist(str_split(tcga$ensembl,"[.]",simplify=T))[,1]
tcga <- tcga[,-3]
cancer <- apply(tcga[,c(46:644)], 1, median)
normal <- apply(tcga[,c(4:45)], 1, median)

library(biomaRt)
mouse <- useDataset("mmusculus_gene_ensembl", mart = useMart("ensembl"))
human <- useMart("ensembl",dataset = "hsapiens_gene_ensembl")
bm <- getLDS(attributes = c("external_gene_name","ensembl_gene_id"),mart = human,filters = "external_gene_name",
             attributesL = c("external_gene_name","ensembl_gene_id"),martL = mouse,values = tcga$symbol)
write.table(x = bm,file = "~/project/colon cancer/TCGA/human_to_mouse_symbol_ensembl.txt",sep = "\t",quote = FALSE,row.names = FALSE)

bm <- read.table("~/project/colon cancer/TCGA/human_to_mouse_symbol_ensembl.txt",sep = "\t",header = TRUE)
tcga_clu <- read.table("~/project/colon cancer/TCGA/COAD_FPKM_add_symbol.txt",sep = "\t")

plotClusterBoxplot <- function(vars_id,limits,breaks,labely) {
  library(dplyr)
  #tcga_clu <- merge(y = tcga,x = bm, by.y= "ensembl",by.x="Gene.stable.ID")
  #tcga_clu <- tcga_clu[,-c(1,3,4,5)]
  #write.table(tcga_clu,"~/project/colon cancer/TCGA/COAD_FPKM_add_symbol.txt",sep = "\t",quote = FALSE)
  gene <- tcga_clu[tcga_clu$Gene.name %in% vars_id,]
  library(reshape2)
  df <- melt(data = gene,id.vars = "Gene.name",variable.name = "sample",value.name = "value")
  df <- cbind(df,type = c(rep("Normal",each=42),rep("Tumor",each=599)))
  df <- na.omit(df)
  
  library(ggpubr)
  library(ggplot2)
  p <- ggplot(data = df, aes(x = factor(type),y = as.numeric(value),fill=factor(type))) +
    stat_boxplot(geom = "errorbar",linetype=1,width=0.3,position = "identity") +
    geom_boxplot(outlier.fill = NA,outlier.shape=NA,width=0.6) +
    stat_compare_means(method = "t.test",paired = FALSE,comparisons = list(c("Normal","Tumor")),label = "p.signif",label.y = c(labely)) +
    theme_bw(base_size = 18,base_family = "sans",base_line_size = 1.1) +
    scale_y_continuous(expand = c(0.02,0.02),limits = limits,breaks = breaks) + 
    scale_x_discrete(expand = c(0.3,0.3),breaks = c("Normal","Tumor"),labels=c("Normal\nn=42","Tumor\nn=599")) +
    scale_fill_manual(values = c("Normal" = "#D9D9D9","Tumor" = "#A50F15")) +
    labs(title = paste(vars_id,"in TCGA",sep = " "),x=NULL,y="Expression (FPKM)",fill=NULL) +
    theme(aspect.ratio = 1/0.75,
      plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,family = "sans"),
      axis.text = element_text(colour = "black",size = 18,family = "sans"),
      axis.text.y = element_text(size = 18,colour = "black",margin = margin(l = 0.1,r = 0.1,unit = "cm")),
      axis.ticks.length = unit(2,'mm'),
      panel.border = element_rect(size = 1.1),
      panel.grid = element_blank()
    ) +
    guides(fill=FALSE)
  ggsave(paste("~/project/colon cancer/TCGA/",vars_id," TCGA FPKM.pdf",sep = ""),p,height = 4.5,width = 3)
  p
}

bm[bm$Gene.name.1 %in% "Trim12c",]

# wnt signal gene
plotClusterBoxplot(vars_id = "AXIN2",limits = c(0,100),breaks = c(0,50,100),labely=90)
plotClusterBoxplot(vars_id = "NLK",limits = c(0,10),breaks = c(0,5,10),labely=9)
plotClusterBoxplot(vars_id = "WIF1",limits = c(0,1),breaks = c(0,0.5,1),labely=0.7)
plotClusterBoxplot(vars_id = "WNT3",limits = c(0,4),breaks = c(0,2,4),labely=3)
plotClusterBoxplot(vars_id = "LEF1",limits = c(0,10),breaks = c(0,5,10),labely=8)
plotClusterBoxplot(vars_id = "WNT7A",limits = c(0,0.1),breaks = c(0,0.1,0.20),labely=0.08)
plotClusterBoxplot(vars_id = "LGR6",limits = c(0,20),breaks = c(0,10,20),labely=15)
plotClusterBoxplot(vars_id = "SOX9",limits = c(0,200),breaks = c(0,100,200),labely=150)
plotClusterBoxplot(vars_id = "SOX4",limits = c(0,100),breaks = c(0,50,100),labely=75)
# H3K27ac cluster3 de novo gene
plotClusterBoxplot(vars_id = "YWHAQ",limits = c(0,210),breaks = c(0,100,200),labely=198)
# nfkb putative gene
plotClusterBoxplot(vars_id = "FZD2",limits = c(0,8),breaks = c(0,4,8),labely=6)
plotClusterBoxplot(vars_id = "MCC",limits = c(0,4),breaks = c(0,2,4),labely=3)
plotClusterBoxplot(vars_id = "SQSTM1",limits = c(0,80),breaks = c(0,40,80),labely=73)
plotClusterBoxplot(vars_id = "LHFPL2",limits = c(0,20),breaks = c(0,10,20),labely=16)
plotClusterBoxplot(vars_id = "PTX3",limits = c(0,6),breaks = c(0,3,6),labely=4)
plotClusterBoxplot(vars_id = "NIPAL1",limits = c(0,20),breaks = c(0,10,20),labely=18)
plotClusterBoxplot(vars_id = "CXCL6",limits = c(0,4),breaks = c(0,2,4),labely=3.5)
plotClusterBoxplot(vars_id = "IL1R1",limits = c(0,20),breaks = c(0,10,20),labely=16)
plotClusterBoxplot(vars_id = "FST",limits = c(0,4),breaks = c(0,2,4),labely=2.7)
plotClusterBoxplot(vars_id = "IFITM1",limits = c(0,1200),breaks = c(0,600,1200),labely=1000)
plotClusterBoxplot(vars_id = "IFITM2",limits = c(0,300),breaks = c(0,100,200,300),labely=280)
plotClusterBoxplot(vars_id = "IFITM3",limits = c(0,1200),breaks = c(0,600,1200),labely=1000)
plotClusterBoxplot(vars_id = "IFITM5",limits = c(0,0.4),breaks = c(0,0.2,0.4),labely=0.3)
plotClusterBoxplot(vars_id = "MTA3",limits = c(0,8),breaks = c(0,4,8),labely=6)
plotClusterBoxplot(vars_id = "CLCN3",limits = c(0,60),breaks = c(0,30,60),labely=45)
plotClusterBoxplot(vars_id = "PTPN22",limits = c(0,6),breaks = c(0,3,6),labely=4)
plotClusterBoxplot(vars_id = "TRIM12c",limits = c(0,15),breaks = c(0,5,10,15),labely=12)
plotClusterBoxplot(vars_id = "TRIM5",limits = c(0,15),breaks = c(0,5,10,15),labely=12)
plotClusterBoxplot(vars_id = "ETS2",limits = c(0,400),breaks = c(0,100,200,300,400),labely=350)
plotClusterBoxplot(vars_id = "ESR1",limits = c(0,0.3),breaks = c(0,0.1,0.2,0.3),labely=0.25)






ifitm <- c("IFITM1","IFITM2","IFITM3","IFITM5","IFITM6")
data <- unique(tcga_clu[tcga_clu$Gene.name %in% ifitm,])
library(reshape2)
df <- melt(data = data,id.vars = "Gene.name",variable.name = "sample",value.name = "value")
df <- cbind(df,type = c(rep("Normal",each=168),rep("Tumor",each=2396)))
df <- na.omit(df)

library(ggpubr)
library(ggplot2)
ggplot(data = df, aes(x = factor(type),y = as.numeric(value),fill=factor(type))) +
  stat_boxplot(geom = "errorbar",linetype=1,width=0.3,position = "identity") +
  geom_boxplot(outlier.fill = NA,outlier.shape=NA,width=0.6) +
  stat_compare_means(method = "t.test",paired = FALSE,comparisons = list(c("Normal","Tumor")),label = "p.signif",label.y = c(730)) +
  theme_bw(base_size = 18,base_family = "sans",base_line_size = 1.1) +
  scale_y_continuous(expand = c(0.02,0.02),limits = c(0,800),breaks = c(0,500,800)) + 
  scale_x_discrete(expand = c(0.3,0.3),breaks = c("Normal","Tumor"),labels=c("Normal\nn=42","Tumor\nn=599")) +
  scale_fill_manual(values = c("Normal" = "#D9D9D9","Tumor" = "#A50F15")) +
  labs(title = "IFITM family in TCGA",x=NULL,y="Expression (FPKM)",fill=NULL) +
  theme(aspect.ratio = 1/0.75,
        plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,family = "sans"),
        axis.text = element_text(colour = "black",size = 18,family = "sans"),
        axis.text.y = element_text(size = 18,colour = "black",margin = margin(l = 0.1,r = 0.1,unit = "cm")),
        axis.ticks.length = unit(2,'mm'),
        panel.border = element_rect(size = 1.1),
        panel.grid = element_blank()
  ) +
  guides(fill=FALSE)
ggsave(paste("~/project/colon cancer/TCGA/","IFITM family"," TCGA FPKM.pdf",sep = ""),height = 4,width = 3.6)

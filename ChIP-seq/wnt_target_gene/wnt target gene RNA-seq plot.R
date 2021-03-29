# wnt target gene
rm(list = ls())
setwd("~/project/colon cancer/RNA-seq/edgeR/masigpro_edgeR_filter/")

# RNA-seq gene expression level -------------------------------------------
cluster1 <- read.delim("RNA-seq_masigpro_BP_cluster1.txt",sep = "\t",check.names = FALSE,header = TRUE)
wnt <- cluster1[c(1,2,4,6,253,422),]
gene <- unique(unlist(strsplit(unlist(wnt$geneID),split = "/")))
write.table(gene,"~/project/colon cancer/wnt target gene/wnt_44_target_gene.txt",quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)


htcount <- read.table("counts_filter.txt",check.names = FALSE)
rownames(htcount) <- substr(rownames(htcount),1,18)
htcount <- as.data.frame(t(apply(htcount,1,scale)))
htcount[,"Week0"] = rowMeans(htcount[,c(1,2,3)]) 
htcount[,"Week2"] = rowMeans(htcount[,c(4,5,6)]) 
htcount[,"Week4"] = rowMeans(htcount[,c(7,8,9)]) 
htcount[,"Week7"] = rowMeans(htcount[,c(10,11,12)])
htcount[,"Week10"] = rowMeans(htcount[,c(13,14,15)]) 

gene_detail = read.delim("~/project/colon cancer/RNA-seq/mart convert gene/mart_export_mouse_geneid_symbol.txt", sep="\t", header=T)
symbol_to_ensembl = function(symbol_list){
  ensembl = gene_detail[match(symbol_list, gene_detail$external_gene_name),c("ensembl_gene_id", "external_gene_name")]
  colnames(ensembl) = c("ensembl", "symbol")
  ensembl = unique(na.omit(ensembl))
  rownames(ensembl) <- ensembl[,1]
  #ensembl = ensembl[,2]
  return(ensembl) 
}

cluster1_count = symbol_to_ensembl(symbol_list = gene)
cluster1_ht <- htcount[rownames(cluster1_count),]
cluster1_ht$ensembl <- rownames(cluster1_ht)
library(dplyr)
data <- inner_join(cluster1_count,cluster1_ht,by="ensembl")
rownames(data) <- data$symbol
df <- as.matrix(data[,-c(1:17)])
colnames(df) <- c("0week","2weeks","4weeks","7weeks","10weeks")

library(ComplexHeatmap)
library(circlize)
Times <- c("0week","2weeks","4weeks","7weeks","10weeks")
top_annotation <- HeatmapAnnotation(
  type = c(rep("Normal",1),rep(c("Inflam","Tumor"),each=2)),
  time = rep(Times,each=1),
  col = list(
    type = c("Normal" = "#DADAEB","Inflam" = "#9E9AC8","Tumor" =  "#6A51A3"),
    time = c("0week" = "#1B9E77","2weeks"="#66A61E","4weeks"="#7570B3","7weeks"="#E7298A","10weeks"="#D95F02")
  ),
  gap = unit(1,'points'),
  simple_anno_size = unit(0.25,"cm"),show_annotation_name = FALSE,show_legend = c(TRUE,TRUE,TRUE),
  annotation_legend_param = list(
    type = list(title = "Type",at = c("Normal","Inflam","Tumor"),labels = c("Normal","Inflam","Tumor")),
    time = list(title = "Time",at = Times,labels = Times)
    ),
  annotation_height = unit(2,"cm"),
  annotation_name_gp = gpar(fontsize = 15,fontfamily="sans")
)

ht <- Heatmap(matrix = df,name = "Z-score",
        cluster_columns = FALSE,cluster_rows = TRUE,show_row_dend = FALSE,
        show_column_names = FALSE,column_names_rot = 45,row_names_side = "right",
        col = colorRamp2(breaks = c(-1,0,1),colors = c("navy","white","firebrick3")),
        top_annotation = top_annotation,
        row_names_gp = gpar(fontsize=12,fontfamily="sans"),
        column_names_gp = gpar(fontsize=16,fontfamily="sans"),
        column_title = "Wnt target gene expression (n=44)",
        column_title_gp = gpar(fontsize=12,fontfamily="sans"),
        height = 5,
        width = 3)

draw(ht,merge_legend=TRUE,heatmap_legend_side="left")
ggsave("~/project/colon cancer/wnt target gene/wnt target 44 gene expression heatmap.pdf")

# single gene expression in different times -------------------------------
rm(list = ls())
library(ggpubr)
library(ggplot2)
library(reshape2)
fpkm <- read.table("~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody2000bp/mm10_genebody_2000bp_merge_rna.csv",
                   sep = ",",header = TRUE,check.names = FALSE)
cluster1_ht <- fpkm[which(fpkm$tracking_id %in% rownames(cluster1_count)),]
colnames(cluster1_ht)[1] <- "ensembl"
library(dplyr)
data <- inner_join(cluster1_count,cluster1_ht,by="ensembl")

plotFPKM <- function(gene_data,name,ymax) {

  df <- melt(data = gene_data,id.vars = "symbol",variable.name = "sample",value.name = "fpkm")
  group <- factor(rep(Times,each=3),levels = Times)
  df <- cbind(df,group)
  ggplot(data = df,aes(x = group,y = log2(fpkm + 1))) + 
    geom_point(shape=21,colour="white",fill="#1B9E77",alpha=0.6,size=3) +
    theme_classic(base_size = 16,base_family = "sans",base_line_size = 1.1) 
    labs(title = name,x=NULL,y="RNA-seq gene expression\n log2 (FPKM + 1)",size=NULL) +
    scale_y_continuous(expand = c(0,0.01),limits = c(0,ymax)) +
    theme(
      plot.margin = margin(t=0,b=0,l=0,r=0,unit = "cm"),
      plot.title = element_text(size = 16,face = "bold",margin = margin(b=0.5,unit = "cm"),hjust = 0.5,vjust = 0.5),
      axis.title.y = element_text(margin = margin(t = 0,b = 0,l = 0,r = 0.5,unit = "cm")),
      axis.text = element_text(colour = "black",size = 14,face = "plain",margin = margin(t = 0.3,b = 0.3,l = 0.3,r = 0.3,unit = "cm")),
      axis.text.x = element_text(angle = 40,hjust = 1)) +
    guides(size=FALSE)
  #stat_compare_means(comparisons = list(c("0week","2weeks"),c("2weeks","4weeks"),c("4weeks","7weeks"),c("7weeks","10weeks")),label = "p.signif",method = "t.test") +
  ggsave(paste("~/project/colon cancer/wnt target gene/",name," RNA-seq fpkm.pdf",sep = ""),p,height = 3.8,width = 3)
  p
}

axin <- data[2,-1]
plotFPKM(gene_data = axin,name = "Axin2")
nlk <- data[7,-1]
plotFPKM(gene_data = nlk,name = "Nlk")
isl1 <- data[34,-c(1,18:22)]
plotFPKM(gene_data = isl1,name = "Isl1",ymax = 10)
shh <- data[5,-1]
plotFPKM(gene_data = shh,name = "Shh",ymax = 4)
wif1 <- data[9,-1]
plotFPKM(gene_data = wif1,name = "Wif1",ymax = 10)
bcl9l <- data[39,-1]
plotFPKM(gene_data = bcl9l,name = "Bcl9l",ymax = 4)


# human HCT116 cell line --------------------------------------------------
library(reshape2)
wang <- read.table("~/project/colon cancer/wnt target gene/genes.fpkm_tracking",sep = "\t",header = TRUE)
isl1 <- wang[which("ISL1" %in% wang$gene_short_name),]
yaojie <- read.table("~/project/colon cancer/wnt target gene/yaojie_HCT116_FPKM.txt",sep = "\t",header = TRUE)
colnames(yaojie)[1] = "symbol"

plotFPKM_human <- function(gene_data,name,ymax) {
  df <- melt(gene_data,id.vars = "symbol",variable.name = "sample",value.name = "fpkm")
  p <- ggplot(data = df,aes(x = sample,y = fpkm)) + 
    geom_point(shape=21,colour="white",fill="#1B9E77",alpha=0.6,size=3) +
    theme_classic(base_size = 16,base_family = "sans",base_line_size = 1.1) +
    labs(title = name,x=NULL,y="RNA-seq gene\nexpression FPKM ",size=NULL) +
    scale_y_continuous(expand = c(0,0.01),limits = c(0,ymax)) +
    theme(
      plot.margin = margin(t=0,b=0,l=0,r=0,unit = "cm"),
      plot.title = element_text(size = 16,face = "bold",margin = margin(b=0.5,unit = "cm"),hjust = 0.5,vjust = 0.5),
      axis.title.y = element_text(margin = margin(t = 0,b = 0,l = 0,r = 0.5,unit = "cm")),
      axis.text = element_text(colour = "black",size = 14,face = "plain",margin = margin(t = 0.3,b = 0.3,l = 0.3,r = 0.3,unit = "cm")),
      axis.text.x = element_text(angle = 40,hjust = 1)) +
    guides(size=FALSE)
  ggsave(paste("~/project/colon cancer/wnt target gene/",name," RNA-seq fpkm yaojie.pdf",sep = ""),p,height = 3.8,width = 3)
  p
}
ISL1 <- yaojie[which(yaojie$symbol %in% "ISL1"),]
plotFPKM_human(gene_data = ISL1,name = "ISL1",ymax = 1)
SHH <- yaojie[which(yaojie$symbol %in% "SHH"),] 
plotFPKM_human(gene_data = SHH,name = "SHH",ymax = 0.1)
WIF1 <- yaojie[which(yaojie$symbol %in% "WIF1"),]
plotFPKM_human(gene_data = WIF1,name = "WIF1",ymax = 0.1)

# class RNA-seq masigpro cluster1 wnt singling pathway gene with chromHMM states --------
wnt_gene <- read.table("~/project/colon cancer/wnt target gene/wnt_44_target_gene.txt")
names(wnt_gene) <- "gene_name"
file <- read.csv("~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody2000bp/mm10_genebody_2000bp_merge_rna.csv")
file = na.omit(file)

df_sub = file[,c("gene_name", "week0", "week2", "week4", "week7", "week10")]
#colnames(df_sub) <- c("ensembl","week0", "week2", "week4", "week7", "week10")
state = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13")
element = c("enhancer","enhancer","enhancer", "promoter", "promoter", "promoter", "promoter", "promoter", "repressed", "repressed",
            "enhancer","quiescent","heterochromatin")
#element = c("heterochromatin", "quiescent", "heterochromatin", "repressed", "repressed", "promoter", "enhancer", "enhancer", "enhancer", "promoter", "promoter", "promoter", "promoter")
dict = setNames(element, state)
df_sub$week0 = unname(dict[df_sub$week0])
df_sub$week2 = unname(dict[df_sub$week2])
df_sub$week4 = unname(dict[df_sub$week4])
df_sub$week7 = unname(dict[df_sub$week7])
df_sub$week10 = unname(dict[df_sub$week10])

library("RColorBrewer")
heterochromatin_color = brewer.pal(n = 11, name = "PRGn")[3]
promoter_color = brewer.pal(n = 9, name = "Greens")[8]
enhancer_color = brewer.pal(n = 9, name = "Reds")[8]
repressed_color = brewer.pal(n = 9, name = "Blues")[7]
quit_color = brewer.pal(n = 9, name = "Greys")[3]

# gene_detail = read.delim("~/project/colon cancer/RNA-seq/mart convert gene/mart_export_mouse_geneid_symbol.txt", sep="\t", header=T)
# symbol_to_ensembl = function(symbol_list){
#   ensembl = gene_detail[match(symbol_list, gene_detail$ensembl_gene_id),c("ensembl_gene_id", "external_gene_name")]
#   colnames(ensembl) = c("ensembl", "symbol")
#   ensembl = unique(na.omit(ensembl))
#   rownames(ensembl) <- ensembl[,1]
#   #ensembl = ensembl[,2]
#   return(ensembl) 
# }

library(dplyr)
# a = symbol_to_ensembl(symbol_list = df_sub$ensembl)
# library(dplyr)
#data_frame <- inner_join(a,wnt_gene,by="symbol")
data_frame <- unique(inner_join(df_sub,wnt_gene,by="gene_name"))
data_frame <- data_frame[-c(36,45),]
rownames(data_frame) <- data_frame$gene_name

enhancer <- data_frame[data_frame$week10 == "enhancer",]
enhancer <- count(enhancer$week0)
enhancer$state <- "enhancer"


promoter <- data_frame[data_frame$`10weeks` == "promoter",]

heterochromatin  <- data_frame[data_frame$`10weeks` == "heterochromatin",]

repressed <- data_frame[data_frame$`10weeks` == "repressed",]
quicsent <- data_frame[data_frame$`10weeks` == "quiescent",]

ht1 <- Heatmap(matrix = as.matrix(enhancer[,-1]),name = "State",
               col = c("enhancer" = enhancer_color,"promoter" = promoter_color,"repressed" = repressed_color,
                       "quiescent" = quit_color,"heterochromatin" = heterochromatin_color),
               show_column_names = FALSE,cluster_columns = FALSE,
               show_row_names = TRUE,cluster_rows = TRUE,
               #top_annotation = top_annotation,
               row_names_gp = gpar(fontsize=16,fontfamily="sans"),
               column_names_gp = gpar(fontsize=16,fontfamily="sans"),
               column_names_rot = 45,
               column_title = "Dominant states",
               column_title_gp = gpar(fontsize=12,fontfamily="sans"),
               height = 5,
               width = 2.2
               )

draw(ht1,merge_legend=TRUE,heatmap_legend_side="left")
library(ggplot2)
pdf("~/project/colon cancer/wnt target gene/wnt_target_gene_heatmap_enhancer.pdf",height = 5,width = 5.3)
draw(ht + ht1,merge_legend=TRUE,heatmap_legend_side="left")
dev.off()

# wnt signal enhancer rna-seq heatmap --------------------------------------------
enhancer_gene <- c("Axin2","Tcf7","Cd44","Nlk","Lgr5","Wif1","Sdc1","Lect2","Csnk1e","Ptk7","Vangl2","Notch1","Ptpro","Nkd1","Tle3","Tmem131l",
                   "Rnf43","Znrf3","Bcl9l","Ccnd1","Apcdd1","Wnt10a","Mdfi","Wnt6","Isl1","Notum","Foxl1")
# write.table(enhancer_gene,file = "~/project/colon cancer/wnt target gene/wnt_enhancer_gene.txt",sep = "\n",
#             quote = FALSE,row.names = FALSE,col.names = FALSE)
new_data <- data[which(data$symbol %in% enhancer_gene),]
rownames(new_data) <- new_data$symbol
new_data <- as.matrix(new_data[,-c(1:17)])

ht <- Heatmap(matrix = new_data,name = "Z-score",
              cluster_columns = FALSE,cluster_rows = FALSE,
              show_column_names = FALSE,column_names_rot = 45,row_names_side = "right",
              col = colorRamp2(breaks = c(-1.3,0,1.3),colors = c("navy","white","firebrick3")),
              top_annotation = top_annotation,
              row_names_gp = gpar(fontsize=16,fontfamily="sans"),
              column_names_gp = gpar(fontsize=16,fontfamily="sans"),
              column_title = "Wnt signal pathway gene",
              column_title_gp = gpar(fontsize=12,fontfamily="sans"),
              height = 5,
              width = 3)
draw(ht,merge_legend=TRUE,heatmap_legend_side="left")


# barplot draw with 0week and 10weeks -------------------------------------
b <- table(data_frame$week0)
c <- table(data_frame$week2)
d <- table(data_frame$week4)
e <- table(data_frame$week7)
f <- table(data_frame$week10)
bar <- cbind(b,c,d,e,f)
write.table(bar,file = "~/project/colon cancer/wnt target gene/barplot.txt",sep = "\t",quote = FALSE,row.names = TRUE)

rm(list = ls())
file <- read.table("~/project/colon cancer/wnt target gene/barplot.txt",sep = "\t")
file$type <- rownames(file)
colnames(file) <- c("week0","week2","week4","week7","week10","type")
library(reshape2)
library(ggplot2)
df <- melt(data = file,id.vars="type",variable.name="time",value.name = "number")
df$type <- factor(c("enhancer","promoter","quiescent","repressed"),levels = c("enhancer","promoter","quiescent","repressed"))
#df$time <- factor(c("week0","week2","week4","week7","week10"),levels = c("week0","week2","week4","week7","week10"))
ggplot(df, aes(x = type,y = number, fill=time)) +
  geom_bar(position = position_dodge(0.9),stat = "identity",width = 0.9) +
  geom_text(aes(label=number),vjust=1.2,position = position_dodge(0.9),size=4,color="white") +
  #geom_line(aes(group=type),stat = "identity",position = position_dodge(0.9)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,25)) +
  theme_classic(base_size = 16,base_family = "sans",base_line_size = 1.1) +
  labs(x=NULL,y="Number of element DSGs",fill=NULL) +
  scale_fill_manual(values = c("#1B9E77","#66A61E","#7570B3","#E7298A","#D95F02"),labels=c("0week","2weeks","4weeks","7weeks","10weeks")) +
  theme(axis.title.x = element_text(margin = margin(t = 0.4,r = 0,b = 0,l = 0,unit = "cm")),
        axis.title.y = element_text(margin = margin(t = 0,r = 0.4,b = 0,l = 0,unit = "cm")),
        axis.text.x = element_text(size = 16,colour = "black",margin = margin(t = 0.1,r = 0,b = 0,l = 0,unit = "cm")),
        axis.text.y = element_text(size = 16,colour = "black",margin = margin(t = 0,r = 0.1,b = 0,l = 0,unit = "cm")),
        axis.ticks.length = unit(0.25,'cm')) +
  theme(legend.position = "right",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 16,family = "sans",lineheight = 1),
        legend.key.height = unit(0.7,"cm"),
        legend.key.size = unit(0.7,"cm")
  )
ggsave(filename = "~/project/colon cancer/wnt target gene/barplot.pdf",height = 3,width = 7)

# wnt gene intensity ------------------------------------------------------
rm(list = ls())
file <- read.csv("~/project/colon cancer/wnt target gene/wnt_target_gene_tss_line.txt",sep = "\t",header = TRUE,check.names = FALSE)
data <- file[,1:6001]
data <- file[,-2]
library(reshape2)
df <- melt(data = data,id.vars = "bins",variable.name = "location",value.name = "number")
df$bins <- factor(c("ctrl-1-H3K27ac","10weeks-1-H3K27ac","7weeks-1-H3K27ac","2weeks-1-H3K27ac","4weeks-1-H3K27ac"),
                  levels = c("ctrl-1-H3K27ac","2weeks-1-H3K27ac","4weeks-1-H3K27ac","7weeks-1-H3K27ac","10weeks-1-H3K27ac"))
df$location <- as.numeric(df$location)
library(ggplot2)
ggplot(df[1:30000,], aes(x = location,y = number,group=bins,colour=bins)) + 
  geom_line(size=0.9) +
  theme_classic(base_family = "sans",base_size = 16,base_line_size = 1.1) +
  scale_y_continuous(expand = c(0,0),limits = c(0,6)) +
  scale_colour_manual(values = c("#1B9E77","#66A61E","#7570B3","#E7298A","#D95F02"),labels=c("0week","2weeks","4weeks","7weeks","10weeks")) +
  scale_x_continuous(expand = c(0,0),limits=c(0,6000),breaks=c(0,3000,6000),labels=c("-3kb","TSS","3kb")) +
  labs(title = "Wnt signal gene (n=44)",x="Distance from TSS (kb)",y="H3K27ac ChIP-seq intensity (RPKM)",colour=NULL) +
  theme(aspect.ratio = 1/0.8,
        plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,margin = margin(b = 0.5,unit = "cm")),
        plot.margin = margin(t = 0,r = 0.2,b = 0,l = 0,unit = "cm"),
        axis.title.x = element_text(margin = margin(t = 0.4,r = 0,b = 0,l = 0,unit = "cm")),
        axis.title.y = element_text(margin = margin(t = 0,r = 0.4,b = 0,l = 0,unit = "cm")),
        axis.text.x = element_text(size = 16,colour = "black",margin = margin(t = 0.1,r = 0,b = 0,l = 0,unit = "cm")),
        axis.text.y = element_text(size = 16,colour = "black",margin = margin(t = 0,r = 0.1,b = 0,l = 0,unit = "cm")),
        axis.ticks.length = unit(0.25,'cm')) +
  theme(legend.position = c(1,0.85),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 16,family = "sans",lineheight = 0.7),
        legend.key.height = unit(0.5,"cm")
  )
ggsave(filename = "~/project/colon cancer/wnt target gene/wnt_target_gene_peak_center.pdf",width = 4.88,height = 4.6)

# wnt enhancer gene n=21 scale a=3000 b=3000 genebody=5000 --------------------------------------------------
file <- read.csv("~/project/colon cancer/wnt target gene/wnt_enhancer_gene_scale_line.txt",sep = "\t",header = TRUE,check.names = FALSE)
data <- file[,1:11001]
data <- data[,-2]
library(reshape2)
df <- melt(data = data,id.vars = "bins",variable.name = "location",value.name = "number")
df$bins <- factor(c("ctrl-1-H3K27ac","10weeks-1-H3K27ac","7weeks-1-H3K27ac","2weeks-1-H3K27ac","4weeks-1-H3K27ac"),
                  levels = c("ctrl-1-H3K27ac","2weeks-1-H3K27ac","4weeks-1-H3K27ac","7weeks-1-H3K27ac","10weeks-1-H3K27ac"))
df$location <- as.numeric(df$location)
library(ggplot2)
ggplot(df, aes(x = location,y = number,group=bins,colour=bins)) + 
  #geom_line() +
  stat_smooth(geom = "smooth",method = "loess",formula = y ~ x,se = FALSE,span=0.01) +
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1) +
  scale_y_continuous(expand = c(0,0),limits = c(0,6)) +
  scale_colour_manual(values = c("#1B9E77","#66A61E","#7570B3","#E7298A","#D95F02"),labels=c("0week","2weeks","4weeks","7weeks","10weeks")) +
  scale_x_continuous(expand = c(0,0),limits=c(0,11000),breaks=c(0,3000,8000,11000),labels=c("-3kb","TSS","TES","3kb")) +
  labs(title = "Wnt signal 10weeks enhancer DSGs (n=21)",x=NULL,y="H3K27ac ChIP-seq\nintensity (RPKM)",colour=NULL) +
  theme(aspect.ratio = 0.6/1,
        plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,margin = margin(b = 0.5,unit = "cm")),
        plot.margin = margin(t = 0,r = 0,b = 0,l = 0,unit = "cm"),
        axis.title.x = element_text(margin = margin(t = 0.4,r = 0,b = 0,l = 0,unit = "cm")),
        axis.title.y = element_text(margin = margin(t = 0,r = 0.4,b = 0,l = 0,unit = "cm")),
        axis.text.x = element_text(size = 18,colour = "black",margin = margin(t = 0.1,r = 0,b = 0,l = 0,unit = "cm")),
        axis.text.y = element_text(size = 18,colour = "black",margin = margin(t = 0,r = 0.1,b = 0,l = 0,unit = "cm")),
        axis.ticks.length = unit(0.25,'cm')) +
  theme(legend.position = c(0.8,0.8),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 18,family = "sans",lineheight = 0.7),
        legend.key.height = unit(0.5,"cm")
  )
ggsave(filename = "~/project/colon cancer/wnt target gene/wnt_target_enhancer_gene_scale.pdf",width = 6,height = 3)

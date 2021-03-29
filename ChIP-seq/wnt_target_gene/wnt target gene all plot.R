rm(list = ls())
setwd("~/project/colon cancer/RNA-seq/edgeR/masigpro_edgeR_filter/")

# RNA-seq gene expression level -------------------------------------------
# wnt 44 target gene 
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

# signal gene expression in different times--------------------------------------------------
rm(list = ls())
library(ggpubr)
library(ggplot2)
library(reshape2)
wnt_target <- read.table("~/project/colon cancer/wnt target gene/wnt_44_target_gene.txt")
colnames(wnt_target) <- "gene_name"
fpkm <- read.table("~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody2000bp/mm10_genebody_2000bp_merge_rna.csv",
                   sep = ",",header = TRUE,check.names = FALSE)
cluster1_ht <- fpkm[which(fpkm$gene_name %in% wnt_target$gene_name),]



plotFPKM <- function(gene_data,name,ymax) {
  
  df <- melt(data = gene_data,id.vars = "gene_name",variable.name = "sample",value.name = "fpkm")
  group <- factor(rep(Times,each=3),levels = Times)
  df <- cbind(df,group)
  p <- ggplot(data = df,aes(x = group,y = log2(fpkm + 1))) + 
    geom_point(shape=21,colour="white",fill="#1B9E77",alpha=0.6,size=3) +
    theme_classic(base_size = 16,base_family = "sans",base_line_size = 1.1) +
    labs(title = name,x="week(s)",y="RNA-seq gene expression\n log2 (FPKM + 1)",size=NULL) +
    scale_y_continuous(expand = c(0,0.01),limits = c(0,ymax)) +
    scale_x_discrete(labels=c(0,2,4,7,10)) +
    theme(aspect.ratio = 1/0.618,
      plot.margin = margin(t=0,b=0,l=0,r=0,unit = "cm"),
      plot.title = element_text(size = 16,face = "bold",margin = margin(b=0.5,unit = "cm"),hjust = 0.5,vjust = 0.5),
      axis.title.y = element_text(margin = margin(t = 0,b = 0,l = 0,r = 0.5,unit = "cm")),
      axis.text = element_text(colour = "black",size = 14,face = "plain",margin = margin(t = 0.3,b = 0.3,l = 0.3,r = 0.3,unit = "cm")),
      axis.text.x = element_text(hjust = 0.5)) +
    guides(size=FALSE)
  #stat_compare_means(comparisons = list(c("0week","2weeks"),c("2weeks","4weeks"),c("4weeks","7weeks"),c("7weeks","10weeks")),label = "p.signif",method = "t.test") +
  ggsave(paste("~/project/colon cancer/wnt target gene/",name," RNA-seq fpkm.pdf",sep = ""),p,height = 3.8,width = 3)
  p
}

Times <- c("0week","2weeks","4weeks","7weeks","10weeks")

isl1 <- cluster1_ht[34,c(2:16,22)]
plotFPKM(gene_data = isl1,name = "Isl1",ymax = 6)


# chromatin states --------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
fpkm <- read.table("~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody2000bp/mm10_genebody_2000bp_merge_rna.csv",
                   sep = ",",header = TRUE,check.names = FALSE)
fpkm = na.omit(fpkm)

df_sub <- fpkm[,c("week0_mean","week2_mean","week4_mean","week7_mean","week10_mean","gene_name", "week0", "week2", "week4", "week7", "week10")]
state = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13")
element = c("enhancer","enhancer","enhancer", "promoter", "promoter", "promoter", "promoter", "promoter", "repressed", "repressed",
            "enhancer","quiescent","heterochromatin")
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

data_frame <- df_sub[df_sub$gene_name %in% gene,]
rownames(data_frame) <- 1:nrow(data_frame)
data_frame <- data_frame[-c(28,36,44,47),]
enhancer <- data_frame[data_frame$week10 == "enhancer",]
rownames(enhancer) <- enhancer$gene_name
state <- as.matrix(enhancer[,-c(1:6)])

ht1 <- Heatmap(matrix = state,name = "State",
               col = c("enhancer" = enhancer_color,"promoter" = promoter_color,"repressed" = repressed_color,
                       "quiescent" = quit_color,"heterochromatin" = heterochromatin_color),
               show_column_names = FALSE,cluster_columns = FALSE,
               show_row_names = TRUE,cluster_rows = TRUE,
               top_annotation = top_annotation,
               row_names_gp = gpar(fontsize=16,fontfamily="sans"),
               column_names_gp = gpar(fontsize=16,fontfamily="sans"),
               column_names_rot = 45,
               column_title = "chromatin states",
               column_title_gp = gpar(fontsize=12,fontfamily="sans"),
               height = 5,
               width = 2.2
)
pdf("~/project/colon cancer/wnt target gene/wnt_24_target_gene_heatmap_enhancer.pdf",height = 5,width = 5.3)
draw(ht1,merge_legend=TRUE,heatmap_legend_side="left")
dev.off()

rna <- df[rownames(enhancer),]

ht2 <- Heatmap(matrix = rna,name = "Z-score",
               cluster_columns = FALSE,cluster_rows = TRUE,show_row_dend = FALSE,
               show_column_names = FALSE,column_names_rot = 45,row_names_side = "right",
               col = colorRamp2(breaks = c(-1,0,1),colors = c("navy","white","firebrick3")),
               top_annotation = top_annotation,
               row_names_gp = gpar(fontsize=12,fontfamily="sans"),
               column_names_gp = gpar(fontsize=16,fontfamily="sans"),
               column_title = "RNA",
               column_title_gp = gpar(fontsize=12,fontfamily="sans"),
               height = 5,
               width = 3)

pdf("~/project/colon cancer/wnt target gene/wnt_24_target_gene_heatmap.pdf",height = 5,width = 5)
draw(ht2 + ht1,merge_legend=TRUE,heatmap_legend_side="left")
dev.off()

pdf("~/project/colon cancer/wnt target gene/wnt_24_target_gene_heatmap_rna.pdf",height = 4,width = 3.5)
draw(ht2,merge_legend=TRUE,heatmap_legend_side="left")
dev.off()

write.table(x = enhancer$gene_name,"~/project/colon cancer/wnt target gene/wnt_enhancer_24_target_gene.txt",sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)

# barplot -----------------------------------------------------------------
b <- table(data_frame$week0)
c <- table(data_frame$week2)
d <- table(data_frame$week4)
e <- table(data_frame$week7)
f <- table(data_frame$week10)
bar <- cbind(b,c,d,e,f)
write.table(bar,file = "~/project/colon cancer/wnt target gene/barplot.txt",sep = "\t",quote = FALSE,row.names = TRUE)

file <- read.delim("~/project/colon cancer/wnt target gene/barplot.txt",sep = "\t")
file$type <- rownames(file)
colnames(file) <- c("week0","week2","week4","week7","week10","type")

df <- melt(data = file,id.vars="type",variable.name="time",value.name = "number")
df$type <- factor(c("enhancer","promoter","quiescent","repressed"),levels = c("enhancer","promoter","quiescent","repressed"))

ggplot(df, aes(x = type,y = number, fill=time)) +
  geom_bar(position = position_dodge(0.9),stat = "identity",width = 0.9) +
  geom_text(aes(label=number),vjust=1.2,position = position_dodge(0.9),size=4,color="white") +
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
plotTSS <- function(filename,limits,number) {
  file <- read.csv(filename,sep = "\t",header = TRUE,check.names = FALSE)
  data <- file[,1:6002]
  data <- file[,-2]
  library(reshape2)
  df <- melt(data = data,id.vars = "bins",variable.name = "location",value.name = "number")
  df$bins <- factor(c("7weeks-1-H3K27ac","10weeks-1-H3K27ac","ctrl-1-H3K27ac","2weeks-1-H3K27ac","4weeks-1-H3K27ac"),
                    levels = c("ctrl-1-H3K27ac","2weeks-1-H3K27ac","4weeks-1-H3K27ac","7weeks-1-H3K27ac","10weeks-1-H3K27ac"))
  df$location <- as.numeric(df$location)
  library(ggplot2)
  p <- ggplot(df[1:30000,], aes(x = location,y = number,group=bins,colour=bins)) + 
    geom_line(size=0.9) +
    theme_classic(base_family = "sans",base_size = 16,base_line_size = 1.1) +
    scale_y_continuous(expand = c(0,0),limits = limits) +
    scale_colour_manual(values = c("#1B9E77","#66A61E","#7570B3","#E7298A","#D95F02"),labels=c("0week","2weeks","4weeks","7weeks","10weeks")) +
    scale_x_continuous(expand = c(0,0),limits=c(0,6000),breaks=c(0,3000,6000),labels=c("-3kb","TSS","3kb")) +
    labs(title = paste("Wnt signal gene (n=",number,")",sep = ""),x="Distance from TSS (kb)",y="H3K27ac ChIP-seq intensity (RPKM)",colour=NULL) +
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
  basename <- basename(filename)
  pdfname <- paste("~/project/colon cancer/wnt target gene/profile/",gsub("_line.txt",".pdf",basename),sep = "")
  pngname <- paste("~/project/colon cancer/wnt target gene/profile/",gsub("_line.txt",".png",basename),sep = "")
  ggsave(filename = pdfname,p,width = 4.88,height = 4.6)
  ggsave(filename = pngname,p,width = 4.88,height = 4.6,dpi = 300)
  p
}

plotTSS(filename = "~/project/colon cancer/wnt target gene/profile/wnt_enhancer_24_target_gene_tss_line.txt",limits = c(0,20),number = 24)
plotTSS(filename = "~/project/colon cancer/wnt target gene/profile/wnt_44_target_gene_tss_line.txt",limits = c(0,15),number = 44)

# wnt gene scale a=3000 b=3000 genebody=5000 --------------------------------------------------
plotLine <- function(filename,limits,number) {
  file <- read.csv(filename,sep = "\t",header = TRUE,check.names = FALSE)
  data <- file[,1:11002]
  data <- data[,-2]
  library(reshape2)
  df <- melt(data = data,id.vars = "bins",variable.name = "location",value.name = "number")
  df$bins <- factor(c("7weeks-H3K27ac","10weeks-H3K27ac","ctrl-H3K27ac","2weeks-H3K27ac","4weeks-H3K27ac"),
                    levels = c("ctrl-H3K27ac","2weeks-H3K27ac","4weeks-H3K27ac","7weeks-H3K27ac","10weeks-H3K27ac"))
  df$location <- as.numeric(df$location)
  library(ggplot2)
  p <- ggplot(df, aes(x = location,y = number,group=bins,colour=bins)) + 
    stat_smooth(geom = "smooth",method = "loess",formula = y ~ x,se = FALSE,span=0.01) +
    theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1) +
    scale_y_continuous(expand = c(0,0),limits = limits,breaks = c(0,5,10)) +
    scale_colour_manual(values = c("#1B9E77","#66A61E","#7570B3","#E7298A","#D95F02"),labels=c("0week","2weeks","4weeks","7weeks","10weeks")) +
    scale_x_continuous(expand = c(0,0),limits=c(0,11000),breaks=c(0,3000,8000,11000),labels=c("-3kb","TSS","TES","3kb")) +
    labs(title = paste("Wnt signal gene (n=",number,")",sep = ""),x=NULL,y = "H3K27ac ChIP-seq\nintensity (RPKM)",colour=NULL) +
    theme(aspect.ratio = 0.6/1,
          plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,margin = margin(b = 0.5,unit = "cm")),
          plot.margin = margin(t = 0,r = 0.3,b = 0,l = 0.3,unit = "cm"),
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
  basename <- basename(filename)
  pdfname <- paste("~/project/colon cancer/wnt target gene/profile/",gsub("_line.txt",".pdf",basename),sep = "")
  pngname <- paste("~/project/colon cancer/wnt target gene/profile/",gsub("_line.txt",".png",basename),sep = "")
  ggsave(filename = pdfname,p,width = 5.2,height = 3)
  ggsave(filename = pngname,p,width = 5.2,height = 3,dpi = 300)
  p
}

plotLine(filename = "~/project/colon cancer/wnt target gene/profile/wnt_enhancer_24_target_gene_scale_line.txt",limits = c(0,15),number = 24)
plotLine(filename = "~/project/colon cancer/wnt target gene/profile/wnt_44_target_gene_scale_line.txt",limits = c(0,10),number = 44)

# wnt gene in TCGA --------------------------------------------------------
rm(list = ls())
library(biomaRt)

file <- read.delim("~/project/colon cancer/wnt target gene/wnt_44_target_gene.txt",sep = "\t",header = FALSE,col.names = "gene_name")
#file <- read.delim("~/project/colon cancer/wnt target gene/wnt_enhancer_24_target_gene.txt",sep = "\t",header = FALSE,col.names = "gene_name")
tcga <- read.delim("~/project/colon cancer/TCGA/COAD_FPKM_add_symbol.txt",sep = "\t",check.names = FALSE)

mouse <- useDataset("mmusculus_gene_ensembl", mart = useMart("ensembl"))
human <- useMart("ensembl",dataset = "hsapiens_gene_ensembl")
bm <- getLDS(attributes = c("external_gene_name","ensembl_gene_id"),mart = mouse,filters = "external_gene_name",
             attributesL = c("external_gene_name","ensembl_gene_id"),martL = human,values = file$gene_name)

data <- tcga[tcga$Gene.name %in% bm$Gene.name.1,]

library(reshape2)
df <- melt(data = data,id.vars = "Gene.name",variable.name = "sample",value.name = "value")
df <- cbind(df,type = c(rep("Normal",each=42*nrow(data)),rep("Tumor",each=599*nrow(data))))
df <- na.omit(df)

library(ggpubr)
library(ggplot2)
p <- ggplot(data = df, aes(x = factor(type),y = as.numeric(log2(value + 1)),fill=factor(type))) +
  stat_boxplot(geom = "errorbar",linetype=1,width=0.3,position = "identity") +
  geom_boxplot(outlier.fill = NA,outlier.shape=NA,width=0.6) +
  stat_compare_means(method = "t.test",paired = FALSE,comparisons = list(c("Normal","Tumor")),label = "p.signif",label.y = c(10)) +
  theme_bw(base_size = 18,base_family = "sans",base_line_size = 1.1) +
  scale_y_continuous(expand = c(0.02,0.02),limits = c(0,12),breaks = c(0,4,8,12)) + 
  scale_x_discrete(expand = c(0.3,0.3),breaks = c("Normal","Tumor"),labels=c("Normal\nn=42","Tumor\nn=599")) +
  scale_fill_manual(values = c("Normal" = "#D9D9D9","Tumor" = "#A50F15")) +
  labs(title = "Wnt signal gene(24) \nin TCGA",x=NULL,y="Expression log2(FPKM+1)",fill=NULL) +
  theme(aspect.ratio = 1/0.75,
        plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,family = "sans"),
        axis.text = element_text(colour = "black",size = 18,family = "sans"),
        axis.text.y = element_text(size = 18,colour = "black",margin = margin(l = 0.1,r = 0.1,unit = "cm")),
        axis.ticks.length = unit(2,'mm'),
        panel.border = element_rect(size = 1.1),
        panel.grid = element_blank()
  ) +
  guides(fill=FALSE)
ggsave(paste("~/project/colon cancer/wnt target gene/FPKM plot/wnt 24 gene TCGA FPKM.pdf",sep = ""),p,height = 4.5,width = 3)
p


# wnt 44 gene fpkm  -------------------------------------------------------
fpkm <- read.csv("~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody2000bp/mm10_genebody_2000bp_merge_rna.csv",check.names = FALSE)
file <- read.delim("~/project/colon cancer/wnt target gene/wnt_44_target_gene.txt",sep = "\t",header = FALSE,col.names = "gene_name")

data <- fpkm[fpkm$gene_name %in% file$gene_name,]
data <- data[,c(22,2:16)]
data <- data[,c("gene_name","control_0","control_1","control_2","2week_0","2week_1","2week_2","4week_0","4week_1","4week_2",
                  "7week_0","7week_1","7week_2","10week_0","10week_1","10week_2")]
                   
library(reshape2)
df <- melt(data = data,id.vars = "gene_name",variable.name = "sample",value.name = "value")
library(stringr)
df["time"] <- str_split_fixed(df$sample,pattern = "_",n = 2)[,1]
df["rep"] <- str_split_fixed(df$sample,pattern = "_",n = 2)[,2]
df <- na.omit(df)

df$time <- factor(rep(c("control","2week","4week","7week","10week"),each=141),levels = c("control","2week","4week","7week","10week"))

library(ggpubr)
library(ggplot2)
p <- ggplot(data = df, aes(x = factor(sample),y = as.numeric(value),fill=time)) +
  stat_boxplot(geom = "errorbar",linetype=1,width=0.45,position = "identity") +
  geom_boxplot(outlier.fill = NA,outlier.shape=NA,width=0.9)  +
  #geom_point(aes(fill=time,colour=time),alpha=0.5,size=2,shape=21,position = position_jitter(0.5)) +
  #stat_compare_means(method = "t.test",paired = FALSE,comparisons = list(c("Normal","Tumor")),label = "p.signif",label.y = c(10)) +
  theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) +
  scale_y_continuous(expand = c(0.01,0.01),limits = c(0,60),breaks = c(0,20,40,60)) +
  scale_x_discrete(expand = c(0.05,0.05),labels=c("rep1","rep2\n0week","rep3","rep1","rep2\n2weeks","rep3","rep1","rep2\n4weeks","rep3",
                   "rep1","rep2\n7weeks","rep3","rep1","rep2\n10weeks","rep3")) +
  scale_fill_manual(values = c("#1B9E77","#66A61E","#7570B3","#E7298A","#D95F02")) +
  labs(title = "Wnt signal gene (44)",x=NULL,y="Expression (FPKM)",fill=NULL) +
  theme(#aspect.ratio = 1/0.75,
        plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,family = "sans"),
        axis.text.x = element_text(colour = "black",size = 18,family = "sans",angle = 45,hjust = 1),
        axis.text.y = element_text(size = 18,colour = "black",margin = margin(l = 0.1,r = 0.1,unit = "cm")),
        axis.ticks.length = unit(2,'mm')
        #panel.border = element_rect(size = 1.1),
        #panel.grid = element_blank()
  ) +
  guides(fill=FALSE)
ggsave(paste("~/project/colon cancer/wnt target gene/FPKM plot/wnt 44 gene AOMDSS FPKM.pdf",sep = ""),p,height = 4.5,width = 6)
p

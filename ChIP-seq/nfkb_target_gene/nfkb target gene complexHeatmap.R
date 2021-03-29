rm(list = ls())
setwd("~/project/colon cancer/nfkb target gene/")
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

refTosymbol <- read.table("~/project/colon cancer/nfkb target gene/nfkb refTosymbol getBM.txt",sep = "\t")

Times <- c("0week","2weeks","4weeks","7weeks","10weeks")
rna <- read.delim(file = "~/project/colon cancer/RNA-seq/DESeq2_maSigPro/merge_htseqCount.txt",sep = "\t",
                  header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
rna <- rna[,-16]
columns <-c()
for (ti in Times) {
  for (i in 1:3) {
    column <- paste(ti,i,sep = "_rep")
    columns <- c(columns,column)
  }
}
colnames(rna) <- columns
rna <- rna[rowSums(rna > 3) > 2,]

rna$ensembl <- substr(rownames(rna),1,18)
# byown
gene_detail <- read.table('~/project/colon cancer/RNA-seq/mart convert gene/mart_export_mouse_geneid_symbol.txt',sep = "\t")
names(gene_detail) <- c("ensembl","symbol")
library(dplyr)
rna <- inner_join(x = rna,y = gene_detail,by="ensembl")
rna_data <- subset(rna,symbol %in% refTosymbol$external_gene_name)
rownames(rna_data) <- rna_data$symbol

data <- data.frame(apply(rna_data[,1:15],2,scale))
rownames(data) <- rna_data$symbol
df <- data.frame(t(apply(data[,1:15],1,scale)))

# get rowMeans()
df[,"week0"] <- rowMeans(df[,1:3]) 
df[,"week2"] <- rowMeans(df[,4:6])
df[,"week4"] <- rowMeans(df[,7:9]) 
df[,"week7"] <- rowMeans(df[,10:12]) 
df[,"week10"] <- rowMeans(df[,13:15]) 

pcor <- read.delim("~/project/colon cancer/nfkb target gene/nfkb distal correlation bigger than 0 gene.txt",sep = "\t")
#row_order <- intersect(row_order,H3K27ac_rpm$gene_name)

H3K27ac <- as.matrix(read.table('H3K27ac_addrpm_100kb_mean.txt',sep = "\t",header = TRUE,check.names = FALSE))
H3K4me1 <- as.matrix(read.table('H3K4me1_addrpm_100kb_mean.txt',sep = "\t",header = TRUE,check.names = FALSE))
H3K4me3 <- as.matrix(read.table('H3K4me3_addrpm_100kb_mean.txt',sep = "\t",header = TRUE,check.names = FALSE))

H3K27ac_data <- H3K27ac[pcor$symbol,]
H3K4me1_data <- H3K4me1[pcor$symbol,]
H3K4me3_data <- H3K4me3[pcor$symbol,]
rna_df <- df[rownames(H3K27ac_data),]
rna_df  <- as.matrix(rna_df[,16:20])

cor_H3K27ac <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], H3K27ac_data[i, ] ))
cor_H3K4me1 <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], H3K4me1_data[i, ] ))
cor_H3K4me3 <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], H3K4me3_data[i, ] ))
cor_list <- as.matrix(cbind(cor_H3K27ac,cor_H3K4me1,cor_H3K4me3))
rownames(cor_list) = rownames(rna_df)

top_annotation <- HeatmapAnnotation(
  type = c(rep("Normal",1),rep(c("Inflam","Tumor"),each=2)),
  time = rep(Times,each=1),
  #Rep = rep(1:3,5),
  col = list(
    type = c("Normal" = "#DADAEB","Inflam" = "#9E9AC8","Tumor" =  "#6A51A3"),
    time = c("0week" = "#1B9E77","2weeks"="#66A61E","4weeks"="#7570B3","7weeks"="#E7298A","10weeks"="#D95F02")
    #Rep = c("1" = "#FEE0D2","2"=  "#FCBBA1","3"= "#FC9272")
  ),
  gap = unit(1,'points'),
  simple_anno_size = unit(0.25,"cm"),show_annotation_name = FALSE,show_legend = c(TRUE,TRUE,TRUE),
  annotation_legend_param = list(
    type = list(title = "Type",at = c("Normal","Inflam","Tumor"),labels = c("Normal","Inflam","Tumor")),
    time = list(title = "Time",at = Times,labels = Times)),
  #Rep = list(title = "Rep",title_position = "topleft"),
  annotation_name_gp = gpar(fontzize = 16,fontfamily="sans")
)
histone_top_annotation <- HeatmapAnnotation(
  Histone = c("H3K27ac","H3K4me1","H3K4me3"), #,"H3K27me3","H3K9me3"
  col = list(
    Histone = c("H3K27ac" = "#A50F15","H3K4me1" = "#FF7F00","H3K4me3" =  "#006D2C")), #,"H3K27me3" = "#2171B5","H3K9me3" = "#9970AB"
  gap = unit(1,'points'),
  simple_anno_size = unit(0.25,"cm"),show_annotation_name = FALSE,show_legend = c(TRUE,TRUE,TRUE),
  annotation_legend_param = list(
    Histone = list(title = "Histone\nmodification",
                   at = c("H3K27ac","H3K4me1","H3K4me3"), #,"H3K27me3","H3K9me3"),
                   labels = c("H3K27ac","H3K4me1","H3K4me3"))), #,"H3K27me3","H3K9me3"))),
  annotation_name_gp = gpar(fontzize = 18,fontfamily="sans"))

labels = c("Il1r1","Ptx3","Vcam1","Aqp1","Maml2","Esr1","Phlda1","Sqstm1","Aoc3","Fzd2","Lifr","Enpp2","Mcc","Pdgfrb","Cd248","Il13ra2")
unknown = c("Gorab","Sdc4","Plpp3","Nipal1","C1s1","Serpina3f","Lhfpl2","Fst","Zfp395","Gdnf","Zfpm2","Lbh")
a <- data.frame(ht@row_order)
b <- a[labels,]
c <- a[unknown,]

a <- row_order(ht)
a1 <- ht@row_names_param$labels[unlist(a[1])]
a2 <- ht@row_names_param$labels[unlist(a[2])]
a3 <- ht@row_names_param$labels[unlist(a[3])]
write.table(x = a,file = "~/project/colon cancer/nfkb target gene/nfkb_target_gene_138.txt",sep = "\t",quote = FALSE,col.names = FALSE)
write.table(x = a1,file = "~/project/colon cancer/nfkb target gene/nfkb_cluster1_gene.txt",sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
write.table(x = a2,file = "~/project/colon cancer/nfkb target gene/nfkb_cluster2_gene.txt",sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
write.table(x = a3,file = "~/project/colon cancer/nfkb target gene/nfkb_cluster3_gene.txt",sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)

ht@row_names_param$labels

left_annotation <- rowAnnotation(
  gene = anno_mark(at=c(1,18,24,45,62,64,71,75,86,87,105,108,130,132,134,138),
                   labels = c("Il1r1","Ptx3","Vcam1","Aqp1","Maml2","Esr1","Phlda1","Sqstm1","Aoc3","Fzd2","Lifr","Enpp2","Mcc",
                              "Pdgfrb","Cd248","Il13ra2"),
                   labels_gp = gpar(col= "navy",fontsize=14),
                   side = "left"),
  unknown = anno_mark(#at = c(133,39,33,71,34,99,92,11,41,74,72,82),
                      at = c(7,14,30,33,48,92,97,99,102,106,107,125),
                      labels = c("Gorab","Sdc4","Plpp3","Nipal1","C1s1","Serpina3f","Lhfpl2","Fst","Zfp395","Gdnf","Zfpm2","Lbh"),
                      labels_gp = gpar(col="pink",fontsize=14),side = "left"))

row_tree <- hclust(dist(rna_df),method = "complete")
row_order <- row_tree$labels[row_tree$order]

ht <- Heatmap(as.matrix(rna_df), column_title = "RNA",border = TRUE,
              cluster_rows = FALSE,cluster_columns = FALSE,row_names_centered = FALSE,
              show_column_names = FALSE,
              show_row_dend = FALSE,
              row_names_gp = gpar(fontsize=4,fontfamily="sans"),
              col = colorRamp2(c(-1,0,1),c("navy","white","firebrick3")),
              row_order = rev(row_order),
              row_km = 3,
              row_title = "cluster%s",
              row_title_gp = gpar(font=2,fontsize=18,fontfamily="sans"),
              column_title_gp = gpar(fontsize=18,fontfamily="sans"),
              top_annotation = top_annotation,
              left_annotation = left_annotation,
              show_row_names = TRUE,
              heatmap_legend_param = list(title = "RNA"),
              height = 5,
              width = unit(3,"cm"))
ht_cor <- Heatmap(as.matrix(cor_list), column_title = "Cor",border = FALSE,
               col = colorRamp2(c(-1,0,1),c("#377EB8","white","#E41A1C")),
               cluster_columns = FALSE,cluster_rows = FALSE,
               show_row_names = FALSE,show_column_names = FALSE,
               column_title_gp = gpar(fontsize=18),
               row_names_gp = gpar(fontsize=4),
               #row_order = rev(b),
               top_annotation = histone_top_annotation,
               height = 5,
               width = unit(2.5,"cm"),show_heatmap_legend = TRUE,
               heatmap_legend_param = list(title = "Correlation"))
pdf("~/project/colon cancer/nfkb target gene/nfkb target gene 138.pdf",height = 7,width = 6)
draw(ht+ ht_cor,merge_legend=TRUE,heatmap_legend_side="right",annotation_legend_side="right",ht_gap=unit(1,"mm"),
     column_title_gp=gpar(fontsize=18,fontfamily="sans"))
dev.off()

# chromatin state heatmap -------------------------------------------------
state <- read.csv("~/project/colon cancer/chip-seq/luo_figure/gene_base_summary_merge_rna.csv",check.names = FALSE)
state <- na.omit(state)

df_sub = state[,c("symbol", "week0", "week2", "week4", "week7", "week10")]
colnames(df_sub) <- c("ensembl","0week","2weeks","4weeks","7weeks","10weeks")
state = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13")
element = c("heterochromatin", "quiescent", "heterochromatin", "repressed", "repressed", "promoter", "enhancer", "enhancer", "enhancer", "promoter", "promoter", "promoter", "promoter")
dict = setNames(element, state)
df_sub$`0week` = unname(dict[df_sub$`0week`])
df_sub$`2weeks` = unname(dict[df_sub$`2weeks`])
df_sub$`4weeks` = unname(dict[df_sub$`4weeks`])
df_sub$`7weeks` = unname(dict[df_sub$`7weeks`])
df_sub$`10weeks` = unname(dict[df_sub$`10weeks`])

library("RColorBrewer")
heterochromatin_color = brewer.pal(n = 11, name = "PRGn")[3]
promoter_color = brewer.pal(n = 9, name = "Greens")[8]
enhancer_color = brewer.pal(n = 9, name = "Reds")[8]
repressed_color = brewer.pal(n = 9, name = "Blues")[7]
quit_color = brewer.pal(n = 9, name = "Greys")[3]

gene_detail = read.delim("~/project/colon cancer/RNA-seq/mart convert gene/mart_export_mouse_geneid_symbol.txt", sep="\t", header=T)
symbol_to_ensembl = function(symbol_list){
  ensembl = gene_detail[match(symbol_list, gene_detail$ensembl_gene_id),c("ensembl_gene_id", "external_gene_name")]
  colnames(ensembl) = c("ensembl", "symbol")
  ensembl = unique(na.omit(ensembl))
  rownames(ensembl) <- ensembl[,1]
  #ensembl = ensembl[,2]
  return(ensembl) 
}

chromstate = symbol_to_ensembl(symbol_list = df_sub$ensembl)
library(dplyr)
data_frame <- chromstate[which(chromstate$symbol %in% rownames(rna_df)),]
data_frame <- unique(inner_join(data_frame,df_sub,by="ensembl"))
data_frame <- data_frame[-c(18),]
rownames(data_frame) <- data_frame$symbol
data_frame <- data_frame[,3:7]
data_frame <- data_frame[rownames(rna_df),]

ht_state <- Heatmap(matrix = as.matrix(data_frame),name = "State",border = FALSE,
                    show_row_names = FALSE,show_column_names = FALSE,
                    #row_order = rev(row_order),
                    top_annotation = top_annotation,
                    col = c("enhancer" = enhancer_color,"promoter" = promoter_color,"repressed" = repressed_color,
                            "quiescent" = quit_color,"heterochromatin" = heterochromatin_color),
                    column_title = "State",
                    column_title_gp = gpar(fontsize=18,fontfamily="sans"),
                    height = 5,
                    width = unit(2.5,"cm")
)
pdf("~/project/colon cancer/nfkb target gene/nfkb target gene 138 add chromhmm.pdf",height = 7,width = 7)
draw(ht+ht_cor+ht_state,merge_legend=TRUE,heatmap_legend_side="right",annotation_legend_side="right",ht_gap=unit(1,"mm"),
     column_title_gp=gpar(fontsize=18,fontfamily="sans"))
dev.off()


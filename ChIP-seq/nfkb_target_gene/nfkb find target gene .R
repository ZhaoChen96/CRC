rm(list = ls())
setwd("~/project/colon cancer/nfkb target gene/")
library(stringr)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)



mouse <- useMart(dataset = "mmusculus_gene_ensembl",biomart = "ensembl",host="http://uswest.ensembl.org")
#### query inculding Refseq gene 
#listFilters(mouse)[grep("Refseq",as.character(listFilters(mouse)$description)),]
#filters <- listFilters(mouse)
#grep("refseq", filters$name, ignore.case=T, value=T)

nfkb <- read.table("~/laboratory/labmeeting/Sup. Table2.txt",sep = "\t",header = TRUE)
nfkb$ACCESSION.NUMBER <- str_split_fixed(nfkb$ACCESSION.NUMBER,"[.]",2)[,1]
#### second method
# library(tidyr)
# nfkb %>% separate(ACCESSION.NUMBER,c("Name","num"),"[.]")

refTosymbol <- getBM(filters = "refseq_mrna",attributes = c("refseq_mrna","ensembl_gene_id","external_gene_name"),
                     mart = mouse,values = nfkb$ACCESSION.NUMBER)
write.table(refTosymbol,file = "~/project/colon cancer/nfkb target gene/nfkb refTosymbol getBM.txt",sep = "\t",quote = FALSE)
refTosymbol <- read.table("~/project/colon cancer/nfkb target gene/nfkb refTosymbol getBM.txt",sep = "\t")

 # RNA-seq heatmap ---------------------------------------------------------
#### use read count 
# first method: total read count apply column scale, then subset 240 gene apply row scale
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
# biomart
ensembl_symbol <- getBM(attributes = c("entrezgene_id","external_gene_name","ensembl_gene_id"),filters = "ensembl_gene_id",
                        values = rna$ensembl,mart = mouse)
index <- ensembl_symbol$external_gene_name[match(rna$ensembl,ensembl_symbol$ensembl_gene_id)]
rna$symbol <- ensembl_symbol$external_gene_name[match(rna$ensembl,ensembl_symbol$ensembl_gene_id)]
rna_data <- subset(rna,symbol %in% refTosymbol$external_gene_name)
rownames(rna_data) <- rna_data$symbol

data <- data.frame(apply(rna_data[,1:15],2,scale))
rownames(data) <- rna_data$symbol
df <- data.frame(t(apply(data[,1:15],1,scale)))
write.table(df,file = "~/project/colon cancer/nfkb target gene/nfkb RNA-seq matrix.txt",sep = "\t",quote = FALSE,row.names = TRUE)
df <- read.table("~/project/colon cancer/nfkb target gene/nfkb RNA-seq matrix.txt",sep = "\t")

# get rowMeans()
df[,"week0"] <- rowMeans(df[,1:3]) 
df[,"week2"] <- rowMeans(df[,4:6])
df[,"week4"] <- rowMeans(df[,7:9]) 
df[,"week7"] <- rowMeans(df[,10:12]) 
df[,"week10"] <- rowMeans(df[,13:15]) 
# get rowMedians()
# library(miscTools)
# df[,"week0"] <- rowMedians(df[,1:3]) 
# df[,"week2"] <- rowMedians(df[,4:6])
# df[,"week4"] <- rowMedians(df[,7:9]) 
# df[,"week7"] <- rowMedians(df[,10:12]) 
# df[,"week10"] <- rowMedians(df[,13:15]) 

median(as.matrix(df[,c(16:20)]))
mean(as.matrix(df[,c(16:20)]))
min(as.matrix(df[,c(16:20)]))

row_tree <- hclust(dist(df),method = "complete")
row_order <- row_tree$labels[row_tree$order]

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
ht <- Heatmap(as.matrix(df[,c(16:20)]), column_title = "RNA",
              col = colorRamp2(c(-1,0,1),c("navy","white","firebrick3")),
              row_order = rev(row_tree$order),
              cluster_rows = FALSE,cluster_columns = FALSE,row_names_centered = FALSE,
              show_row_names = FALSE,show_column_names = FALSE,row_names_gp = gpar(fontsize=4),
              top_annotation = top_annotation,
              #left_annotation = ha,
              heatmap_legend_param = list(title = "RNA"),
              height = 5,
              width = unit(3,"cm"))
draw(ht,merge_legend=TRUE)

# GO analysis -------------------------------------------------------------
library(clusterProfiler)
gene <- df[which(rowSums(df[,c(17,18)]) > rowSums(df[,c(19,20)])),]
gene <- gene[which(gene[,16] < rowMeans(gene[,c(17,18)])),]
b <- row_tree$labels[row_tree$order][95:170] #2,4weeks 
c <- row_tree$labels[row_tree$order][1:62]

bp <- enrichGO(gene = c,
               keyType = "SYMBOL",
               OrgDb = "org.Mm.eg.db",
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05
)
bp_0.05 <- subset(bp@result,bp@result$pvalue < 0.05,)[c(1,2,3,5,6,7,8,9),]
# b [c(1,7:13),]
# total gene[c(1,3,4,5,6,8,9,12),]
plotbp <- function(data,tf) {
  p <- ggplot(data,aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill=-log10(pvalue))) +
    geom_bar(stat = "identity",width = 0.8,position = position_dodge(width = 0.8)) +
    theme_classic(base_size = 18,base_family = "sans") +
    scale_x_continuous(expand = c(0,0),limits = c(0,max(ceiling(-log10(bp_0.05$pvalue))))) +
    scale_fill_gradient(low = "#FEE6CE",high = "#F16913") +
    labs(title = paste("Biological process of",tf,"putative target gene ",sep = " "),x="-log10(pValue)",y=NULL) +
    guides(fill=FALSE) +
    theme(plot.title = element_text(size = 20,hjust = 1.2,vjust = 0.5),
          plot.margin = margin(t=3,r=3,b=3,l=3,unit = "pt"),
          axis.text = element_text(colour = 'black',size = 18))
  p
}
plotbp(data = bp_0.05,tf = "nfkb")
ggsave(filename = paste("~/project/colon cancer/nfkb target gene/nfkb_down_bp.pdf",sep = ""),height = 3,width = 7.5)

entrezgene <- getBM(mart = mouse,attributes = c('external_gene_name','entrezgene_id'),filters = "external_gene_name",values = c)    
kegg <- enrichKEGG(gene = entrezgene[,2],
                   keyType = "kegg",
                   organism = "mmu",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05
)
library(scales)
kegg_0.05 <- kegg@result[kegg@result$pvalue < 0.05,][1:10,]
plotkegg <- function(data,tf){
  p <- ggplot(data,aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill=-log10(pvalue))) +
    geom_bar(stat = "identity",width = 0.8,position = position_dodge(width = 0.8)) +
    theme_classic(base_size = 14,base_family = "sans") +
    scale_x_continuous(expand = c(0,0),limits = c(0,16)) +#limits = c(0,max(ceiling(-log10(kegg@result$pvalue))))) +
    scale_fill_gradient(low = "#C6DBEF",high = "#08519C") +
    labs(title = paste("KEGG pathway of",tf,"putative target gene ",sep = " "),x="-log10(pValue)",y=NULL) +
    guides(fill=FALSE) +
    theme(plot.title = element_text(size = 14,hjust = 1.3,vjust = 0.5),
          plot.margin = margin(t=0,r=0.3,b=0,l=0,unit = "cm"),
          axis.text = element_text(colour = 'black',size = 14))
  p
}
plotkegg(data = kegg_0.05,tf = "nfkb")

# correlation Complexheatmap ----------------------------------------------
# histone modification between TSS upstream and downstream 1500kb
H3K27ac_rpm <- read.table("~/project/colon cancer/chip-seq/tf taget gene/three_regions/H3K27ac_rpm_1500bp.txt",sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
row_order <- intersect(row_order,H3K27ac_rpm$gene_name)
write.table(x = row_order,file = "~/project/colon cancer/nfkb target gene/nfkb target gene 241.txt",sep = "\t",row.names = FALSE,quote = FALSE)
H3K27ac_data <- subset(H3K27ac_rpm,H3K27ac_rpm$gene_name %in% row_order)
H3K27ac_df <- data.frame(H3K27ac_data[,2:16])
rownames(H3K27ac_df) <- H3K27ac_data$gene_name

H3K27ac_df[,"week0"] <- rowMeans(H3K27ac_df[,1:3])
H3K27ac_df[,"week2"] <- rowMeans(H3K27ac_df[,4:6])
H3K27ac_df[,"week4"] <- rowMeans(H3K27ac_df[,7:9])
H3K27ac_df[,"week7"] <- rowMeans(H3K27ac_df[,10:12])
H3K27ac_df[,"week10"] <- rowMeans(H3K27ac_df[,13:15])

H3K4me1_rpm <- read.table("~/project/colon cancer/chip-seq/tf taget gene/three_regions/H3K4me1_rpm_1500bp.txt",sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
H3K4me1_data <- subset(H3K4me1_rpm, H3K4me1_rpm$gene_name %in% row_order)
#H3K4me1_data <- H3K4me1_rpm[which(b %in% H3K4me1_rpm$`gene name` ),]
H3K4me1_df <- data.frame(H3K4me1_data[,2:16])
rownames(H3K4me1_df) <- H3K4me1_data$gene_name

H3K4me1_df[,"week0"] <- rowMeans(H3K4me1_df[,1:3])
H3K4me1_df[,"week2"] <- rowMeans(H3K4me1_df[,4:6])
H3K4me1_df[,"week4"] <- rowMeans(H3K4me1_df[,7:9])
H3K4me1_df[,"week7"] <- rowMeans(H3K4me1_df[,10:12])
H3K4me1_df[,"week10"] <- rowMeans(H3K4me1_df[,13:15])

H3K4me3_rpm <- read.table("~/project/colon cancer/chip-seq/tf taget gene/three_regions/H3K4me3_rpm_1500bp.txt",sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
H3K4me3_data <- subset(H3K4me3_rpm, H3K4me3_rpm$gene_name %in% row_order)
#H3K4me3_data <- H3K4me3_rpm[which(b %in% H3K4me3_rpm$`gene name` ),]
H3K4me3_df <- data.frame(H3K4me3_data[,2:16])
rownames(H3K4me3_df) <- H3K4me3_data$gene_name

H3K4me3_df[,"week0"] <- rowMeans(H3K4me3_df[,1:3])
H3K4me3_df[,"week2"] <- rowMeans(H3K4me3_df[,4:6])
H3K4me3_df[,"week4"] <- rowMeans(H3K4me3_df[,7:9])
H3K4me3_df[,"week7"] <- rowMeans(H3K4me3_df[,10:12])
H3K4me3_df[,"week10"] <- rowMeans(H3K4me3_df[,13:15])

H3K27me3_rpm <- read.table("~/project/colon cancer/chip-seq/tf taget gene/three_regions/H3K27me3_rpm_1500bp.txt",sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
H3K27me3_data <- subset(H3K27me3_rpm, H3K27me3_rpm$gene_name %in% row_order)
#H3K27me3_data <- H3K27me3_rpm[which(b %in% H3K27me3_rpm$`gene name` ),]
H3K27me3_df <- data.frame(H3K27me3_data[,2:16])
rownames(H3K27me3_df) <- H3K27me3_data$gene_name

H3K27me3_df[,"week0"] <- rowMeans(H3K27me3_df[,1:3])
H3K27me3_df[,"week2"] <- rowMeans(H3K27me3_df[,4:6])
H3K27me3_df[,"week4"] <- rowMeans(H3K27me3_df[,7:9])
H3K27me3_df[,"week7"] <- rowMeans(H3K27me3_df[,10:12])
H3K27me3_df[,"week10"] <- rowMeans(H3K27me3_df[,13:15])

H3K9me3_rpm <- read.table("~/project/colon cancer/chip-seq/tf taget gene/three_regions/H3K9me3_rpm_1500bp.txt",sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
H3K9me3_data <- subset(H3K9me3_rpm, H3K9me3_rpm$gene_name %in% row_order)
#H3K9me3_data <- H3K9me3_rpm[which(b %in% H3K9me3_rpm$`gene name` ),]
H3K9me3_df <- data.frame(H3K9me3_data[,2:16])
rownames(H3K9me3_df) <- H3K9me3_data$gene_name

H3K9me3_df[,"week0"] <- rowMeans(H3K9me3_df[,1:3])
H3K9me3_df[,"week2"] <- rowMeans(H3K9me3_df[,4:6])
H3K9me3_df[,"week4"] <- rowMeans(H3K9me3_df[,7:9])
H3K9me3_df[,"week7"] <- rowMeans(H3K9me3_df[,10:12])
H3K9me3_df[,"week10"] <- rowMeans(H3K9me3_df[,13:15])

rna_df <- df[rownames(H3K27ac_df),]
rna_df  <- as.matrix(rna_df[,16:20])
cor_H3K27ac <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], as.matrix(H3K27ac_df[,16:20])[i, ] ))
cor_H3K4me1 <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], as.matrix(H3K4me1_df[,16:20])[i, ] ))
cor_H3K4me3 <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], as.matrix(H3K4me3_df[,16:20])[i, ] ))
cor_H3K27me3 <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], as.matrix(H3K27me3_df[,16:20])[i, ] ))
cor_H3K9me3 <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], as.matrix(H3K9me3_df[,16:20])[i, ] ))
cor_list <- as.matrix(cbind(cor_H3K27ac,cor_H3K4me1,cor_H3K4me3,cor_H3K27me3,cor_H3K9me3))
rownames(cor_list) = rownames(rna_df)

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
left_annotation <- rowAnnotation(
  gene = anno_mark(at=c(3,37,45,84,118,126,135,142,158,159,192,196,229,232,234,241),
                   labels = c("Il1r1","Ptx3","Vcam1","Aqp1","Maml2","Esr1","Phlda1","Sqstm1","Aoc3","Fzd2","Lifr","Enpp2","Mcc",
                   "Pdgfrb","Cd248","Il13ra2"),
                   labels_gp = gpar(col= "navy",fontsize=14),
                   side = "left"),
  unknown = anno_mark(at = c(16,32,55,67,90,168,179,183,187,193,194,223),
                      labels = c("Gorab","Sdc4","Plpp3","Nipal1","C1s1","Serpina3f","Lhfpl2","Fst","Zfp395","Gdnf","Zfpm2","Lbh"),
                      labels_gp = gpar(col="pink",fontsize=14),side = "left"))
  #row = anno_text(rownames(rna_df),location = 1,just="right",gp = gpar(fontsize=4,fontfamily="sans",col="navy")))


ht <- Heatmap(as.matrix(rna_df), column_title = "RNA",border = TRUE,
              cluster_rows = FALSE,cluster_columns = FALSE,row_names_centered = FALSE,
              show_row_names = FALSE,show_column_names = FALSE,
              show_row_dend = TRUE,
              row_names_gp = gpar(fontsize=4),
              col = colorRamp2(c(-1,0,1),c("navy","white","firebrick3")),
              row_order = rev(row_order),
              row_km = 3,
              row_title = "cluster%s",
              row_title_gp = gpar(font=2,fontsize=18),
              column_title_gp = gpar(fontsize=18),
              top_annotation = top_annotation,
              left_annotation = left_annotation,
              heatmap_legend_param = list(title = "RNA"),
              height = 5,
              width = unit(3,"cm"))
ht1 <- Heatmap(cor_list,column_title = "Proximal",
               col = colorRamp2(c(-1,0,1),c("#377EB8","white","#E41A1C")),
               cluster_columns = FALSE,cluster_rows = FALSE,
               show_row_names = FALSE,show_column_names = FALSE,
               row_names_gp = gpar(fontsize=6),
               #row_order = rev(b),
               top_annotation = histone_top_annotation,
               width = unit(3,"cm"),
               heatmap_legend_param = list(title = "Correlation"))
draw(ht ,row_km=3, merge_legend=TRUE)
+ ht1
# histone modifiction upstream and downstream 1.5kb - 10kb  ---------------
Markers = c("H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me3")
for (marker in Markers){
  file <- read.table(paste("~/project/colon cancer/chip-seq/tf taget gene/three_regions/",marker,"_addrpm_10kb.txt",sep = ""),sep = "\t",header = TRUE,check.names = FALSE)
  #b <- intersect(b,file$gene_name)
  file <- unique(subset(file,file$gene_name %in% row_order))
  df <- as.matrix(file[,2:16])
  row_mean_c <- apply(df[,1:3],1,mean)
  row_mean_2 <- apply(df[,4:6],1,mean)
  row_mean_4 <- apply(df[,7:9],1,mean)
  row_mean_7 <- apply(df[,10:12],1,mean)
  row_mean_10 <- apply(df[,13:15],1,mean)
  data <- cbind(row_mean_c,row_mean_2,row_mean_4,row_mean_7,row_mean_10)
  rownames(data) <- unique(file$gene_name)
  write.table(data,file = paste(marker,'_addrpm_10kb_mean.txt',sep = ""),sep = "\t",quote = FALSE,row.names = TRUE)
}

H3K27ac <- as.matrix(read.table('H3K27ac_addrpm_10kb_mean.txt',sep = "\t",header = TRUE,check.names = FALSE))
H3K4me1 <- as.matrix(read.table('H3K4me1_addrpm_10kb_mean.txt',sep = "\t",header = TRUE,check.names = FALSE))
H3K4me3 <- as.matrix(read.table('H3K4me3_addrpm_10kb_mean.txt',sep = "\t",header = TRUE,check.names = FALSE))
H3K27me3 <- as.matrix(read.table('H3K27me3_addrpm_10kb_mean.txt',sep = "\t",header = TRUE,check.names = FALSE))
H3K9me3 <- as.matrix(read.table('H3K9me3_addrpm_10kb_mean.txt',sep = "\t",header = TRUE,check.names = FALSE))

cor_H3K27ac <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], H3K27ac[i, ] ))
cor_H3K4me1 <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], H3K4me1[i, ] ))
cor_H3K4me3 <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], H3K4me3[i, ] ))
cor_H3K27me3 <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], H3K27me3[i, ] ))
cor_H3K9me3 <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], H3K9me3[i, ] ))
cor_list <- as.matrix(cbind(cor_H3K27ac,cor_H3K4me1,cor_H3K4me3,cor_H3K27me3,cor_H3K9me3))
names(cor_list) = rownames(rna_df)

ht2 <- Heatmap(cor_list, column_title = "Flanking",
               col = colorRamp2(c(-1,0,1),c("#377EB8","white","#E41A1C")),
               cluster_columns = FALSE,cluster_rows = FALSE,
               show_row_names = FALSE,show_column_names = FALSE,
               row_names_gp = gpar(fontsize=4),
               #row_order = rev(b),
               top_annotation = histone_top_annotation,
               height = 5,
               width = unit(2.5,"cm"),show_heatmap_legend = FALSE,
               heatmap_legend_param = list(title = "Correlation"))
draw(ht+ht1+ht2,merge_legend=TRUE)


# histone modification upstream and downstream 10kb - 100kb ---------------
for (marker in Markers){
  file <- read.table(paste("~/project/colon cancer/chip-seq/tf taget gene/three_regions/",marker,"_addrpm_100kb.txt",sep = ""),sep = "\t",header = TRUE,check.names = FALSE)
  file <- unique(subset(file,file$gene_name %in% row_order))
  df <- as.matrix(file[,2:16])
  row_mean_c <- apply(df[,1:3],1,mean)
  row_mean_2 <- apply(df[,4:6],1,mean)
  row_mean_4 <- apply(df[,7:9],1,mean)
  row_mean_7 <- apply(df[,10:12],1,mean)
  row_mean_10 <- apply(df[,13:15],1,mean)
  data <- cbind(row_mean_c,row_mean_2,row_mean_4,row_mean_7,row_mean_10)
  rownames(data) <- unique(file$gene_name)
  write.table(data,file = paste(marker,'_addrpm_100kb_mean.txt',sep = ""),sep = "\t",quote = FALSE,row.names = TRUE)
}

H3K27ac <- as.matrix(read.table('H3K27ac_addrpm_100kb_mean.txt',sep = "\t",header = TRUE,check.names = FALSE))
H3K4me1 <- as.matrix(read.table('H3K4me1_addrpm_100kb_mean.txt',sep = "\t",header = TRUE,check.names = FALSE))
H3K4me3 <- as.matrix(read.table('H3K4me3_addrpm_100kb_mean.txt',sep = "\t",header = TRUE,check.names = FALSE))
H3K27me3 <- as.matrix(read.table('H3K27me3_addrpm_100kb_mean.txt',sep = "\t",header = TRUE,check.names = FALSE))
H3K9me3 <- as.matrix(read.table('H3K9me3_addrpm_100kb_mean.txt',sep = "\t",header = TRUE,check.names = FALSE))

cor_H3K27ac <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], H3K27ac[i, ] ))
cor_H3K4me1 <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], H3K4me1[i, ] ))
cor_H3K4me3 <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], H3K4me3[i, ] ))
cor_H3K27me3 <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], H3K27me3[i, ] ))
cor_H3K9me3 <- sapply(seq.int(dim(rna_df)[1]), function(i) cor(rna_df[i, ], H3K9me3[i, ] ))
#cor_list <- as.matrix(cbind(cor_H3K27ac,cor_H3K4me1,cor_H3K4me3,cor_H3K27me3,cor_H3K9me3))
cor_list <- as.matrix(cbind(cor_H3K27ac,cor_H3K4me1,cor_H3K4me3))
rownames(cor_list) = rownames(rna_df)

ht3 <- Heatmap(as.matrix(cor_list), column_title = "Cor",border = FALSE,
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
pdf("tf target gene.pdf",height = 10,width = 10)
draw(ht+ht3,merge_legend=TRUE)
dev.off()


# chromHMM state ----------------------------------------------------------
#state <- read.csv("~/project/colon cancer/chip-seq/luo_figure/gene_base_summary_merge_rna.csv",check.names = FALSE)
state <- read.csv("~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody2000bp/mm10_genebody_2000bp_merge_rna.csv",check.names = FALSE)
state <- na.omit(state)

df_sub = state[,c("gene_name", "week0", "week2", "week4", "week7", "week10")]
state = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13")
element = c("enhancer", "enhancer", "enhancer", "promoter", "promoter", "promoter", "promoter", "promoter", "repressed", "repressed",
            "enhancer", "quiescent", "heterochromatin")
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
data_frame <- chromstate[which(chromstate$symbol %in% row_order),]
data_frame <- unique(inner_join(data_frame,df_sub,by="ensembl"))
data_frame <- data_frame[-c(45,77),]
rownames(data_frame) <- data_frame$symbol
data_frame <- data_frame[,3:7]
data_frame <- data_frame[rownames(rna_df),]

ht_state <- Heatmap(matrix = as.matrix(data_frame),name = "State",border = FALSE,
                    show_row_names = FALSE,show_column_names = FALSE,
                    row_order = rev(row_order),
                    top_annotation = top_annotation,
                    col = c("enhancer" = enhancer_color,"promoter" = promoter_color,"repressed" = repressed_color,
                            "quiescent" = quit_color,"heterochromatin" = heterochromatin_color),
                    column_title = "State",
                    column_title_gp = gpar(fontsize=18,fontfamily="sans"),
                    height = 5,
                    width = unit(2.5,"cm")
)
pdf("~/project/colon cancer/nfkb target gene/nfkb target gene with distal region add genename2.pdf",height = 7,width = 8)
draw(ht+ht3+ht_state,merge_legend=TRUE,heatmap_legend_side="right",annotation_legend_side="right",column_title_gp=gpar(fontsize=18,fontfamily="sans"))
dev.off()

draw(ht_state)
a <- row_order(ht)
b <- data.frame(ht@row_order)
a1 <- ht@row_names_param$labels[unlist(a[1])]
a2 <- ht@row_names_param$labels[unlist(a[2])]
a3 <- ht@row_names_param$labels[unlist(a[3])]
write.table(x = a1,file = "~/project/colon cancer/nfkb target gene/nfkb_cluster1_gene.txt",sep = "\t",quote = FALSE)
write.table(x = a2,file = "~/project/colon cancer/nfkb target gene/nfkb_cluster2_gene.txt",sep = "\t",quote = FALSE)
write.table(x = a3,file = "~/project/colon cancer/nfkb target gene/nfkb_cluster3_gene.txt",sep = "\t",quote = FALSE)

# barplot -----------------------------------------------------------------
nfkb_target_gene <- read.delim("~/project/colon cancer/nfkb target gene/nfkb target gene 241.txt",sep = "\t",col.names = "gene_name")
data_frame <- df_sub[df_sub$gene_name %in% nfkb_target_gene$gene_name,]

infla <- data_frame[data_frame$week2 == "enhancer" & data_frame$week4 == "enhancer",]
rownames(infla) <- infla$gene_name
infla <- as.matrix(infla[-1])
ht_state <- Heatmap(matrix = infla,name = "State",
                    show_row_names = TRUE,show_column_names = FALSE,
                    #top_annotation = top_annotation,
                    col = c("enhancer" = enhancer_color,"promoter" = promoter_color,"repressed" = repressed_color,
                            "quiescent" = quit_color,"heterochromatin" = heterochromatin_color),
                    column_title = "chromHMM states",
                    column_title_gp = gpar(fontsize=12,fontfamily="sans"),
                    height = 5,
                    #width = unit(2.5,"cm")
)
draw(ht_state,merge_legend=TRUE)

library(dplyr)
rownames(data_frame) <- 1:nrow(data_frame)
data_frame <- data_frame[-157,]
rownames(data_frame) <- data_frame$gene_name
data_frame <- data_frame[,-1]
colnames(data_frame) <- c("week0","week2","week4","week7","week10")
table_0 <- data_frame %>% count(week0)
colnames(table_0) <- c("state","num")
table_2 <- data_frame %>% count(week2)
colnames(table_2) <- c("state","num")
table_4 <- data_frame %>% count(week4)
colnames(table_4) <- c("state","num")
table_7 <- data_frame %>% count(week7)
colnames(table_7) <- c("state","num")
table_10 <- data_frame %>% count(week10)
colnames(table_10) <- c("state","num")
#table_10 <- table_10[-2,]
table <- rbind(table_0,table_2,table_4,table_7,table_10)
table <- cbind(table,time=rep(Times,each=4))

table$time <- factor(rep(Times,each=4),levels = Times)
ggplot(table, aes(x = state,y = num, fill=time)) +
  geom_bar(position = position_dodge(0.9),stat = "identity",width = 0.9) +
  geom_text(aes(label=num),vjust=1.2,position = position_dodge(0.9),size=4,color="white") +
  #geom_line(aes(group=type),stat = "identity",position = position_dodge(0.9)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,150)) +
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
ggsave("~/project/colon cancer/nfkb target gene/nfkb target gene 241 barplot.pdf",height = 3,width = 8)

cluster4 <- read.table("~/project/colon cancer/chip-seq/maSigPro_lily/control/H3K27ac/H3K27ac_cluster4.txt")
cluster2 <- read.table("~/project/colon cancer/RNA-seq/edgeR/masigpro_edgeR_filter/RNA-seq_masigpro_genelist_cluster2.txt")
intersect(cluster4$cluster4,rownames(data_frame))
intersect(cluster2$cluster2,rownames(infla))


# hypergenometric test ----------------------------------------------------
library(stats)
p1 <- phyper(q = 298-1,m = 16402,n =35815,k = 663,lower.tail = F)
p2 <- phyper(q = 300-1,n = 16402,m = 35815,k = 480,lower.tail = F)
p3 <- phyper(q = 54-1,n = 16402,m = 35815,k = 161,lower.tail = F)
p4 <- phyper(q = 300-1,n = 16402,m = 35815,k = 614,lower.tail = F)

pnfkb <- phyper(q = 132-1,m = 16402,n = 35815,k = 241,lower.tail = F)

# correlation between histone modification ---------------------------------
summary(cor_list)
apply(cor_list,2,mean)
apply(cor_list,2,sd)
apply(cor_list,2, var)

cor_data <- data.frame(cor_list)
cor_data$symbol <- rownames(cor_data) 
library(reshape2)
cor_df <- melt(cor_data,id.vars = "symbol",variable.name = "sample",value.name = "cor")

library(dplyr)
cor_df1 <- cor_df %>% group_by(sample) %>% mutate(med=median(cor))

library(ggsignif)
compaired <- list(c("cor_H3K27ac","cor_H3K4me1"),c("cor_H3K27ac","cor_H3K4me3"),
                  c("cor_H3K27ac","cor_H3K27me3"), c("cor_H3K27ac","cor_H3K9me3"))
my_compaired <- list(c("cor_H3K27ac","cor_H3K4me1"),c("cor_H3K4me1","cor_H3K4me3"),c("cor_H3K27ac","cor_H3K4me3"))
p + geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = T,test = t.test) 

library(ggpubr)

ggplot(cor_df1,aes(x = sample,y = cor,fill=sample)) +
  geom_violin() +
  geom_boxplot(width=.1,fill="black",outlier.colour = NA) +
  stat_summary(fun = median,geom = "point",fill="white",shape=21,size=2.5) +
  scale_fill_manual(values = c("#E41A1C","#FF7F00","#4DAF4A"),guide=FALSE) +
  scale_x_discrete(expand = c(0.2,0.2),
                   breaks = c("cor_H3K27ac","cor_H3K4me1","cor_H3K4me3"),
                   labels = c("H3K27ac","H3K4me1","H3K4me3")) +
  theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) +
  labs(x=NULL,y="Correlation",title = "Distal region (n=241)") +
  theme(
    plot.background= element_rect(size = 3),
    plot.margin = margin(t=0.5,r = 0.5,b = 0.5,l = 0.5,unit = "cm"),
    plot.title = element_text(colour = "black",size = 18,hjust = 0.5),
    axis.title.y = element_text(colour = "black",size = 18),
    axis.text = element_text(colour = 'black',size = 18)
    #axis.text.x = element_text(angle = 45,hjust = 1)
  ) +
  stat_compare_means(comparisons = my_compaired,method = "t.test",step.increase = 0.07,label = "p.signif",paired = TRUE)

ggplot(cor_df1,aes(x = sample,y = cor,color=sample)) +
  stat_boxplot(geom = "errorbar",linetype=1,width=0.4,position = "identity") +
  geom_boxplot(outlier.fill = NA,outlier.shape=NA,width=0.6) +
  geom_point(alpha=0.5,position = position_jitter(width = .25,height = 0.06),size=2,shape=19) + #position = position_jitter(width = 0.3,height = 0)) +
  #geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = T,test = t.test) +
  scale_color_manual(values = c("#E41A1C","#FF7F00","#4DAF4A"),guide=FALSE) +
  scale_y_continuous(limits = c(-1,1.5),breaks = c(-1,0,1)) +
  scale_x_discrete(expand = c(0.2,0.2),
                   breaks = c("cor_H3K27ac","cor_H3K4me1","cor_H3K4me3"),
                   labels = c("H3K27ac","H3K4me1","H3K4me3")) +
  theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) +
  labs(x=NULL,y="Correlation",title = "Distal region (n=241)") +
  theme(
    plot.background= element_rect(size = 3),
    plot.margin = margin(t=0.5,r = 0.5,b = 0.5,l = 0.5,unit = "cm"),
    plot.title = element_text(colour = "black",size = 18,hjust = 0.5),
    axis.title.y = element_text(colour = "black",size = 18),
    axis.text = element_text(colour = 'black',size = 18)
    #axis.text.x = element_text(angle = 45,hjust = 1)
  ) +
  stat_compare_means(comparisons = my_compaired,method = "t.test",step.increase = 0.07,label = "p.signif",paired = TRUE) 
ggsave(filename = "~/project/colon cancer/nfkb target gene/nfkb target gene distal region cor.pdf",height = 4,width = 4.7)

ggplot(cor_df1,aes(x = sample,y = cor,color=sample)) +
  stat_boxplot(geom = "errorbar",linetype=1,width=0.4,position = "identity") +
  geom_boxplot(outlier.fill = NA,outlier.shape=NA,width=0.6) +
  geom_point(alpha=0.5,position = position_jitter(width = .25,height = 0.06),size=2,shape=19) + #position = position_jitter(width = 0.3,height = 0)) +
  #geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = T,test = t.test) +
  stat_compare_means(comparisons = compaired,method = "t.test",step.increase = 0.07,label = "p.signif",paired = TRUE) +
  scale_color_manual(values = c("#E41A1C","#FF7F00","#4DAF4A","#377EB8","#984EA3"),guide=FALSE) +
  #scale_y_continuous(expand = c(0.01,0.1)) +
  scale_x_discrete(expand = c(0.2,0.2),
                   breaks = c("cor_H3K27ac","cor_H3K4me1","cor_H3K4me3","cor_H3K27me3","cor_H3K9me3"),
                   labels = c("H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me3")) +
  theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) +
  labs(x=NULL,y="Correlation",title = "Flanking region (n=241)") +
  theme(
    #plot.margin = margin(t=0.5,r = 0.1,b = 0.1,l = 0.1,unit = "cm"),
    plot.title = element_text(colour = "black",size = 18,hjust = 0.5),
    axis.title.y = element_text(colour = "black",size = 18),
    axis.text = element_text(colour = 'black',size = 18),
    axis.text.x = element_text(angle = 45,hjust = 1)
    )
ggsave(filename = "~/project/colon cancer/nfkb target gene/nfkb target gene flanking region cor.pdf",height = 4.6,width = 6)

dnorm(cor_data$cor_H3K27ac)
plot(x = cor_data$cor_H3K27ac,y = dnorm(cor_data$cor_H3K27ac))
y <- quantile(cor_data$cor_H3K27ac,c(.75,.5,.25,0))

hist(cor_data$cor_H3K27ac,breaks = 10,col = "red")
rug(jitter(cor_data$cor_H3K27ac),amount=0.01)
lines(density(cor_data$cor_H3K27ac),col="blue",lwd=2)

plot(density(cor_data$cor_H3K27ac))
rug(cor_data$cor_H3K27ac,col = "brown")
boxplot(cor_data,main="Distal region")
t.test(cor_data$cor_H3K27ac,cor_data$cor_H3K4me1,paired = FALSE)

# select correlation bigger than 0.6 --------------------------------------
corright <- cor_df[cor_df$cor > 0 & cor_df$sample == "cor_H3K27ac",]
write.table(x = corright,file = "~/project/colon cancer/nfkb target gene/nfkb distal correlation bigger than 0 gene.txt",sep = "\t",quote = FALSE)

corright <- read.table("~/project/colon cancer/nfkb target gene/nfkb distal correlation bigger than 0.6 gene.txt",sep = "\t")
corright <- corright[corright$cor > 0.85 & corright$sample == "cor_H3K27ac",]
report <- corright$symbol[-c(2,3,6,7,9,16:19,21,22,24)]
report <- ht@row_names_param$labels[report]
ht@row_names_param$labels[81]
ht@row_order
row_order(ht)
unknown <- corright$symbol[c(2,3,6,7,9,16:19,21,22,24)] 

library(plotrix)
pdf('~/project/colon cancer/nfkb target gene/nfkb 138 gene cor distribution.pdf',height = 4,width=4)
pie(x = c(16,12,41,69), 
    col = c("#3182BD","#6BAED6","#BDD7E7","#DEEBF7"),border = "white")
legend('topright',c("cor<=0.6","cor>0.6","no report","report"),cex = 0.7,fill = rev(c("#3182BD","#6BAED6","#BDD7E7","#DEEBF7")),border = FALSE)
dev.off()

# nfkb 241 gene FPKM ------------------------------------------------------
rm(list = ls())
plotBox <- function(filename,number) {
  fpkm <- read.csv("~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody2000bp/mm10_genebody_2000bp_merge_rna.csv",check.names = FALSE)
  file <- read.delim(filename,sep = "\t",header = FALSE,col.names = "gene_name")
  
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
  
  df$time <- factor(rep(c("control","2week","4week","7week","10week"),each=nrow(data)*3),levels = c("control","2week","4week","7week","10week"))
  
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
    labs(title = paste("nfkb signal gene (",number,")",sep = ""),x=NULL,y="Expression (FPKM)",fill=NULL) +
    theme(#aspect.ratio = 1/0.75,
      plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,family = "sans"),
      axis.text.x = element_text(colour = "black",size = 18,family = "sans",angle = 45,hjust = 1),
      axis.text.y = element_text(size = 18,colour = "black",margin = margin(l = 0.1,r = 0.1,unit = "cm")),
      axis.ticks.length = unit(2,'mm')
      #panel.border = element_rect(size = 1.1),
      #panel.grid = element_blank()
    ) +
    guides(fill=FALSE)
  ggsave(paste("~/project/colon cancer/nfkb target gene/nfkb ",number," gene AOMDSS FPKM.pdf",sep = ""),p,height = 4.5,width = 6)
  p
}

plotBox(filename = "~/project/colon cancer/nfkb target gene/nfkb_target_gene_138.txt",number = 138)
plotBox(filename = "~/project/colon cancer/nfkb target gene/nfkb target gene 241.txt",number = 241)

write.table(file$gene_name,"~/project/colon cancer/nfkb target gene/nfkb_target_gene_138.txt",col.names = FALSE,row.names = FALSE,quote = FALSE)


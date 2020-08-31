rm(list = ls())
setwd('~/project/colon cancer/RNA-seq/DESeq2_maSigPro/')
file_path <- "~/project/colon cancer/RNA-seq/DESeq2_maSigPro/"
library(DESeq2)
require(BiocParallel)
library(ggplot2)

inputf <- read.table('merge_htseqCount.txt',header = T,row.names = 1,check.names = FALSE,stringsAsFactors = FALSE)
data <- inputf[rowSums(inputf > 3) >2,]
Times <- c("control","2weeks","4weeks","7weeks","10weeks")
sample_total <- c()
for (time in c("control","2weeks","4weeks","7weeks","10weeks")){
  for (rep in c(1:3)){
    sample <- paste(time,rep,sep = "_rep")
    sample_total <- c(sample_total,sample)
  }
}
names(data) <- sample_total
condition <- factor(rep(Times,each=3),levels = Times)
cData <- data.frame(row.names = colnames(data),condition)
dds <- DESeqDataSetFromMatrix(countData = data,colData = cData,design = ~ condition)
dds <- DESeq(dds)
dds

library(biomaRt)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "mmusculus_gene_ensembl")

# every weeks vs control --------------------------------------------------
Degs <- function(p) {
  for (time in Times[2:5]) {
    res <- results(dds,contrast = c("condition",time,"control"))
    res <- res[order(res$pvalue),]
    summary(res)
    sig <- res[!is.na(res$pvalue) & res$pvalue < p,]
    sig.deseq <- rownames(sig)
    gene_list = paste(file_path,time,"/res_control_",time,".txt",sep = "")
    #write.table(x = sig,file = gene_list,sep = '\t',row.names = TRUE,quote = FALSE)
    resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized = TRUE)),by='row.names',sort=FALSE)
    names(resdata)[1] <- "Gene_id"
    resdata_name <- paste(file_path,time,"/res_control_",time,".txt",sep = "")
    write.table(resdata,resdata_name, quote = F,col.names = T,row.names = F,sep = "\t")
  
    table(res$padj < p)
    resdata$change <- as.factor(
      ifelse(
        resdata$pvalue < p & abs(resdata$log2FoldChange) > 2,
        ifelse(resdata$log2FoldChange > 2,'Up','Down'),
        'Nd'
      )
    )
    a <- as.data.frame(table(resdata$change))
    
    high <- resdata[which(resdata$pvalue < p & resdata$change=='Up'),][,1]
    low <- resdata[which(resdata$pvalue < p & resdata$change=='Down'),][,1]
    high <- substr(high,1,18)
    low <- substr(low,1,18)
    high_gene_name <- getBM(mart = mart,attributes = c("external_gene_name","entrezgene_id","ensembl_gene_id"),filters = "ensembl_gene_id",values = high)
    low_gene_name <- getBM(mart = mart,attributes = c("external_gene_name","entrezgene_id","ensembl_gene_id"),filters = "ensembl_gene_id",values = low)
    high_gene_list <- paste(file_path,time,"/high_genelist_",time,p,".txt",sep = "")
    low_gene_list <- paste(file_path,time,"/low_genelist_",time,p,".txt",sep = "")
    write.table(x = high_gene_name,file = high_gene_list,quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(x = low_gene_name,file = low_gene_list,quote = FALSE,row.names = FALSE,col.names = FALSE)
    
    if (time %in% c("2weeks","4weeks")) {
      ymax <- 30
    pdf_name <- paste(file_path,time,"/volcano_plot_",p,".pdf",sep = "")
    attach(resdata)
    ggplot(data = resdata,aes(x = log2FoldChange, y = -log10(pvalue),color=change)) +
      geom_point(alpha=0.5,size =2,shape=19) +
      xlim(-15,15) +
      theme_bw() +
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",lineend = "square")
      ) +
      scale_y_continuous(expand = c(0,0),limits = c(0,ymax)) +
      ggtitle(paste("Different expression genes between",time,"and control",sep = " ")) +
      labs(y="-log10(Pvalue)",color=NULL) +
      theme(
        plot.title = element_text(size = 14,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust =5),
        axis.title.x = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust = -1,color = "black"),
        axis.title.y = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust = 3,color = "black"),
        axis.text = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust = 0.5,color = "black")
      ) + 
      theme(legend.text = element_text(family = 'Helvetica',face = 'plain',size = 12)) + 
      theme(legend.position = c(1,1),legend.justification = c(1,1),legend.background = element_blank()) +
      #scale_color_brewer(palette = "Set1",limits=c('Up','Down','Nd')) +
      scale_colour_manual(name='',values = c("#E41A1C","#377EB8","grey90"),limits=c('Up','Down','Nd'),
                          labels=c(paste('Up:',a[3,2],sep = " "),paste('Down:',a[1,2],sep = " "),paste('Nd:',a[2,2],sep = " "))) + 
      geom_vline(xintercept = c(-2,2),lty=2,col='grey',lwd=0.5) +
      geom_hline(yintercept = -log10(0.01),lty=2,col='grey',lwd=0.5) +
      theme(plot.margin = unit(c(1,1.5,1,1.5),"cm"))
    ggsave(pdf_name,height = 7,width = 7)
    }
    else {
      ymax <- 40
      pdf_name <- paste(file_path,time,"/volcano_plot_",p,".pdf",sep = "")
      attach(resdata)
      ggplot(data = resdata,aes(x = log2FoldChange, y = -log10(pvalue),color=change)) +
        geom_point(alpha=0.5,size =2,shape=19) +
        xlim(-15,15) +
        theme_bw() +
        theme(
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black",lineend = "square")
        ) +
        scale_y_continuous(expand = c(0,0),limits = c(0,ymax)) +
        ggtitle(paste("Different expression genes between",time,"and control",sep = " ")) +
        labs(y="-log10(Pvalue)",color=NULL) +
        theme(
          plot.title = element_text(size = 14,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust =5),
          axis.title.x = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust = -1,color = "black"),
          axis.title.y = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust = 3,color = "black"),
          axis.text = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust = 0.5,color = "black")
        ) + 
        theme(legend.text = element_text(family = 'Helvetica',face = 'plain',size = 12)) + 
        theme(legend.position = c(1,1),legend.justification = c(1,1),legend.background = element_blank()) +
        #scale_color_brewer(palette = "Set1",limits=c('Up','Down','Nd')) +
        scale_colour_manual(name='',values = c("#E41A1C","#377EB8","grey90"),limits=c('Up','Down','Nd'),
                            labels=c(paste('Up:',a[3,2],sep = " "),paste('Down:',a[1,2],sep = " "),paste('Nd:',a[2,2],sep = " "))) + 
        geom_vline(xintercept = c(-2,2),lty=2,col='grey',lwd=0.5) +
        geom_hline(yintercept = -log10(0.01),lty=2,col='grey',lwd=0.5) +
        theme(plot.margin = unit(c(1,1.5,1,1.5),"cm"))
      ggsave(pdf_name,height = 7,width = 7)
    }
  }
} 


Degs(p = 0.01)
Degs(p = 0.05)



# cancer vs inflammation --------------------------------------------------
library(biomaRt)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "mmusculus_gene_ensembl")

levels <- c(rep("control",3),rep(c("inflammation","cancer"),each=6))
sample_total <- c()
for (level in c("inflammation","cancer")) {
  for (i in 1:6) {
    sample <- paste(level,i,sep = "_rep")
    sample_total <- c(sample_total,sample)
  }
}
sample_total <- c("control_rep1","control_rep2","control_rep3",sample_total)
names(data) <- sample_total
condition <- factor(levels,levels = c("control","inflammation","cancer"))
cData <- data.frame(condition,row.names = colnames(inputf))
dds <- DESeqDataSetFromMatrix(countData = data,colData = cData,design = ~ condition)
dds <- DESeq(dds)
dds

degs_process <- function(level,control,p,ymax,dir) {
  res <- results(dds,contrast = c("condition",level,control))
  res <- res[order(res$pvalue),]
  summary(res)
  sig <- res[!is.na(res$pvalue) & res$pvalue < 0.01,]
  sig.deseq <- rownames(sig)
  gene_list = paste(file_path,dir,"/res_control_",level,".txt",sep = "")
  write.table(x = sig,file = gene_list,sep = '\t',row.names = TRUE,quote = FALSE)
  resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized = TRUE)),by='row.names',sort=FALSE)
  names(resdata)[1] <- "Gene_id"
  resdata_name <- paste(file_path,dir,"/res_control_",level,".txt",sep = "")
  write.table(resdata,resdata_name, quote = F,col.names = T,row.names = F,sep = "\t")
  
  table(res$padj < p)
  resdata$change <- as.factor(
    ifelse(
      resdata$pvalue < p & abs(resdata$log2FoldChange) > 2,
      ifelse(resdata$log2FoldChange > 2,'Up','Down'),
      'Nd'
    )
  )
  a <- as.data.frame(table(resdata$change))
  
  high <- resdata[which(resdata$pvalue < p & resdata$change=='Up'),][,1]
  low <- resdata[which(resdata$pvalue < p & resdata$change=='Down'),][,1]
  high <- substr(high,1,18)
  low <- substr(low,1,18)
  high_gene_name <- getBM(mart = mart,attributes = c("external_gene_name","entrezgene_id","ensembl_gene_id"),filters = "ensembl_gene_id",values = high)
  low_gene_name <- getBM(mart = mart,attributes = c("external_gene_name","entrezgene_id","ensembl_gene_id"),filters = "ensembl_gene_id",values = low)
  high_gene_list <- paste(file_path,dir,"/high_genelist_",level,p,".txt",sep = "")
  low_gene_list <- paste(file_path,dir,"/low_genelist_",level,p,".txt",sep = "")
  write.table(x = high_gene_name,file = high_gene_list,quote = FALSE,row.names = FALSE,col.names = FALSE)
  write.table(x = low_gene_name,file = low_gene_list,quote = FALSE,row.names = FALSE,col.names = FALSE)
  
  pdf_name <- paste(file_path,dir,"/volcano_plot_",level,p,".pdf",sep = "")
  attach(resdata)
  ggplot(data = resdata,aes(x = log2FoldChange, y = -log10(pvalue),color=change)) +
    geom_point(alpha=0.5,size =2,shape=19) +
    xlim(-15,15) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black",lineend = "square")
    ) +
    scale_y_continuous(expand = c(0,0),limits = c(0,ymax)) +
    ggtitle(paste("Different expression genes between",level,"and ",control,sep = " ")) +
    labs(y="-log10(Pvalue)",color=NULL) +
    theme(
      plot.title = element_text(size = 14,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust =5),
      axis.title.x = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust = -1,color = "black"),
      axis.title.y = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust = 3,color = "black"),
      axis.text = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust = 0.5,color = "black")
    ) + 
    theme(legend.text = element_text(family = 'Helvetica',face = 'plain',size = 12)) + 
    theme(legend.position = c(1,1),legend.justification = c(1,1),legend.background = element_blank()) +
    #scale_color_brewer(palette = "Set1",limits=c('Up','Down','Nd')) +
    scale_colour_manual(name='',values = c("#E41A1C","#377EB8","grey90"),limits=c('Up','Down','Nd'),
                        labels=c(paste('Up:',a[3,2],sep = " "),paste('Down:',a[1,2],sep = " "),paste('Nd:',a[2,2],sep = " "))) + 
    geom_vline(xintercept = c(-2,2),lty=2,col='grey',lwd=0.5) +
    geom_hline(yintercept = -log10(0.01),lty=2,col='grey',lwd=0.5) +
    theme(plot.margin = unit(c(1,1.5,1,1.5),"cm"))
  ggsave(pdf_name,height = 7,width = 7)
}

model <- c("inflammation","cancer")
degs_process(level = "inflammation",p = 0.01,control = "control",ymax = 30,dir = level)
degs_process(level = "cancer",p = 0.01,control = "control",ymax = 40,dir = level)
degs_process(level = "cancer",control = "inflammation",p = 0.01,ymax = 50,dir = "inflammation vs cancer")


# Visualization -----------------------------------------------------------
# heatmap correlation
vsd <- getVarianceStabilizedData(dds)
heatmap(cor(vsd),cexCol=0.75,cexRow=0.75,margins = c(5,5))

# heatmap gene expression
library(pheatmap)
df <- as.data.frame(resdata[which(resdata$pvalue < 0.01 & resdata$change=='Up'),])
row.names(df) <- df$Gene_id
df <- df[,-c(2:7,23)]
df <- df[,-1]
df <- df[rowSums(df)>0,]
pheatmap(df,scale = "row",cluster_cols = F,cluster_rows = T,show_rownames = F,angle_col = 45,
         color = colorRampPalette(rev(brewer.pal(n=3,name = "RdYlBu")))(100)) 

library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(c(-4,0,4),c("lightblue","white","red"))
Heatmap(df,cluster_rows = FALSE,cluster_columns = FALSE,show_column_dend = FALSE,col = col_fun)


# PCA analysis
library(RColorBrewer)
pr <- prcomp(t(vsd))
# first
plot(pr$x,col="white",main="PC plot",xlim=c(-80,80),tck=-0.01,las=2,mgp=c(1.5,0.5,0))
points(pr$x[,1],pr$x[,2],pch=19,type = "p",col=rep(brewer.pal(5,"Dark2"),each=3),
       cex=1.5)
# second
plot(pr$x,col="white",main="PC plot",xlim=c(-100,100),ylim=c(-60,60),tck=-0.01,mar=c(2,2.5,2,2.5),axes=FALSE,mgp=c(2,0.5,0))
axis(side=1,las=0,tck=-0.01,pos =-61,at = seq(-100,100,100),labels = seq(-100,100,100),mgp=c(1,0.5,0))
axis(side=2,las=2,tck=-0.01,pos =-102,at = seq(-60,60,60),labels = seq(-60,60,60),mgp=c(1,0.5,0))
points(pr$x[,1],pr$x[,2],pch=19,type = "p",col=rep(brewer.pal(5,"Dark2"),each=3),
     cex=2,alpha=0.5)
legend("topright",inset = 0.03,Times,pch = 19,col = brewer.pal(5,"Dark2"),cex = 0.8,adj = c(0,0.5),text.col = brewer.pal(5,"Dark2"))


biplot(pr,cex=c(1,0.5),main="Biplot",
       col=c("black","grey"))

# volcano plot shows the p-value against the fold change in each gene
# first
plot(res$log2FoldChange,-log(res$padj),pch=15)
points(sig$log2FoldChange,-log(sig$padj),
       col="grey",pch=15)
library("calibrate")# if not installed,run 'install.packages("calibrate")
#textxy(sig$log2FoldChange,-log(sig$padj),rownames(sig),cex=0.9)

# second
pdf('volcano plot 0.01.pdf',width = 7,height = 7)
attach(resdata)
ggplot(data = resdata,aes(x = log2FoldChange, y = -log10(pvalue),color=change)) +
  geom_point(alpha=0.5,size =2,shape=19) +
  xlim(-15,15) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black",lineend = "square")
  ) +
  scale_y_continuous(expand = c(0,0),limits = c(0,50)) +
  ggtitle('Different expression genes between 10weeks and control') +
  labs(y="-log10(Pvalue)",color=NULL) +
  theme(
    plot.title = element_text(size = 14,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust =5),
    axis.title.x = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust = -1,color = "black"),
    axis.title.y = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust = 3,color = "black"),
    axis.text = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust = 0.5,color = "black")
    ) + 
  theme(legend.text = element_text(family = 'Helvetica',face = 'plain',size = 12)) + 
  theme(legend.position = c(1,1),legend.justification = c(1,1),legend.background = element_blank()) +
  #scale_color_brewer(palette = "Set1",limits=c('Up','Down','Nd')) +
  scale_colour_manual(name='',values = c("#E41A1C","#377EB8","grey90"),limits=c('Up','Down','Nd'),labels=c('Up: 930','Down: 1257','Nd: 17586')) + 
  geom_vline(xintercept = c(-2,2),lty=2,col='grey',lwd=0.5) +
  geom_hline(yintercept = -log10(0.01),lty=2,col='grey',lwd=0.5) +
  theme(plot.margin = unit(c(1,1.5,1,1.5),"cm"))
dev.off()

# gene ontology analysis


# view Rcolorbrewer color RGB
library(RColorBrewer)
brewer.pal(9,"Set1")

# venn cancer vs inflammation
library(grid)
library(VennDiagram)

high_infla <- read.table("~/project/colon cancer/RNA-seq/DESeq2_maSigPro/inflammation/high_genelist_inflammation0.01.txt",
                         sep = " ",col.names = c("gene_name","entrez","ensembl"))
high_cancer <- read.table("~/project/colon cancer/RNA-seq/DESeq2_maSigPro/cancer/high_genelist_cancer0.01.txt",
                          sep = " ",col.names = c("gene_name","entrez","ensembl"))
a <- intersect(high_infla$gene_name,high_cancer$gene_name)
low_infla <- read.table("~/project/colon cancer/RNA-seq/DESeq2_maSigPro/inflammation/low_genelist_inflammation0.01.txt",
                        sep = " ",col.names = c("gene_name","entrez","ensembl"))
low_cancer <- read.table("~/project/colon cancer/RNA-seq/DESeq2_maSigPro/cancer/low_genelist_cancer0.01.txt",
                          sep = " ",col.names = c("gene_name","entrez","ensembl"))
b <- intersect(low_infla$gene_name,low_cancer$gene_name)
pdf()
venn.plot <- draw.pairwise.venn(area1 = 984,area2 = 457,cross.area = 231,scaled = F,
                                category = c("Cancer\nn=984","Inflammation\nn=457"),
                                col = c("#E41A1C","#377EB8"),
                                #fill = c('#00BFC5','#F8765D'),
                                #fill=c("#E41A1C","#377EB8"),
                                #lty = 'blank',
                                lwd = 6,
                                cex = 2,
                                ext.pos = 30,
                                ext.line.lwd = 2,
                                fontfamily = "serif",
                                fontface = "plain",
                                cat.cex = 2,
                                cat.col = c("#E41A1C","#377EB8"),
                                cat.pos = c(185,170),
                                cat.dist = c(0.08,0.08),
                                cat.fontfamily = 'serif',
                                cat.fontface = 'plain',
                                label.col = c("#E41A1C","black","#377EB8"),
                                margin=0.05,
                                height=3000,wight=3000,
                                resolution=300,
                                alpha = c(0.4,0.4)
                                )
grid.draw(venn.plot)
dev.off()


# select top genes plot expression of each gene 
sig.ordered <- sig[order(sig$padj),]
for(gene in head(rownames(sig.ordered))) {
  boxplot(vsd[gene,which(condition=="control")],vsd[gene, which(condition=="2weeks")],
          vsd[gene, which(condition=="4weeks")],vsd[gene, which(condition=="7weeks")],
          vsd[gene, which(condition=="10weeks")],
          main=paste(gene,signif(sig[gene,"padj"],2)),names=Times)
  readline()
}

for(gene in head(rownames(sig.ordered))) {
  barplot(vsd[gene,],las=2,col=as.numeric(as.factor(condition)),main=gene,cex.names=0.9)
  readline()
}



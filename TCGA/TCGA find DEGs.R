rm(list = ls())
setwd('~/project/colon cancer/TCGA/')
library(ggplot2)
library(BiocParallel)
library(DESeq2)

# data perpare ------------------------------------------------------------
inputf <- read.csv('htseq_counts.txt',sep='\t',header = TRUE,stringsAsFactors = FALSE)
row.names(inputf) <- inputf[,'gene_name']
inputf <- inputf[,c(-1)]
fpkm()
condition <- factor(c(rep('Normal',44),rep('Tumor',571)),
                    levels = c('Normal','Tumor'))
cData <- data.frame(row.names = colnames(inputf),condition)
dds <- DESeqDataSetFromMatrix(inputf,colData = cData,design = ~ condition)
dds <- DESeq(dds)
dds
res <- results(dds,contrast = c('condition','Tumor','Normal'))
res <- res[order(res$pvalue),]
summary(res)
write.table(x = res,file = 'res_Tumor_Normal.txt',sep = '\t',row.names = TRUE,quote = FALSE)
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized = TRUE)),by='row.names',sort=FALSE)
names(resdata)[1] <- "gene_name"
write.table(resdata,'resdata_Tumor_Normal.txt', quote = F,col.names = T,row.names = F,sep = "\t") 
table(res$padj < 0.01)
resdata$change <- as.factor(
  ifelse(
    resdata$pvalue < 0.01 & abs(resdata$log2FoldChange) > 2,
    ifelse(resdata$log2FoldChange > 2,'Tumor','Normal'),
    'Nd'
  )
)
table(resdata$change)
attach(resdata)
ggplot(data = resdata,aes(x = log2FoldChange, y = -log10(pvalue),color=change)) +
  geom_point(alpha=0.8,size =1) +
#  geom_text(aes(label=label_1),vjust=0.4,hjust=-0.2,alpha=.9,size=3) +
#  geom_text(aes(label=label_2),vjust=0.4,hjust=1.3,alpha=.9,size=3) +
  ylim(0,305) +
  xlim(-15,15) +
  theme_bw(base_size = 16,base_family = 'Times') +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank() 
  ) +
  ggtitle('Different expression genes of colorcal cancer in TCGA') + #ͼ????
  theme(title = element_text(size = rel(1.2),lineheight = .9,family = 'Times',face = 'plain')) + #???⼰??????????
  theme(axis.text = element_text(size = rel(1.0),lineheight = .9,family = 'Times',face = 'plain')) + #?̶?????ǩ????
  #scale_fill_discrete(limits=c(labels=c('Tumor:4324','Normal:1607','Nd:50178')) + #ͼ????ǩlimits?̶?λ?ã?labels?????ı???ǩ
  theme(legend.text = element_text(family = 'Times',face = 'plain',size = rel(1.0))) + #ͼ??????
  theme(legend.position = c(.9,.9),legend.justification = c(1,1),legend.background = element_blank()) +#ͼ??λ??
  scale_colour_manual(name='',values = c('red','blue','grey'),limits=c('Tumor','Normal','Nd')) + #??ɫ????
  geom_vline(xintercept = c(-2,2),lty=2,col='grey',lwd=0.5) +
  geom_hline(yintercept = -log10(0.01),lty=2,col='grey',lwd=0.5)
##缺失值是NA即没有counts富集的数据

tcga <- read.csv("~/project/colon cancer/TCGA/resdata_Tumor_Normal.txt",sep = "\t")
cancer <- apply(tcga[,c(52:622)], 1, median)
normal <- apply(tcga[,c(8:51)], 1, median)

tcga <- read.csv("~/project/colon cancer/TCGA/htseq_counts.txt",sep = "\t")
tcga_scale <- t(apply(tcga[,c(2:616)],1, scale))
rownames(tcga_scale) <- tcga$gene_name
cancer <- apply(tcga_scale[,c(46:615)], 1, mean)
normal <- apply(tcga_scale[,c(2:45)], 1, mean)

data <- as.data.frame(cbind(tcga$gene_name,normal,cancer))
names(data) <- c("gene_name","normal","tumor")

library(biomaRt)
mouse <- useDataset("mmusculus_gene_ensembl", mart = useMart("ensembl"))
human <- useMart("ensembl",dataset = "hsapiens_gene_ensembl")
bm <- getLDS(attributes = c("hgnc_symbol"),mart = human,filters = "hgnc_symbol",
             attributesL = c("mgi_symbol"),martL = mouse,values = tcga$gene_name)
library(dplyr)
tcga_clu2 <- merge(x = cluster1,y = bm, by.x= "cluster1",by.y="MGI.symbol")
tcga_clu2 <- merge(x = tcga_clu2, y = data,by.x="HGNC.symbol",by.y="gene_name",all.x=TRUE)
tcga_clu2 <- tcga_clu2[,-1]
library(reshape2)
df <- melt(data = tcga_clu2,id.vars = "cluster1",variable.name = "type",value.name = "value")
df <- na.omit(df)

library(ggplot2)
ggplot(data = df, aes(x = factor(type),y = as.numeric(value),fill=factor(type))) +
  stat_boxplot(geom = "errorbar",linetype=1,width=0.7,position = "identity") +
  geom_boxplot(outlier.fill = NA,outlier.shape=NA,width=0.7) +
  scale_y_continuous(limits = c(0,2000))













# Visualization ------------------------------------------------------------------
#样本相关性 correlation heatmap of samples
vsd <- getVarianceStabilizedData(dds)
heatmap(cor(vsd),cexCol = 0.75,cexRow = 0.75)
#PCAplot
pr <- prcomp(t(vsd))
plot(pr$x,col='red',main="PC plot")
#text(pr$x[,1],pr$x[,2],labels = colnames(vsd),cex = 0.7)
biplot(pr,main="Biplot",col=c("blue","red"))

# DEGs --------------------------------------------------------------------
res <- read.csv("res_Tumor_Normal.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
resdata <- read.csv("resdata_Tumor_Normal.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
colnames(res) <- c('gene_name','baseMean','log2FoldChange','IfcSE','stat','pvalue','padj')
plot(x)


###计算normal 和cancer的均值
normal <- rowMeans(resdata[8:51])
cancer <- rowMeans(resdata[52:622])
res <- cbind(res,normal)
res <- cbind(res,cancer)
attach(res)
high <- res[which(padj < 0.01 & log2FoldChange > 2),][,c(1,3,6,7)]
low <- resdata[which(pvalue < 0.01 & log2FoldChange < -2),][,c(1,3,6,7)]
write.table(x = res,file = 'res_cancer_normal.txt',sep = '\t',quote = FALSE,row.names = FALSE,col.names = TRUE)
write.table(x = high,file = 'new_high_genelist.txt',sep = '\t',quote = FALSE,row.names = FALSE,col.names = TRUE)
write.table(x = low,file = 'new_low_genelist.txt',sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)

g <- function(){
  t <- function(x) return(x^2)
  return(t)
  }



library(car)
x = log2(normal)
y = log2(cancer)
c <- cbind(x,y)
c <- as.data.frame(c)
n <- grep("-Inf",c$x)
m <- grep("-Inf",c$y)
newdata <- c[-c(n,m),]
smoothScatter(newdata,colramp = colorRampPalette(c("white", "yellow")),xlim = c(-5,20),ylim = c(-5,20),
              xlab = "log2FC of Normal tissue",ylab = "log2FC of Tumor tissue")
abline(lm(newdata$y ~ newdata$x))

fit <- lm(newdata$y ~ newdata$x,newdata)
summary(fit)
cor(newdata$y,newdata$x,method = "spearman")
scatterplot(newdata$y ~ newdata$x,newdata,col=carPalette()[-1])


 # GO analysis -------------------------------------------------------------
library(clusterProfiler)
library(biomaRt)
bp <- read.table("~/colon cancer/TCGA/TCGA_BP_low.txt",header = T,sep = "\t")
pathway <- read.csv("~/colon cancer/TCGA/TCGA_pathway_KEGG.txt",header = T,sep = "\t")
min(bp$Fold.Enrichment)
max(bp$Fold.Enrichment)
min(pathway$Fold.Enrichment)
max(pathway$Fold.Enrichment)
bp <- subset(bp,bp$PValue < 0.01 & bp$Fold.Enrichment > 2)
pathway <- subset(pathway,pathway$PValue < 0.01 & pathway$Fold.Enrichment > 2)

library(stringr)
bp <- head(bp[which(bp$PValue < 0.01),],10)
bp$Term <- str_split_fixed(bp$Term,'~',2)[,2]
kegg <- head(pathway[which(pathway$PValue < 0.01),],10)
kegg$Term <- str_split_fixed(kegg$Term,':',2)[,2]
data <- rbind(bp,kegg)
data$PValue <- -log10(data$PValue)
levels(data$Category)[levels(data$Category)=='GOTERM_BP_DIRECT'] <- "Biological Process"
levels(data$Category)[levels(data$Category)=='KEGG_PATHWAY'] <- 'KEGG Pathway'
nameorder <- data$Term[order(data$Category,data$PValue,decreasing = TRUE)]
data$Term <- factor(data$Term,levels = nameorder)
num <- c(1:20)
data <- cbind(data,num)

attach(data)
library(ggplot2)
ggplot(data = data,aes(y = order(num,PValue),x = Term,fill=Category)) +
  geom_bar(position = 'dodge',stat = 'identity',width = 0.8) +
  #scale_color_brewer(palette = "Set1",limits=c("Biological Process","KEGG Pathway")) +
  labs(y='-log10(Pvalue)',x=NULL) +
  ylim(0,20) +
#  scale_x_discrete(limits=rev(levels(Term))) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "grey80"),
    panel.grid.major.y = element_blank(), 
    panel.grid.major.x = element_line(colour='grey90',linetype = 'dashed'),
    panel.grid.minor.x = element_blank()#element_line(colour='grey80',linetype = 'dashed')
  ) +
  ggtitle('GO analysis of TCGA database corectal caner DEGs') +
  theme(plot.title = element_text(size = rel(1.5),lineheight = .9,family = 'Times',face = 'bold',hjust = 0.5,vjust = 0),
        axis.title = element_text(size = rel(1.2),lineheight = .9,family = 'Times',face = 'bold'), 
        axis.text = element_text(size = rel(1.2),lineheight = .9,family = 'Times',face = 'bold',colour = 'black'),
        legend.position = 'top', 
        legend.title = element_text(size = rel(1.1),family = "Times",face = "bold"),
        legend.text = element_text(size = rel(1.1),family = "Times",face = "bold")) +
  scale_fill_discrete(limits=c("Biological Process","KEGG Pathway")) 



library(psych)
library(GPArotation)
options(digits = 2)
covariances <- ability.cov$cov
correlations <- cov2cor(covariances)
fa.parallel(correlations,n.obs = 112,fa="both",n.iter = 100)

fa <- fa(correlations,nfactors = 2,rotate = "none",fm = "pa")
fa.varimax <- fa(correlations,nfactors=2,rotate="varimax",fm="pa")
fa.promax <- fa(correlations,nfactors=2,rotate="promax",fm="pa")
fsm <- function(oblique) {
  if (class(oblique)[2]== "fa" & is.null(oblique$Phi)){
    warning("Object dosen't look like oblique EFA")
  } else {
    P <- unclass(oblique$loading)
    F <- P %*% oblique$Phi
    colnames(F) <- c("PA1","PA2")
    return(F)
  }
}
factor.plot(fa.promax,labels = rownames(fa.promax$loadings))
fa.diagram(fa.promax,simple = FALSE)
fa.tests <- fa(Harman74.cor$cov,nfactors=4,rotate="promax")



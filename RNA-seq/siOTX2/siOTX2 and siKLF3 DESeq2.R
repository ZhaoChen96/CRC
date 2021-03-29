rm(list = ls())
setwd("~/project/colon cancer/siOTX2/")

library(DESeq2)
library(ggplot2)


#inputf <- read.delim("tophat_featureCounts.txt",comment.char = "#")
inputf <- read.delim("star_featureCounts.txt",sep = "\t",comment.char = "#")
names(inputf) <- c("Geneid","Chr","Start","End","Starnd","Length","siKLF3-1","siKLF3-2","siNC-1","siNC-2","siOTX2-1","siOTX2-2")
rownames(inputf) <- inputf$Geneid

# OTX2
inputf <- inputf[,c(9:12)]
data <- inputf[rowSums(inputf > 3) > 2,]
cData <- data.frame(row.names = colnames(data),
                    condition = rep(factor(c("control","OTX2")),each=2))
# KLF3
inputf <- inputf[,c(7:10)]
data <- inputf[rowSums(inputf > 3) > 2,]
cData <- data.frame(row.names = colnames(data),
                    condition = rep(factor(c("KLF3","control")),each=2))

dds <- DESeqDataSetFromMatrix(countData = data,colData = cData,design = ~ condition)
dds <- DESeq(dds)
dds

res <- results(dds,contrast = c("condition","KLF3","control"))
res <- results(dds,contrast = c("condition","OTX2","control"))
res <- res[order(res$pvalue),]
summary(res)

resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized = TRUE)),by='row.names',sort=FALSE)
sig <- res[!is.na(res$pvalue) & res$pvalue < 0.05,]
sig.deseq <- rownames(sig)

table(res$padj < 0.05)
resdata$change <- as.factor(
  ifelse(
    resdata$pvalue < 0.05 & abs(resdata$log2FoldChange) > 0.5,
    ifelse(resdata$log2FoldChange > 0.5,'Up','Down'), # 2^0.5 = 1.414
    'Nd'
  )
)

up <- resdata[resdata$change == "Up",]
down <- resdata[resdata$change == "Down",]
colnames(up)
up <- up[,c("Row.names","log2FoldChange","pvalue","siKLF3-1","siKLF3-2","siNC-1","siNC-2","change")]
down <- down[,c("Row.names","log2FoldChange","pvalue","siKLF3-1","siKLF3-2","siNC-1","siNC-2","change")]
write.table(up,file = "~/project/colon cancer/siOTX2/siKLF3_up_pvalue0.05_gene.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(down,file = "~/project/colon cancer/siOTX2/siKLF3_down_pvalue0.05_gene.txt",sep = "\t",quote = FALSE,row.names = FALSE)

up <- up[,c("Row.names","log2FoldChange","pvalue","siOTX2-1","siOTX2-2","siNC-1","siNC-2","change")]
down <- down[,c("Row.names","log2FoldChange","pvalue","siOTX2-1","siOTX2-2","siNC-1","siNC-2","change")]
write.table(up,file = "~/project/colon cancer/siOTX2/siOTX2_up_pvalue0.05_gene.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(down,file = "~/project/colon cancer/siOTX2/siOTX2-1_down_pvalue0.05_gene.txt",sep = "\t",quote = FALSE,row.names = FALSE)


klfdown <- read.table("~/project/colon cancer/siOTX2/siKLF3_down_pvalue0.05_gene.txt",header = TRUE)
otxdown <- read.table("~/project/colon cancer/siOTX2/siOTX2_down_pvalue0.05_gene.txt",header = TRUE)

a <- intersect(klfdown$Row.names,otxdown$Row.names)
b <- sort(a)
write.table(b,file = "~/project/colon cancer/siOTX2/siKLF3_siOTX2_down_overlap.txt",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

library(clusterProfiler)
bp <- enrichGO(gene = down$Row.names,
               keyType = "SYMBOL",
               OrgDb = "org.Hs.eg.db",
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05
)
bp_0.05 <- subset(bp@result,bp@result$pvalue < 0.05,)[c(1,3,4,6,7,8,9,10,11,12),]

plotbp <- function(data,tf) {
  p <- ggplot(data,aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill=-log10(pvalue))) +
    geom_bar(stat = "identity",width = 0.8,position = position_dodge(width = 0.8)) +
    theme_classic(base_size = 18,base_family = "sans") +
    scale_x_continuous(expand = c(0,0),limits = c(0,max(ceiling(-log10(bp_0.05$pvalue))))) +
    scale_fill_gradient(low = "black",high = "black") +
    labs(title = paste("BP of",tf,"down regulation gene",sep = " "),x="-log10(pValue)",y=NULL) +
    guides(fill=FALSE) +
    theme(plot.title = element_text(size = 20,hjust = 1.2,vjust = 0.5),
          plot.margin = margin(t=0.3,r=0.3,b=0.3,l=0.3,unit = "cm"),
          axis.text = element_text(colour = 'black',size = 18))
  p
}
plotbp(data = bp_0.05,tf = "siKLF3")
ggsave(filename = paste("~/project/colon cancer/siOTX2/siKLF3_down_bp.pdf",sep = ""),height = 3.9,width = 7)

library(biomaRt)
human <- useMart(dataset = "hsapiens_gene_ensembl",biomart = "ensembl",host="http://uswest.ensembl.org")
entrezgene <- getBM(mart = human,attributes = c('external_gene_name','entrezgene_id'),filters = "external_gene_name",values = klfdown$Row.names)    
kegg <- enrichKEGG(gene = entrezgene[,2],
                   keyType = "kegg",
                   organism = "hsa",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05
)
library(scales)
kegg_0.05 <- kegg@result[kegg@result$pvalue < 0.05,]#[1:10,]
plotkegg <- function(data,tf){
  p <- ggplot(data,aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill=-log10(pvalue))) +
    geom_bar(stat = "identity",width = 0.8,position = position_dodge(width = 0.8)) +
    theme_classic(base_size = 14,base_family = "sans") +
    scale_x_continuous(expand = c(0,0),limits = c(0,3)) +#limits = c(0,max(ceiling(-log10(kegg@result$pvalue))))) +
    scale_fill_gradient(low = "#08519C",high = "#08519C") +
    labs(title = paste("KEGG pathway of",tf,"down regulation gene",sep = " "),x="-log10(pValue)",y=NULL) +
    guides(fill=FALSE) +
    theme(plot.title = element_text(size = 14,hjust = 1.1,vjust = 0.5),
          plot.margin = margin(t=0,r=0.3,b=0,l=0,unit = "cm"),
          axis.text = element_text(colour = 'black',size = 14))
  p
}
plotkegg(data = kegg_0.05,tf = "siKLF3")
ggsave(filename = paste("~/project/colon cancer/siOTX2/siKLF3_down_kegg.pdf",sep = ""),height = 3.9,width = 4.8)



library(clusterProfiler)
bp <- enrichGO(gene = otxdown$Row.names,
               keyType = "SYMBOL",
               OrgDb = "org.Hs.eg.db",
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05
)
bp_0.05 <- subset(bp@result,bp@result$pvalue < 0.05,)[c(1,4,5,6,7,8,9,10,11,12),]

plotbp <- function(data,tf) {
  p <- ggplot(data,aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill=-log10(pvalue))) +
    geom_bar(stat = "identity",width = 0.8,position = position_dodge(width = 0.8)) +
    theme_classic(base_size = 18,base_family = "sans") +
    scale_x_continuous(expand = c(0,0),limits = c(0,max(ceiling(-log10(bp_0.05$pvalue))))) +
    scale_fill_gradient(low = "black",high = "black") +
    labs(title = paste("BP of",tf,"down regulation gene",sep = " "),x="-log10(pValue)",y=NULL) +
    guides(fill=FALSE) +
    theme(plot.title = element_text(size = 20,hjust = 1.2,vjust = 0.5),
          plot.margin = margin(t=0.3,r=0.3,b=0.3,l=0.3,unit = "cm"),
          axis.text = element_text(colour = 'black',size = 18))
  p
}
plotbp(data = bp_0.05,tf = "siOTX2")
ggsave(filename = paste("~/project/colon cancer/siOTX2/siOTX2_down_bp.pdf",sep = ""),height = 5,width = 9)

library(biomaRt)
human <- useMart(dataset = "hsapiens_gene_ensembl",biomart = "ensembl",host="http://uswest.ensembl.org")
entrezgene <- getBM(mart = human,attributes = c('external_gene_name','entrezgene_id'),filters = "external_gene_name",values = otxdown$Row.names)    
kegg <- enrichKEGG(gene = entrezgene[,2],
                   keyType = "kegg",
                   organism = "hsa",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05
)
library(scales)
kegg_0.05 <- kegg@result[kegg@result$pvalue < 0.05,]#[1:10,]
plotkegg <- function(data,tf){
  p <- ggplot(data,aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill=-log10(pvalue))) +
    geom_bar(stat = "identity",width = 0.8,position = position_dodge(width = 0.8)) +
    theme_classic(base_size = 14,base_family = "sans") +
    scale_x_continuous(expand = c(0,0),limits = c(0,3)) +#limits = c(0,max(ceiling(-log10(kegg@result$pvalue))))) +
    scale_fill_gradient(low = "#08519C",high = "#08519C") +
    labs(title = paste("KEGG pathway of",tf,"down regulation gene",sep = " "),x="-log10(pValue)",y=NULL) +
    guides(fill=FALSE) +
    theme(plot.title = element_text(size = 14,hjust = 1.1,vjust = 0.5),
          plot.margin = margin(t=0.3,r=0.3,b=0.3,l=0.3,unit = "cm"),
          axis.text = element_text(colour = 'black',size = 14))
  p
}
plotkegg(data = kegg_0.05,tf = "siOTX2")
ggsave(filename = paste("~/project/colon cancer/siOTX2/siOTX2_down_kegg.pdf",sep = ""),height = 3,width = 4.8)


vsd <- getVarianceStabilizedData(dds)
heatmap(cor(vsd),cexCol = 1,cexRow = 1)
pr <- prcomp(t(vsd))
plot(pr$x,col="white",main="PC plot")
text(pr$x[,1],pr$x[,2],labels=colnames(vsd),
     cex=0.7)

biplot(pr,cex=c(1,0.5),main="KLF3",
       col=c("black","grey"))
pvalue 



# HCT116 FPKM -------------------------------------------------------------
file <- read.delim("genes.fpkm_tracking",sep = "\t")
data <- file[,c(5,10,14,18,22,26,30)]
write.table(data,file = "siKLF3 and siOTX2 FPKM.txt",sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)

df <- data[data$gene_short_name == "KLF3",]
df <- data[data$gene_short_name == "OTX2",]

library(reshape2)
df <- melt(df[,1:5],id.vars="gene_short_name",variable.name = "sample",value.name = "fpkm")
df <- cbind(df,type=rep(c("siKLF3","siNC"),each=2))

ggplot(df,aes(x = factor(type),y = fpkm,color=type)) +
  geom_point(shape=19,alpha=0.8) +
  #geom_point(shape=21,colour="white",fill="#1B9E77",alpha=0.6) + #,position = position_jitter(width = 0.25,height = 0)
  theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) +
  labs(title = "KLF3",x=NULL,y="Gene expression (FPKM)",color=NULL) +
  scale_y_continuous(expand = c(0,0.01),limits = c(0,9),breaks = c(0,3,6,9)) +
  #scale_x_discrete(breaks=times,labels=c(0,2,4,7,10)) +
  theme(aspect.ratio = 1/0.618,
        plot.margin = margin(t=2,b=2,l=0.5,r=0.5,unit = "cm"),
        plot.title = element_text(size = 18,face = "bold",margin = margin(b=0,unit = "cm"),hjust = 0.5,vjust = 0.5),
        axis.title.y = element_text(margin = margin(t = 0,b = 0,l = 0,r = 0.2,unit = "cm")),
        axis.text = element_text(colour = "black",size = 18,margin = margin(t = 0.1,b = 0.1,l = 0.1,r = 0.1,unit = "cm")),
        axis.text.x = element_text(angle = 0,hjust = 0.5)) +
  guides(color=FALSE)


# overlap with chip-seq data ----------------------------------------------
klfdown <- read.table("~/project/colon cancer/siOTX2/siKLF3_down_pvalue0.05_gene.txt",header = TRUE)
otxdown <- read.table("~/project/colon cancer/siOTX2/siOTX2_down_pvalue0.05_gene.txt",header = TRUE)
klfup <- read.delim("~/project/colon cancer/siOTX2/siKLF3_up_pvalue0.05_gene.txt",header = TRUE)
otxup <- read.delim("~/project/colon cancer/siOTX2/siOTX2_up_pvalue0.05_gene.txt",header = TRUE)

###### Cistrome Data: human KLF3
file <- read.delim("~/project/colon cancer/siOTX2/klf3 target gene /38790_gene_score_5fold.txt",sep = "\t",comment.char = "#",header = FALSE)
colnames(file) <- c("chr","start","end","refseq","score","strand","symbol")

a <- unique(file$symbol)[1:500]
a <- unique(file$symbol)[500:1000]
a <- unique(file$symbol)[1000:1500]
a <- unique(file$symbol)[1500:2000]
b <- intersect(klfdown$Row.names,a)
b
c <- intersect(klfup$Row.names,a)
c

###### Cistrome Data: human OTX2
# 6540 cell human OTX2 epithelium
inputf <- read.delim("~/project/colon cancer/siOTX2/50896_gene_score_5fold.txt",comment.char = "#",header = FALSE)
colnames(inputf) <- c("chr","start","end","refseq","score","strand","symbol")

# D341 cell 人髓母细胞瘤
file1 <- read.delim("~/project/colon cancer/siOTX2/OTX2 human D341 cell/74350_gene_score_5fold.txt",comment.char = "#",header = FALSE)
colnames(file1) <- c("chr","start","end","refseq","score","strand","symbol")
file2 <- read.delim("~/project/colon cancer/siOTX2/OTX2 human D341 cell/74352_gene_score_5fold.txt",comment.char = "#",header = FALSE)
colnames(file2) <- c("chr","start","end","refseq","score","strand","symbol")
file3 <- read.delim("~/project/colon cancer/siOTX2/OTX2 human D341 cell/74358_gene_score_5fold.txt",comment.char = "#",header = FALSE)
colnames(file3) <- c("chr","start","end","refseq","score","strand","symbol")

list <- c(unique(file1$symbol)[1500:2000],unique(file2$symbol)[1500:2000],unique(file3$symbol)[1500:2000])
list <- c(unique(file1$symbol)[1:500],unique(file2$symbol)[1:500],unique(file3$symbol)[1:500])
list <- unique(list)

a <- unique(list)[1:500]
a <- unique(list)[500:1000]
a <- unique(list)[1000:1500]
a <- unique(list)[1500:2000]

b <- intersect(otxdown$Row.names,list)
b
c <- intersect(otxup$Row.names,list)
c













killDbConnections = function () {
  all_cons = dbListConnections(MySQL())
  for(con in all_cons){dbDisconnect(con)}
}
args=commandArgs(T)
signature = args[1]       ## signature genelist
dataset = args[2]         ## selected datasets
parameter = args[3]       ## parameters
signatures = strsplit(signature,",")[[1]]       ## split the signature genelist into array
if(is.na(signatures[4])){signatures[4] = ""}
datasets = strsplit(dataset,",")[[1]]           ## split the datasets info into array
datasets = datasets[datasets != ""]
parameters = strsplit(parameter,",")[[1]]
dbs = "GE_SF"             ## database
symbol1 = parameters[1]
symbol2 = parameters[2]
symbol3 = parameters[3]
symbol4 = parameters[4]
method = parameters[5]          ## "pearson" or "spearman" or "kendall"
outputdir = parameters[6]       ## outputdir


.libPaths(c("--------",.libPaths()))      ## add the specific directory of RMySQL R package
suppressPackageStartupMessages(library("RMySQL"))
killDbConnections()
#mydb = dbConnect(MySQL(), user='--------', password='--------', dbname=dbs)     ## connet to mysql database
mydb = dbConnect(MySQL(), user='--------', password='--------', dbname=dbs,unix.socket="--------") 

## Gene A:
if(signatures[3] != ""){
  symbol1 = paste(symbol1,symbol3,sep = "/")
  signatures_a = c(signatures[1],signatures[3])
}else{signatures_a = signatures[1]}
df_a = as.matrix(array(data = 0,dim = c(0,length(signatures_a))))
for(i in 1:length(datasets)){
  table = datasets[i]
  df_t=t(dbGetQuery(mydb,paste("SELECT * FROM ",table," WHERE geneid IN ('",paste(signatures_a, collapse = "','"),"')",sep="")))
  colnames(df_t) = df_t[1,]
  df_t = df_t[-1,,drop = F]
  df_a = rbind(df_a,df_t)
}
storage.mode(df_a) = "numeric"
df_a = df_a[,signatures_a,drop = F]
if(signatures[3] != ""){
  df_a = log2(df_a + 0.001)
  df_a = df_a[,1,drop = F] - df_a[,2,drop = F]
  df_a = 2^df_a
  colnames(df_a) = signatures_a[1]
}
## Gene B:
if(signatures[4] != ""){
  symbol2 = paste(symbol2,symbol4,sep = "/")
  signatures_b = c(signatures[2],signatures[4])
}else{signatures_b = signatures[2]}
df_b = as.matrix(array(data = 0,dim = c(0,length(signatures_b))))
for(i in 1:length(datasets)){
  table = datasets[i]
  df_t=t(dbGetQuery(mydb,paste("SELECT * FROM ",table," WHERE geneid IN ('",paste(signatures_b, collapse = "','"),"')",sep="")))
  colnames(df_t) = df_t[1,]
  df_t = df_t[-1,,drop = F]
  df_b = rbind(df_b,df_t)
}
storage.mode(df_b) = "numeric"
df_b = df_b[,signatures_b,drop = F]
if(signatures[4] != "" ){
  df_b = log2(df_b + 0.001)
  df_b = df_b[,1,drop = F] - df_b[,2,drop = F]
  df_b = 2^df_b
  colnames(df_b) = signatures_b[1]
}

df = cbind(df_a,df_b)



### build the integrated table
### single gene
colnames(df)[colnames(df) == signatures[1]] = symbol1
colnames(df)[colnames(df) == signatures[2]] = symbol2

pdf(file = outputdir,title="Result Display",width = 6,height = 5.5)
par(mar=c(4.5, 5.1, 1.1, 2.1))
options(warn=-1)
rvalue = signif(cor(x = df[,1], y = df[,2],method = method),2)
cpvalue = signif(as.numeric(cor.test(x = df[,1], y = df[,2],method = method)[3]),2)
plot(x = log2(df[,1] + 1),y = log2(df[,2] + 1),main = NULL,cex.lab = 1.5,
     xlab = paste("log2(",colnames(df)[1]," TPM)",sep = ""),
     ylab = paste("log2(",colnames(df)[2]," TPM)",sep = ""), pch = 19,cex=0.5)
range_y = range(log2(df[,2] + 1))
text(x = max(log2(df[,1] + 1)) * 1.03,y = max(log2(df[,2] + 1)) - (range_y[2] - range_y[1])/15,labels = paste("p-value = ",as.character(cpvalue),"\nR = ",as.character(rvalue),sep=""),cex = 1.3,col = "black",pos = 2)
a = dev.off()









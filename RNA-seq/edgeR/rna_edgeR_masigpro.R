rm(list = ls())


# edgeR : colon cancer RNA-seq htcount-------------------------------------------------------------------
library(limma)
library(edgeR)

setwd("~/project/colon cancer/RNA-seq/edgeR/masigpro_chen/")

inputf <- read.table(file = "~/project/colon cancer/RNA-seq/DESeq2_maSigPro/merge_htseqCount.txt",header = TRUE,row.names = 1,check.names = FALSE)
sample_total <- c()
for (time in c("control","2weeks","4weeks","7weeks","10weeks")){
  for (rep in c(1:3)){
    sample <- paste(time,rep,sep = "_rep")
    sample_total <- c(sample_total,sample)
       }
}
names(inputf) <- sample_total
# 2.make matrix
countMatrix <- as.matrix(inputf)
# 3.make DEGlist
times <- c("control","2weeks","4weeks","7weeks","10weeks")
group <- factor(rep(times,each=3),levels = times)
y <- DGEList(counts = countMatrix,group = group)
# 4.filter the data: the data=0
keep <- filterByExpr(y)
keep <- rowSums(cpm(y)>1) >=2
y <- y[keep, , keep.lib.sizes=FALSE]
# 5.TMM normalization
y <- calcNormFactors(y,method = "TMM")
#write.table(y$counts,file = "~/project/colon cancer/RNA-seq/edgeR/masigpro_chen/counts_filter.txt",sep = "\t")
write.table(y$counts,file = "~/project/colon cancer/RNA-seq/edgeR/nonfilter/counts_notfilter.txt",sep = "\t")
# 6.design a matrix
design <- model.matrix(~group)
rownames(design) <- colnames(y)
# estimate disper
y <- estimateDisp(y,design,robust = TRUE)
y$common.dispersion
#dispersion=0.124,bcv = 0.35
plotBCV(y)
plotMD(cpm(y,log = TRUE),column = 1)
abline(h=0,col="red",lty=2,lwd=2)
#differential expression
# To perform likelihood ratio tests:
fit <- glmFit(y,design)
lrt <- glmLRT(fit)
topTags(lrt)
colnames(design)
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]
summary(decideTests(lrt))
plotMD(lrt)
abline(h=c(-1,1), col="blue")

# To perform quasi-likelihood F-tests:
fit <- glmQLFit(y,design,robust = TRUE)
qlf <- glmQLFTest(fit,coef = 2)
topTags(qlf)
plotQLDisp(fit)

qlf <- glmQLFTest(fit)
topTags(qlf)
top <- rownames(topTags(qlf))
cpm(y)[top,]

logcpm <- cpm(y,log = TRUE)
plotMDS(logcpm)

#pick the DEGs name
summary(dt <- decideTestsDGE(qlf))
isDE <- as.logical(dt)
DEnames <- rownames(y)[isDE]
head(DEnames)
plotSmear(qlf,de.tags = DEnames)
abline(h=c(-1,1),col="blue")


# RNA-seq maSigPro --------------------------------------------------------
rm(list = ls())
library(maSigPro)
library(MASS)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)


Times <- rep(c(0,2,4,7,10),each=3)
Replicates <- rep(1:5,each=3) 
RNA <- rep(1,15)
ss.edesign <- cbind(Times,Replicates,RNA)
sample_vector <- c()
for (time in c("control","2weeks","4weeks","7weeks","10weeks")){
  for (rep in c(1:3)){
    sample <- paste(time,rep,sep = "_rep")
    sample_vector <- c(sample_vector,sample)
  }
}
rownames(ss.edesign) <- sample_vector
design <- make.design.matrix(ss.edesign,degree = 4)

file <- read.table(file = "~/project/colon cancer/RNA-seq/edgeR/filter/counts_filter.txt",sep = "\t",check.names = FALSE)
rownames(file) <- substr(rownames(file),1,18)
fit <- p.vector(file,design,counts = TRUE,theta = 10)
# p.vector() returns a list of values:
fit$G #total number of input gene 16099
fit$i #returns the number of significant genes: 6113
pie(c(14330,40984),labels = c("significant genes:14330 ","non significant genes: 40984"),col = c("pink","white")) 
fit$SELEC # a matrix with the significant genes and their expression values
# selects the best regression model for each gene using stepwise regression
tstep <- T.fit(fit,step.method = "backward",alfa = 0.05)
save(tstep,file = "~/project/colon cancer/RNA-seq/edgeR/filter/RNA-seq_masigpro_20200709.Rdata")
# creates lists of significant genes for a set of variables 
sigs <- get.siggenes(tstep,rsq = 0.7,vars = "all")

pdf("~/project/colon cancer/RNA-seq/edgeR/RNA-seq_masigpro_nonfilter_k8.pdf")
# visualisation tools for gene expression values in a time course experiment
result <- see.genes(sigs$sig.genes,k=4,newX11 = FALSE)
dev.off()


# This plot can be made for specic genes or for groups of genes where the median will be computed.
# a specic gene 
Gpr18 <- file[rownames(file)=="Gpr18",] 
PlotGroups(Gpr18,edesign = ss.edesign,show.fit = T,dis = design$dis,groups.vector = design$groups.vector)

# a gene expression profiles:check the homogeneity of the clusters
PlotProfiles(tstep$sig.profiles,cond = colnames(tstep$sig.profiles),repvect = rep(c(1:5),each=3))


# genes in each clusters --------------------------------------------------
library(biomaRt)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "mmusculus_gene_ensembl")

mytheme <- theme_bw() + 
  theme(panel.grid = element_blank(), 
        panel.grid.major.x = element_line(linetype = "dashed",colour = "grey80"),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size = rel(0.9))) +
  theme(
    plot.title = element_text(size = 14,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 1,vjust = 0.5),
    axis.title = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain'),  
    axis.text = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain',colour = 'black')
  )
  

for (i in 1:4){
  cluster <- paste("cluster",i,sep = "")
  clu <- result$cut[result$cut == i]
  gene_table <- data.frame(names(clu))
  names(gene_table) <- cluster
  table_name <- paste("~/project/colon cancer/RNA-seq/edgeR/filter/RNA-seq_masigpro_genelist_",cluster,".txt",sep = "")
  write.table(gene_table,file = table_name,sep = "\t",quote = F)
  
  #GO analysis: biology process
  ego <- enrichGO(gene = names(clu),
                  keyType = "ENSEMBL",
                  OrgDb = "org.Mm.eg.db",
                  ont = "BP",
                 pAdjustMethod = "BH",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01)
  GO_list <- paste("~/project/colon cancer/RNA-seq/edgeR/filter/RNA-seq_masigpro_BP_",cluster,".txt",sep = "")
  write.table(ego@result, GO_list,sep = "\t",quote = F,row.names = F,col.names = T,na = "NA",eol = "\n")
  BP_pdf <- paste("~/project/colon cancer/RNA-seq/edgeR/filter/RNA-seq_masigpro_BP_",cluster,".pdf",sep = "")
  
  ###barplot
  #pdf(BP_pdf)
  #p <- barplot(ego,showCategory=15)
  #print(p)
  #dev.off()
  
  ###ggplot
  plottitle <- paste("RNA-seq maSigPro biological process:",cluster,sep = " ")
  ggplot(ego, aes(-log10(pvalue),Description),showCategory=10) + 
    geom_bar(stat = "identity",width = 0.7,position = position_dodge(0.7),fill="#bdd7e7") + 
    ggtitle(plottitle) +
    labs(y="") +
    xlim(0,ceiling(max(-log10(ego@result$pvalue)))) +
    mytheme 
  ggsave(BP_pdf,width = 10,height = 5)
  
  
  #GO analysis: kegg
  entrezgene <- getBM(mart = mart,attributes = c("entrezgene_id","external_gene_name"),filters = "ensembl_gene_id",values = names(clu))
  ego <- enrichKEGG(gene = entrezgene[,1],
                    keyType = "kegg",
                    organism = "mmu",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01)
  kegg_list <- paste("~/project/colon cancer/RNA-seq/edgeR/nonfilter/RNA-seq_masigpro_KEGG_",cluster,".txt",sep = "")
  write.table(ego@result, kegg_list,sep = "\t",quote = F,row.names = F,col.names = T,na = "NA",eol = "\n")
  kegg_pdf <- paste("~/project/colon cancer/RNA-seq/edgeR/nonfilter/RNA-seq_masigpro_KEGG_",cluster,".pdf",sep = "")
  #pdf(kegg_pdf,width = 10,height = 10)
  #p <- barplot(ego,showCategory = 15)
  #print(p)
  #dev.off()
  
  # ggplot
  KEGG_pdf <- paste("~/project/colon cancer/RNA-seq/edgeR/filter/RNA-seq_masigpro_KEGG_",cluster,".pdf",sep = "")
  plottitle <- paste("RNA-seq maSigPro KEGG pathway:",cluster,sep = " ")
  ggplot(ego, aes(-log10(pvalue),Description),showCategory=10) + 
    geom_bar(stat = "identity",width = 0.7,position = position_dodge(0.7),fill="skyblue") + 
    ggtitle(plottitle) +
    labs(y="") +
    xlim(0,ceiling(max(-log10(ego@result$pvalue)))) +
    mytheme 
  
  #ggarrange(p1,p2,ncol = 1,nrow = 2,labels = c("A","B"))
  ggsave(KEGG_pdf,width = 10,height = 5)
}




mart_mouse_symbol <- read.delim("~/project/colon cancer/RNA-seq/mart_export_mouse_geneid_symbol.txt")
a <- intersect(entrezgene$ensembl_gene_id,mart_mouse_symbol$ensembl_gene_id)

ensembl_to_entrez <- function(entrezgene) {
  entrez_map <- as.data.frame(unlist(as.list(org.Mm.egENSEMBL2EG)))
  index <- match(entrezgene$ensembl_gene_id,rownames(entrez_map))
  entrez_id <- mart_mouse_symbol[match(entrezgene$ensembl_gene_id,mart_mouse_symbol$ensembl_gene_id),"ensembl_gene_id"]
}


a <- c("12006","21414","170758","12125","16842","21802","12443")
b <- getBM(mart = mart,attributes = c("entrezgene_id","external_gene_name"),filters = "entrezgene_id",values = a)
c <- c("21926","16176","26414","12702")
d <- getBM(mart = mart,attributes = c("entrezgene_id","external_gene_name"),filters = "entrezgene_id",values = c)

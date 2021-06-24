rm(list=ls())
library(edgeR)
library(ggplot2)
library("maSigPro")
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(enrichplot)


Markers <- c("H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me3")
marker <- Markers[1]

Times <- c(rep(c(0,2,4,7,10),each=3),rep(c(0,2,4,7,10),each=3))
Replications <- rep(1:10,each=3)
Control <- c(rep(1,15),rep(0,15))
marker <- c(rep(0,15),rep(1,15))
ss.edesign <- cbind(Times,Replications,Control,marker)
print(ss.edesign)



sample_vector <- c()
for (mark in c(marker,"ctrl")) {
  for (time in c("ctrl","2weeks","4weeks","7weeks","10weeks")){
    for (rep in c(1:3)){
      sample <- paste(time,rep,mark,sep = "-")
      sample_vector <- c(sample_vector,sample)
    }
  }
}

print(sample_vector)
rownames(ss.edesign) <- sample_vector
design <- make.design.matrix(ss.edesign,degree = 4)
print(design)

# scale
readCounts = paste("/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/",marker,"/",marker,"_readCounts.txt",sep = "")
file <- read.delim(file = readCounts,sep = "\t",check.names = FALSE,header = T,row.names = 1)

df = data.matrix(file)
df = scale(df, center = TRUE, scale = TRUE)

fit <- p.vector(df,design,Q = 0.05, MT.adjust = "BH", min.obs = 5)
fit$G
fit$i
tstep <- T.fit(fit,step.method = "backward",alfa = 0.05)
tstep_name <- paste("/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/",marker,"/",marker,"_20200724.Rdata",sep = "")
save(tstep,file = tstep_name)

sigs <- get.siggenes(tstep,rsq = 0.7,vars = "all")

pdf_name <- paste("/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/",marker,"/",marker,"_maSigPro_k9.pdf",sep = "")
pdf(pdf_name,height = 10,width = 10)
result <- see.genes(sigs$sig.genes,k=9,newX11 = FALSE)
dev.off()

# GO
library(biomaRt)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "mmusculus_gene_ensembl")

mytheme <- theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed",colour = "grey80"),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size = rel(0.9))) +
  theme(
    plot.title = element_text(size = 14,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 1.1,vjust = 2.6),
    axis.title = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain'),
    axis.text = element_text(size = 12,lineheight = .9,family = 'Helvetica',face = 'plain',colour = 'black')
  )

outputDir = paste("/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/",marker,"/",marker,"_",sep = "")

for (i in 1:9){
  cluster <- paste("cluster",i,sep = "")
  clu <- result$cut[result$cut == i]
  print(length(names(clu)))
  entrezgene <- getBM(mart = mart,attributes = c("entrezgene_id","external_gene_name"),filters = "ensembl_gene_id",values = names(clu))
  gene_table <- data.frame(entrezgene[,2])
  print(length(gene_table))
  names(gene_table) <- cluster
  table_name <- paste(outputDir,cluster,".txt",sep = "")
  write.table(gene_table,file = table_name,sep = "\t",quote = F)


  #GO analysis: biology process
  ego <- enrichGO(gene = entrezgene[,2],
                  keyType = "SYMBOL",
                  OrgDb = "org.Mm.eg.db",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.01)
  GO_list <- paste(outputDir,cluster,"_BP.txt",sep = "")
  write.table(ego@result, GO_list,sep = "\t",quote = F,row.names = F,col.names = T,na = "NA",eol = "\n")
  BP_pdf <- paste(outputDir,cluster,"_BP.pdf",sep = "")

  ###ggplot
  plottitle <- paste(marker,"maSigPro biological process:",cluster,sep = " ")
  ggplot(ego, aes(-log10(pvalue),Description),showCategory=15) +
    geom_bar(stat = "identity",width = 0.7,position = position_dodge(0.7),fill="#bdd7e7") +
    ggtitle(plottitle) +
    labs(y="") +
    xlim(0,ceiling(max(-log10(ego@result$pvalue)))) +
    mytheme
    #theme(axis.text.y = element_text(size = 14, family= 'Helvetica',face = 'plain',colour = "#bdd7e7"))
  ggsave(BP_pdf,width = 10,height = 5)


  #GO analysis: kegg
  ego <- enrichKEGG(gene = entrezgene[,1],
                    keyType = "kegg",
                    organism = "mmu",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01)
  kegg_list <- paste(outputDir,cluster,"_KEGG.txt",sep = "")
  write.table(ego@result, kegg_list,sep = "\t",quote = F,row.names = F,col.names = T,na = "NA",eol = "\n")

  # ggplot
  KEGG_pdf <- paste(outputDir,cluster,"_KEGG.pdf",sep = "")
  plottitle <- paste(marker,"maSigPro KEGG pathway:",cluster,sep = " ")
  ggplot(ego, aes(-log10(pvalue),Description),showCategory=15) +
    geom_bar(stat = "identity",width = 0.7,position = position_dodge(0.7),fill="skyblue") +
    ggtitle(plottitle) +
    labs(y="") +
    xlim(0,ceiling(max(-log10(ego@result$pvalue)))) +
    mytheme

  #ggarrange(p1,p2,ncol = 1,nrow = 2,labels = c("A","B"))
  ggsave(KEGG_pdf,width = 10,height = 5)
}

rm(list = ls())
library(ggplot2)

klftarget <- c("ADD3","CALCOCO1","ACBD4","C4orf3","CCNG2","AHNAK2","EHF","HRK","PDCD4","ELF3","BHLHE40","KDM3A")

otx2target <- c("CYP24A1","DUSP6","F2R","IL17RD","S100A10" ,"CPA4", "LAMA3","KRT7","H2BC12","NCOA7","EHF","UNC13D")
otx2targetd341 <- c("ENO2","PPFIA4","BNIP3L","TLE6","GRAMD1B","NDRG1","BACH1","HCFC2","GDF11","PDCD4","KIAA1958","HAS3","SYT1","TMEM198")

tcga <- read.delim("~/project/colon cancer/TCGA/COAD_FPKM_add_symbol.txt",check.names = FALSE)
#row.names(tcga) <- tcga$Gene.name

##### KLF3
plotPoint <- function(i,limits) {
  gene <- tcga[tcga$Gene.name %in% c("KLF3",klftarget[i]),]
  rownames(gene) <- gene$Gene.name
  gene <- gene[,-1]
  df <- as.data.frame(t(gene))
  cor <- cor(df[,1],df[,2],method="spearman")
  corlist <- c(corlist,cor)
  print(cor)
  p <- ggplot(df,aes(x = log2(KLF3),y = log2(df[,2]))) +
    geom_point() +
    # scale_x_continuous(limits = limits) +
    # scale_y_continuous(limits = limits) +
    labs(title = paste("Cor:",round(cor,2),sep = " "),x="log2(KLF3 FPKM)",y=paste("log2(",colnames(df)[2]," FPKM)",sep = "")) +
    theme_bw(base_size = 18,base_family = "sans",base_line_size = 1.1) +
    theme(aspect.ratio = 1/1,
      panel.grid = element_blank(),
          plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
          axis.title = element_text(size = 18,colour = "black"),
          axis.text = element_text(size = 18,colour = "black"))
  ggsave(paste("TCGA_KLF3_",colnames(df)[2],"cor.png",sep = ""),p,width = 10,height = 10,
         path = "~/project/colon cancer/siOTX2/klf3 target gene /",dpi = 300,units = "cm")
  p
}

plotPoint(i = 7,limits = c(0,8))
plotPoint(i = 9,limits = c(0,8))
plotPoint(i = 10,limits = c(0,8))
plotPoint(i = 12,limits = c(0,8))

# OTX2 --------------------------------------------------------------------
# 6540 cell
plotPoint <- function(i,limits) {
  gene <- tcga[tcga$Gene.name %in% c("OTX2",otx2target[i]),]
  rownames(gene) <- gene$Gene.name
  gene <- gene[,-1]
  df <- as.data.frame(t(gene))
  cor <- cor(df[,1],df[,2],method="spearman")
  p <- ggplot(df,aes(x = log2(OTX2),y = log2(df[,1]))) +
    geom_point() +
    scale_x_continuous(limits = limits) +
    scale_y_continuous(limits = limits) +
    labs(title = paste("Cor:",round(cor,2),sep = " "),x="log2(OTX2 FPKM)",y=paste("log2(",colnames(df)[1]," FPKM)",sep = "")) +
    theme_bw(base_size = 18,base_family = "sans",base_line_size = 1.1) +
    theme(aspect.ratio = 1/1,
          panel.grid = element_blank(),
          plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
          axis.title = element_text(size = 18,colour = "black"),
          axis.text = element_text(size = 18,colour = "black"))
  ggsave(paste("TCGA_OTX2_",colnames(df)[1],"cor.png",sep = ""),p,width = 10,height = 10,
         path = "~/project/colon cancer/siOTX2/OTX2 human cell/6540 cell/",dpi = 300,units = "cm")
  p
}

plotPoint(i = 4,limits = c(-10,10))
plotPoint(i = 8,limits = c(-10,10))
plotPoint(i = 11,limits = c(-10,10))
plotPoint(i = 12,limits = c(-10,10))

# D341 cell
i <- 14
plotPoint <- function(i,limits) {
  gene <- tcga[tcga$Gene.name %in% c("OTX2",otx2targetd341[i]),]
  rownames(gene) <- gene$Gene.name
  gene <- gene[,-1]
  df <- as.data.frame(t(gene))
  cor <- cor(df[,1],df[,2],method="spearman")
  p <- ggplot(df,aes(x = log2(OTX2),y = log2(df[,1]))) +
    geom_point() +
    scale_x_continuous(limits = limits) +
    scale_y_continuous(limits = limits) +
    labs(title = paste("Cor:",round(cor,2),sep = " "),x="log2(OTX2 FPKM)",y=paste("log2(",colnames(df)[1]," FPKM)",sep = "")) +
    theme_bw(base_size = 18,base_family = "sans",base_line_size = 1.1) +
    theme(aspect.ratio = 1/1,
          panel.grid = element_blank(),
          plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
          axis.title = element_text(size = 18,colour = "black"),
          axis.text = element_text(size = 18,colour = "black"))
  ggsave(paste("TCGA_OTX2_",colnames(df)[1],"cor.png",sep = ""),p,width = 10,height = 10,
         path = "~/project/colon cancer/siOTX2/OTX2 human cell/D341 cell/",dpi = 300,units = "cm")
  p
}

plotPoint(i = 2,limits = c(-10,10))
plotPoint(i = 4,limits = c(-10,10))
plotPoint(i = 9,limits = c(-10,10))
plotPoint(i = 11,limits = c(-10,10))
plotPoint(i = 12,limits = c(-10,10))
plotPoint(i = 13,limits = c(-10,10))
plotPoint(i = 14,limits = c(-10,10))
ggplot(df,aes(x = log2(OTX2),y = log2(df[,2]))) +
  geom_point() +
  scale_x_continuous(limits = limits) +
  scale_y_continuous(limits = limits) +
  labs(title = paste("Cor:",round(cor,2),sep = " "),x="log2(OTX2 FPKM)",y=paste("log2(",colnames(df)[2]," FPKM)",sep = "")) +
  theme_bw(base_size = 18,base_family = "sans",base_line_size = 1.1) +
  theme(aspect.ratio = 1/1,
        panel.grid = element_blank(),
        plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
        axis.title = element_text(size = 18,colour = "black"),
        axis.text = element_text(size = 18,colour = "black"))
ggsave(paste("TCGA_OTX2_",colnames(df)[2],"cor.png",sep = ""),width = 10,height = 10,
       path = "~/project/colon cancer/siOTX2/OTX2 human cell/D341 cell/",dpi = 300,units = "cm")

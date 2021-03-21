rm(list = ls())
setwd("~/project/colon cancer/chip-seq/chromHMM/geneCount/")
library(ggplot2)
library(reshape2)
library(RColorBrewer)
geneCountDir <- "~/project/colon cancer/chip-seq/chromHMM/geneCount/"


plotBox <- function(filename) {
  region <- paste(unlist(strsplit(basename(filename),split = "_"))[2],unlist(strsplit(basename(filename),split = "_"))[3],sep = "")
  file <- read.csv(filename,check.names = FALSE)
  state = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13")
  element = c("enhancer", "enhancer", "enhancer", "promoter", "promoter", "promoter", "promoter", "promoter", "repressed", "repressed",
              "enhancer", "quiescent", "heterochromatin")
  dict = setNames(element, state)
  file$week0 = unname(dict[file$week0])
  file$week2 = unname(dict[file$week2])
  file$week4 = unname(dict[file$week4])
  file$week7 = unname(dict[file$week7])
  file$week10 = unname(dict[file$week10])
  
  colnames(file)
  #### 0week quiescent to 10week enhancer
  data <- file[file$week0 == "quiescent" & file$week10 == "enhancer",]
  data <- data[,c("gene_name","control_0","control_1","control_2","2week_0","2week_1","2week_2","4week_0","4week_1","4week_2",
                  "7week_0","7week_1","7week_2","10week_0","10week_1","10week_2")]
  df <- melt(data = data,id.vars="gene_name",variable.name = "sample",value.name = "fpkm")
  library(stringr)
  df$time <- unlist(str_split_fixed(df$sample,pattern = "_",2))[,1]
  df$time <- factor(df$time,levels = c("control","2week","4week","7week","10week"))
  
  library(ggpubr)
  my_comparisons <- list(c("control","10week"),c("2week","10week"),c("4week","10week"),c("7week",'10week'))
  
  p1 <- ggplot(data = df,aes(x = time,y = log2(fpkm+1),fill=time)) +
    stat_boxplot(geom = "errorbar",linetype=1,width=0.6,position = "identity") +
    geom_boxplot(outlier.fill = NA,outlier.shape=NA,width=0.8) +
    stat_compare_means(comparisons = my_comparisons,label = "p.signif",label.y = c(14,13,12,11)) +
    theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1) +
    labs(title = "0week quiescent to 10weeks enhancer",y="Expresion log2(FPKM + 1)",x="week(s)",fill=NULL) +
    scale_y_continuous(expand = c(0,0),limits = c(0,16)) +
    scale_x_discrete(limits=c("control","2week","4week","7week","10week"),labels=c(0,2,4,7,10)) +
    scale_fill_manual(values = c("#1B9E77","#66A61E","#7570B3","#E7298A","#D95F02"),
                      limits=c("control","2week","4week","7week","10week"),
                      labels = c("0week","2weeks","4weeks","7weeks","10weeks")) +
    theme(aspect.ratio = 0.618/1,
          plot.margin = margin(t = 0,r = 0.3,b = 0,l = 0,unit = "cm"),
          plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
          axis.text.y = element_text(size = 18,colour = 'black'),
          axis.text.x = element_text(size = 20,colour = "black",angle = 0,hjust = 0.5),
          axis.ticks.length = unit(2,"mm")) +
    theme(legend.position = "top",
          legend.text = element_text(size = 18))
  ggsave(filename = paste(geneCountDir,region,"/",region," 0week quiescent to 10week enhancer boxplot.pdf",sep = ""),p1,height = 5,width = 7.3)
  p1
  
  #### 0week promoter to 10week enhancer
  data <- file[file$week0 == "promoter" & file$week10 == "enhancer",]
  data <- data[,c("gene_name","control_0","control_1","control_2","2week_0","2week_1","2week_2","4week_0","4week_1","4week_2",
                  "7week_0","7week_1","7week_2","10week_0","10week_1","10week_2")]
  df <- melt(data = data,id.vars="gene_name",variable.name = "sample",value.name = "fpkm")
  library(stringr)
  df$time <- unlist(str_split_fixed(df$sample,pattern = "_",2))[,1]
  df$time <- factor(df$time,levels = c("control","2week","4week","7week","10week"))
  
  library(ggpubr)
  my_comparisons <- list(c("control","10week"),c("2week","10week"),c("4week","10week"),c("7week",'10week'))
  
  p2 <- ggplot(data = df,aes(x = time,y = log2(fpkm+1),fill=time)) +
    stat_boxplot(geom = "errorbar",linetype=1,width=0.6,position = "identity") +
    geom_boxplot(outlier.fill = NA,outlier.shape=NA,width=0.8) +
    stat_compare_means(comparisons = my_comparisons,label = "p.signif",label.y = c(14,13,12,11)) +
    theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1) +
    labs(title = "0week promoter to 10weeks enhancer",y="Expresion log2(FPKM + 1)",x="week(s)",fill=NULL) +
    scale_y_continuous(expand = c(0,0),limits = c(0,16)) +
    scale_x_discrete(limits=c("control","2week","4week","7week","10week"),labels=c(0,2,4,7,10)) +
    scale_fill_manual(values = c("#1B9E77","#66A61E","#7570B3","#E7298A","#D95F02"),
                      limits=c("control","2week","4week","7week","10week"),
                      labels = c("0week","2weeks","4weeks","7weeks","10weeks")) +
    theme(aspect.ratio = 0.618/1,
          plot.margin = margin(t = 0,r = 0.3,b = 0,l = 0,unit = "cm"),
          plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
          axis.text.y = element_text(size = 18,colour = 'black'),
          axis.text.x = element_text(size = 20,colour = "black",angle = 0,hjust = 0.5),
          axis.ticks.length = unit(2,"mm")) +
    theme(legend.position = "top",
          legend.text = element_text(size = 18))
  ggsave(filename = paste(geneCountDir,region,"/",region," 0week promoter to 10week enhancer boxplot.pdf",sep = ""),p2,height = 5,width = 7.3)
  p2
}

plotBox(filename = "~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody0bp/mm10_genebody_0bp_merge_rna.csv")  
plotBox(filename = "~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody2000bp/mm10_genebody_2000bp_merge_rna.csv") 
plotBox(filename = "~/project/colon cancer/chip-seq/chromHMM/geneCount/genetss4000bp/mm10_genetss_4000bp_merge_rna.csv") 


# genebody 0bp ------------------------------------------------------------
library(clusterProfiler)
library(AnnotationDbi)
file <- read.csv("genebody0bp/mm10_genebody_0bp_merge_rna.csv",check.names = FALSE)
state = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13")
element = c("enhancer", "enhancer", "enhancer", "promoter", "promoter", "promoter", "promoter", "promoter", "repressed", "repressed",
            "enhancer", "quiescent", "heterochromatin")
dict = setNames(element, state)
file$week0 = unname(dict[file$week0])
file$week2 = unname(dict[file$week2])
file$week4 = unname(dict[file$week4])
file$week7 = unname(dict[file$week7])
file$week10 = unname(dict[file$week10])

# quiescent to enhancer
data <- file[file$week0 == "quiescent" & file$week10 == "enhancer",]
bp <- enrichGO(gene = data$gene_name,
               keyType = "SYMBOL",
               OrgDb = "org.Mm.eg.db",
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05
)
bp_0.05 <- bp@result[which(bp@result$pvalue < 0.05),][c(1:3,7,10,14,19),]
ggplot(data = bp_0.05,aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill=-log10(pvalue))) +
  geom_bar(stat = "identity",width = 0.8,position = position_dodge(width = 0.8),fill="#6A51A3") +
  theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) + 
  scale_x_continuous(expand = c(0,0),limits = c(0,30),breaks = c(0,10,20,30)) +
  #scale_fill_gradient(low = "#FEE6CE",high = "#F16913") +
  labs(title = "Biological process of 0week quiescent to 10weeks enhancer ",x="-log10(pValue)",y=NULL) +
  guides(fill=FALSE) +
  theme(aspect.ratio = 1/0.618,
        plot.title = element_text(size = 18,hjust = 0.95,vjust = 0.5),
        plot.margin = margin(t=0,r=0.3,b=0,l=0.3,unit = "cm"),
        axis.text = element_text(colour = 'black',size = 18),
        axis.ticks.y = element_blank())
ggsave(filename = "genebody0bp/genebody0bp bp of 0week quiescent to 10weeks enhancer.pdf",height = 3,width = 7.7)

# promoter to enhancer
promoter <- file[file$week0 == "promoter" & file$week10 == "enhancer",]
promoter_bp <- enrichGO(gene = promoter$gene_name,
               keyType = "SYMBOL",
               OrgDb = "org.Mm.eg.db",
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05
)
promoter_bp_0.05 <- promoter_bp@result[which(promoter_bp@result$pvalue < 0.05),][c(1:4,10,18,34),]
promoter_bp_0.05$Description[6] <- "apoptotic signaling pathway in response to DNA damage"
ggplot(data = promoter_bp_0.05,aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill=-log10(pvalue))) +
  geom_bar(stat = "identity",width = 0.8,position = position_dodge(width = 0.8),fill="#6A51A3") +
  theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) +
  scale_x_continuous(expand = c(0,0),limits = c(0,10),breaks = c(0,5,10)) +
  #scale_fill_gradient(low = "#FEE6CE",high = "#F16913") +
  labs(title = "Biological process of 0week promoter to 10weeks enhancer",x="-log10(pValue)",y=NULL) +
  guides(fill=FALSE) +
  theme(aspect.ratio = 1/0.618,
        plot.title = element_text(size = 18,hjust = 1,vjust = 0.5),
        plot.margin = margin(t=0,r=0.4,b=0,l=0.3,unit = "cm"),
        axis.text = element_text(colour = 'black',size = 18),
        axis.ticks.y = element_blank())
ggsave(filename = "genebody0bp/bp of 0week promoter to 10weeks enhancer.pdf",height = 3,width = 8)


# genebody 2000bp ---------------------------------------------------------
file <- read.csv("genebody2000bp/mm10_genebody_2000bp_merge_rna.csv",check.names = FALSE)
state = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13")
element = c("enhancer", "enhancer", "enhancer", "promoter", "promoter", "promoter", "promoter", "promoter", "repressed", "repressed",
            "enhancer", "quiescent", "heterochromatin")
dict = setNames(element, state)
file$week0 = unname(dict[file$week0])
file$week2 = unname(dict[file$week2])
file$week4 = unname(dict[file$week4])
file$week7 = unname(dict[file$week7])
file$week10 = unname(dict[file$week10])

# quiescent to enhancer
data <- file[file$week0 == "quiescent" & file$week10 == "enhancer",]
bp <- enrichGO(gene = data$gene_name,
               keyType = "SYMBOL",
               OrgDb = "org.Mm.eg.db",
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05
)
bp_0.05 <- bp@result[which(bp@result$pvalue < 0.05),][c(1:4,12,22,26),]
ggplot(data = bp_0.05,aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill=-log10(pvalue))) +
  geom_bar(stat = "identity",width = 0.8,position = position_dodge(width = 0.8),fill="#6A51A3") +
  theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) + 
  scale_x_continuous(expand = c(0,0),limits = c(0,30),breaks = c(0,10,20,30)) +
  #scale_fill_gradient(low = "#FEE6CE",high = "#F16913") +
  labs(title = "Biological process of 0week quiescent to 10weeks enhancer ",x="-log10(pValue)",y=NULL) +
  guides(fill=FALSE) +
  theme(aspect.ratio = 1/0.618,
        plot.title = element_text(size = 18,hjust = 0.95,vjust = 0.5),
        plot.margin = margin(t=0,r=0.3,b=0,l=0.3,unit = "cm"),
        axis.text = element_text(colour = 'black',size = 18),
        axis.ticks.y = element_blank())
ggsave(filename = "genebody2000bp/genebody2000bp bp of 0week quiescent to 10weeks enhancer.pdf",height = 3,width = 7.7)

# promoter to enhancer
promoter <- file[file$week0 == "promoter" & file$week10 == "enhancer",]
promoter_bp <- enrichGO(gene = promoter$gene_name,
                        keyType = "SYMBOL",
                        OrgDb = "org.Mm.eg.db",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05
)
promoter_bp_0.05 <- promoter_bp@result[which(promoter_bp@result$pvalue < 0.05),][c(1:7),]

ggplot(data = promoter_bp_0.05,aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill=-log10(pvalue))) +
  geom_bar(stat = "identity",width = 0.8,position = position_dodge(width = 0.8),fill="#6A51A3") +
  theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) +
  scale_x_continuous(expand = c(0,0),limits = c(0,10),breaks = c(0,5,10)) +
  #scale_fill_gradient(low = "#FEE6CE",high = "#F16913") +
  labs(title = "Biological process of 0week promoter to 10weeks enhancer",x="-log10(pValue)",y=NULL) +
  guides(fill=FALSE) +
  theme(aspect.ratio = 1/0.618,
        plot.title = element_text(size = 18,hjust = 1,vjust = 0.5),
        plot.margin = margin(t=0,r=0.4,b=0,l=0.3,unit = "cm"),
        axis.text = element_text(colour = 'black',size = 18),
        axis.ticks.y = element_blank())
ggsave(filename = "genebody2000bp/genebody2000bp bp of 0week promoter to 10weeks enhancer.pdf",height = 3,width = 8)


# genetss 4000bp ----------------------------------------------------------
file <- read.csv("genetss4000bp/mm10_genetss_4000bp_merge_rna.csv",check.names = FALSE)
state = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13")
element = c("enhancer", "enhancer", "enhancer", "promoter", "promoter", "promoter", "promoter", "promoter", "repressed", "repressed",
            "enhancer", "quiescent", "heterochromatin")
dict = setNames(element, state)
file$week0 = unname(dict[file$week0])
file$week2 = unname(dict[file$week2])
file$week4 = unname(dict[file$week4])
file$week7 = unname(dict[file$week7])
file$week10 = unname(dict[file$week10])

# quiescent to enhancer
data <- file[file$week0 == "quiescent" & file$week10 == "enhancer",]
bp <- enrichGO(gene = data$gene_name,
               keyType = "SYMBOL",
               OrgDb = "org.Mm.eg.db",
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05
)
bp_0.05 <- bp@result[which(bp@result$pvalue < 0.05),][c(1:6,14),]
ggplot(data = bp_0.05,aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill=-log10(pvalue))) +
  geom_bar(stat = "identity",width = 0.8,position = position_dodge(width = 0.8),fill="#6A51A3") +
  theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) + 
  scale_x_continuous(expand = c(0,0),limits = c(0,8),breaks = c(0,4,8)) +
  #scale_fill_gradient(low = "#FEE6CE",high = "#F16913") +
  labs(title = "Biological process of 0week quiescent to 10weeks enhancer ",x="-log10(pValue)",y=NULL) +
  guides(fill=FALSE) +
  theme(aspect.ratio = 1/0.618,
        plot.title = element_text(size = 18,hjust = 0.95,vjust = 0.5),
        plot.margin = margin(t=0,r=0.3,b=0,l=0.3,unit = "cm"),
        axis.text = element_text(colour = 'black',size = 18),
        axis.ticks.y = element_blank())
ggsave(filename = "genetss4000bp/genetss4000bp bp of 0week quiescent to 10weeks enhancer.pdf",height = 3,width = 8.5)

# promoter to enhancer
promoter <- file[file$week0 == "promoter" & file$week10 == "enhancer",]
promoter_bp <- enrichGO(gene = promoter$gene_name,
                        keyType = "SYMBOL",
                        OrgDb = "org.Mm.eg.db",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05
)
promoter_bp_0.05 <- promoter_bp@result[which(promoter_bp@result$pvalue < 0.05),][c(1,2,5:7,14,26),]

ggplot(data = promoter_bp_0.05,aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill=-log10(pvalue))) +
  geom_bar(stat = "identity",width = 0.8,position = position_dodge(width = 0.8),fill="#6A51A3") +
  theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) + 
  scale_x_continuous(expand = c(0,0),limits = c(0,20),breaks = c(0,10,20)) +
  #scale_fill_gradient(low = "#FEE6CE",high = "#F16913") +
  labs(title = "Biological process of 0week promoter to 10weeks enhancer",x="-log10(pValue)",y=NULL) +
  guides(fill=FALSE) +
  theme(aspect.ratio = 1/0.618,
        plot.title = element_text(size = 18,hjust = 1,vjust = 0.5),
        plot.margin = margin(t=0,r=0.4,b=0,l=0.3,unit = "cm"),
        axis.text = element_text(colour = 'black',size = 18),
        axis.ticks.y = element_blank())
ggsave(filename = "genetss4000bp/genetss4000bp bp of 0week promoter to 10weeks enhancer.pdf",height = 3,width = 8)





rm(list = ls())

plotScale <- function(cluster,ymax,marker) {
  file <- read.csv(paste("/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/masigpro_rna/masigpro_rna_",cluster,
                         "/RNA-seq_",cluster,"_",marker,"_scale_line.txt",sep = ""),
                   sep = "\t",header = TRUE,check.names = FALSE)
  data <- file[,1:30002]
  data <- data[,-2]
  library(reshape2)
  df <- melt(data = data,id.vars = "bins",variable.name = "location",value.name = "number")
  print(df[1:5,])
  Times <- c("ctrl","2weeks","4weeks","7weeks","10weeks")
  sample_vector <-c()
  for (time in Times) {
    sample <- paste(time,marker,sep = "-")
    sample_vector <- c(sample_vector,sample)
  }

  # H3K27me3
  df$bins <- factor(c(paste("4weeks-",marker,sep = ""),paste("ctrl-",marker,sep = ""),paste("2weeks-",marker,sep = ""),
                      paste("7weeks-",marker,sep = ""),paste("10weeks-",marker,sep = "")

  ),
  levels = sample_vector)
  df$location <- as.numeric(df$location)

  library(ggplot2)
  p <- ggplot(df, aes(x = location,y = number,group=bins,colour=bins)) +
    geom_line(size=0.9) +
    #stat_smooth(geom = "smooth",method = "loess",formula = y ~ x,se = FALSE,span=0.01) +
    theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1) +
    scale_y_continuous(expand = c(0,0),limits = c(0,ymax)) +
    scale_colour_manual(values = c("#1B9E77","#66A61E","#7570B3","#E7298A","#D95F02"),labels=c("0week","2weeks","4weeks","7weeks","10weeks")) +
    scale_x_continuous(expand = c(0,0),limits=c(0,30000),breaks=c(0,10000,20000,30000),labels=c("-10kb","TSS","TES","10kb")) +
    labs(title = paste(cluster,sep = ""),x=NULL,y="H3K27me3 ChIP-seq\nintensity (RPKM)",colour=NULL) +
    theme(aspect.ratio = 0.6/1,
          plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,margin = margin(b = 0.5,unit = "cm")),
          plot.margin = margin(t = 0,r = 0.3,b = 0,l = 0,unit = "cm"),
          axis.title.x = element_text(margin = margin(t = 0.4,r = 0,b = 0,l = 0,unit = "cm")),
          axis.title.y = element_text(margin = margin(t = 0,r = 0.4,b = 0,l = 0,unit = "cm")),
          axis.text.x = element_text(size = 18,colour = "black",margin = margin(t = 0.1,r = 0,b = 0,l = 0,unit = "cm")),
          axis.text.y = element_text(size = 18,colour = "black",margin = margin(t = 0,r = 0.1,b = 0,l = 0,unit = "cm")),
          axis.ticks.length = unit(0.25,'cm')) +
    theme(legend.position = "none")
    # theme(legend.position = c(0.8,0.8),
    #       legend.background = element_blank(),
    #       legend.key = element_blank(),
    #       legend.text = element_text(size = 18,family = "sans",lineheight = 0.7),
    #       legend.key.height = unit(0.5,"cm")
    #       )
  ggsave(p,filename = paste("/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/masigpro_rna/masigpro_rna_",cluster,"/",
                         marker,"_RNA-seq_",cluster,"_scale.pdf",sep = ""),width = 5,height = 3.5)
  ggsave(filename = paste("/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/masigpro_rna/masigpro_rna_",cluster,"/",
                          marker,"_RNA-seq_",cluster,"_scale.png",sep = ""),width = 4.60,height = 2.60,dpi =300)
}

plotScale(cluster = "cluster1",ymax = 5,marker = "H3K27me3")
plotScale(cluster = "cluster2",ymax = 4,marker = "H3K27me3")
plotScale(cluster = "cluster3",ymax = 6,marker = "H3K27me3")
plotScale(cluster = "cluster4",ymax = 5,marker = "H3K27me3")

plotScale2 <- function(cluster,ymax,marker) {
  file <- read.csv(paste("/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/masigpro_rna/masigpro_rna_",cluster,
                         "/RNA-seq_",cluster,"_",marker,"_scale_line.txt",sep = ""),
                   sep = "\t",header = TRUE,check.names = FALSE)
  data <- file[,1:30002]
  data <- data[,-2]
  library(reshape2)
  df <- melt(data = data,id.vars = "bins",variable.name = "location",value.name = "number")
  print(df[1:5,])
  Times <- c("ctrl","2weeks","4weeks","7weeks","10weeks")
  sample_vector <-c()
  for (time in Times) {
    sample <- paste(time,marker,sep = "-")
    sample_vector <- c(sample_vector,sample)
  }

  # H3K9me3
  df$bins <- factor(c(paste("2weeks-",marker,sep = ""),paste("ctrl-",marker,sep = ""),paste("10weeks-",marker,sep = ""),
                      paste("7weeks-",marker,sep = ""),paste("4weeks-",marker,sep = "")
  ),
  levels = sample_vector)
  df$location <- as.numeric(df$location)

  library(ggplot2)
  p <- ggplot(df, aes(x = location,y = number,group=bins,colour=bins)) +
    geom_line(size=0.9) +
    #stat_smooth(geom = "smooth",method = "loess",formula = y ~ x,se = FALSE,span=0.01) +
    theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1) +
    #scale_y_continuous(expand = c(0,0),limits = c(0,ymax)) +
    scale_colour_manual(values = c("#1B9E77","#66A61E","#7570B3","#E7298A","#D95F02"),labels=c("0week","2weeks","4weeks","7weeks","10weeks")) +
    scale_x_continuous(expand = c(0,0),limits=c(0,30000),breaks=c(0,10000,20000,30000),labels=c("-10kb","TSS","TES","10kb")) +
    labs(title = paste(cluster,sep = ""),x=NULL,y=paste(marker," ChIP-seq\nintensity (RPKM)",sep = ""),colour=NULL) +
    theme(aspect.ratio = 0.6/1,
          plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,margin = margin(b = 0.5,unit = "cm")),
          plot.margin = margin(t = 0,r = 0.3,b = 0,l = 0,unit = "cm"),
          axis.title.x = element_text(margin = margin(t = 0.4,r = 0,b = 0,l = 0,unit = "cm")),
          axis.title.y = element_text(margin = margin(t = 0,r = 0.4,b = 0,l = 0,unit = "cm")),
          axis.text.x = element_text(size = 18,colour = "black",margin = margin(t = 0.1,r = 0,b = 0,l = 0,unit = "cm")),
          axis.text.y = element_text(size = 18,colour = "black",margin = margin(t = 0,r = 0.1,b = 0,l = 0,unit = "cm")),
          axis.ticks.length = unit(0.25,'cm')) +
    theme(legend.position = "none")
    # theme(legend.position = c(0.8,0.8),
    #       legend.background = element_blank(),
    #       legend.key = element_blank(),
    #       legend.text = element_text(size = 18,family = "sans",lineheight = 0.7),
    #       legend.key.height = unit(0.5,"cm")
    #       )
  ggsave(p,filename = paste("/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/masigpro_rna/masigpro_rna_",cluster,"/",
                         marker,"_RNA-seq_",cluster,"_scale.pdf",sep = ""),width = 5,height = 3.5)
  ggsave(filename = paste("/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/masigpro_rna/masigpro_rna_",cluster,"/",
                          marker,"_RNA-seq_",cluster,"_scale.png",sep = ""),width = 4.60,height = 2.60,dpi =300)
}
plotScale2(cluster = "cluster1",ymax = 6,marker = "H3K9me3")
plotScale2(cluster = "cluster2",ymax = 6,marker = "H3K9me3")
plotScale2(cluster = "cluster3",ymax = 6,marker = "H3K9me3")
plotScale2(cluster = "cluster4",ymax = 6,marker = "H3K9me3")
rm(list = ls())
setwd("~/project/colon cancer/nfkb target gene/")

library(ggplot2)
file <- read.delim("nfkb_target_138_gene_tss_line.txt",sep = "\t")

plotProfile <- function(file,cluster,ymax) {
  data <- file[,1:6002]
  data <- data[,-2]
  library(reshape2)
  df <- melt(data = data,id.vars = "bins",variable.name = "location",value.name = "number")
  df$bins <- factor(c("7weeks-1-H3K27ac","10weeks-1-H3K27ac","ctrl-1-H3K27ac","2weeks-1-H3K27ac","4weeks-1-H3K27ac"),
                    levels = c("ctrl-1-H3K27ac","2weeks-1-H3K27ac","4weeks-1-H3K27ac","7weeks-1-H3K27ac","10weeks-1-H3K27ac"))
  df$location <- as.numeric(df$location)
  p <- ggplot(df, aes(x = location,y = number,group=bins,colour=bins)) + 
    #geom_line(size=0.9) +
    geom_smooth(stat = "smooth",method = "loess",formula = y ~ x,se = FALSE,span=0.01) +
    theme_classic(base_family = "sans",base_size = 16,base_line_size = 1.1) +
    scale_y_continuous(expand = c(0,0),limits = c(0,ymax)) +
    scale_colour_manual(values = c("#1B9E77","#66A61E","#7570B3","#E7298A","#D95F02"),labels=c("0week","2weeks","4weeks","7weeks","10weeks")) +
    scale_x_continuous(expand = c(0,0),limits=c(0,6000),breaks=c(0,3000,6000),labels=c("-3kb","TSS","3kb")) +
    labs(title = "",x="Distance from TSS (kb)",y="H3K27ac ChIP-seq intensity (RPKM)",colour=NULL) +
    theme(aspect.ratio = 1/0.8,
          plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,margin = margin(b = 0.5,unit = "cm")),
          plot.margin = margin(t = 0,r = 0.2,b = 0,l = 0,unit = "cm"),
          axis.title.x = element_text(margin = margin(t = 0.4,r = 0,b = 0,l = 0,unit = "cm")),
          axis.title.y = element_text(margin = margin(t = 0,r = 0.4,b = 0,l = 0,unit = "cm")),
          axis.text.x = element_text(size = 16,colour = "black",margin = margin(t = 0.1,r = 0,b = 0,l = 0,unit = "cm")),
          axis.text.y = element_text(size = 16,colour = "black",margin = margin(t = 0,r = 0.1,b = 0,l = 0,unit = "cm")),
          axis.ticks.length = unit(0.25,'cm')) +
    theme(legend.position = c(0.9,0.85),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size = 16,family = "sans",lineheight = 0.7),
          legend.key.height = unit(0.5,"cm")
    )
  ggsave(paste(filename = "~/project/colon cancer/nfkb target gene/nfkb_",cluster,"_gene_TSS_line.pdf"),p,width = 4.88,height = 4.6)
  p
}

cluster1 <- read.delim("~/project/colon cancer/nfkb target gene/nfkb_cluster1_gene_tss_line.txt",sep = "\t")
cluster2 <- read.delim("~/project/colon cancer/nfkb target gene/nfkb_cluster2_gene_tss_line.txt",sep = "\t")
cluster3 <- read.delim("~/project/colon cancer/nfkb target gene/nfkb_cluster3_gene_tss_line.txt",sep = "\t")
plotProfile(file = cluster1,cluster = "cluster1",ymax = 20)
plotProfile(file = cluster2,cluster = "cluster2",ymax = 12)
plotProfile(file = cluster3,cluster = "cluster3",ymax = 10)


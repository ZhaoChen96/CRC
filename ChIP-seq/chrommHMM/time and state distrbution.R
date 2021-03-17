rm(list = ls())
setwd("~/project/colon cancer/chip-seq/chromHMM/geneCount/")
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(plyr)

plotPercent <- function(filename) {
  generegion <- unlist(strsplit(filename,split = "/"))[1]
  file <- read.csv(filename)
  list <- c(file$week0,file$week2,file$week4,file$week7,file$week10)
  state_list <- c()
  for (i in 1:13) {
    state <- paste("E",i,sep = "") 
    state_list <- c(state_list,state)
  }
  state_list <- sort(state_list)
  week0 <- table(file$week0)
  week2 <- table(file$week2)
  week4 <- table(file$week4)
  week7 <- table(file$week7)
  week10 <- table(file$week10)
  data <- data.frame(as.numeric(week0),as.numeric(week2),as.numeric(week4),as.numeric(week7),as.numeric(week10))
  df <- cbind(state_list,data)
  names(df) <- c("state","gene_week0","gene_week2","gene_week4","gene_week7","gene_week10")
  
  library(reshape2)
  df <- melt(df, variable.name = "time",id.vars = "state",value.name = "number")
  df$state <- factor(df$state,levels = c("E12","E13","E9","E10","E1","E11","E3","E2","E6","E7","E4","E5","E8"))

  ce <- ddply(df, "state",transform, percent= number/sum(number) * 100)
  cd <- ddply(df, "time",transform, percent= number/sum(number) * 100)
  ce$state = factor(ce$state, levels = c("E12","E13","E9","E10","E1","E11","E3","E2","E6","E7","E4","E5","E8"))
  cd$time = factor(c("gene_week0","gene_week2","gene_week4","gene_week7","gene_week10"),
                   levels = c("gene_week0","gene_week2","gene_week4","gene_week7","gene_week10"))
  

  
  ggplot(ce, aes(x = state,y = percent,fill=time)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("#1B9E77","#66A61E","#7570B3","#E7298A","#D95F02"),labels=c("0week","2weeks","4weeks","7weeks","10weeks")) +
    scale_y_continuous(expand = c(0,0),limits=rev(levels(state_list))) +
    theme_classic(base_family = "sans",base_size = 18) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black")) +
    labs(y="Percent (%)",x=NULL,fill=NULL) +
    theme(legend.position = "bottom")
  ggsave(paste(generegion,"/different time in each state.pdf",sep = ""),width = 7,height = 3.7)
  
  heterochromatin_color = brewer.pal(n = 11, name = "PRGn")[3]
  promoter_color = brewer.pal(n = 9, name = "Greens")[8]
  enhancer_color = brewer.pal(n = 9, name = "Reds")[8]
  repressed_color = brewer.pal(n = 9, name = "Blues")[7]
  quit_color = brewer.pal(n = 9, name = "Greys")[3]
  colour_list <- c(quit_color,rep(c(heterochromatin_color,repressed_color),each=2),rep(enhancer_color,3),rep(promoter_color,5))
  
  ggplot(df,aes(y = time, x = cd$percent, fill=state)) +
    geom_bar(stat = "identity") +
    scale_y_discrete(limits=rev(levels(cd$time)),labels=rev(c("0week","2weeks","4weeks","7weeks","10weeks"))) +
    scale_x_continuous(expand = c(0,0)) +
    labs(x="Percent (%)",y=NULL,fill=NULL) +
    #scale_fill_manual(values = colour_list) +
    # set3
    #scale_fill_manual(values = c(quit_color,"#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69",
    #                             "#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F")) +
    # paired
    #scale_fill_manual(values = c(quit_color,"#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F",
    #                             "#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928")) +
    scale_fill_manual(values = c(quit_color,heterochromatin_color,"#6BAED6","#08519C",brewer.pal(n=4,name = "Reds"),
                                 "#74C476","#41AB5D","#238B45","#006D2C","#00441B")) +
    theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1) +
    theme(plot.margin = margin(t = 0.5,r = 0.5,b = 0.5,l = 0.5,unit = "cm"),
          axis.text = element_text(colour = "black",size = 18),
          axis.ticks.length = unit(2,"mm")) 
  guides(fill=FALSE)
  ggsave(filename = paste(generegion,"/different states in each time.pdf",sep = ""),width = 7,height = 3.7)
}

filename <- "genebody0bp/mm10_genebody_0bp_merge_rna.csv"
plotPercent(filename = filename)
filename <- "genebody2000bp/mm10_genebody_2000bp_merge_rna.csv"
plotPercent(filename = filename)
filename <- "genetss4000bp/mm10_genetss_4000bp_merge_rna.csv"
plotPercent(filename = filename)

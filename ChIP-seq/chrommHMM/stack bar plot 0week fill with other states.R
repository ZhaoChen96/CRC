rm(list = ls())
setwd("~/project/colon cancer/chip-seq/chromHMM/geneCount/")
library(dplyr)
library(ggplot2)
geneCountDir <- "~/project/colon cancer/chip-seq/chromHMM/geneCount/"


library("RColorBrewer")
heterochromatin_color = brewer.pal(n = 11, name = "PRGn")[3]
promoter_color = brewer.pal(n = 9, name = "Greens")[8]
enhancer_color = brewer.pal(n = 9, name = "Reds")[8]
repressed_color = brewer.pal(n = 9, name = "Blues")[7]
quit_color = brewer.pal(n = 9, name = "Greys")[3]

plotBar <- function(filename) {
  file <- read.csv(filename,check.names = FALSE)
  region <- paste(unlist(strsplit(basename(filename),split = "_"))[2],unlist(strsplit(basename(filename),split = "_"))[3],sep = "")
  df = na.omit(file)
  
  df_sub = df[,c("gene_name", "week0", "week2", "week4", "week7", "week10")]
  state = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13")
  element = c("enhancer", "enhancer", "enhancer","promoter",  "promoter", "promoter", "promoter", "promoter",
              "repressed", "repressed","enhancer", "quiescent","heterochromatin")
  dict = setNames(element, state)
  df_sub$week0 = unname(dict[df_sub$week0])
  df_sub$week2 = unname(dict[df_sub$week2])
  df_sub$week4 = unname(dict[df_sub$week4])
  df_sub$week7 = unname(dict[df_sub$week7])
  df_sub$week10 = unname(dict[df_sub$week10])
  
  # 0week to 10weeks each state including -----------------------------------
  enhancer <- df_sub[df_sub$week10 == "enhancer",]
  enhancer <- count(enhancer$week0)
  enhancer$state <- "enhancer"
  
  promoter <- df_sub[df_sub$week10 == "promoter",]
  promoter <- count(promoter$week0)
  promoter$state <- "promoter"
  
  repressed <- df_sub[df_sub$week10 == "repressed",]
  repressed <- count(repressed$week0)
  repressed$state <- "repressed"
  
  heterochromatin <- df_sub[df_sub$week10 == "heterochromatin",]
  heterochromatin <- count(heterochromatin$week0)
  heterochromatin$state <- "heterochromatin"
  
  quiescent <- df_sub[df_sub$week10 == "quiescent",]
  quiescent <- count(quiescent$week0)
  quiescent$state <- "quiescent"
  
  data <- rbind(enhancer,promoter,repressed,heterochromatin,quiescent)
  data$x<- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                       levels = c("quiescent","heterochromatin","repressed","enhancer","promoter"))
  data$freq <- as.numeric(data$freq)
  
  library(plyr)
  library(ggplot2)
  ce <- ddply(data,"state",transform,percent_state = freq/sum(freq) * 100)
  ggplot(ce,aes(x = state,y = percent_state,fill=x)) +
    geom_bar(stat = "identity",width = 0.9) +
    scale_x_discrete(limits=c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c(enhancer = enhancer_color,promoter = promoter_color, repressed = repressed_color,
                                 heterochromatin = heterochromatin_color,quiescent = quit_color),
                      breaks = c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) +
    labs(title = "10weeks",x=NULL,y="Percentage (%)",fill=NULL) +
    theme(plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,face = "bold"),
          axis.title = element_text(colour = "black",size = 18),
          axis.text.y = element_text(size = 18,colour = "black"),
          axis.text.x = element_text(size = 18,colour = "black",angle = 45,hjust = 1),
          axis.ticks.length = unit(2,units = "mm")) +
    #guides(fill=FALSE)
    theme(legend.position = "right")
  ggsave(filename = paste(geneCountDir,region,"/",region,"_10weeks state contribution_legend.pdf",sep = ""),height = 4,width = 5)
  
  # 0week to 7weeks ecah states distrbution  -------------------------------------
  enhancer <- df_sub[df_sub$week7 == "enhancer",]
  enhancer <- count(enhancer$week0)
  enhancer$state <- "enhancer"
  promoter <- df_sub[df_sub$week7 == "promoter",]
  promoter <- count(promoter$week0)
  promoter$state <- "promoter"
  repressed <- df_sub[df_sub$week7 == "repressed",]
  repressed <- count(repressed$week0)
  repressed$state <- "repressed"
  heterochromatin <- df_sub[df_sub$week7 == "heterochromatin",]
  heterochromatin <- count(heterochromatin$week0)
  if (region == "genetss4000bp") {
    heterochromatin <- rbind(c("enhancer",0),heterochromatin[1:4,])
  }
  heterochromatin$state <- "heterochromatin"
  quiescent <- df_sub[df_sub$week7 == "quiescent",]
  quiescent <- count(quiescent$week0)
  quiescent$state <- "quiescent"
  
  data <- rbind(enhancer,promoter,repressed,heterochromatin,quiescent)
  data$x <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                   levels = c("quiescent","heterochromatin","repressed","enhancer","promoter"))
  data$freq <- as.numeric(data$freq)
  
  library(plyr)
  library(ggplot2)
  ce <- ddply(data,"state",transform,percent_state = freq/sum(freq) * 100)
  ggplot(ce,aes(x = state,y = percent_state,fill=x)) +
    geom_bar(stat = "identity",width = 0.9) +
    scale_x_discrete(limits=c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c(enhancer = enhancer_color,promoter = promoter_color, repressed = repressed_color,
                                 heterochromatin = heterochromatin_color,quiescent = quit_color),
                      breaks = c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) +
    labs(title = "7weeks",x=NULL,y="Percentage (%)",fill=NULL) +
    theme(plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,face = "bold"),
          axis.title = element_text(colour = "black",size = 18),
          axis.text.y = element_text(size = 18,colour = "black"),
          axis.text.x = element_text(size = 18,colour = "black",angle = 45,hjust = 1),
          axis.ticks.length = unit(2,units = "mm")) +
    #guides(fill=FALSE)
    theme(legend.position = "right")
  ggsave(filename = paste(geneCountDir,region,"/",region,"_7weeks state contribution_legend.pdf",sep = ""),height = 4,width = 5)
  
  # 0week to 4weeks ---------------------------------------------------------
  enhancer <- df_sub[df_sub$week4 == "enhancer",]
  enhancer <- count(enhancer$week0)
  enhancer$state <- "enhancer"
  promoter <- df_sub[df_sub$week4 == "promoter",]
  promoter <- count(promoter$week0)
  promoter$state <- "promoter"
  repressed <- df_sub[df_sub$week4 == "repressed",]
  repressed <- count(repressed$week0)
  repressed$state <- "repressed"
  heterochromatin <- df_sub[df_sub$week4 == "heterochromatin",]
  heterochromatin <- count(heterochromatin$week0)
  heterochromatin$state <- "heterochromatin"
  quiescent <- df_sub[df_sub$week4 == "quiescent",]
  quiescent <- count(quiescent$week0)
  quiescent$state <- "quiescent"
  
  data <- rbind(enhancer,promoter,repressed,heterochromatin,quiescent)
  data$x <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                   levels = c("quiescent","heterochromatin","repressed","enhancer","promoter"))
  data$freq <- as.numeric(data$freq)
  
  library(plyr)
  library(ggplot2)
  ce <- ddply(data,"state",transform,percent_state = freq/sum(freq) * 100)
  p4 <- ggplot(ce,aes(x = state,y = percent_state,fill=x)) +
    geom_bar(stat = "identity",width = 0.9) +
    scale_x_discrete(limits=c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c(enhancer = enhancer_color,promoter = promoter_color, repressed = repressed_color,
                                 heterochromatin = heterochromatin_color,quiescent = quit_color),
                      breaks = c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    theme_classic(base_size = 18,base_family = "sans") +
    labs(title = "4weeks",x=NULL,y="Percentage (%)",fill=NULL) +
    theme(plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,face = "bold"),
          axis.title = element_text(colour = "black",size = 18),
          axis.text.y = element_text(size = 18,colour = "black"),
          axis.text.x = element_text(size = 18,colour = "black",angle = 45,hjust = 1),
          axis.ticks.length = unit(2,units = "mm")) +
    #guides(fill=FALSE)
    theme(legend.position = "right")
  ggsave(filename = paste(geneCountDir,region,"/",region,"_4weeks state contribution_legend.pdf",sep = ""),height = 4,width = 5)
  
  # 0week to 2weeks ---------------------------------------------------------
  enhancer <- df_sub[df_sub$week2 == "enhancer",]
  enhancer <- count(enhancer$week0)
  enhancer$state <- "enhancer"
  promoter <- df_sub[df_sub$week2 == "promoter",]
  promoter <- count(promoter$week0)
  #promoter <- rbind(promoter[1,],c("heterochromatin",0),promoter[2:4,])
  promoter$state <- "promoter"
  repressed <- df_sub[df_sub$week2 == "repressed",]
  repressed <- count(repressed$week0)
  repressed$state <- "repressed"
  heterochromatin <- df_sub[df_sub$week2 == "heterochromatin",]
  heterochromatin <- count(heterochromatin$week0)
  #heterochromatin <- rbind(heterochromatin[1:2,],c("promoter",0),heterochromatin[3:4,])
  heterochromatin$state <- "heterochromatin"
  quiescent <- df_sub[df_sub$week2 == "quiescent",]
  quiescent <- count(quiescent$week0)
  quiescent$state <- "quiescent"
  
  data <- rbind(enhancer,promoter,repressed,heterochromatin,quiescent)
  data$x <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                   levels = c("quiescent","heterochromatin","repressed","enhancer","promoter"))
  data$freq <- as.numeric(data$freq)
  
  library(plyr)
  library(ggplot2)
  ce <- ddply(data,"state",transform,percent_state = freq/sum(freq) * 100)
  p2 <- ggplot(ce,aes(x = state,y = percent_state,fill=x)) +
    geom_bar(stat = "identity",width = 0.9) +
    scale_x_discrete(limits=c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c(enhancer = enhancer_color,promoter = promoter_color, repressed = repressed_color,
                                 heterochromatin = heterochromatin_color,quiescent = quit_color),
                      breaks = c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    theme_classic(base_size = 18,base_family = "sans") +
    labs(title = "2weeks",x=NULL,y="Percentage (%)",fill=NULL) +
    theme(plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,face = "bold"),
          axis.title = element_text(colour = "black",size = 18),
          axis.text.y = element_text(size = 18,colour = "black"),
          axis.text.x = element_text(size = 18,colour = "black",angle = 45,hjust = 1),
          axis.ticks.length = unit(2,units = "mm")) +
    #guides(fill=FALSE)
    theme(legend.position = "right")
  p2
  ggsave(filename = paste(geneCountDir,region,"/",region,"_2weeks state contribution_legend.pdf",sep = ""),height = 4,width = 5)
}

plotBar(filename = "~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody0bp/mm10_genebody_0bp_merge_rna.csv")
plotBar(filename = "~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody2000bp/mm10_genebody_2000bp_merge_rna.csv")
plotBar(filename = "~/project/colon cancer/chip-seq/chromHMM/geneCount/genetss4000bp/mm10_genetss_4000bp_merge_rna.csv")


plotStack <- function(filename) {
  file <- read.csv(filename,check.names = FALSE)
  region <- paste(unlist(strsplit(basename(filename),split = "_"))[2],unlist(strsplit(basename(filename),split = "_"))[3],sep = "")
  df = na.omit(file)
  
  df_sub = df[,c("gene_name", "week0", "week2", "week4", "week7", "week10")]
  state = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13")
  element = c("enhancer", "enhancer", "enhancer","promoter",  "promoter", "promoter", "promoter", "promoter",
              "repressed", "repressed","enhancer", "quiescent","heterochromatin")
  dict = setNames(element, state)
  df_sub$week0 = unname(dict[df_sub$week0])
  df_sub$week2 = unname(dict[df_sub$week2])
  df_sub$week4 = unname(dict[df_sub$week4])
  df_sub$week7 = unname(dict[df_sub$week7])
  df_sub$week10 = unname(dict[df_sub$week10])
  
  # 2weeks fill with 0week state --------------------------------------------
  enhancer <- df_sub[df_sub$week0 == "enhancer",]
  enhancer <- count(enhancer$week2)
  enhancer$state <- "enhancer"
  promoter <- df_sub[df_sub$week0 == "promoter",]
  promoter <- count(promoter$week2)
  #promoter <- rbind(promoter[1,],c("heterochromatin",0),promoter[2:4,])
  promoter$state <- "promoter"
  repressed <- df_sub[df_sub$week0 == "repressed",]
  repressed <- count(repressed$week2)
  repressed$state <- "repressed"
  heterochromatin <- df_sub[df_sub$week0 == "heterochromatin",]
  heterochromatin <- count(heterochromatin$week2)
  #heterochromatin <- rbind(heterochromatin[1:2,],c("promoter",0),heterochromatin[3:4,])
  heterochromatin$state <- "heterochromatin"
  quiescent <- df_sub[df_sub$week0 == "quiescent",]
  quiescent <- count(quiescent$week2)
  quiescent$state <- "quiescent"
  
  data <- rbind(enhancer,promoter,repressed,heterochromatin,quiescent)
  data$x <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                   levels = c("quiescent","heterochromatin","repressed","enhancer","promoter"))
  data$freq <- as.numeric(data$freq)
  
  library(plyr)
  library(ggplot2)
  ce <- ddply(data,"state",transform,percent_state = freq/sum(freq) * 100)
  ggplot(ce,aes(x = state,y = percent_state,fill=x)) +
    geom_bar(stat = "identity",width = 0.9) +
    scale_x_discrete(limits=c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c(enhancer = enhancer_color,promoter = promoter_color, repressed = repressed_color,
                                 heterochromatin = heterochromatin_color,quiescent = quit_color),
                      breaks = c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    theme_classic(base_size = 18,base_family = "sans") +
    labs(title = "2weeks",x=NULL,y="Percentage (%)",fill=NULL) +
    theme(plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,face = "bold"),
          axis.title = element_text(colour = "black",size = 18),
          axis.text.y = element_text(size = 18,colour = "black"),
          axis.text.x = element_text(size = 18,colour = "black",angle = 45,hjust = 1),
          axis.ticks.length = unit(2,units = "mm")) +
    #guides(fill=FALSE)
    theme(legend.position = "right")
  ggsave(filename = paste(geneCountDir,region,"/",region,"_0week state fill with 2weeks state.pdf",sep = ""),height = 4,width = 5)
  
  # 4weeks fill with 0week state --------------------------------------------
  enhancer <- df_sub[df_sub$week0 == "enhancer",]
  enhancer <- count(enhancer$week4)
  enhancer$state <- "enhancer"
  promoter <- df_sub[df_sub$week0 == "promoter",]
  promoter <- count(promoter$week4)
  #promoter <- rbind(promoter[1,],c("heterochromatin",0),promoter[2:4,])
  promoter$state <- "promoter"
  repressed <- df_sub[df_sub$week0 == "repressed",]
  repressed <- count(repressed$week4)
  repressed$state <- "repressed"
  heterochromatin <- df_sub[df_sub$week0 == "heterochromatin",]
  heterochromatin <- count(heterochromatin$week4)
  #heterochromatin <- rbind(heterochromatin[1:2,],c("promoter",0),heterochromatin[3:4,])
  heterochromatin$state <- "heterochromatin"
  quiescent <- df_sub[df_sub$week0 == "quiescent",]
  quiescent <- count(quiescent$week4)
  quiescent$state <- "quiescent"
  
  data <- rbind(enhancer,promoter,repressed,heterochromatin,quiescent)
  data$x <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                   levels = c("quiescent","heterochromatin","repressed","enhancer","promoter"))
  data$freq <- as.numeric(data$freq)
  
  library(plyr)
  library(ggplot2)
  ce <- ddply(data,"state",transform,percent_state = freq/sum(freq) * 100)
  ggplot(ce,aes(x = state,y = percent_state,fill=x)) +
    geom_bar(stat = "identity",width = 0.9) +
    scale_x_discrete(limits=c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c(enhancer = enhancer_color,promoter = promoter_color, repressed = repressed_color,
                                 heterochromatin = heterochromatin_color,quiescent = quit_color),
                      breaks = c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    theme_classic(base_size = 18,base_family = "sans") +
    labs(title = "4weeks",x=NULL,y="Percentage (%)",fill=NULL) +
    theme(plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,face = "bold"),
          axis.title = element_text(colour = "black",size = 18),
          axis.text.y = element_text(size = 18,colour = "black"),
          axis.text.x = element_text(size = 18,colour = "black",angle = 45,hjust = 1),
          axis.ticks.length = unit(2,units = "mm")) +
    #guides(fill=FALSE)
    theme(legend.position = "right")
  ggsave(filename = paste(geneCountDir,region,"/",region,"_0week state fill with 4weeks state.pdf",sep = ""),height = 4,width = 5)
  
  # 7weeks fill with 0week state --------------------------------------------
  enhancer <- df_sub[df_sub$week0 == "enhancer",]
  enhancer <- count(enhancer$week7)
  if (region == "genetss4000bp") {
    enhancer <- rbind(enhancer[1,],c("heterochromatin",0),enhancer[2:4,])
  }
  enhancer$state <- "enhancer"
  promoter <- df_sub[df_sub$week0 == "promoter",]
  promoter <- count(promoter$week7)
  #promoter <- rbind(promoter[1,],c("heterochromatin",0),promoter[2:4,])
  promoter$state <- "promoter"
  repressed <- df_sub[df_sub$week0 == "repressed",]
  repressed <- count(repressed$week7)
  repressed$state <- "repressed"
  heterochromatin <- df_sub[df_sub$week0 == "heterochromatin",]
  heterochromatin <- count(heterochromatin$week7)
  #heterochromatin <- rbind(heterochromatin[1:2,],c("promoter",0),heterochromatin[3:4,])
  heterochromatin$state <- "heterochromatin"
  quiescent <- df_sub[df_sub$week0 == "quiescent",]
  quiescent <- count(quiescent$week7)
  quiescent$state <- "quiescent"
  
  data <- rbind(enhancer,promoter,repressed,heterochromatin,quiescent)
  data$x <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                   levels = c("quiescent","heterochromatin","repressed","enhancer","promoter"))
  data$freq <- as.numeric(data$freq)
  
  library(plyr)
  library(ggplot2)
  ce <- ddply(data,"state",transform,percent_state = freq/sum(freq) * 100)
  ggplot(ce,aes(x = state,y = percent_state,fill=x)) +
    geom_bar(stat = "identity",width = 0.9) +
    scale_x_discrete(limits=c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c(enhancer = enhancer_color,promoter = promoter_color, repressed = repressed_color,
                                 heterochromatin = heterochromatin_color,quiescent = quit_color),
                      breaks = c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    theme_classic(base_size = 18,base_family = "sans") +
    labs(title = "7weeks",x=NULL,y="Percentage (%)",fill=NULL) +
    theme(plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,face = "bold"),
          axis.title = element_text(colour = "black",size = 18),
          axis.text.y = element_text(size = 18,colour = "black"),
          axis.text.x = element_text(size = 18,colour = "black",angle = 45,hjust = 1),
          axis.ticks.length = unit(2,units = "mm")) +
    #guides(fill=FALSE)
    theme(legend.position = "right")
  ggsave(filename = paste(geneCountDir,region,"/",region,"_0week state fill with 7weeks state.pdf",sep = ""),height = 4,width = 5)
  
  # 10weeks fill with 0week state --------------------------------------------
  enhancer <- df_sub[df_sub$week0 == "enhancer",]
  enhancer <- count(enhancer$week10)
  enhancer$state <- "enhancer"
  promoter <- df_sub[df_sub$week0 == "promoter",]
  promoter <- count(promoter$week10)
  #promoter <- rbind(promoter[1,],c("heterochromatin",0),promoter[2:4,])
  promoter$state <- "promoter"
  repressed <- df_sub[df_sub$week0 == "repressed",]
  repressed <- count(repressed$week10)
  repressed$state <- "repressed"
  heterochromatin <- df_sub[df_sub$week0 == "heterochromatin",]
  heterochromatin <- count(heterochromatin$week10)
  #heterochromatin <- rbind(heterochromatin[1:2,],c("promoter",0),heterochromatin[3:4,])
  heterochromatin$state <- "heterochromatin"
  quiescent <- df_sub[df_sub$week0 == "quiescent",]
  quiescent <- count(quiescent$week10)
  quiescent$state <- "quiescent"
  
  data <- rbind(enhancer,promoter,repressed,heterochromatin,quiescent)
  data$x <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                   levels = c("quiescent","heterochromatin","repressed","enhancer","promoter"))
  data$freq <- as.numeric(data$freq)
  
  library(plyr)
  library(ggplot2)
  ce <- ddply(data,"state",transform,percent_state = freq/sum(freq) * 100)
  ggplot(ce,aes(x = state,y = percent_state,fill=x)) +
    geom_bar(stat = "identity",width = 0.9) +
    scale_x_discrete(limits=c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c(enhancer = enhancer_color,promoter = promoter_color, repressed = repressed_color,
                                 heterochromatin = heterochromatin_color,quiescent = quit_color),
                      breaks = c("quiescent","heterochromatin","repressed","enhancer","promoter")) +
    theme_classic(base_size = 18,base_family = "sans") +
    labs(title = "10weeks",x=NULL,y="Percentage (%)",fill=NULL) +
    theme(plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,face = "bold"),
          axis.title = element_text(colour = "black",size = 18),
          axis.text.y = element_text(size = 18,colour = "black"),
          axis.text.x = element_text(size = 18,colour = "black",angle = 45,hjust = 1),
          axis.ticks.length = unit(2,units = "mm")) +
    #guides(fill=FALSE)
    theme(legend.position = "right")
  ggsave(filename = paste(geneCountDir,region,"/",region,"_0week state fill with 10weeks state.pdf",sep = ""),height = 4,width = 5)
}

plotStack(filename = "~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody0bp/mm10_genebody_0bp_merge_rna.csv")
plotStack(filename = "~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody2000bp/mm10_genebody_2000bp_merge_rna.csv")
plotStack(filename = "~/project/colon cancer/chip-seq/chromHMM/geneCount/genetss4000bp/mm10_genetss_4000bp_merge_rna.csv")

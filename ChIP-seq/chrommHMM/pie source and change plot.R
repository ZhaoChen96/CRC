rm(list = ls())
setwd("~/project/colon cancer/chip-seq/chromHMM/geneCount/")
library(ggplot2)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
geneCountDir <- "~/project/colon cancer/chip-seq/chromHMM/geneCount/"

plotPies(filename = "~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody0bp/mm10_genebody_0bp_merge_rna.csv")
plotPies(filename = "~/project/colon cancer/chip-seq/chromHMM/geneCount/genebody2000bp/mm10_genebody_2000bp_merge_rna.csv")

plotPies <- function(filename) {
  region <- paste(unlist(strsplit(basename(filename),split = "_"))[2],unlist(strsplit(basename(filename),split = "_"))[3],sep = "")
  df <- read.csv(filename)
  df <- df[,c("gene_name","week0","week2","week4","week7","week10")]
  df <- na.omit(df)
  
  #whole genome------------
  state = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13")
  element = c("enhancer", "enhancer", "enhancer", "promoter", "promoter", "promoter", "promoter", "promoter", "repressed", "repressed",
              "enhancer", "quiescent", "heterochromatin")
  dict = setNames(element, state)
  df$week0 = unname(dict[df$week0])
  df$week2 = unname(dict[df$week2])
  df$week4 = unname(dict[df$week4])
  df$week7 = unname(dict[df$week7])
  df$week10 = unname(dict[df$week10])
  
  heterochromatin_color = brewer.pal(n = 11, name = "PRGn")[3]
  promoter_color = brewer.pal(n = 9, name = "Greens")[8]
  enhancer_color = brewer.pal(n = 9, name = "Reds")[8]
  repressed_color = brewer.pal(n = 9, name = "Blues")[7]
  quit_color = brewer.pal(n = 9, name = "Greys")[3]
  
  quit_transform = df[(df$week10 == "enhancer") & (df$week0 == "quitscent"),]
  #write.csv(quit_transform, paste(geneCountDir,region,"/",region,"_quitscent_to_enhancer_gene.csv",sep = ""),quote = FALSE)
  
  prom_transform_enh = df[(df$week10 == "enhancer") & (df$week0 == "promoter"),]
  #write.csv(prom_transform_enh, paste(geneCountDir,region,"/",region,"_promoter_to_enhancer_gene.csv",sep = ""),quote = FALSE)
  
  # 0-10weeks enhancer total source  ----------------------------------------
  # 2th
  df_sub_old_enhancer_source = df[(df$week10 == "enhancer"),]
  df_sub_old_enhancer_source = df_sub_old_enhancer_source %>% count(week0)
  colnames(df_sub_old_enhancer_source) = c("group", "count")
  # # 1th
  # df_sub_enhancer_source= df_sub[(df_sub$week10 == "enhancer") & (df_sub$week0 != "enhancer"),]
  # df_sub_enhancer_source = df_sub_enhancer_source %>% count(week0)
  # colnames(df_sub_enhancer_source) = c("group", "count")
  
  library(ggplot2)
  library(scales)
  label_value <- paste(round(df_sub_old_enhancer_source$count/sum(df_sub_old_enhancer_source$count) * 100, 1), '%', sep = '')
  label <- paste(df_sub_old_enhancer_source$group,label_value)
  label_value
  df_sub_old_enhancer_source$group <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                                             levels = c("enhancer","quiescent","promoter","repressed","heterochromatin"))
  ggplot(df_sub_old_enhancer_source, aes(x="", y=count, fill=group)) +
    geom_bar(width = 1, stat = "identity") + 
    geom_text(aes(y = c(12000,0,1800,5000,300),x=c(1.1,1.4,1.1,1.2,1.1),label = label_value),size=5) +
    coord_polar("y", start=0,direction = 1) +
    labs(title="10weeks enhancer from 0week states",x=NULL,y=NULL,fill=NULL) +
    scale_fill_manual(values=c(enhancer_color,quit_color,promoter_color,repressed_color,heterochromatin_color))+
    theme_minimal(base_size = 18,base_family = "sans") +
    theme(panel.grid=element_blank()) +
    theme(plot.title = element_text(hjust = -0.2,vjust = 0.5,size = 18),
          axis.text.x=element_blank()) 
  ggsave(filename = paste(geneCountDir,region,"/",region,"_10weeks enhancer from 0week states pie chart.pdf",sep = ""),height = 4,width = 6)
  
  # 0week to 10weeks  states change into  heterochromatin -------------------
  heterochromatin_change <- df[df$week10 == "heterochromatin",]
  heterochromatin_change <- heterochromatin_change %>% count(week0)
  colnames(heterochromatin_change) <- c("group","count")
  heterochromatin_change$group <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                                         levels = c("quiescent","heterochromatin","repressed","enhancer","promoter"))
  label_value <- paste(round(heterochromatin_change$count/sum(heterochromatin_change$count) * 100, 1), '%', sep = '')
  label <- paste(heterochromatin_change$group,label_value)
  
  ggplot(heterochromatin_change, aes(x = "",y = count,fill=group)) +
    geom_bar(stat = "identity",width = 1) +
    geom_text(aes(y = c(8,500,19,1400,44),x=c(1.3,1.1,1.4,1.1,1.2),label = label_value),size=5) +
    theme_minimal(base_size = 18,base_family = "sans") +
    labs(title="10weeks heterochromatin from 0week states",x=NULL,y=NULL,fill=NULL) + 
    coord_polar(theta = 'y', start = 0,direction = 1) +
    scale_fill_manual(values = c(quit_color,heterochromatin_color,repressed_color,enhancer_color,promoter_color)) +
    theme(panel.grid=element_blank()) +
    theme(plot.title = element_text(hjust = -0.2,vjust = 0.5,size = 18),
          axis.text.x=element_blank()) +
    theme(legend.position = "right")
  ggsave(filename = paste(geneCountDir,region,"/",region,"_10weeks heterochromatin from 0week states.pdf",sep = ""),height = 4,width = 6)
  
  #quiescent
  df_sub_quiescent_source= df[(df$week0 == "quiescent"),]
  df_sub_quiescent_source = df_sub_quiescent_source %>% count(week10)
  colnames(df_sub_quiescent_source) = c("group", "count")
  
  label_value <- paste(round(df_sub_quiescent_source$count/sum(df_sub_quiescent_source$count) * 100, 1), '%', sep = '')
  label <- paste(df_sub_quiescent_source$group,label_value)
  label
  df_sub_quiescent_source$group <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                                          levels = c("quiescent","enhancer","heterochromatin","repressed","promoter"))
  ggplot(df_sub_quiescent_source, aes(x="", y=count, fill=group)) +
    geom_bar(width = 1, stat = "identity") + 
    geom_text(aes(y = c(3500,1200,9,14000,383),x=c(1.2,1.1,1.4,1.1,1.3),label = label_value),size=5) +
    coord_polar("y", start=0,direction = 1) +
    labs(title="0week quiescent to 10weeks states",x=NULL,y=NULL,fill=NULL) +
    scale_fill_manual(values=c(quit_color,enhancer_color,heterochromatin_color,repressed_color,promoter_color))+
    theme_minimal(base_size = 18,base_family = "sans") +
    theme(panel.grid=element_blank()) +
    theme(plot.title = element_text(hjust = -0.2,vjust = 0.5,size = 18),
          axis.text.x=element_blank()) 
  ggsave(filename = paste(geneCountDir,region,"/",region,"_0week quiescent to 10weeks states pie chart.pdf",sep = ""),height = 4,width = 6)
  
  #repressed
  df_sub_repressed_source= df[(df$week0 == "repressed"),]
  df_sub_repressed_source = df_sub_repressed_source %>% count(week10)
  colnames(df_sub_repressed_source) = c("group", "count")
  
  label_value <- paste(round(df_sub_repressed_source$count/sum(df_sub_repressed_source$count) * 100, 1), '%', sep = '')
  label <- paste(df_sub_repressed_source$group,label_value)
  label_value
  df_sub_repressed_source$group <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                                             levels = c("repressed","quiescent","enhancer","heterochromatin","promoter"))
  ggplot(df_sub_repressed_source, aes(x="", y=count, fill=group)) +
    geom_bar(width = 1, stat = "identity") + 
    geom_text(aes(y = c(400,21,7,1800,6000),x=c(1.2,1.3,1.4,1.1,1),label = label_value),size=5) +
    coord_polar("y", start=0,direction = 1) +
    labs(title="0week repressed to 10weeks states",x=NULL,y=NULL,fill=NULL) +
    scale_fill_manual(values=c(repressed_color,quit_color,enhancer_color,heterochromatin_color,promoter_color))+
    theme_minimal(base_size = 18,base_family = "sans") +
    theme(panel.grid=element_blank()) +
    theme(plot.title = element_text(hjust = -0.2,vjust = 0.5,size = 18),
          axis.text.x=element_blank()) 
  ggsave(filename = paste(geneCountDir,region,"/",region,"_0week repressed to 10weeks states pie chart.pdf",sep = ""),height = 4,width = 6)
  
  #promoter
  df_sub_promoter_source= df[(df$week0 == "promoter"),]
  df_sub_promoter_source = df_sub_promoter_source %>% count(week10)
  colnames(df_sub_promoter_source) = c("group", "count")
  
  label_value <- paste(round(df_sub_promoter_source$count/sum(df_sub_promoter_source$count) * 100, 1), '%', sep = '')
  label <- paste(df_sub_promoter_source$group,label_value)
  label
  df_sub_promoter_source$group <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                                          levels = c("enhancer","promoter","repressed","quiescent","heterochromatin"))
  ggplot(df_sub_promoter_source, aes(x="", y=count, fill=group)) +
    geom_bar(width = 1, stat = "identity") + 
    geom_text(aes(y = c(3800,8,1200,100,250),x=c(1,1.4,1.1,1.3,1.2),label = label_value),size=5) +
    coord_polar("y", start=0,direction = 1) +
    labs(title="0week promoter to 10weeks states",x=NULL,y=NULL,fill=NULL) +
    scale_fill_manual(values=c(enhancer_color,promoter_color,repressed_color,quit_color,heterochromatin_color))+
    theme_minimal(base_size = 18,base_family = "sans") +
    theme(panel.grid=element_blank()) +
    theme(plot.title = element_text(hjust = -0.2,vjust = 0.5,size = 18),
          axis.text.x=element_blank()) 
  ggsave(filename = paste(geneCountDir,region,"/",region,"_0week promoter to 10weeks states pie chart.pdf",sep = ""),height = 4,width = 6)
}






get_enhancer_source(filename = "~/project/colon cancer/chip-seq/chromHMM/geneCount/genetss4000bp/mm10_genetss_4000bp_merge_rna.csv")
get_enhancer_source <- function(filename) {
  region <- paste(unlist(strsplit(basename(filename),split = "_"))[2],unlist(strsplit(basename(filename),split = "_"))[3],sep = "")
  df <- read.csv(filename)
  df <- df[,c("gene_name","week0","week2","week4","week7","week10")]
  df <- na.omit(df)
  
  #whole genome------------
  state = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13")
  element = c("enhancer", "enhancer", "enhancer", "promoter", "promoter", "promoter", "promoter", "promoter", "repressed", "repressed",
              "enhancer", "quiescent", "heterochromatin")
  dict = setNames(element, state)
  df$week0 = unname(dict[df$week0])
  df$week2 = unname(dict[df$week2])
  df$week4 = unname(dict[df$week4])
  df$week7 = unname(dict[df$week7])
  df$week10 = unname(dict[df$week10])
  
  heterochromatin_color = brewer.pal(n = 11, name = "PRGn")[3]
  promoter_color = brewer.pal(n = 9, name = "Greens")[8]
  enhancer_color = brewer.pal(n = 9, name = "Reds")[8]
  repressed_color = brewer.pal(n = 9, name = "Blues")[7]
  quit_color = brewer.pal(n = 9, name = "Greys")[3]
  
  quit_transform = df[(df$week10 == "enhancer") & (df$week0 == "quitscent"),]
  #write.csv(quit_transform, paste(geneCountDir,region,"/",region,"_quitscent_to_enhancer_gene.csv",sep = ""),quote = FALSE)
  
  prom_transform_enh = df[(df$week10 == "enhancer") & (df$week0 == "promoter"),]
  #write.csv(prom_transform_enh, paste(geneCountDir,region,"/",region,"_promoter_to_enhancer_gene.csv",sep = ""),quote = FALSE)
  
  # 0-10weeks enhancer total source  ----------------------------------------
  # 2th
  df_sub_old_enhancer_source = df[(df$week10 == "enhancer"),]
  df_sub_old_enhancer_source = df_sub_old_enhancer_source %>% count(week0)
  colnames(df_sub_old_enhancer_source) = c("group", "count")
  # # 1th
  # df_sub_enhancer_source= df_sub[(df_sub$week10 == "enhancer") & (df_sub$week0 != "enhancer"),]
  # df_sub_enhancer_source = df_sub_enhancer_source %>% count(week0)
  # colnames(df_sub_enhancer_source) = c("group", "count")
  
  library(ggplot2)
  library(scales)
  label_value <- paste(round(df_sub_old_enhancer_source$count/sum(df_sub_old_enhancer_source$count) * 100, 1), '%', sep = '')
  label <- paste(df_sub_old_enhancer_source$group,label_value)
  label
  df_sub_old_enhancer_source$group <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                                             levels = c("enhancer","promoter","quiescent","repressed","heterochromatin"))
  ggplot(df_sub_old_enhancer_source, aes(x="", y=count, fill=group)) +
    geom_bar(width = 1, stat = "identity") + 
    geom_text(aes(y = c(15000,11,7500,2500,650),x=c(1.1,1.4,1.1,1.2,1.1),label = label_value),size=5) +
    coord_polar("y", start=0,direction = 1) +
    labs(title="10weeks enhancer from 0week states",x=NULL,y=NULL,fill=NULL) +
    scale_fill_manual(values=c(enhancer_color,promoter_color,quit_color,repressed_color,heterochromatin_color))+
    theme_minimal(base_size = 18,base_family = "sans") +
    theme(panel.grid=element_blank()) +
    theme(plot.title = element_text(hjust = -0.2,vjust = 0.5,size = 18),
          axis.text.x=element_blank()) 
  ggsave(filename = paste(geneCountDir,region,"/",region,"_10weeks enhancer from 0week states pie chart.pdf",sep = ""),height = 4,width = 6)
  
  # 0week to 10weeks  states change into  heterochromatin -------------------
  heterochromatin_change <- df[df$week10 == "heterochromatin",]
  heterochromatin_change <- heterochromatin_change %>% count(week0)
  colnames(heterochromatin_change) <- c("group","count")
  heterochromatin_change$group <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                                         levels = c("quiescent","heterochromatin","repressed","promoter","enhancer"))
  label_value <- paste(round(heterochromatin_change$count/sum(heterochromatin_change$count) * 100, 1), '%', sep = '')
  label <- paste(heterochromatin_change$group,label_value)
  label
  
  ggplot(heterochromatin_change, aes(x = "",y = count,fill=group)) +
    geom_bar(stat = "identity",width = 1) +
    geom_text(aes(y = c(6,400,40,1200,80),x=c(1.4,1.1,1.3,1,1.2),label = label_value),size=5) +
    theme_minimal(base_size = 18,base_family = "sans") +
    labs(title="10weeks heterochromatin from 0week states",x=NULL,y=NULL,fill=NULL) + 
    coord_polar(theta = 'y', start = 0,direction = 1) +
    scale_fill_manual(values = c(quit_color,heterochromatin_color,repressed_color,enhancer_color,promoter_color)) +
    theme(panel.grid=element_blank()) +
    theme(plot.title = element_text(hjust = -0.2,vjust = 0.5,size = 18),
          axis.text.x=element_blank()) +
    theme(legend.position = "right")
  ggsave(filename = paste(geneCountDir,region,"/",region,"_10weeks heterochromatin from 0week states.pdf",sep = ""),height = 4,width = 6)
  
  #quiescent
  df_sub_quiescent_source= df[(df$week0 == "quiescent"),]
  df_sub_quiescent_source = df_sub_quiescent_source %>% count(week10)
  colnames(df_sub_quiescent_source) = c("group", "count")
  
  label_value <- paste(round(df_sub_quiescent_source$count/sum(df_sub_quiescent_source$count) * 100, 1), '%', sep = '')
  label <- paste(df_sub_quiescent_source$group,label_value)
  label
  df_sub_quiescent_source$group <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                                          levels = c("quiescent","enhancer","heterochromatin","repressed","promoter"))
  ggplot(df_sub_quiescent_source, aes(x="", y=count, fill=group)) +
    geom_bar(width = 1, stat = "identity") + 
    geom_text(aes(y = c(3300,1200,9,14000,383),x=c(1.1,1.2,1.4,1,1.3),label = label_value),size=5) +
    coord_polar("y", start=0,direction = 1) +
    labs(title="0week quiescent to 10weeks states",x=NULL,y=NULL,fill=NULL) +
    scale_fill_manual(values=c(quit_color,enhancer_color,heterochromatin_color,repressed_color,promoter_color))+
    theme_minimal(base_size = 18,base_family = "sans") +
    theme(panel.grid=element_blank()) +
    theme(plot.title = element_text(hjust = -0.2,vjust = 0.5,size = 18),
          axis.text.x=element_blank()) 
  ggsave(filename = paste(geneCountDir,region,"/",region,"_0week quiescent to 10weeks states pie chart.pdf",sep = ""),height = 4,width = 6)
  
  #repressed
  df_sub_repressed_source= df[(df$week0 == "repressed"),]
  df_sub_repressed_source = df_sub_repressed_source %>% count(week10)
  colnames(df_sub_repressed_source) = c("group", "count")
  
  label_value <- paste(round(df_sub_repressed_source$count/sum(df_sub_repressed_source$count) * 100, 1), '%', sep = '')
  label <- paste(df_sub_repressed_source$group,label_value)
  label
  df_sub_repressed_source$group <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                                          levels = c("repressed","quiescent","enhancer","heterochromatin","promoter"))
  ggplot(df_sub_repressed_source, aes(x="", y=count, fill=group)) +
    geom_bar(width = 1, stat = "identity") + 
    geom_text(aes(y = c(700,21,65,2200,7000),x=c(1.2,1.4,1.3,1.1,1),label = label_value),size=5) +
    coord_polar("y", start=0,direction = 1) +
    labs(title="0week repressed to 10weeks states",x=NULL,y=NULL,fill=NULL) +
    scale_fill_manual(values=c(repressed_color,quit_color,enhancer_color,heterochromatin_color,promoter_color))+
    theme_minimal(base_size = 18,base_family = "sans") +
    theme(panel.grid=element_blank()) +
    theme(plot.title = element_text(hjust = -0.2,vjust = 0.5,size = 18),
          axis.text.x=element_blank()) 
  ggsave(filename = paste(geneCountDir,region,"/",region,"_0week repressed to 10weeks states pie chart.pdf",sep = ""),height = 4,width = 6)
  
  #promoter
  df_sub_promoter_source= df[(df$week0 == "promoter"),]
  df_sub_promoter_source = df_sub_promoter_source %>% count(week10)
  colnames(df_sub_promoter_source) = c("group", "count")
  
  label_value <- paste(round(df_sub_promoter_source$count/sum(df_sub_promoter_source$count) * 100, 1), '%', sep = '')
  label <- paste(df_sub_promoter_source$group,label_value)
  label
  df_sub_promoter_source$group <- factor(c("enhancer","heterochromatin","promoter","quiescent","repressed"),
                                         levels = c("enhancer","promoter","repressed","quiescent","heterochromatin"))
  ggplot(df_sub_promoter_source, aes(x="", y=count, fill=group)) +
    geom_bar(width = 1, stat = "identity") + 
    geom_text(aes(y = c(6000,8,2500,300,1000),x=c(1,1.4,1.1,1.3,1.2),label = label_value),size=5) + 
    coord_polar("y", start=0,direction = 1) +
    labs(title="0week promoter to 10weeks states",x=NULL,y=NULL,fill=NULL) +
    scale_fill_manual(values=c(enhancer_color,promoter_color,repressed_color,quit_color,heterochromatin_color))+
    theme_minimal(base_size = 18,base_family = "sans") +
    theme(panel.grid=element_blank()) +
    theme(plot.title = element_text(hjust = -0.2,vjust = 0.5,size = 18),
          axis.text.x=element_blank()) 
  ggsave(filename = paste(geneCountDir,region,"/",region,"_0week promoter to 10weeks states pie chart.pdf",sep = ""),height = 4,width = 6)
}




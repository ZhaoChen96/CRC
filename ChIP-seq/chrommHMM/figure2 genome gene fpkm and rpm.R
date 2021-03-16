rm(list = ls())
setwd("~/project/colon cancer/chip-seq/chromHMM/geneCount/")
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(circlize)


#df <- read.delim("mm10_genebody_0bp_merge_rna.csv",sep = ",",check.names = FALSE)
#df <- read.delim("mm10_genebody_2000bp_merge_rna.csv",sep = ",",check.names = FALSE)
df <- read.delim("mm10_genetss_4000bp_merge_rna.csv",sep = ",",check.names = FALSE)
colnames(df)

heterochromatin_color = brewer.pal(n = 11, name = "PRGn")[3]
promoter_color = brewer.pal(n = 9, name = "Greens")[8]
enhancer_color = brewer.pal(n = 9, name = "Reds")[8]
repressed_color = brewer.pal(n = 9, name = "Blues")[7]
quit_color = brewer.pal(n = 9, name = "Greys")[4]
colour_list = c(quit_color,heterochromatin_color,rep(repressed_color,2),rep(enhancer_color,4),rep(promoter_color,5))

plotstate <- function(time,xdata) {
  p <- ggplot(data, aes(y=state, x=xdata,fill=state)) +
    stat_boxplot(geom = "errorbar",linetype=1,width=0.8,position = "identity") +
    geom_boxplot(outlier.fill = NA,outlier.shape=NA,width=0.8) +
    theme_bw(base_family = "sans",base_size = 18) +
    labs(title = time,x="Gene Expression log2(FPKM+1)",y=NULL) +
    scale_x_continuous(limits = c(0,13),breaks = c(0,3,6,9,12)) +
    scale_y_discrete(limits=rev(levels(data$state))) +
    scale_fill_manual(values = colour_list) +
    theme(#aspect.ratio = 2/1,
      plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 17,colour = "black"),
      axis.ticks.y = element_blank(),
      axis.ticks.length.x = unit(2,"mm"),
      panel.grid.major = element_line(size = 0.3),
      panel.grid.minor = element_line(size = 0.3)) +
    guides(fill=FALSE)
  ggsave(paste("genetss_4000bp_hist_state_expression_",time,".pdf",sep = ""), p,height = 6.9,width = 3.8)
  p
}

# week0
df_sub <- df[,c("week0","control_0","control_1","control_2")]
data <- melt(df_sub,id.vars = "week0",variable.name = "rep",value.name = "fpkm")
colnames(data) <- c("state","rep","fpkm")
data$fpkm <- log2(data$fpkm + 1)
data$state <- factor(data$state,levels = c("E12","E13","E9","E10","E1","E11","E3","E2","E6","E7","E4","E5","E8"))
plotstate(time = "0week",xdata = data$fpkm)
# week2
df_sub <- df[,c("week2","2week_1","2week_0","2week_2")]
data <- melt(df_sub,id.vars = "week2",variable.name = "rep",value.name = "fpkm")
colnames(data) <- c("state","rep","fpkm")
data$fpkm <- log2(data$fpkm + 1)
data$state <- factor(data$state,levels = c("E12","E13","E9","E10","E1","E11","E3","E2","E6","E7","E4","E5","E8"))
plotstate(time = "2weeks",xdata = data$fpkm)
# week4
df_sub <- df[,c("week4","4week_0","4week_1","4week_2")]
data <- melt(df_sub,id.vars = "week4",variable.name = "rep",value.name = "fpkm")
colnames(data) <- c("state","rep","fpkm")
data$fpkm <- log2(data$fpkm + 1)
data$state <- factor(data$state,levels = c("E12","E13","E9","E10","E1","E11","E3","E2","E6","E7","E4","E5","E8"))
plotstate(time = "4weeks",xdata = data$fpkm)
# week7
df_sub <- df[,c("week7","7week_0","7week_1","7week_2")]
data <- melt(df_sub,id.vars = "week7",variable.name = "rep",value.name = "fpkm")
colnames(data) <- c("state","rep","fpkm")
data$fpkm <- log2(data$fpkm + 1)
data$state <- factor(data$state,levels = c("E12","E13","E9","E10","E1","E11","E3","E2","E6","E7","E4","E5","E8"))
plotstate(time = "7weeks",xdata = data$fpkm)
# week10
df_sub <- df[,c("week10","10week_0","10week_1","10week_2")]
data <- melt(df_sub,id.vars = "week10",variable.name = "rep",value.name = "fpkm")
colnames(data) <- c("state","rep","fpkm")
data$fpkm <- log2(data$fpkm + 1)
data$state <- factor(data$state,levels = c("E12","E13","E9","E10","E1","E11","E3","E2","E6","E7","E4","E5","E8"))
plotstate(time = "10weeks",xdata = data$fpkm)



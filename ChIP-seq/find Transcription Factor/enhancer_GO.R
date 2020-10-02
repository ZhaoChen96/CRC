rm(list = ls())
setwd("~/project/colon cancer/chip-seq/peak_length_distribution/")
library(ggplot2)
library(stringr)


# bed read number ---------------------------------------------------------
data <- read.csv('number.txt',sep = "",header = FALSE)
names(data) <- c("read_number","sample")
data$sample <- str_split_fixed(data$sample,"_",2)[,1]
data$group <- str_split_fixed(data$sample,"-",2)[,2]
data$read_number <- data$read_number/1000000
data <- data[c(-36,-5,-12,-19,-26,-33),]
times <- c("ctrl","2weeks","4weeks","7weeks","10weeks")
new_times <- c("0week","2weeks","4weeks","7weeks","10weeks")
markers <- c("H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me3","Input")
sample_vector <- c()
sample_list <- c()
for (marker in markers){
  for (time in times){
    labels <- paste(time,marker,sep = "-")
    sample_vector <- c(sample_vector,labels)
  }
}
for (marker in markers){
  for (time in new_times){
    labels <- paste(time,marker,sep = "-")
    sample_list <- c(sample_list,labels)
  }
}

data$sample <- factor(data$sample,levels = sample_vector)
data$group <- factor(data$group,levels = markers)

ggplot(data, aes(y = sample,x = read_number,colour=group)) + 
  geom_segment(aes(yend=sample),xend=0) +
  geom_point(size=3) +
  theme_classic(base_size = 16,base_family = "sans") +
  scale_x_continuous(expand = c(0,0),limits = c(0,150)) +
  scale_y_discrete(limits=rev(sample_vector),labels=rev(sample_list)) +
  scale_colour_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","grey50"),labels=factor(markers)) +
  labs(title = "ChIP-seq histone modifications total reads number",x="Number of reads (million)",y = NULL,colour=NULL) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90",size = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  theme(plot.margin = margin(t = 0.5,r = 0.5,b = 0.5,l = 0.5,unit = "cm"),
        aspect.ratio = 1/0.618,
        axis.line.y = element_blank(),
        axis.line.x = element_line(size=1.2),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 1.2),
        axis.ticks.length.x = unit(0.25,"cm"),
        plot.title = element_text(hjust = 0.8,vjust = 0.5),
        axis.text.y = element_text(colour=rep(rev(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","grey50")),each=5),size = 16),
        axis.text.x = element_text(colour="black",size = 16,margin = margin(t = 0.2,b = 0.2,unit = "cm")),
        axis.title.x = element_text(size = 20)) +
  theme(legend.position = "top",legend.text = element_text(size = 16)) +
  guides(colour = guide_legend(nrow = 2,byrow = TRUE))
ggsave(filename = "chip-seq reads number bedfile.pdf",height = 8,width = 8)

attach(data)
pdf('new reads number.pdf',height = 8,width = 8,"cm")
ggplot(data = data,aes(y = sample,x = read_number,colour=group)) +
  geom_segment(aes(yend=sample),xend=0,colour="grey40") +
  geom_point(size=3) +
  scale_y_discrete(limits=rev(sample_vector)) +
  scale_x_continuous(expand = c(0,0),limits = c(0,150)) +
  scale_colour_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","grey50")) +
  theme_bw(base_family = "sans",base_size = 16) +
  theme(aspect.ratio = 1/0.618,
        plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
    panel.grid.minor.x = element_line(colour = "grey90",size = 0.3),
    panel.grid.major.x = element_line(colour = "grey90",size = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(colour="black",size = 16,margin = margin(t = 0.2,b = 0.2,unit = "cm")),
  ) +
  theme(legend.position = "top") +
  guides(colour = guide_legend(nrow = 2,byrow = TRUE)) +
  labs(title = "ChIP-seq total reads number",x="Number of reads (million)",y=NULL,colour=NULL) 
dev.off()


# pre-data ----------------------------------------------------------------
data <- read.csv('peak_number.txt',header = FALSE,sep = '')
names(data) <- c("peak_number","sample")
library(stringr)
data$file <- str_split_fixed(data$sample,'/',9)[,9]
data$sample <- str_split_fixed(data$file,'_',2)[,1]
data$time <- str_split_fixed(data$file,'-',2)[,1]
data$group <- str_split_fixed(data$sample,"-",2)[,2]
data <- data[c(-1,-4,-11,-14,-20,),]

label_vector = c()
sample_vector <- c()
times <- c("ctrl","2weeks","4weeks","7weeks","10weeks")
markers = c("H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me3")
for (marker in markers){
  for (time in times){
    label <- paste(time,marker,sep = "-")
    label_vector <- c(label_vector,label)
  }
}
labels <- factor(x = label_vector)
for (marker in markers){
  for (time in new_times){
    label <- paste(time,marker,sep = "-")
    sample_vector <- c(sample_vector,label)
  }
}


# first plot --------------------------------------------------------------------
attach(data)
pdf('peak number.pdf',width = 8,height = 8)
ggplot(data = data,aes(y = factor(sample),x = peak_number,fill=group)) +
  geom_bar(stat = "identity",width = 0.7) +
  geom_text(aes(label=peak_number),vjust=0.3,colour="black",hjust=0)  +
  #guides(x=guide_legend(reverse = FALSE)) 
  scale_y_discrete(limits=rev(label_vector),labels=rev(sample_vector)) +
  scale_fill_brewer(palette = "Pastel1",limits= markers) +
  scale_x_continuous(expand = c(0,0),limits = c(0,50000)) +
  theme_bw(base_size = 16,base_family = "sans") +
  theme(legend.position = "top") +
  theme(aspect.ratio = 1/0.618,
        plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
    panel.grid.major = element_blank(),
    axis.text = element_text(colour = "black")
  ) +
  labs(title = "Macs2 peak calling reads number",x="Number of Macs2 peak calling",y = NULL,fill=NULL) +
  guides(fill=guide_legend(nrow = 2,byrow = TRUE))
  
dev.off()

# second plot -------------------------------------------------------------
ggplot(data,aes(x = peak_number,y = sample,fill=group)) +
  geom_bar(stat = "identity",width = 0.8) +
  scale_y_discrete(limits=rev(label_vector),labels=rev(sample_vector)) +
  scale_fill_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00"),limits=markers) +
  scale_x_continuous(expand = c(0,0),limits = c(0,40000),breaks = c(0,10000,20000,30000,40000)) +
  theme_classic(base_family = "sans",base_size = 16) +
  labs(title = "Macs2 peak calling reads number",x="Number of Macs2 peak calling",y = NULL,fill=NULL) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90",size = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  theme(plot.margin = margin(t = 0.5,r = 0.5,b = 0.5,l = 0.5,unit = "cm"),
        plot.title = element_text(size =18,hjust = 0.5,vjust = 0.5),
        aspect.ratio = 1/0.618,
        axis.line.y = element_blank(),
        axis.line.x = element_line(size=1.2),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 1.2),
        axis.ticks.length.x = unit(0.25,"cm"),
        axis.text.y = element_text(colour=rep(rev(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")),each=5),size = 16),
        axis.text.x = element_text(colour="black",size = 16,margin = margin(t = 0.2,b = 0.2,unit = "cm")),
        axis.title.x = element_text(size = 18)) +
  theme(legend.position = "top",legend.text = element_text(size = 12)) +
  guides(fill=guide_legend(nrow = 2,byrow = TRUE))
ggsave(filename = "chip-seq peak number.pdf",width = 8,height = 8)


# peak length  ------------------------------------------------------------
rm(list = ls())
setwd("~/colon cancer/chip-seq/peak_length_distribution/peakInfo/")
library(ggplot2)
library(stringr)
library(beanplot)
pdf_name <- paste(getwd(),"peak_length.pdf",sep = "/")
pdf(pdf_name,width = 20,height = 15)
par(mfrow = c(2,3))
filelist <- list.files(path = getwd(),pattern = "*_peak_length.sort.txt$")
for (file in filelist){
  data <- read.csv(file = file,sep = "\t",header = FALSE)
  names(data) <- c("peak_length","sample")
  data$time <- str_split_fixed(data$sample,"-",2)[,1]
#  data$sample <- factor(data$sample,levels = label_vector)
  attach(data)
  model <- str_split_fixed(file,"_",3)[,1]
  beanplot(peak_length ~ time,ll=0.003,border = NA,ylim = ylim,col = c("black", "white"),
           overallline = "median",xlab=NULL,ylab="Peak length (Kb)",main=model)
}
dev.off()

a <- data[which(data$peak_length > 10000),]
attach(data)
pdf('H3K27ac_violin.pdf',height = 6.18,width = 6.18)
ggplot(data = data,aes(x=factor(sample),y = peak_length),color=sample) +
  geom_violin(
    trim = FALSE, #trim 保留提琴尾部
    scale = 'count', #图面积和观测值数目成正比
  ) +
  ylim(0,5000) +
  geom_boxplot(width=.06,fill='grey80',outlier.colour = NA) +
  stat_summary(fun.y=median,geom = "point",fill='white',shape=21,size=2) +
  theme_bw() +
  ggtitle('The stats of sample peak length') +
  labs(x=NULL,y="Peak length(kb)") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank()
    ) +
  theme(
    axis.line.x.bottom = element_line(colour = "black",linetype = "solid"),
    axis.line.y = element_line(colour = "black",linetype = "solid"),
    plot.title = element_text(size = rel(1.3),lineheight = .9,family = 'Times',face = 'plain',hjust = 0.5,vjust = 0),
    axis.title = element_text(size = rel(1.2),lineheight = .9,family = 'Times',face = 'plain'),  
    axis.text.y = element_text(size = rel(1.2),lineheight = .9,family = 'Times',face = 'plain',colour = 'black'),
    axis.text.x = element_text(size = rel(1.2),lineheight = .9,family = 'Times',face = 'plain',colour = 'black',angle = 30,hjust = 1)
  ) +
  scale_x_discrete(limits=c("ctrl-H3K27ac","2weeks-H3K27ac","4weeks-H3K27ac","7weeks-H3K27ac","10weeks-H3K27ac"))
dev.off()

# peak distribution -------------------------------------------------------
rm(list = ls())
setwd("~/colon cancer/chip-seq/peak_length_distribution/peakInfo/")
library(ggplot2)
library(stringr)

filelist <- list.files(path = getwd(),pattern = "*_peak_distribution.txt$",full.names = TRUE,recursive = FALSE)
for (file in filelist){
  data <- read.csv(file = file,sep = "\t",header = FALSE)
  names(data) <- c("term","number","sample")
  outputDir <- "~/colon cancer/chip-seq/peak_length_distribution/peakInfo/"
  model <- str_split_fixed(file,"/",8)[,8]
  model <- str_split_fixed(model,"_",4)[,2]
  times <- c("ctrl","2weeks","4weeks","7weeks","10weeks")
  labels_vector <- c()
  for (time in times){
    labels <- paste(time,model,sep = "-")
    labels_vector <- c(labels_vector,labels)
  }
  data$term <- factor(data$term,levels = c("non-coding","TTS","3' UTR","5' UTR","exon",
                                           "promoter-TSS","Intergenic","intron"))
  pdf_name <- paste(outputDir,model,"_peak_distribution.pdf",sep = "")
  pdf(pdf_name,width = 8,height = 6)
  attach(data)
  p <- ggplot(data = data,aes(x = sample,y = number,fill=term)) +
    geom_bar(stat = "identity",width = 0.7) +
    scale_x_discrete(limits=rev(labels_vector)) +
    scale_fill_brewer(palette ="Set3",direction = -1) +
    coord_flip() +
    theme_bw() +
    ggtitle('The stats of sample peak distribution on genebody') +
    labs(x=NULL,y=NULL,fill=NULL) +
    guides(fill=guide_legend(reverse = TRUE)) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_blank()
    ) +
    theme(
      axis.line = element_line(colour = "black",linetype = "solid",size = rel(0.5)),
      legend.position = "bottom",
      plot.title = element_text(size = rel(1.3),lineheight = .9,family = 'Times',face = 'plain',hjust = 0.5,vjust = 0),
      axis.title = element_text(size = rel(1.2),lineheight = .9,family = 'Times',face = 'plain'),  
      axis.text.y = element_text(size = rel(1.2),lineheight = .9,family = 'Times',face = 'plain',colour = 'black'),
      axis.text.x = element_text(size = rel(1.2),lineheight = .9,family = 'Times',face = 'plain',colour = 'black')
    )
  print(p)
  dev.off()
}


# Super enhancer to gene --------------------------------------------------------
rm(list = ls())
setwd("~/colon cancer/enhancer/")
library(ggplot2)
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
file <- read.csv(file = '10weeks-H3K27ac_Input_peaks_AllEnhancers_GENE_TO_ENHANCER.txt',sep = '',header = T) 
data <- file[which(file$IS_SUPER=='1'),] #挑出superenhancer
gene_name <- as.character(data[,1])
library(stringr)
ego_BP <- enrichGO(gene_name,
                   OrgDb = org.Mm.eg.db,
                   keyType = 'SYMBOL',
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01)
dotplot(ego_BP,font.size=14)
plotGOgraph(ego_BP)
barplot(ego_BP)
barplot(ego_BP,showCategory = 10,width=0.5,title="The GO_BP enrichment analysis of all Superenhancers related genes in 10weeks H3K27ac")+ 
  scale_size(range=c(-1,1))+
  scale_x_discrete(labels=function(ego_BP) str_wrap(ego_BP,width = 30))
write.table(ego_BP@result,'10weeks-H3K27ac-Upenhancer.txt',sep = "\t",quote=FALSE,col.names=TRUE)



rm(list = ls())

library(devtools)
library(ggord)
library(ggbiplot)
library(ggplot2)
library(ggforce)

data("wine")
wine.pca <- prcomp(wine,scale. = TRUE)
ggbiplot(wine.pca,obs.scale = 1,var.scale = 1,groups = wine.class,ellipse = TRUE) +
  scale_color_discrete(name=" ") +
  theme_bw() +
  theme(legend.direction = "horizontal",legend.position = "top") 


# PCA analysis ------------------------------------------------------------
library(RColorBrewer)
vsd <- getVarianceStabilizedData(dds)
pr <- prcomp(t(vsd))

# first
plot(pr$x,col="white",main="PC plot",xlim=c(-80,80),tck=-0.01,las=2,mgp=c(1.5,0.5,0))
points(pr$x[,1],pr$x[,2],pch=19,type = "p",col=rep(brewer.pal(5,"Dark2"),each=3),
       cex=1.5)

# second
plot(pr$x,col="white",main="PC plot",xlim=c(-100,100),ylim=c(-60,60),tck=-0.01,mar=c(2,2.5,2,2.5),axes=FALSE,mgp=c(2,0.5,0))
axis(side=1,las=0,tck=-0.01,pos =-61,at = seq(-100,100,100),labels = seq(-100,100,100),mgp=c(1,0.5,0))
axis(side=2,las=2,tck=-0.01,pos =-102,at = seq(-60,60,60),labels = seq(-60,60,60),mgp=c(1,0.5,0))
points(pr$x[,1],pr$x[,2],pch=19,type = "p",col=rep(brewer.pal(5,"Dark2"),each=3),
       cex=2,alpha=0.5)
legend("topright",inset = 0.03,Times,pch = 19,col = brewer.pal(5,"Dark2"),cex = 0.8,adj = c(0,0.5),text.col = brewer.pal(5,"Dark2"))

# third
biplot(pr,cex=c(1,0.5),main="Biplot",
       col=c("black","grey"))

# forth 
# ellipse.prob=0.68 置信区间, var.axes draw arrows for the variable
mytheme <- theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(
    plot.margin = margin(t = 0.5,l = 0.5,b = 0.5,r = 0.5,unit = "cm"),
    plot.title = element_text(size = 16,lineheight = .9,family = 'Helvetica',face = 'bold',hjust = 0.5,vjust = 0.5,margin = margin(b=0.5,unit = "cm")),
    axis.title = element_text(size = 14,lineheight = .9,family = 'Helvetica',face = 'plain',hjust = 0.5,vjust = 0.5),  
    axis.text.x = element_text(size = 14,lineheight = .9,family = 'Helvetica',face = 'plain',color="black",hjust = 0.5,vjust = 0.5,
                             margin = margin(t=0.15,b=0.3,unit = "cm")),
    axis.text.y = element_text(size = 14,lineheight = .9,family = 'Helvetica',face = 'plain',color="black",hjust = 0.5,vjust = 0.5,
                               margin = margin(l=0.3,r=0.15,unit = "cm")),
    axis.ticks.length = unit(.15,"cm")
) 

pdf("~/project/colon cancer/RNA-seq/DESeq2_maSigPro/PCA plot.pdf",height = 6,width = 6)
ggbiplot(pr,groups = condition,ellipse = TRUE,ellipse.prob = 0.75,labels = NULL,
         var.axes = FALSE,alpha = 0.01) + 
  geom_point(aes(colour=condition),shape=19,size=3,alpha=0.8) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(expand = c(0,0),limits=c(-2.5,2),breaks = c(-2,0,2)) +
  scale_y_continuous(expand = c(0,0),limits=c(-2.5,2),breaks = c(-2,0,2)) + 
  #coord_fixed(ratio = 1/1) +
  mytheme +
  labs(title = "RNA-seq PC plot",x="PC1 (46.1%)",y="PC2 (16.4%)",colour=NULL) +
  theme(aspect.ratio = 1,legend.position = c(0.9,0.2),legend.direction = "vertical",legend.text.align = 1,
        legend.text = element_text(size = 12,lineheight = .7,family = 'Helvetica',face = 'plain',color="black",hjust = 0.5,vjust = 0.5))
dev.off()


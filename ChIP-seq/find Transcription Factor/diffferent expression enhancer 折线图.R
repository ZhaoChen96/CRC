# 折线图orange data set ------------------------------------------------------
opar <- par(no.readonly = TRUE)
par(mfrow=c(1,2))
t1 <- subset(Orange,Tree==1)
plot(t1$age,t1$circumference,xlab="Age (days)",ylab="Circumference (mm)",main="Orange Tree 1 Growth")

Orange$Tree <- as.numeric(Orange$Tree)
ntrees <- max(Orange$Tree)

xrange <- range(Orange$age)
yrange <- range(Orange$circumference)

plot(xrange,yrange,type="n",xlab="Age (days)",ylab="Circumference (mm)")
colors <- rainbow(ntrees)
linetype <- c(1:ntrees)
plotchar <- seq(18,18+ntrees,1)

for (i in 1:ntrees) {
  tree <- subset(Orange,Tree==i)
  lines(tree$age, tree$circumference,
        type="b",
        lwd=2,
        lty=linetype[i],
        col=colors[i],
        pch=plotchar[i]
        )
}

title("Tree Growth","example of line plot")
legend(xrange[1],yrange[2],
       1:ntrees,
       cex=0.8,
       col=colors,
       pch = plotchar,
       lty = linetype,
       title = "Tree")
par(opar)

# CRC mouse ChIP-seq H3K27ac find transcription factor DEGs 折线图 ------------------------------
time <- c(2,4,7,10)
times <- c("2weeks","4weeks","7weeks","10weeks")
Up <- c(417,493,1398,4641)
Down <- c(5008,1853,2497,1440)
data <- data.frame(time,Up,Down,stringsAsFactors = FALSE)

#ggplot2
library(reshape2)
df <- melt(data,id.vars = "time",variable.name = "condition",value.name = "number")

library(ggplot2)
ggplot(data = df,aes(x=time,y=number,color=condition)) + 
  geom_line() +
  geom_point() +
  ylim(0,max(df$number)) +
  scale_color_manual(values = c("skyblue","#bdd7e7")) +
  labs(x="Time (weeks)",y = "Different expression enhancer number") +
  #scale_color_brewer(palette = "Blues") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size = rel(0.6))) +
  theme(legend.position = c(0.8,0.2)) +
  guides(fill=guide_legend(title = NULL))
  
  








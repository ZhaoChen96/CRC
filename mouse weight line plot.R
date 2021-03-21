rm(list = ls())
file <- read.table("~/project/colon cancer/chip-seq/figure/mouse weight.txt",sep = "\t",header = TRUE,check.names = FALSE)

file <- file[,-c(1,29,30)]
file[is.na(file)] <- 0
rownames(file) <- file$Day
data <- data.frame(lapply(file[,-1],function(x) x/x[1]))
data <- cbind(Day = file$Day,data)
library(reshape2)
data <- melt(data = data,id.vars = "Day",variable.name = "sample",value.name = "weight")
data <- cbind(data,time = c(rep(c("week0","week2","week4"),each=92),rep("week7",each=138),rep("week10",each=184)))
#data[is.na(data)] <- 0
data$Day <- as.numeric(data$Day)
data$time <- factor(data$time,levels = c("week0","week2","week4","week7","week10"))
#data <- na.omit(data)

df <- data.frame(tapply(data$weight,list(data$Day,data$time),mean))
df$Day <- rownames(df)
df$Day <- as.numeric(df$Day)
df <- melt(df,id.vars = "Day",variable.name = "time",value.name = "weight")
df <- inner_join(df,a,by=c("Day","time","weight"))
df <- df[df$weight !=0,]


ggplot(data = df,aes(x = Day,y = weight,colour=time)) +
  geom_line(aes(group=time),size=1.1) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=weight-se, ymax=weight+se),position=position_dodge(0.05),size=0.8) +
  theme_classic(base_size = 18,base_family = "sans",base_line_size = 1.1) +
  scale_y_continuous(expand = c(0,0),limits = c(0.8,1.3)) +
  scale_x_continuous(expand = c(0.01,0),breaks = c(0,7,14,21,28,35,42,49,56,63,70)) +
  scale_colour_manual(values = c("week0" = "#1B9E77","week2"="#66A61E","week4"="#7570B3","week7"="#E7298A","week10"="#D95F02"),
                      breaks = c("week0","week2","week4","week7","week10"),
                      labels = c( "0week", "2weeks", "4weeks","7weeks", "10weeks")) +
  labs(y="Relative weight",colour=NULL,title = "Weight of AOM/DSS mouse") +
  theme(#aspect.ratio = 0.5/1,
    plot.title = element_text(colour = "black",size = 18,hjust = 0.5),
    axis.text = element_text(colour = "black",size = 18),
    axis.title = element_text(colour = "black",size = 18)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 18,family = "sans")) 
  guides(colour=guide_legend(nrow = 2,byrow = TRUE))
ggsave("~/project/colon cancer/chip-seq/figure/mouse weight legend right.pdf",height = 3,width = 7)

exp_se <- describe(data[,c(3:4)])

summary(df)

library(pastecs)
stat.desc(data) 

library(Rmisc)
a <- summarySE(data = data,measurevar = "weight",groupvars = c("Day","time"))

library(gcookbook)
ce <- subset(cabbage_exp,Cultivar == "c39")











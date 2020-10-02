rm(list = ls())
setwd("/data3/zhaochen/project/colon_cancer/colon_chip/DiffBind")
bamDir <- "/data3/zhaochen/project/colon_cancer/colon_chip/samtools/"
bams <- dir(bamDir)
bams <- bams[-grep(".bai",bams)]
bams <- bams[-grep("H3K9me2",bams)]
bams <- bams[grep(".bam",bams)]

SampleID <- sapply(strsplit(bams,split = "_"),"[",1)
names(SampleID) <- bams



ctrl <- which(substring(bams,1,1)=="c")
two <- which(substring(bams,1,1)==2)
four <- which(substring(bams,1,1)==4)
seven <- which(substring(bams,1,1)==7)
ten <- which(substring(bams,1,1)==1)

controls <- grep("Input",bams)

ControlID <- bams
names(ControlID) <- bams
ControlID[ctrl] <- rep(bams[intersect(ctrl,controls)],each=6)
ControlID[two] <- rep(bams[intersect(two,controls)],each=6)
ControlID[four] <- rep(bams[intersect(four,controls)],each=6)
ControlID[seven] <- rep(bams[intersect(seven,controls)],each=6)
ControlID[ten] <- rep(bams[intersect(ten,controls)],each=6)

bamReads <- paste(bamDir,bams,sep = "")
bamControl <- paste(bamDir,ControlID,sep = "")

Tissue <- bams
names(Tissue) <- bams

Tissue[c(ctrl)] <- "Normal"
Tissue[c(two,four)] <-"Inflammation"
Tissue[c(seven,ten)] <- "Tumor"
Markers <- c("H3K27ac","H3K27me3","H3K4me1","H3K4me3","H3K9me3","Input")
Factor <- rep(Markers,15)
Condition <- rep(c("10weeks","2weeks","4weeks","7weeks","control"),each=18)
Replicate <- rep(rep(1:3,each=6),5)
comptab <- cbind(SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl)
comptab <- comptab[-controls,]
comptab

macs2Dir <- "/data3/zhaochen/project/colon_cancer/colon_chip/macs2/"
Peaks <- dir(macs2Dir,recursive=TRUE)[grep(".broadPeak",dir(macs2Dir,recursive=TRUE))]
Peaks <- Peaks[-grep("H3K9me2",Peaks)]
Peaks <- paste(macs2Dir,Peaks,sep="")
Peaks <- c(rep(Peaks[1:5],3),rep(Peaks[6:10],3),rep(Peaks[11:15],3),rep(Peaks[16:20],3),rep(Peaks[21:25],3))
PeakCaller <- Peaks
PeakCaller[1:length(PeakCaller)] <- "bed"

comptab <- cbind(comptab,Peaks,PeakCaller)
comptab <- as.data.frame(comptab)


for (i in 1:ncol(comptab)) {
comptab[,i] <- as.character(comptab[,i])
}
rownames(comptab) <- NULL

comptab
write.table(comptab,file = "sampleSheet.csv",sep=",",quote = FALSE,row.names = FALSE)

library(DiffBind)
# read peakset
DiffDir <- "/data3/zhaochen/project/colon_cancer/colon_chip/DiffBind"
file <- "/data3/zhaochen/project/colon_cancer/colon_chip/DiffBind/sampleSheet.csv"
file <- dba(sampleSheet = file)
names(file)
file
# counting reads
data <- dba.count(file,bUseSummarizeOverlaps = TRUE)
dba.save(DBA = data,file = "allmarkers_counts",pre = "dba",dir = DiffDir,ext = "Rdata")
out <- as.data.frame(data$binding)
write.table(out,file = "allmarkers_counts_binding.txt",sep = "\t",quote = FALSE,row.names = FALSE)
# plot correlation heatmap
pdf("allCalledPeaksHeatmap.pdf")
dba.plotHeatmap(data)
dev.off()
# plotPCA
pdf('allPCA.pdf')
dba.plotPCA(data,attributes = DBA_FACTOR,label = NULL)
dev.off()

# establishing a contrast
data <- dba.contrast(data,categories = c(DBA_FACTOR,DBA_CONDITION))

# performing the different analysis
data <- dba.analyze(data,method = DBA_ALL_METHODS)
dba.show(data,bContrasts = TRUE)
a <- dba.show(data,bContrasts = TRUE)
write.table(a,file = "allmarkers-contrast.txt",sep = "\t",quote = FALSE,row.names = TRUE)
pdf('Heatmap.pdf')
dba.plotHeatmap(data,correlations = FALSE,scale='row',margin = 10,ColAttributes = DBA_FACTOR)
dev.off()

# retrieve differentially bound sites
for (i in 1:nrow(a)) {
  filenames <- paste('contrast',i,sep = "")
  data.DB <- dba.report(data,th = 0.05,bUsePval = TRUE,method = DBA_ALL_METHODS,contrast = i,file = filenames,
                        DataType = DBA_DATA_FRAME,initString = "allMehtod",dir=DiffDir)
}



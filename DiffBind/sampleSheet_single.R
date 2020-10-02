rm(list=ls())
library(DiffBind)
library(ggplot2)

Markers <- c("H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me3")
for (marker in Markers[3:5]) {
  setwd(paste("/data3/zhaochen/project/colon_cancer/colon_chip/DiffBind/diff_",marker,sep = ""))
  bamDir <- "/data3/zhaochen/project/colon_cancer/colon_chip/samtools/"
  bams <- dir(bamDir)
  bams <- bams[-grep(".bai",bams)]
  bams <- bams[grep(".bam",bams)]
  histone <- bams[grep(marker,bams)]
  histone

  SampleID <- sapply(strsplit(histone,split = "_"),"[",1)
  names(SampleID) <- histone
  ctrl <- which(substring(histone,1,1)=="c")
  two <- which(substring(histone,1,1)==2)
  four <- which(substring(histone,1,1)==4)
  seven <- which(substring(histone,1,1)==7)
  ten <- which(substring(histone,1,1)==1)

  controls <- bams[grep("Input",bams)]
  ControlID <- controls

  bamReads <- paste(bamDir,histone,sep = "")
  bamControl <- paste(bamDir,ControlID,sep = "")

  Tissue <- histone
  names(Tissue) <- histone

  Tissue[c(ctrl)] <- "Normal"
  Tissue[c(two,four)] <-"Inflammation"
  Tissue[c(seven,ten)] <- "Tumor"

  Factor <- rep(marker,15)
  Condition <- rep(c("10weeks","2weeks","4weeks","7weeks","control"),each=3)
  Replicate <- rep(1:3,5)
  Treatment <- rep("Full-Media",15)
  comptab <- cbind(SampleID,Tissue,Factor,Condition,Replicate,Treatment,bamReads,ControlID,bamControl)
  comptab

  macs2Dir <- "/data3/zhaochen/project/colon_cancer/colon_chip/macs2/"
  Peaks <- dir(macs2Dir,recursive=TRUE)[grep(".broadPeak",dir(macs2Dir,recursive=TRUE))]
  Peaks <- Peaks[grep(marker,Peaks)]
  Peaks <- rep(paste(macs2Dir,Peaks,sep=""),each=3)
  PeakCaller <- Peaks
  PeakCaller[1:length(PeakCaller)] <- "bed"

  comptab <- cbind(comptab,Peaks,PeakCaller)
  comptab <- as.data.frame(comptab)

  for (i in 1:ncol(comptab)) {
    comptab[,i] <- as.character(comptab[,i])
  }
  rownames(comptab) <- NULL

  comptab
  write.table(comptab,file = paste(marker,"_sampleSheet.csv",sep = ""),row.names=FALSE,sep = ",")


  markerDir = paste("/data3/zhaochen/project/colon_cancer/colon_chip/DiffBind/diff_",marker,"/",sep = "")
  # read peakset
  #file <- paste(markerDir,marker,"_sampleSheet.csv",sep = "")
  file <- dba(sampleSheet = comptab)
  names(file)
  file
  pdf("CalledPeaksHeatmap.pdf")
  dba.plotHeatmap(file)
  dev.off()
  # counting reads
  data <- dba.count(file,bUseSummarizeOverlaps = TRUE)
  dba.save(DBA = data,file = paste(marker,"_count",sep = ""),pre = "dba_",dir = markerDir,ext = "RData")
  out <- as.data.frame(data$binding)
  write.table(out,paste(marker,"_binding.txt",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE)
  # plot correlation heatmap
  pdf("CalledPeaksHeatmap.pdf")
  dba.plotHeatmap(data)
  dev.off()
  # plot PCA
  pdf("plotPCA.pdf")
  dba.plotPCA(data,attributes = DBA_CONDITION,label = NULL)
  dev.off()

  # establishing a contrast
  data <- dba.contrast(data,categories = DBA_CONDITION)
  #load("dba_H3K27ac_count.RData")
  #data <- dba.contrast(DBAobject,categories = DBA_CONDITION)

  # performing the different analysis
  data <- dba.analyze(data,method = DBA_ALL_METHODS)
  dba.show(data,bContrasts = TRUE)
  a <- dba.show(data,bContrasts = TRUE)
  write.table(a,file = "allMethods_contrast.txt",sep = "\t",quote = FALSE,row.names = TRUE)
  pdf('Heatmap.pdf')
  dba.plotHeatmap(data,correlations = FALSE,scale='row',margin = 10,ColAttributes = DBA_FACTOR)
  dev.off()

  for (i in 1:10) {
    pdf(paste(markerDir,'Venn_contrast',i,".pdf",sep = ""))
    dba.plotVenn(data,contrast = i,method = DBA_ALL_METHODS)
    dev.off()
  }

  pdf('Heatmap.pdf')
  dba.plotHeatmap(data,correlations = FALSE,scale='row',margin = 12,
                  ColAttributes = DBA_CONDITION)
  dev.off()

  # retrieve differentially bound sites
  for (i in 1:10) {
    filename <- paste('contrast',i,sep = "")
    # if do not specify method, default DESeq
    data.DB <- dba.report(data,th = 0.05,contrast = i,method = DBA_ALL_METHODS,bUsePval=TRUE)
    out <- as.data.frame(data.DB)
    filename_all <- paste("allMethod_contrast",i,".txt",sep = "")
    write.table(out,sep = "\t",file = filename_all,quote = FALSE,col.names = NA)
    # edgeR, filename need change
    # data_edgeR.DB <- dba.report(data,th = 0.05,contrast = i,method = DBA_EDGER,bUsePval = TRUE,file = filename,
    #                             pre = "dba_EDGER",DataType=DBA_DATA_FRAME,initString = "edgeR")
    # out_egder <- as.data.frame(data_edgeR.DB)
    # filename_edger <- paste('edgeR_contrast',i,".txt",sep = "")
    # write.table(out_egder,sep = "\t",file = filename_edger,quote = FALSE,col.names = NA)
    # # DESeq2, filename need change
    # data_deseq2.DB <- dba.report(data,th = 0.05,contrast = i,method = DBA_DESEQ2,bUsePval = TRUE,file = filename,
    #                             pre = "dba_DESeq2",DataType=DBA_DATA_FRAME,initString = "deSeq2")
    # out_deseq2 <- as.data.frame(data_deseq2.DB)
    # filename_deseq2 <- paste('deseq2_contrast',i,'.txt',sep = "")
    # write.table(out_deseq2,sep = "\t",file = filename_deseq2,quote = FALSE,col.names = NA)
  }

}







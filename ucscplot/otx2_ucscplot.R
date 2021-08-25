rm(list = ls())
library(karyoploteR)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BiocFileCache)
library(rtracklayer)
library(ggplot2)
library("org.Mm.eg.db")

hmm.file0 <- import("/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/bed/state/ctrl_13_dense_new.bed")
hmm.file2 <- import("/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/bed/state/2weeks_13_dense_new.bed")
hmm.file4 <- import("/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/bed/state/4weeks_13_dense_new.bed")
hmm.file7 <- import("/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/bed/state/7weeks_13_dense_new.bed")
hmm.file10 <- import("/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/bed/state/10weeks_13_dense_new.bed")

hmm0 <- toGRanges(hmm.file0)
hmm2 <- toGRanges(hmm.file2)
hmm4 <- toGRanges(hmm.file4)
hmm7 <- toGRanges(hmm.file7)
hmm10 <- toGRanges(hmm.file10)

#Otx2 <- "/data3/zhaochen/project/colon_cancer/colon_chip/peakUCSCplot/otx2/Otx2-1-75821.bw"
Otx2 <- "/data3/zhaochen/project/colon_cancer/colon_chip/peakUCSCplot/otx2/Otx2-2-48516.bw"
#Otx2 <- "/data3/zhaochen/project/colon_cancer/colon_chip/peakUCSCplot/otx2/Otx2-3-48514.bw"

H3K27ac <- c(week10 <- "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input/10weeks-H3K27ac_bs1bp_ratio.bw",
             week7 <- "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input/7weeks-H3K27ac_bs1bp_ratio.bw",
             week4 <- "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input/4weeks-H3K27ac_bs1bp_ratio.bw",
             week2 <- "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input/2weeks-H3K27ac_bs1bp_ratio.bw",
             week0 <- "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input/ctrl-H3K27ac_bs1bp_ratio.bw")

H3K4me1 <- c(week10 <- "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input/10weeks-H3K4me1_bs1bp_ratio.bw",
             week7 <- "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input/7weeks-H3K4me1_bs1bp_ratio.bw",
             week4 <- "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input/4weeks-H3K4me1_bs1bp_ratio.bw",
             week2 <- "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input/2weeks-H3K4me1_bs1bp_ratio.bw",
             week0 <- "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input/ctrl-H3K4me1_bs1bp_ratio.bw")

H3K4me3 <- c(week10 <- "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input/10weeks-H3K4me3_bs1bp_ratio.bw",
             week7 <- "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input/7weeks-H3K4me3_bs1bp_ratio.bw",
             week4 <- "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input/4weeks-H3K4me3_bs1bp_ratio.bw",
             week2 <- "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input/2weeks-H3K4me3_bs1bp_ratio.bw",
             week0 <- "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input/ctrl-H3K4me3_bs1bp_ratio.bw")


plotUCSC <- function (region,yH3K4me3,yH3K4me1,yH3K27ac) {
  gene.region <- toGRanges(region)
  pp <- getDefaultPlotParams(plot.type = 1)
  pp$leftmargin <- 0.3
  pp$topmargin <- 0.3
  pp$bottommargin <- 8
  pp$ideogramheight <- 0
  pp$data1height <- 200
  pp$data1inmargin <- 1

  kp <- plotKaryotype(zoom = gene.region, cex=1,genome = "mm10",cytobands = NULL,plot.params = pp)
  gene.data <- makeGenesDataFromTxDb(TxDb.Mmusculus.UCSC.mm10.knownGene,
                                    karyoplot = kp,
                                    plot.transcripts = TRUE,
                                    plot.transcripts.structure = TRUE)

  gene.data <- addGeneNames(gene.data)
  gene.data <- mergeTranscripts(gene.data)
  kpPlotGenes(kp, data = gene.data,r0 = 0,r1 = 0.05,gene.name.cex = 1)
  kpAddBaseNumbers(kp, tick.dist = 50000,minor.tick.dist = 10000,cex = 1,tick.len = 2,digits = 2,add.units = TRUE)
  kpPlotRegions(kp,data = hmm10,col = hmm10$itemRgb,r0 = 0.08,r1 = 0.1)
  kpPlotRegions(kp,data = hmm7,col = hmm7$itemRgb,r0 = 0.11,r1 = 0.13)
  kpPlotRegions(kp,data = hmm4,col = hmm4$itemRgb,r0 = 0.14,r1 = 0.16)
  kpPlotRegions(kp,data = hmm2,col = hmm2$itemRgb,r0 = 0.17,r1 = 0.19)
  kpPlotRegions(kp,data = hmm0,col = hmm0$itemRgb,r0 = 0.20,r1 = 0.22)
  kpAddLabels(kp,labels = "10weeks",r0 = 0.08,r1 = 0.1,cex=1,label.margin = 0.01)
  kpAddLabels(kp,labels = "7weeks",r0 = 0.11,r1 = 0.13,cex=1,label.margin = 0.01)
  kpAddLabels(kp,labels = "4weeks",r0 = 0.14,r1 = 0.16,cex=1,label.margin = 0.01)
  kpAddLabels(kp,labels = "2weeks",r0 = 0.17,r1 = 0.19,cex=1,label.margin = 0.01)
  kpAddLabels(kp,labels = "0week",r0 = 0.20,r1 = 0.22,cex=1,label.margin = 0.01)

  total.tracks <- length(Otx2) + length(H3K27ac) + length(H3K4me1) + length(H3K4me3)

  out.at <- autotrack(1:length(H3K4me3), total.tracks,margin = 0.3,r0 = 0.23)
  kpAddLabels(kp,labels = "H3K4me3",r0 = out.at$r0,r1 = out.at$r1,cex=1,srt=0,pos = 1,label.margin = 0.23)

  for (i in seq_len(length(H3K4me3))) {
    at <- autotrack(i,length(H3K4me3),r0 = out.at$r0,r1 = out.at$r1,margin = 0.1)
    kp <- kpPlotBigWig(kp,data = H3K4me3[i],ymax = yH3K4me3,r0 = at$r0,r1 = at$r1,col = c("#4DAF4A"),border = NA) #"visible.region"
    computed.max <- ceiling(kp$latest.plot$computed.values$ymax)
    kpAxis(kp,ymin = 0,ymax = yH3K4me3,numticks = 2,r0 = at$r0,r1 = at$r1,col = "#4DAF4A",text.col = "#4DAF4A",cex=0.8)
    label <-  unlist(strsplit(basename(H3K4me3[i]),split = "-"))[1]
    if (label == "ctrl") {
      label <- "0week"
    }
    kpAddLabels(kp,labels = label,r0 = at$r0,r1 = at$r1,cex=1,label.margin = 0.065)
  }

  out.at <- autotrack((length(H3K4me3) +1):(length(H3K4me3) + length(H3K4me1)),total.tracks,margin = 0.1,r0 = 0.23)
  kpAddLabels(kp,labels = "H3K4me1",r0 = out.at$r0,r1 = out.at$r1,cex=1,srt=0,pos = 1,label.margin = 0.23)

  for (i in seq_len(length(H3K4me1))) {
    at <- autotrack(i,length(H3K4me1),r0 = out.at$r0,r1 = out.at$r1,margin = 0.1)
    kp <- kpPlotBigWig(kp,data = H3K4me1[i],ymax = yH3K4me1,r0 = at$r0,r1 = at$r1,col = c("#FF7F00"),border = NA)
    computed.max <- ceiling(kp$latest.plot$computed.values$ymax)
    kpAxis(kp,ymin = 0,ymax = yH3K4me1,numticks = 2,r0 = at$r0,r1 = at$r1,col = "#FF7F00",text.col = "#FF7F00",cex=0.8)
    label <- unlist(strsplit(basename(H3K4me1[i]),split = "-"))[1]
    if (label == "ctrl") {
      label <- "0week"
    }
    kpAddLabels(kp,labels = label,r0 = at$r0,r1 = at$r1,cex=1,label.margin = 0.065)
  }

  out.at <- autotrack((length(H3K4me1) + length(H3K4me3) +1):(length(H3K4me1) + length(H3K4me3) + length(H3K27ac)),total.tracks,margin = 0.1,r0 = 0.23)
  kpAddLabels(kp,labels = "H3K27ac",r0 = out.at$r0,r1 = out.at$r1,cex=1,srt=0,pos = 1,label.margin = 0.23)

  for (i in seq_len(length(H3K27ac))) {
    at <- autotrack(i,length(H3K27ac),r0 = out.at$r0,r1 = out.at$r1,margin = 0.1)
    kp <- kpPlotBigWig(kp,data = H3K27ac[i],ymax = yH3K27ac,r0 = at$r0,r1 = at$r1,col = c("#E41A1C"),border = NA)
    computed.max <- ceiling(kp$latest.plot$computed.values$ymax)
    kpAxis(kp,ymin = 0,ymax = yH3K27ac,numticks = 2,r0 = at$r0,r1 = at$r1,col = "#E41A1C",text.col = "#E41A1C",cex=0.8)
    label <- unlist(strsplit(basename(H3K27ac[i]),split = "-"))[1]
    if (label == "ctrl") {
      label <- "0week"
    }
    kpAddLabels(kp,labels = label,r0 = at$r0,r1 = at$r1,cex=1,label.margin = 0.065)
  }

  out.at <- autotrack((length(H3K4me1) + length(H3K4me3)+length(H3K27ac) +1):total.tracks,total.tracks,margin = 0.3,r0 = 0.23)
  kp <- kpPlotBigWig(kp,data = Otx2,ymax = "visible.region",r0 = out.at$r0,r1 = out.at$r1,col = c("navy"),border = NA)
  computed.max <- ceiling(kp$latest.plot$computed.values$ymax)
  kpAxis(kp,ymin = 0,ymax = computed.max,numticks = 2,r0 = out.at$r0,r1 = out.at$r1,col = "navy",text.col = "navy",cex=0.8)
  kpAddLabels(kp,labels = "Otx2",r0 = out.at$r0,r1 = out.at$r1,cex=1,srt=0,pos = 1,label.margin = 0.09)
}




# name <- "Axin2"
#plotUCSC(region = "chr11:108907828-108960783",yH3K4me3 = 200,yH3K4me1 = 60,yH3K27ac = 80)

# name <- "Thap6"
# plotUCSC(region = "chr5:91952389-91982066",yH3K4me3 = 200,yH3K4me1 = 60,yH3K27ac = 80)

# name <- "Cxcl5"
# plotUCSC(region = "chr5:90756360-90764624",yH3K4me3 = 100,yH3K4me1 = 60,yH3K27ac = 60)

# name <- "Lhfpl2"
# plotUCSC(region = "chr13:93920183-94333022",yH3K4me3 = 100,yH3K4me1 = 60,yH3K27ac = 60)

# name <- "Nipal1"
# plotUCSC(region = "chr5:72617795-72701078",yH3K4me3 = 200,yH3K4me1 = 60,yH3K27ac = 80)


plot <- function(name,region,yH3K4me3,yH3K4me1,yH3K27ac) {
  plotUCSC(region = region,yH3K4me3 = yH3K4me3,yH3K4me1 = yH3K4me1,yH3K27ac = yH3K27ac)
  filename <- paste("/data3/zhaochen/project/colon_cancer/colon_chip/peakUCSCplot/otx2/new_",name,"_ucscplot_otx2.pdf",sep = "")
  shell_cmd <- paste("mv Rplots.pdf ",filename,sep = "")
  out <- system(shell_cmd,intern = TRUE)
  cat(out)
}

#plot(name="Nipal1",region = "chr5:72624512-72694361",yH3K4me3 = 100,yH3K4me1 = 60,yH3K27ac = 60)
#plot(name="Lhfpl2",region = "chr13:94007796-94245409",yH3K4me3 = 100,yH3K4me1 = 60,yH3K27ac = 60)
#plot(name="Pabpc1",region = "chr15:36582349-36622285",yH3K4me3 = 200,yH3K4me1 = 60,yH3K27ac = 80)
plot(name="Mta3",region = "chr17:83650600-83824858",yH3K4me3 = 100,yH3K4me1 = 60,yH3K27ac = 60)
#plot(name="Cd52",region = "chr4:134069814-134107716",yH3K4me3 = 100,yH3K4me1 = 60,yH3K27ac = 60)
#plot(name="Ifitm_family",region = "chr7:140,934,161-141,033,177",yH3K4me3 = 60,yH3K4me1 = 60,yH3K27ac = 60)
#plot(name="Phlda1",region = "chr10:111,504,942-111,509,963",yH3K4me3 = 200,yH3K4me1 = 60,yH3K27ac = 100)
#plot(name="Cxcl_family",region = "chr5:90,739,130-90,923,838",yH3K4me3 = 60,yH3K4me1 = 60,yH3K27ac = 60)
#plot(name="Ywhaq",region="chr12:21,375,358-21,458,565",yH3K4me3 = 60,yH3K4me1 = 60,yH3K27ac = 60)
#plot(name = "Ahsa1",region = "chr12:87266388-87277700",yH3K4me3 = 200,yH3K4me1 = 60,yH3K27ac = 80)
#plot(name = "Ascl4",region = "chr10:85927227-85930952",yH3K4me3 = 90,yH3K4me1 = 30,yH3K27ac = 40)
#plot(name = "Pih1d2",region = "chr9:50617817-50627007",yH3K4me3 = 90,yH3K4me1 = 30,yH3K27ac = 30)

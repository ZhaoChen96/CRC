rm(list=ls())
library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BiocFileCache)
library(rtracklayer)
library(ggplot2)
library("org.Hs.eg.db")

H3K27ac <- c(
  #t78 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/T78-H3K27ac_bs1bp_ratio.bw",
  t53 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/T53-H3K27ac_bs1bp_ratio.bw",
  #t48 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/T48-H3K27ac_bs1bp_ratio.bw",
  #t45 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/T45-H3K27ac_bs1bp_ratio.bw",
  #t43 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/T43-H3K27ac_bs1bp_ratio.bw",
  #t37 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/T37-H3K27ac_bs1bp_ratio.bw",
  #t28 <- "/data1/chenjidong2/project/H3K27ac/bam/T28-H3K27ac_bwa_hg19_sorted_rmdup.bam.bw",
  t26 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/T26-H3K27ac_bs1bp_ratio.bw",
  #t21 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/T21-H3K27ac_bs1bp_ratio.bw",
  t9 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/T9-H3K27ac_bs1bp_ratio.bw",
  #t5 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/T5-H3K27ac_bs1bp_ratio.bw",
  #n78 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/N78-H3K27ac_bs1bp_ratio.bw",
  n53 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/N53-H3K27ac_bs1bp_ratio.bw",
  #n48 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/N48-H3K27ac_bs1bp_ratio.bw",
  #n45 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/N45-H3K27ac_bs1bp_ratio.bw",
  #n43 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/N43-H3K27ac_bs1bp_ratio.bw"
  #n37 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/N37-H3K27ac_bs1bp_ratio.bw",
  #n28 <- "/data1/chenjidong2/project/H3K27ac/bam/N28-H3K27ac_bwa_hg19_sorted_rmdup.bam.bw",
  n26 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/N26-H3K27ac_bs1bp_ratio.bw",
  #n21 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/N21-H3K27ac_bs1bp_ratio.bw"
  n9 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/N9-H3K27ac_bs1bp_ratio.bw"
  #n5 <- "/data3/zhaochen/project/colon_cancer/patient/chip-seq/bigwig/N5-H3K27ac_bs1bp_ratio.bw"
)

plotUCSC <- function (region,ymax) {
  gene.region <- toGRanges(region)
  pp <- getDefaultPlotParams(plot.type = 1)
  pp$leftmargin <- 0.3
  pp$topmargin <- 0.3
  pp$bottommargin <- 8
  pp$ideogramheight <- 0
  pp$data1height <- 200
  pp$data1inmargin <- 1

  kp <- plotKaryotype(zoom = gene.region, cex=1,genome = "hg19",cytobands = NULL,plot.params = pp)
  gene.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                    karyoplot = kp,
                                    plot.transcripts = TRUE,
                                    plot.transcripts.structure = TRUE)

  gene.data <- addGeneNames(gene.data)
  gene.data <- mergeTranscripts(gene.data)
  kpPlotGenes(kp, data = gene.data,r0 = 0,r1 = 0.05,gene.name.cex = 1)
  #kpAddBaseNumbers(kp, tick.dist = 50000,minor.tick.dist = 10000,cex = 1,tick.len = 2,digits = 2,add.units = TRUE)
  kpAddBaseNumbers(kp, tick.dist = 50000,minor.tick.dist = 10000,cex = 1,tick.len = 2,digits = 2,add.units = TRUE)


  for (i in seq_len(length(H3K27ac))) {
    at <- autotrack(i,length(H3K27ac),r0 = 0.08,r1 = 0.5,margin = 0.1)
    label <- unlist(strsplit(basename(H3K27ac[i]),split = "-"))[1]
    type <- substr(label,1,1)
    if (type == "N") {
      kp <- kpPlotBigWig(kp,data = H3K27ac[i],ymax = ymax,r0 = at$r0,r1 = at$r1,col = c("black"),border = NA)
      kpAxis(kp,ymin = 0,ymax = ymax,numticks = 2,r0 = at$r0,r1 = at$r1,col = "black",text.col = "black",cex=0.8)
      if (i == 5) {
        kpAddLabels(kp,labels = "Normal",r0 = at$r0,r1 = at$r1,cex=1,label.margin = 0.08,col="black")
      }
    } else if (type == "T") {
      kpAddLabels(kp,labels = "H3K27ac",r0 = at$r0,r1 = at$r1,cex=1,label.margin = 0.17,col="black")
      kp <- kpPlotBigWig(kp,data = H3K27ac[i],ymax = ymax,r0 = at$r0,r1 = at$r1,col = c("red"),border = NA)
      kpAxis(kp,ymin = 0,ymax = ymax,numticks = 2,r0 = at$r0,r1 = at$r1,col = "red",text.col = "red",cex=0.8)
      if (i== 2) {
        kpAddLabels(kp,labels = "Tumor",r0 = at$r0,r1 = at$r1,cex=1,label.margin = 0.08,col="red")
      }
    }
    computed.max <- ceiling(kp$latest.plot$computed.values$ymax)
  }
}

plot <- function(name,region,ymax) {
  plotUCSC(region = region,ymax=ymax)
  filename <- paste("/data3/zhaochen/project/colon_cancer/patient/chip-seq/pyGenomeTrack/cellline_otx2_chip/new_",name,"_ucscplot.pdf",sep = "")
  shell_cmd <- paste("mv Rplots.pdf ",filename,sep = "")
  out <- system(shell_cmd,intern = TRUE)
  cat(out)
}

#plot(name = "AXIN2",region = "chr17:63,465,572-63,700,233",ymax = 100) #21,28,53
#plot(name = "YWHAQ",region = "chr2:9700475-9955955",ymax = 100) #78,53,26
#plot(name = "AHSA1",region = "chr14:77,916,473-77,940,490",ymax = 100)

plot(name = "PABPC1",region = "chr8:101,711,234-101,743,980",ymax = 200) #9,26,53
#plot(name = "MTA3",region = "chr2:42,750,000-42,985,000",ymax = 100)
#plot(name = "PTX3",region = "chr8:128,743,058-128,759,430",ymax = 100)
#plot(name = "CD52",region = "chr1:26638000-26648203",ymax = 100) # 45,53,78
#plot(name = "IFITM_family",region = "chr11:296315-341800",ymax = 200)
#plot(name = "PHLDA1",region = "chr12:76,417,582-76,427,868",ymax = 200)
#plot(name = "Cxcl_family",region = "chr4:74554878-74987960",ymax = 100) #78,28,9
# plot(name = "CXCL1",region = "chr4:74,733,284-74,738,629",ymax = 100)
#plot(name = "MYC",region = "chr8:128,743,058-128,759,430",ymax = 100)




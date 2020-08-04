#!/usr/bin/env python3
#######################################################################
# File Name: htseq.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Wed 29 Jul 2020 04:22:40 PM CST
# Description: 
# History: 
#######################################################################
import os
import HTSeq
import numpy as np
from matplotlib import pyplot
import pandas as pd


class peakUCSC():
    def __init__(self):
        # bed
        #self.H3K27ac_dict = {"ctrl-H3K27ac":106.890540,"2weeks-H3K27ac":121.182308,"4weeks-H3K27ac":113.026669,
        #                    "7weeks-H3K27ac":88.780726,"10weeks-H3K27ac":110.733206}
        # bam file
        self.H3K27ac_bam = {"ctrl-H3K27ac":32.206495,"2weeks-H3K27ac":35.186986,"4weeks-H3K27ac":34.270315,
                            "7weeks-H3K27ac":27.069006,"10weeks-H3K27ac":31.946400}
        self.H3K4me1_dict = {"ctrl-H3K4me1":100.507372,"2weeks-H3K4me1":91.275744,"4weeks-H3K4me1":110.464083,
                             "7weeks-H3K4me1":103.184888,"10weeks-H3K4me1":93.454727}
        self.H3K4me3_dict = {"ctrl-H3K4me3":112.068098,"2weeks-H3K4me3":143.140564,"4weeks-H3K4me3":111.670897,
                            "7weeks-H3K4me3":125.350389,"10weeks-H3K4me3":62.792181}
        self.H3K27me3_dict = {"ctrl-H3K27me3":101.469167,"2weeks-H3K27me3":114.354223,"4weeks-H3K27me3":105.623027,
                              "7weeks-H3K27me3":95.733744,"10weeks-H3K27me3":81.325928}
        self.H3K9me3_dict = {"ctrl-H3K9me3": 63.783317, "2weeks-H3K9me3": 99.063894, "4weeks-H3K9me3": 86.891546,
                              "7weeks-H3K9me3": 53.068023, "10weeks-H3K9me3": 113.493940}
        self.pool_bamDir = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bam"

    def HTSeq(self,bamlist):
        # Axin2
        #window = HTSeq.GenomicInterval("chr11", 108914532, 108954079, "+")
        # Elf3
        window = HTSeq.GenomicInterval("chr1",135253574,135258472,"-")
        coverage = HTSeq.GenomicArray("auto",stranded=True,typecode="i")
        a = []
        samplelist = []
        for bamfile in bamlist:
            sample = os.path.basename(bamfile).split("_")[0]
            marker = sample.split("-")[0]
            samplelist.append(sample)
            bamfile = HTSeq.BAM_Reader(bamfile)
            for almnt in bamfile:
                if almnt.aligned:
                    almnt.iv.length = 1
                    coverage[almnt.iv] += 1

            normalization = np.fromiter(coverage[window],dtype=float)/p.H3K27ac_bam[sample]
            a.append(normalization)
        b = np.array(a)
        df = pd.DataFrame(b.T)
        df.columns = samplelist
        data = df[["ctrl-H3K27ac","2weeks-H3K27ac","4weeks-H3K27ac","7weeks-H3K27ac","10weeks-H3K27ac"]]
        data.to_csv("/data3/zhaochen/project/colon_cancer/colon_chip/peakUCSCplot/H3K27ac_Elf3.txt",sep="\t",index=False)




if __name__ == "__main__":
    step = 0

    p = peakUCSC()

    if step < 1:
        bamlist = []
        for root,dirs,files in os.walk(p.pool_bamDir):
            for file in files:
                if os.path.splitext(file)[1] == ".bam":
                    if "H3K27ac" in file:
                        bamlist.append(os.path.join(root,file))


        p.HTSeq(bamlist)

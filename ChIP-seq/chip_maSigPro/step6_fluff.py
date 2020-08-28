#!/usr/bin/env python3
#######################################################################
# File Name: step5_fluff.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Sun 16 Aug 2020 08:21:10 PM CST
# Description: 
# History: 
#######################################################################
import os,sys
import pandas as pd

class fluff():

    def __init__(self):
        self.chipmSPDir = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro"
        self.poolbamDir = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bam"
        self.regionDict = {
            "H3K27ac": "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K27ac/mm10_TSS_10kb.txt",
            "H3K4me1" : "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K27ac/mm10_TSS_10kb.txt",
            "H3K4me3" : "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K4me3/mm10_TSS_2kb.txt",
            "H3K27me3" : "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K27me3/mm10_genebody_10000.txt",
            "H3K9me3" : "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/mm10_annotation_genebody.bed"}

    def heatmap(self,marker,region,bamlist):
        output = os.path.join(f.chipmSPDir,marker,marker + "_rpkm_heatmap")
        cmd = "fluff heatmap -f %s -d %s -r -P 21 -o %s" % (region,bamlist,output)
        os.system(cmd)
        return cmd




if __name__ == "__main__":
    f = fluff()
    step = 0

    if step < 1:
        Markers = ["H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me3"]
        Times = ["ctrl","2weeks","4weeks","7weeks","10weeks"]
        for marker in Markers[1:4]:
            bamlist = []
            for time in Times:
                bamfile = os.path.join(f.poolbamDir,time + "-" + marker + "_merge.bam")
                bamlist.append(bamfile)

            bamlist = " ".join(bamlist)
            #print(bedlist)
            #print(f.regionDict[marker])
            f.heatmap(marker=marker,region=f.regionDict[marker],bamlist=bamlist)




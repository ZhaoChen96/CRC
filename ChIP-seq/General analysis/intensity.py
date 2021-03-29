#!/usr/bin/env python3
#######################################################################
# File Name: intensity.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Tue 07 Jan 2020 05:09:44 PM CST
# Description: 
# History: 
#######################################################################
import os,re
import sys
import pandas as pd
import numpy as np
from subprocess import *

class signal_intentity():

    def __init__(self):
        self.macs2 = "/data3/zhaochen/project/colon_cancer/colon_chip/macs2"
        self.intentityDir = "/data3/zhaochen/project/colon_cancer/colon_chip/intensity"
        self.bw_Input = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input"
        self.bw_without = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_withoutInput"
        self.sicer = "/data3/zhaochen/project/colon_cancer/colon_chip/SICER/H3K9me3"

    def run(self,cmd):
        p = Popen(cmd,shell=True)
        p.wait()
        return p

    def merge(self,broadPeaklist,merge_peak):
        tmp = merge_peak.replace("merged","tmp")
        sortbed = tmp.replace("tmp","sort")
        os.system("cat %s > %s " % (broadPeaklist,tmp))
        os.system("sort -k1,1 -k2,2n %s > %s" % (tmp,sortbed))
        cmd = "bedtools merge -i %s > %s" % (sortbed,merge_peak)
        return cmd

    def mutliBigwigsummary(self):
        cmd = "multi"
        return cmd

    def computeMatrix(self,bigwiglist,broadPeaklist,after,before,outName,outfilenamedata):
        cmd = "computeMatrix reference-point -S %s -R %s --referencePoint center -a %s -b %s -p 14 --skipZeros -bs 100 " \
              "-o %s --outFileSortedRegions %s" % (bigwiglist,broadPeaklist,after,before,outName,outfilenamedata)
        return cmd

    def computeMatrix_sclae(self,bigwig,broadPeak,after,before,outName,outfilenamedata):
        cmd = "computeMatrix scale-regions -S %s -R %s --startLabel TSS --endLabel TES --regionBodyLength 5000 -a %s -b %s -p 21 " \
              "--skipZeros -bs 10 -o %s --outFileNameMatrix %s" % (bigwig,broadPeak,after,before,outName,outfilenamedata)
        return cmd

    def plotProfile(self,matrix,name,title):
        colour = "'#1B9E77' '#D95F02' '#7570B3' '#E7298A' '#66A61E'"
        outfilenamedata = name.replace(".pdf",".txt")
        cmd = "plotProfile -m %s --colors %s -o %s --outFileNameData %s --plotWidth 7 --plotHeight 8 --legendLocation best --plotTitle %s " \
              "--samplesLabel ctrl 2weeks 4weeks 7weeks 10weeks" % (matrix,colour,name,outfilenamedata,title)
        return cmd

    def profile(self,matrix,name,title):
        colour = "Red"
        outfilenamedata = name.replace(".pdf", ".txt")
        cmd = "plotProfile -m %s --colors %s -o %s --outFileNameData %s --plotWidth 7 --plotHeight 8 --legendLocation best --plotTitle %s"\
              % (matrix,colour,name,outfilenamedata,title)
        return cmd


if __name__=="__main__":
    step = 5
    s = signal_intentity()
    times = ["ctrl","2weeks","4weeks","7weeks","10weeks"]
    markers = ["H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me3"]

    if step < 1:
        for marker in markers[3:5]:
            print(marker)
            term = []
            b = []
            for time in times:
                sample = "%s-%s" % (time,marker)
                broadPeak = os.path.join(s.macs2, sample + "/" + sample + "_Input_peaks.broadPeak")
                bigwig = os.path.join(s.bw_Input, sample + "_ratio.bw")
                outName = os.path.join(s.intentityDir, marker, sample + "_center_matrix.mat.gz")
                outfilenamedata = os.path.join(s.intentityDir, marker, sample + "_center_matrix.bed")
                cmd = s.computeMatrix(bigwiglist=bigwig,broadPeaklist=broadPeak,after=50000,before=50000,outName=outName,
                                      outfilenamedata=outfilenamedata)
                print(cmd)
                #s.run(cmd)

    if step < 2:
        for root,dirs,files in os.walk(s.intentityDir):
            for file in files:
                if "_center_matrix.mat.gz" in file:
                    sample = file.split(".")[0]
                    file = os.path.join(root,file)
                    name = os.path.join(root, sample + ".pdf")
                    title = sample.split("_")[0]
                    #cmd = s.plotProfile(matrix=file,name=name,title=title)
                    cmd1 = s.profile(matrix=file,name=name,title=title)
                    s.run(cmd1)

    # bigwig compare with input
    if step < 3:
        for marker in markers:
            term = []
            b = []
            for time in times:
                sample = "%s-%s" % (time,marker)
                broadPeak = os.path.join(s.macs2, sample + "/" + sample + "_Input_peaks.broadPeak")
                bigwig = os.path.join(s.bw_Input, sample + "_ratio.bw")
                outName = os.path.join(s.intentityDir, marker, sample + "_scale_matrix.mat.gz")
                outfilenamedata = os.path.join(s.intentityDir, marker, sample + "_scale_matrix.tab")
                cmd = s.computeMatrix_sclae(bigwig=bigwig,broadPeak=broadPeak,after=50000,before=50000,outName=outName,
                                      outfilenamedata=outfilenamedata)
                #s.run(cmd)
                name = os.path.join(s.intentityDir, marker, sample + "_scale.pdf")
                cmd1 = s.profile(matrix=outName,name=name,title=sample)
                # s.run(cmd1)
                outfilename = name.replace(".pdf", ".txt")
                os.system("sed -i '1d' %s" % outfilename)

    # bigwig not compare with input
    if step > 4:
        for marker in markers:
            for time in times:
                sample = "%s-%s" % (time,marker)
                broadPeak = os.path.join(s.macs2, sample + "/" + sample + "_Input_peaks.broadPeak")
                bigwig = os.path.join(s.bw_without, sample + ".bw")
                outName = os.path.join(s.intentityDir, marker, sample + "_center_matrix_noInput.mat.gz")
                outfilenamedata = os.path.join(s.intentityDir, marker, sample + "_center_matrix_noInput.tab")
                cmd = s.computeMatrix(bigwiglist=bigwig, broadPeaklist=broadPeak, after=5000, before=5000,
                                            outName=outName,
                                            outfilenamedata=outfilenamedata)
                #print(cmd)
                s.run(cmd)
                name = os.path.join(s.intentityDir, marker, sample + "_center_noInput.pdf")
                cmd1 = s.profile(matrix=outName, name=name, title=sample)
                #print(cmd1)
                s.run(cmd1)
                outfilename = name.replace(".pdf", ".txt")
                os.system("sed -i '1d' %s" % outfilename)



    if step < 5:
        for root,dirs,files in os.walk(s.intentityDir):
            for file in files:
                if "_center_matrix.txt" in file:
                    sample = file.split("_")[0]
                    file = os.path.join(root,file)
                    os.system("sed -i '1d' %s" % file)

    # SICER with H3K9me3
    if step < 0:
        for time in times:
            filelist = []
            for i in range(1, 4):
                file = os.path.join(s.sicer, "%s-R%s.rmdup-W1000-G10000-FDR0.01-island.bed" % (time,i))
                filelist.append(file)
            fi = " ".join(filelist)
            merge_peak = os.path.join(s.intentityDir, time + "_H3K9me3_peaks_merged.broadPeak")
            # cmd = s.merge(broadPeaklist=fi, merge_peak=merge_peak)
            # s.run(cmd)

            for marker in markers:
                b = []
                for time in times:
                    sample = "%s-%s" % (time, marker)
                    bigwig = os.path.join(s.bw_Input, sample + "_ratio.bw")
                    broadPeak = os.path.join(s.macs2, sample + "/" + sample + "_Input_peaks.broadPeak")
                #     b.append(bigwig)
                # bigwiglist = " ".join(b)
                outName = os.path.join(s.intentityDir, sample, "_sicer_center_matrix.mat.gz")
                outfilenamedata = os.path.join(s.intentityDir, marker, sample + "_sicer_center_matrix.tab")
                cmd1 = s.computeMatrix(bigwiglist=bigwiglist, broadPeaklist=broadPeak, after=50000, brfore=50000,outName=outName,
                                       outfilenamedata=outfilenamedata)
                s.run(cmd1)













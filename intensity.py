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

    def computeMatrix(self,bigwiglist,broadPeaklist,after,brfore,outName):
        cmd = "computeMatrix reference-point -S %s -R %s --referencePoint center -a %s -b %s -p 14 --skipZeros -bs 100 " \
              "-o %s" % (bigwiglist,broadPeaklist,after,brfore,outName)
        return cmd

    def plotProfile(self,matrix,name):
        cmd = "plotProfile -m %s --colors red orange yellow green blue -o %s" % (matrix,name)
        return cmd


if __name__=="__main__":
    step = 1
    s = signal_intentity()

    if step < 1:
        times = ["ctrl","2weeks","4weeks","7weeks","10weeks"]
        markers = ["H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me2","H3K9me3"]
        for marker in markers:
            term = []
            b = []
            for time in times:
                sample = "%s-%s" % (time,marker)
                broadPeak = os.path.join(s.macs2, sample + "/" + sample + "_Input_peaks.broadPeak")
                bigwig = os.path.join(s.bw_Input, sample + "_ratio.bw")
                term.append(broadPeak)
                b.append(bigwig)
            broadPeaklist = " ".join(term)
            bigwiglist = " ".join(b)
            if marker in markers[0:3]:
                merge_peak = os.path.join(s.intentityDir, marker + "_peaks_merged.broadPeak")
                cmd = s.merge(broadPeaklist=broadPeaklist,merge_peak=merge_peak)
                s.run(cmd)
                outName = os.path.join(s.intentityDir, marker + "_center_matrix.mat.gz")
                cmd1 = s.computeMatrix(bigwiglist=bigwiglist,broadPeaklist=merge_peak,after=5000,brfore=5000,outName=outName)
                s.run(cmd1)



    if step < 2:
        for root,dirs,files in os.walk(s.intentityDir):
            for file in files:
                if "_center_matrix.mat.gz" in file:
                    sample = file.split(".")[0]
                    file = os.path.join(root,file)
                    name = os.path.join(root, sample + ".pdf")
                    cmd = s.plotProfile(matrix=file,name=name)
                    s.run(cmd)











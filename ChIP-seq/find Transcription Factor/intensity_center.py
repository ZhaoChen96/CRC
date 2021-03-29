#!/usr/bin/env python3
#######################################################################
# File Name: intensity_center.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Sun 15 Nov 2020 08:49:18 PM CST
# Description: 
# History: 
#######################################################################
import os,sys
import pandas as pd
import numpy as np

rpmDir = "/data3/zhaochen/project/colon_cancer/colon_chip/enhancer/find_TF/rpm"

def computeMatrix_center(sample,bedfile,bigwig):
    out = os.path.join(rpmDir, "%s/%s_center_matrix.gz" % (sample,sample))
    outFileName = os.path.join(rpmDir, "%s/%s_center_matrix.txt" % (sample,sample))
    outFileSorted = os.path.join(rpmDir, "%s/%s_center_region.bed" % (sample,sample))
    cmd = "computeMatrix reference-point --referencePoint center -p 21 -S %s -R %s -a 5000 -b 5000 -bs 1 --skipZeros -o %s " \
          "--outFileNameMatrix %s --outFileSortedRegions %s" % (bigwig, bedfile, out, outFileName, outFileSorted)
    os.system(cmd)
    #print(cmd)
    return cmd

def plotProfile(matrix,sample):
    pdf = matrix.replace("center_matrix.gz","center.pdf")
    outFileName = matrix.replace("center_matrix.gz","center_line.txt")
    cmd = "plotProfile -m %s -o %s --samplesLabel %s --outFileNameData %s" % (matrix,pdf,sample,outFileName)
    os.system(cmd)
    return cmd

Times = ["ctrl","2weeks","4weeks","7weeks","10weeks"]

step = 0

if step < 1:
    #bigwigDir = "/data3/zhaochen/project/colon_cancer/chip/pooled_bigwig"
    bigwigDir = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input"
    bigwig_list = []
    sample_list = []
    for time in Times:
        sample = "%s-H3K27ac" % time
        #bigwig = os.path.join(bigwigDir,"%s-1-H3K27ac_combined_1_paired.merged.nodup.pooled_x_%s-1-Input_combined_1_paired.merged.nodup.pooled.fc.signal.bigwig" %(time,time))
        bigwig = os.path.join(bigwigDir,"%s-H3K27ac_bs1bp_ratio.bw" % time)
        file = os.path.join(rpmDir, "%s-H3K27ac/%s-H3K27ac_rm2kb.bed" % (time, time))
        computeMatrix_center(sample=sample, bedfile=file, bigwig=bigwig)
        matrix = os.path.join(rpmDir, "%s-H3K27ac/%s-H3K27ac_center_matrix.gz" % (time, time))
        plotProfile(matrix=matrix, sample=sample)

    #     bigwig_list.append(file)
    #     sample_list.append(sample)
    #
    # sample_list = " ".join(sample_list)
    # bigwig_list = " ".join(bigwig_list)





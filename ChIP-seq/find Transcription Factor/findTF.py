#!/usr/bin/env python3
#######################################################################
# File Name: findTF.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Wed 26 Aug 2020 03:05:14 PM CST
# Description: 
# History: 
#######################################################################
import os,sys
import pandas as pd
import numpy as np

rpmDir = "/data3/zhaochen/project/colon_cancer/colon_chip/enhancer/find_TF/rpm"
timepointDir = "/data3/zhaochen/project/colon_cancer/colon_chip/enhancer/find_TF/timepointTF"
poolDir = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bed"

def merge_peak(treat,ctrl):
    treatment = os.path.join(rpmDir,treat,treat +'_rm2kb.bed')
    control = os.path.join(rpmDir,ctrl,ctrl + '_rm2kb.bed')
    poolbed = os.path.join(timepointDir,treat + "_" + ctrl + "_pool.bed")
    sortbed = os.path.join(timepointDir,treat + "_" + ctrl + "_sort.bed")
    mergebed = os.path.join(timepointDir,treat + "_" + ctrl + "_merge.bed")
    treatbed = os.path.join(poolDir,treat + '_pool.bed')
    treatcount = os.path.join(timepointDir,treat + '_count.txt')
    controlbed = os.path.join(poolDir, ctrl + '_pool.bed')
    controlcount = os.path.join(timepointDir, ctrl + '_' +treat + '_count.txt')
    cmd = 'cat %s %s > %s' % (treatment,control,poolbed)
    #os.system(cmd)
    cmd1 = 'sort -k1,1 -k2,2n %s > %s' % (poolbed,sortbed)
    #os.system(cmd1)
    cmd2 = 'bedtools merge -i %s > %s' % (sortbed,mergebed)
    #os.system(cmd2)
    #os.remove(poolbed)
    #os.remove(sortbed)
    cmd3 = "bedtools intersect -c -wa -a %s -b %s > %s" % (mergebed,treatbed,treatcount)
    #os.system(cmd3)
    cmd4 = "bedtools intersect -c -wa -a %s -b %s > %s" % (mergebed,controlbed,controlcount)
    os.system(cmd4)
    return cmd3

def rpm(treat,ctrl):
    treatcount = os.path.join(timepointDir, treat + '_count.txt')
    controlcount = os.path.join(timepointDir, ctrl + '_' + treat + '_count.txt')
    dict = {"ctrl-H3K27ac": 106.890540, "2weeks-H3K27ac": 121.182308, "4weeks-H3K27ac": 113.026669,
            "7weeks-H3K27ac": 88.780726, "10weeks-H3K27ac": 110.733206}
    df1 = pd.read_csv(controlcount, sep="\t", header=None, names=["chr", "start", "end", "short"])
    df2 = pd.read_csv(treatcount, sep="\t", header=None, names=["chr", "start", "end", "long"])
    df = pd.merge(df1, df2, on=["chr", "start", "end"])
    columns = ["chr", "start", "end", "PeakID", "Not used", "strand", "long", "short", "RPM","down_RPM"]
    if treat in dict:
        df["RPM"] = df.loc[:,"long"] * dict[ctrl] / (df.loc[:,"short"] * dict[treat])
        df["down_RPM"] = df.loc[:, "short"] * dict[treat] / (df.loc[:, "long"] * dict[ctrl])
        df["PeakID"] = [treat + "_Input_peak_" + str(i) for i in range(1, len(df["chr"]) + 1)]
        df["Not used"] = "."
        df["strand"] = "+"
        dfname = os.path.join(timepointDir, "total_" + treat + "_DEEnhancers.bed")
        df.to_csv(dfname, sep="\t", header=0, index=False, columns=columns)

        data = df[df["RPM"] > 2]
        data["RPM"] = data.loc[:,"long"] * dict[ctrl] / (data.loc[:,"short"] * dict[treat])
        data["PeakID"] = [treat + "_Input_peak_" + str(i) for i in range(1, len(data["chr"]) + 1)]
        data["Not used"] = "."
        data["strand"] = "+"
        filename = os.path.join(timepointDir, "Up_" + treat + "_DEEnhancers.bed")
        data.to_csv(filename, sep="\t", header=0, index=False, columns=columns)

        down = df[df["down_RPM"] < 2]
        down["PeakID"] = [treat + "_Input_peak_" + str(i) for i in range(1, len(down["chr"]) + 1)]
        down["Not used"] = "."
        down["strand"] = "+"
        downname = os.path.join(timepointDir, "Down_" + treat + "_DEEnhancers.bed")
        down.to_csv(downname, sep="\t", header=0, index=False, columns=columns)
    return df

def findMotifsGenome(file, outputDir, size):
    cmd = "findMotifsGenome.pl %s mm10 %s -size %s -p 21" % (file, outputDir, size)
    os.system(cmd)
    return cmd

if __name__ == "__main__":
    step = 1
    marker_list = []
    Times = ['10weeks','7weeks','4weeks','2weeks']
    for time in Times:
        marker = "%s-H3K27ac" % time
        marker_list.append(marker)

    if step < 1:
        for i in range(0,3):
            treat = marker_list[i]
            ctrl = marker_list[i+1]
            #merge_peak(treat=treat,ctrl=ctrl)
            rpm(treat=treat,ctrl=ctrl)

    if step < 2:
        models = ["Up_", "total_"]
        for model in models:
            for time in Times[0:3]:
                file = os.path.join(timepointDir, model + time + "-H3K27ac_DEEnhancers.bed")
                outputDir = os.path.join(timepointDir, model.split("_")[0] + "-" + time)
                if not os.path.exists(outputDir):
                    os.makedirs(outputDir)
                #print(file)
                findMotifsGenome(file=file, outputDir=outputDir, size=600)







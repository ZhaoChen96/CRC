#!/usr/bin/env python3
#######################################################################
# File Name: peak_length.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Sat 21 Dec 2019 10:42:45 AM CST
# Description: 
# History: 
#######################################################################
import os,sys
import re
import pandas as pd
import numpy as np

class peak_length():

    def __init__(self):
        self.macs2 = "/data3/zhaochen/project/colon_cancer/colon_chip/macs2"

class peak_distribution():

    def __init__(self):
        self.peakInfo = "/data3/zhaochen/project/colon_cancer/colon_chip/macs2/peakInfo"

if __name__ == '__main__':
    pl = peak_length()
    step = 3
    if step < 1:
        modifications = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "H3K9me2", "H3K9me3"]
        for modification in modifications:
            txt = os.path.join(pl.macs2,modification + "_peak_length.txt")
            sort_txt = os.path.join(pl.macs2, modification + "_peak_length.sort.txt")
            for root,dirs,files in os.walk(pl.macs2):
                for file in files:
                    if os.path.splitext(file)[1] == ".broadPeak":
                        if modification in file:
                            sample = file.split("_")[0]
                            file = os.path.join(root,file)
                            df = pd.read_csv(file, header=None, sep="\t", names=["chr","start","end","peak","value","none",
                                                                     "fold_change","-log10Pvalue","-log10Qvalue","summit_position"])
                            peak_length = df["end"] - df["start"]
                            df1 = pd.DataFrame(peak_length,columns=["peak_length"])
                            df1["sample"] = sample
                            df1.sort_values(by=["sample"],inplace=True)
                            df1.to_csv(txt,mode="a",header=None,sep="\t",index=False)
            os.system("sort -n -k2,2 %s > %s" % (txt,sort_txt))
            os.remove(txt)

    p = peak_distribution()
    if step > 2:
        for root,dirs,files in os.walk(pl.macs2):
            for file in files:
                if "annotation" in file:
                    sample = file.split("_")[0]
                    file = os.path.join(root,file)
                    df = pd.read_csv(file,sep="\t",header=0)
                    df["A"] = [item.split("(")[0].strip() for item in df.iloc[:,7]]
#                    df["A"] = df.iloc[:,7].apply(lambda x:x[:4])
                    a = df["A"].value_counts()
                    a = pd.DataFrame(a)
                    a["sample"] = str(sample)
                    distribution = os.path.join(p.peakInfo, sample + "_peak_distribution.txt")
                    a.to_csv(distribution,header=False,index=True,sep="\t")

    if step < 3:
        for root, dirs, files in os.walk(pl.macs2):
            for file in files:
                if "broadPeak" in file:
                    sample = file.split("_")[0]
                    file = os.path.join(root,file)
                    number = os.path.join(p.peakInfo, "peak_number.txt")
                    os.system("wc -l %s >> %s" % (file, number))
                    
        total = os.path.join(p.peakInfo, "total_peak_number.txt")
        df = pd.read_csv(number, header=None,
                         sep=" ", names=["peak_number", "filename"])
        df["sample"] = df["filename"].str.split("/", expand=True)[7]
        df.drop(["filename"], axis=1, inplace=True)
        df.to_csv(total, sep="\t", header=True, index=False)
        os.remove(number)



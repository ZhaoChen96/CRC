#!/usr/bin/env python3
#######################################################################
# File Name: step5_element_percent_line.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Thu 18 Mar 2021 11:02:42 AM CST
# Description: 
# History: 
#######################################################################
import os,sys
import pandas as pd
import numpy as np 

chromDir = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM"

def state(filename,region):
    df =  pd.DataFrame()
    file = pd.read_csv(filename,sep="\t",header=None,index_col=False,
                       names=["chr","start","end","gene_name","score","strand","ensembl","chra","starta","enda","week0","week2","week4","week7","week10","num"])
    file = file.drop(columns=["chra","starta","enda"])
    Times = ["week0","week2","week4","week7","week10"]
    for time in Times:
        data = file.groupby([time]).sum()
        data["percent"] = data["num"]/sum(data["num"]) * 100
        data = data.drop(columns=["start","end","num"])
        data["time"] = time
        df = pd.concat([df,data],axis=0)

    df["state"] = df.index
    df.to_csv(os.path.join(chromDir,region + "_freq.txt"),sep="\t",index=False,header=True)
        #data.to_csv(os.path.join(chromDir,region + "_%s_freq.txt") % time,sep="\t",index=True,header=True)
    #print(sum(data["percent"]))
    #file = file[:,[]]



if __name__ == "__main__":
    step = 0

    if step < 1:
        regions = ["mm10_genebody_0bp","mm10_genebody_2000bp","mm10_genetss_4000bp"]
        for region in regions:
            filename = os.path.join(chromDir,"%s_segments_wo.bed" % region)
            #print(filename)
            state(filename=filename,region=region)
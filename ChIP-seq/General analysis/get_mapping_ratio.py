#!/usr/bin/env python3
#######################################################################
# File Name: get_mapping_ratio.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Fri 02 Oct 2020 04:28:21 PM CST
# Description: 
# History: 
#######################################################################
import pandas as pd
import numpy as np
import linecache
import os

samtoolsDir = "/data3/zhaochen/project/colon_cancer/colon_chip/samtools"

# paired
#line_number = 9
def get_line_context_paired(filename,linenumber):
    return linecache.getline(filename,linenumber).strip().split(" ")[5].strip("(").strip("%")

# overall
#line_number = 5
def get_line_context_all(filename,linenumber):
    return linecache.getline(filename,linenumber).strip().split(" ")[4].strip("(").strip("%")

dict = {}
times = ["ctrl","2weeks","4weeks","7weeks","10weeks"]
markers = ["H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me3","Input"]

step = 0
if step < 1:
    for marker in markers:
        for time in times:
            for i in range(1,4):
                sample = "%s-%s-%s" % (time,i,marker)
                filename = os.path.join(samtoolsDir,sample + "_flagstat.txt")
                dict[sample] = [get_line_context_paired(filename, linenumber=9)]

    data = pd.DataFrame(dict).T
    data.to_csv("chip-seq_mapping_ratio_paired.txt",sep="\t",index=True)


if step > 2 :
    for root,dirs,files in os.walk(samtoolsDir):
        for file in files:
            if "flagstat.txt" in file:
                sample = file.split("_")[0]
                filename = os.path.join(root,file)
                dict[sample] = [get_line_context(filename, line_number)]

    data = pd.DataFrame(dict).T
    print(data)

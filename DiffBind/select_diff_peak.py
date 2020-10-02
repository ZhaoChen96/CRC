#!/usr/bin/env python3
#######################################################################
# File Name: select_diff_peak.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Tue 08 Sep 2020 03:51:30 PM CST
# Description: 
# History: 
#######################################################################
import os,sys
import pandas as pd
import numpy as np

def select(marker):
    num = [7,9,10,4]
    dict = {}
    for i in num:
        contrast = "contrast%s" % i
        filename = "/data3/zhaochen/project/colon_cancer/colon_chip/DiffBind/diff_%s/allMethod_contrast%s.txt" % (marker,i)
        data = pd.read_csv(filename, sep="\t")
        up = data[data["Fold"] > 2]
        down = data[data["Fold"] < -2]
        dict[contrast] = {'up':len(up),'down':len(down)}

    df = pd.DataFrame(dict)
    dfname = "/data3/zhaochen/project/colon_cancer/colon_chip/DiffBind/diff_%s/%s_diff_peak_number.txt" % (marker,marker)
    df.to_csv(dfname,sep="\t",index=True)

select(marker= "H3K4me3")
select(marker= "H3K27me3")
select(marker= "H3K9me3")




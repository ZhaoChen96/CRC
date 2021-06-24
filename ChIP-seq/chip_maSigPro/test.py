#!/usr/bin/env python3
#######################################################################
# File Name: test.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Wed 22 Jul 2020 03:58:40 PM CST
# Description: 
# History: 
#######################################################################
import os,sys
import pandas as pd


mm10_annotation_genebody = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/mm10_annotation_genebody.bed"
df = pd.read_csv(mm10_annotation_genebody,sep="\t", index_col=False, names=["chr","start","end","gene_name","score","strand","ensembl"],nrows=5)

def f(df,distance1=10000,distance2=1500):
    if df["strand"] == "+":
        df["new_start"] = df["start"] - distance1
        df["new_end"] = df["start"] - distance2
    else:
        df["new_start"] = df["end"] + distance2
        df["new_end"] = df["end"] + distance1
    return df

data = df.apply(f, axis=1)

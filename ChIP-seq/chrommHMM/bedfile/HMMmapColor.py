#!/usr/bin/env python3
#######################################################################
# File Name: HMMmapColor.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Tue 06 Apr 2021 10:30:57 PM CST
# Description: 
# History: 
#######################################################################
import os,sys
import pandas as pd
import numpy as ny

# dict = {1:"#FEE5D9",2:"#CB181D",3:"#FB6A4A",4:"#238B45",5:"#006D2C",6:"#74C476",7:"#41AB5D",8:"#00441B",
#         9:"#6BAED6",10:"#08519C",11:"#FCAE91",12:"#BDBDBD",13:"#9970AB"}
dict = {1:"254,229,217",2:"203,24,29",3:"251,106,74",4:"35,139,69",5:"0,109,44",6:"116,196,118",7:"65,171,93",8:"0,68,27",
        9:"107,174,214",10:"8,81,156",11:"252,174,145",12:"189,189,189",13:"153,112,171"}

def mapColor(filename):
    bedfile = pd.read_csv(filename,sep="\t",comment="t",header=None,
                           names=["chr","start","end","name","score","strand","thickstart","thickend","itemRgb"])
    #print(type(bedfile[0:3]))
    bedfile["itemRgb"] = bedfile["name"].map(dict)
    #print(bedfile[0:3])
    newname = filename.replace(".bed","_new.bed")
    bedfile.to_csv(newname,index=False,header=False,sep="\t")

for root,dirs,files in os.walk("/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/bed/state"):
    for file in files:
        if "dense.bed" in file:
            file = os.path.join(root,file)
            mapColor(filename = file)





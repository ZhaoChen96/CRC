#!/usr/bin/env python3
#######################################################################
# File Name: single_marker_masigpro.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Tue 21 Jul 2020 07:47:23 PM CST
# Description: 
# History: 
#######################################################################
import os,sys
import pandas as pd
import numpy as np
import collections

class add_ensembl():

    def __init__(self):
        self.chip_maSigPro = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro"
        self.H3K27Dir = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K27ac"
        self.mm10_TSS_10kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K27ac/mm10_TSS_10kb.txt"
        self.H3K27me3Dir = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K27me3"
        self.genebody_10kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K27me3/mm10_genebody_10000.txt"
        self.mm10_TSS_2kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K4me3/mm10_TSS_2kb.txt"
        self.genebody = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/mm10_annotation_genebody.bed"

    def add_for_eachmarker(self,BEDfile,RawCounts):
        Dir = os.path.dirname(RawCounts)
        marker = os.path.basename(RawCounts).split("_")[0]
        df = pd.read_csv(RawCounts,sep="\t", index_col=False, header=0,quoting=3)
        df.columns = [item.strip("#'") for item in df]

        if marker == "H3K9me3":
            data = pd.read_csv(BEDfile,sep="\t",names=["chr","start","end","gene_name","score","strand","ensembl"])
            merge_df = pd.merge(df,data,how="left", left_on=["chr","start","end"],right_on=["chr","start","end"])
            merge_df = merge_df.drop_duplicates()
            merge_df = merge_df.dropna()
            merge_df.index = merge_df["ensembl"]
            merge_df = merge_df.drop(["chr","start", "end", "gene_name","score","strand","ensembl"], axis=1)
        else:
            data = pd.read_csv(BEDfile,sep="\t", names=["chr","new_start","new_end","ensembl"])
            merge_df = pd.merge(df,data,how="left", left_on=["chr","start","end"],right_on=["chr","new_start","new_end"])
            merge_df = merge_df.drop_duplicates()
            merge_df = merge_df.dropna()
            merge_df.index = merge_df["ensembl"]
            merge_df = merge_df.drop(["chr","new_start","new_end","start","end","ensembl"],axis=1)

        Times = ["ctrl","2weeks","4weeks","7weeks","10weeks"]
        labels = []
        for mark in [marker,"Input"]:
            for time in Times:
                for i in range(1,4):
                    label = "%s-%s-%s" % (time,i,mark)
                    labels.append(label)

        data_frame = merge_df.reindex(columns=labels)
        data_frame.to_csv(os.path.join(Dir,"%s_readCounts.txt" % marker),sep="\t",header=True,index=True)

    def markervsCtrl(self,readCounts):
        Dir = os.path.dirname(readCounts)
        marker = os.path.basename(readCounts).split("_")[0]
        file = pd.read_csv(readCounts,sep="\t",header=0,index_col=0)

        df = file.iloc[:,range(0,15)]
        df1 = file.iloc[:, :3]
        for i in range(0,5):
            data = pd.concat([df1,df],axis=1,verify_integrity=False)
            df = data
        data.to_csv(os.path.join(Dir,"%s_readCounts1.txt" % marker), sep="\t", header=True,index=True)







if __name__ == "__main__":
    ae = add_ensembl()
    step = 5

    if step < 1:
        H3K27ac_counts = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K27ac/H3K27ac_outRawCounts.tab"
        ae.add_for_eachmarker(BEDfile=ae.mm10_TSS_10kb,RawCounts=H3K27ac_counts)

    if step < 2:
        H3K27me3_counts = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K27me3/H3K27me3_outRawCounts.tab"
        ae.add_for_eachmarker(BEDfile=ae.genebody_10kb,RawCounts=H3K27me3_counts)

    if step < 3:
        H3K4me1_counts = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K4me1/H3K4me1_outRawCounts.tab"
        ae.add_for_eachmarker(BEDfile=ae.mm10_TSS_10kb,RawCounts=H3K4me1_counts)

    if step < 4:
        H3K4me3_counts = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K4me3/H3K4me3_outRawCounts.tab"
        ae.add_for_eachmarker(BEDfile=ae.mm10_TSS_2kb,RawCounts=H3K4me3_counts)

    if step < 5:
        H3K9me3_counts = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K9me3/H3K9me3_outRawCounts.tab"
        ae.add_for_eachmarker(BEDfile=ae.genebody, RawCounts=H3K9me3_counts)

    if step < 6:
        for root,dirs,files in os.walk(ae.chip_maSigPro):
            for file in files:
                if "readCounts.txt" in file:
                    readCounts = os.path.join(root,file)
                    ae.markervsCtrl(readCounts=readCounts)


    


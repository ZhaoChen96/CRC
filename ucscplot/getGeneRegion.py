#!/usr/bin/env python3
#######################################################################
# File Name: getGeneRegion.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Wed 07 Apr 2021 09:57:50 PM CST
# Description: 
# History: 
#######################################################################
import pandas as pd
import numpy as np
import os,sys
import csv

mm10_annotation_genebody = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/mm10_annotation_genebody.bed"

def get_genebody(upstream,downstream):
    df = pd.read_csv(mm10_annotation_genebody, sep="\t", index_col=False,
                     names=["chr", "start", "end", "gene_name", "score", "strand", "ensembl"])
    df["length"] = df["end"] - df["start"]

    df["start"] = [item - upstream for item in df["start"]]
    df["end"] = [item + downstream for item in df["end"]]

    df["value"] = df["end"] - df["start"]
    if df["value"].all() > 0:
        print("right")
    else:
        print("wrong")
    # df = df[["chr", "start", "end", "ensembl"]]
    # df["combine"] = df.apply(lambda x:'%s:%s-%s' % (x["chr"],x["start"],x["end"]),axis=0)
    # df.to_csv(os.path.join("mm10_genebody_%s.txt" % upstream), sep="\t", header=False, index=False)
    df["combine"] = df["chr"]+":"+df["start"].astype(str)+"-"+df["end"].astype(str)
    # data = pd.Series({"region" : [df["chr"]+":"+df["start"].astype(str)+"-"+df["start"].astype(str)+"    "+ df["gene_name"]]
    #                     })
    data = df[["combine","gene_name"]]
    data.to_csv("mm10_genebody_%s.txt" % upstream,sep="\t", header=False, index=False,)
    
def get_genebody_length():
    df = pd.read_csv(mm10_annotation_genebody, sep="\t", index_col=False,
                     names=["chr", "start", "end", "gene_name", "score", "strand", "ensembl"])
    df["length"] = df["end"] - df["start"]
    df["start"] = df["start"] - df["length"]
    df["end"] = df["end"] + df["length"]

    df["value"] = df["end"] - df["start"]
    if df["value"].all() > 0:
        print("right")
    else:
        print("wrong")
    # df = df[["chr", "start", "end", "ensembl"]]
    # df["combine"] = df.apply(lambda x:'%s:%s-%s' % (x["chr"],x["start"],x["end"]),axis=0)
    # df.to_csv(os.path.join("mm10_genebody_%s.txt" % upstream), sep="\t", header=False, index=False)
    df["combine"] = df["chr"]+":"+df["start"].astype(str)+"-"+df["end"].astype(str)
    # data = pd.Series({"region" : [df["chr"]+":"+df["start"].astype(str)+"-"+df["start"].astype(str)+"    "+ df["gene_name"]]
    #                     })
    data = df[["combine","gene_name"]]
    data.to_csv("mm10_genebody_add_length.txt", sep="\t", header=False, index=False,)

#get_genebody(upstream=10000,downstream=10000)
get_genebody_length()
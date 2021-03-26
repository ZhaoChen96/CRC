#!/usr/bin/env python3
#######################################################################
# File Name: wnt_target_gene_profile.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Fri 26 Mar 2021 05:09:19 PM CST
# Description: 
# History: 
#######################################################################
import os,sys
import pandas as pd
import numpy as np

wntDir = "/data3/zhaochen/project/colon_cancer/colon_chip/wnt_target_gene"
bigwigDir = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input"

def get_gene_region(filename):
    wnt_gene = pd.read_csv(filename,header=None,names=["gene_name"])
    annoation = pd.read_csv("/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/mm10_annotation_genebody.bed",
                            names=["chr","start","end","gene_name","score","strand","ensembl"],
                            header=None,sep="\t")
    data = pd.merge(wnt_gene, annoation, on="gene_name")
    output = filename.replace(".txt",".bed")
    data.to_csv(output, sep="\t",columns=["chr", "start", "end", "gene_name", "score", "strand", "ensembl"],index=False,header=False)

def computeMatrix_scale(bigwig,bedfile):
    out = bedfile.replace(".bed","_scale_matrix.gz")
    outFileName = out.replace("scale_matrix.gz","scaled.txt")
    outFileSorted = out.replace("scale_matrix.gz","gene_region.bed")
    cmd = "computeMatrix scale-regions -p 21 -R %s -S %s -b 3000 -a 3000 --regionBodyLength 5000 -bs 1 --skipZeros -o %s " \
          "--outFileNameMatrix %s --outFileSortedRegions %s" % (bedfile,bigwig,out,outFileName,outFileSorted)
    os.system(cmd)
    return cmd

def plotProfile(matrix,sample_list):
    pdf = matrix.replace("_matrix.gz", ".pdf")
    outFileName = matrix.replace("matrix.gz", "line.txt")
    cmd = "plotProfile -m %s -o %s --samplesLabel %s --outFileNameData %s" % (matrix,pdf,sample_list,outFileName)
    os.system(cmd)
    os.system("sed -i '1d' %s" % outFileName)
    return cmd

def computeMatrix_tss(bigwig,bedfile):
    out = bedfile.replace(".bed","_tss_matrix.gz")
    outFileName = out.replace("tss_matrix.gz","tss.txt")
    cmd = "computeMatrix reference-point --referencePoint TSS -p 14 -S %s -R %s -a 3000 -b 3000 -bs 1 --skipZeros -o %s " \
          "--outFileNameMatrix %s " % (bigwig,bedfile,out,outFileName)
    os.system(cmd)
    return cmd

if __name__ == "__main__":
    step = 8

    if step < 0:
        get_gene_region(filename="/data3/zhaochen/project/colon_cancer/colon_chip/wnt_target_gene/wnt_44_target_gene.txt")
        get_gene_region(filename="/data3/zhaochen/project/colon_cancer/colon_chip/wnt_target_gene/wnt_enhancer_24_target_gene.txt")

    if step > 1:
        bigwig_list = []
        sample_list = []
        for root, dirs, files in os.walk(bigwigDir):
            for file in files:
                if "H3K27ac" in file:
                    sample = file.split("_")[0]
                    file = os.path.join(root, file)
                    bigwig_list.append(file)
                    sample_list.append(sample)

        sample_list = " ".join(sample_list)
        bigwig_list = " ".join(bigwig_list)

        ### wnt 44 gene
        bedfile = "/data3/zhaochen/project/colon_cancer/colon_chip/wnt_target_gene/wnt_44_target_gene.bed"
        computeMatrix_scale(bigwig=bigwig_list, bedfile=bedfile)
        matrix = "/data3/zhaochen/project/colon_cancer/colon_chip/wnt_target_gene/wnt_44_target_gene_scale_matrix.gz"
        plotProfile(matrix=matrix, sample_list=sample_list)

        # ### wnt 24 enhancer gene
        # bedfile = "/data3/zhaochen/project/colon_cancer/colon_chip/wnt_target_gene/wnt_enhancer_24_target_gene.bed"
        # computeMatrix_scale(bigwig=bigwig_list, bedfile=bedfile)
        # matrix = "/data3/zhaochen/project/colon_cancer/colon_chip/wnt_target_gene/wnt_enhancer_24_target_gene_scale_matrix.gz"
        # plotProfile(matrix=matrix, sample_list=sample_list)

    if step > 2:
        bigwig_list = []
        sample_list = []
        for root, dirs, files in os.walk(bigwigDir):
            for file in files:
                if "H3K27ac" in file:
                    sample = file.split("_")[0]
                    file = os.path.join(root, file)
                    bigwig_list.append(file)
                    sample_list.append(sample)

        #sample_list = sorted(sample_list)
        #bigwig_list = sorted(bigwig_list)
        sample_list = " ".join(sample_list)
        bigwig_list = " ".join(bigwig_list)

        print(sample_list)
        print(bigwig_list)

        ### wnt 44 gene
        bedfile = "/data3/zhaochen/project/colon_cancer/colon_chip/wnt_target_gene/wnt_44_target_gene.bed"
        computeMatrix_tss(bigwig=bigwig_list, bedfile=bedfile)
        matrix = "/data3/zhaochen/project/colon_cancer/colon_chip/wnt_target_gene/wnt_44_target_gene_tss_matrix.gz"
        plotProfile(matrix=matrix, sample_list=sample_list)

        ### wnt 24 enhancer gene
        # bedfile = "/data3/zhaochen/project/colon_cancer/colon_chip/wnt_target_gene/wnt_enhancer_24_target_gene.bed"
        # computeMatrix_tss(bigwig=bigwig_list, bedfile=bedfile)
        # matrix = "/data3/zhaochen/project/colon_cancer/colon_chip/wnt_target_gene/wnt_enhancer_24_target_gene_tss_matrix.gz"
        # plotProfile(matrix=matrix, sample_list=sample_list)

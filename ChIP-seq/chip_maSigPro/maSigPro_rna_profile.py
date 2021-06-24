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

masigDir = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/masigpro_rna"
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
    out = bedfile.replace("_genebody.bed","_%s_scale_matrix.gz" % marker)
    outFileName = out.replace("scale_matrix.gz","scaled.txt")
    outFileSorted = out.replace("scale_matrix.gz","gene_region.bed")
    cmd = "computeMatrix scale-regions -p 21 -R %s -S %s -b 10000 -a 10000 --regionBodyLength 10000 -bs 1 --skipZeros -o %s " \
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
    out = bedfile.replace("_genebody.bed","_%s_tss_matrix.gz" % marker)
    outFileName = out.replace("tss_matrix.gz","tss.txt")
    cmd = "computeMatrix reference-point --referencePoint TSS -p 14 -S %s -R %s -a 3000 -b 3000 -bs 1 --skipZeros -o %s " \
          "--outFileNameMatrix %s " % (bigwig,bedfile,out,outFileName)
    os.system(cmd)
    return cmd

if __name__ == "__main__":
    step = 8
    Markers = ["H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me3"]

    if step < 0:
        get_gene_region(filename="/data3/zhaochen/project/colon_cancer/colon_chip/wnt_target_gene/wnt_44_target_gene.txt")
        get_gene_region(filename="/data3/zhaochen/project/colon_cancer/colon_chip/wnt_target_gene/wnt_enhancer_24_target_gene.txt")

    if step > 1:
        bigwig_list = []
        sample_list = []
        marker = Markers[3]
        for root, dirs, files in os.walk(bigwigDir):
            for file in files:
                if marker in file:
                    sample = file.split("_")[0]
                    file = os.path.join(root, file)
                    bigwig_list.append(file)
                    sample_list.append(sample)

        sample_list = " ".join(sample_list)
        bigwig_list = " ".join(bigwig_list)

        ### masigpro 4 cluster
        for root, dirs, files in os.walk(masigDir):
            for file in files:
                if "_genebody.bed" in file:
                    file = os.path.join(root, file)
                    computeMatrix_scale(bigwig=bigwig_list, bedfile=file)
                    matrix = file.replace("_genebody.bed", "_%s_scale_matrix.gz" % marker)
                    plotProfile(matrix=matrix, sample_list=sample_list)

    if step < 2:
        bigwig_list = []
        sample_list = []
        marker = Markers[2]
        for root, dirs, files in os.walk(bigwigDir):
            for file in files:
                if marker in file:
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

        ### masigpro 4 cluster
        for root,dirs,files in os.walk(masigDir):
            for file in files:
                if "_genebody.bed" in file:
                    file = os.path.join(root,file)
                    computeMatrix_tss(bigwig=bigwig_list, bedfile=file)
                    matrix = file.replace("_genebody.bed","_%s_tss_matrix.gz" % marker)
                    plotProfile(matrix=matrix, sample_list=sample_list)
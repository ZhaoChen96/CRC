#!/usr/bin/env python3
#######################################################################
# File Name: plotProfile.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Tue 05 Jan 2021 04:35:13 PM CST
# Description: 
# History: 
#######################################################################
import os,sys
import pandas as pd
import numpy as np

nfkbdir = "/data3/zhaochen/project/colon_cancer/colon_chip/inflamma/nfkb_target_gene"
bigwig_Input = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input"

def get_gene_region(file,outfile):
    genelist = pd.read_csv(file,header=None, names=["gene_name"])
    annoation = pd.read_csv(
        "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/mm10_annotation_genebody.bed",
        names=["chr", "start", "end", "gene_name", "score", "strand", "ensembl"],
        header=None, sep="\t")
    data = pd.merge(genelist, annoation, on="gene_name")
    data.to_csv(outfile,sep="\t",columns=["chr", "start", "end", "gene_name", "score", "strand", "ensembl"],
        index=False, header=False)

def computeMatrix_tss(bigwig,bedfile,out):
    #bedfile = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/wnt_target_gene/wnt_target_genebody.bed"
    #out = os.path.join(wntDir,"wnt_target_gene_tss_matrix.gz")
    outFileName = out.replace("tss_matrix.gz","tss.txt")
    #outFileSorted = out.replace("tss_matrix.gz","gene_tss_region.bed")
    cmd = "computeMatrix reference-point --referencePoint TSS -p 14 -S %s -R %s -a 3000 -b 3000 -bs 1 --skipZeros -o %s " \
          "--outFileNameMatrix %s " % (bigwig,bedfile,out,outFileName)
    #print(cmd)
    os.system(cmd)
    return cmd

def plotProfile(matrix,sample_list):
    pdf = matrix.replace("tss_matrix.gz", "tss.pdf")
    #pdf = os.path.join(wntDir, "wnt_enhancer_gene/wnt_enhancer_gene_scale_profile.pdf")
    outFileName = matrix.replace("matrix.gz", "line.txt")
    cmd = "plotProfile -m %s -o %s --samplesLabel %s --outFileNameData %s" % (matrix,pdf,sample_list,outFileName)
    os.system(cmd)
    os.system("sed -i '1d' %s" % outFileName)
    return cmd

if __name__ == "__main__":
    step = 1
    if step < 1:
        for root,dirs,files in os.walk(nfkbdir):
            for file in files:
                if "gene.txt" in file:
                    sample = file.split(".")[0]
                    file = os.path.join(root,file)
                    outputfile = file.replace("txt","bed")
                    bedlist.append(outputfile)
                    #get_gene_region(file=file,outfile=outputfile)
                    #total genelist row number from 0, other file row number from 1, so I add three file form a total gene bed file
        print(bedlist)

    if step < 2:
        bigwig_list = []
        sample_list = []
        for root, dirs, files in os.walk(bigwig_Input):
            for file in files:
                if "H3K27ac_bs1bp_ratio" in file:
                    sample = file.split("_")[0]
                    file = os.path.join(root, file)
                    bigwig_list.append(file)
                    sample_list.append(sample)

        sample_list = " ".join(sample_list)
        bigwig_list = " ".join(bigwig_list)
        #print(bigwig_list)
        #print(sample_list)

        for root,dirs,files in os.walk(nfkbdir):
            for file in files:
                if "gene.bed" in file:
                    bedfile = os.path.join(root,file)
                    out = bedfile.replace(".bed","_tss_matrix.gz")
                    computeMatrix_tss(bigwig=bigwig_list,bedfile=bedfile,out=out)
                    plotProfile(matrix=out, sample_list=sample_list)
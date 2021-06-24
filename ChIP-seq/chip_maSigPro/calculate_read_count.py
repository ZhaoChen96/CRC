#!/usr/bin/env python3
#######################################################################
# File Name: calculate_read_count.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Mon 20 Jul 2020 11:16:24 AM CST
# Description: 
# History: 
#######################################################################
import os,sys
import pandas as pd
import numpy as np
import collections
import itertools

class get_region():
    def __init__(self):
        self.path = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro"
        self.mm10_annotation_genebody = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/mm10_annotation_genebody.bed"  # have removed chrM
        self.mm10_TSS_10kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K27ac/mm10_TSS_10kb.txt"
        self.mm10_TSS_2kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K4me3/mm10_TSS_2kb.txt"
        self.mm10_TSS_1500bp = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/tf_region/mm10_TSS_1500bp.txt"
        self.mm10_down_10kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/tf_region/mm10_TSS_down10kb.txt"
        self.mm10_up_10kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/tf_region/mm10_TSS_up10kb.txt"
        self.mm10_up_100kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/tf_region/mm10_TSS_up100kb.txt"
        self.mm10_down_100kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/tf_region/mm10_TSS_down100kb.txt"

    # for H3K27ac and H3K4me1 get upstream 1500-10000kb
    def get_TSS(self):
        df = pd.read_csv(gr.mm10_annotation_genebody, sep="\t", index_col=False,
                         names=["chr","start","end","gene_name","score","strand","ensembl"])

        def f(df, distance1=100000, distance2=10000):
            if df["strand"] == "+":
                df["new_start"] = df["start"] - distance1
                df["new_end"] = df["start"] - distance2
            else:
                df["new_start"] = df["end"] + distance2
                df["new_end"] = df["end"] + distance1
            return df

        data = df.apply(f, axis=1)
        new_data = data[["chr","new_start","new_end","gene_name","ensembl"]]
        new_data.to_csv(gr.mm10_up_100kb,sep="\t",header=False,index=False)

    # for all marker get downstream 1500-10000kb
    def get_downstream(self):
        df = pd.read_csv(gr.mm10_annotation_genebody, sep="\t", index_col=False,
                         names=["chr","start","end","gene_name","score","strand","ensembl"])

        def f(df,distance1=1500,distance2=10000):
            if df["strand"] == "+":
                df["new_start"] = df["start"] + distance1
                df["new_end"] = df["start"] + distance2
            else:
                df["new_start"] = df["end"] - distance2
                df["new_end"] = df["end"] - distance1
            return df

        data = df.apply(f,axis=1)
        new_data = data[["chr","new_start","new_end","gene_name","ensembl"]]
        new_data.to_csv(gr.mm10_down_10kb, sep="\t", header=False, index=False)

    # for H3K27me3 get genebody upstream and downstream 10kb
    def get_genebody(self, marker, distance):
        df = pd.read_csv(gr.mm10_annotation_genebody, sep="\t", index_col=False,
                         names=["chr","start","end","gene_name","score","strand","ensembl"])

        df["start"] = [item - distance for item in df["start"]]
        df["end"] = [item + distance for item in df["end"]]

        df["value"] = df["end"] - df["start"]
        if df["value"].all() > 0:
            print("right")
        else:
            print("wrong")
        df = df[["chr","start","end","gene_name","ensembl"]]
        df.to_csv(os.path.join(gr.path,marker,"mm10_genebody_%s.txt" % distance),sep="\t", header=False, index=False)

    # for H3K4me3
    def get_2kb(self):
        df = pd.read_csv(gr.mm10_annotation_genebody, sep="\t", index_col=False,
                         names=["chr", "start", "end", "gene_name", "score", "strand", "ensembl"])

        def f(df,distance=1500):
            if df["strand"] == "+":
                df["new_start"] = df["start"] - distance
                df["new_end"] = df["start"] + distance
            else:
                df["new_start"] = df["end"] - distance
                df["new_end"] = df["end"] + distance
            return df

        data = df.apply(f, axis=1)
        new_data = data[["chr", "new_start", "new_end", "gene_name","ensembl"]]
        new_data.to_csv(gr.mm10_TSS_1500bp, sep="\t", header=False, index=False)


class calculate_reads_count():
    def __init__(self):
        self.bamDir = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bam"
        self.poolDir = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/pool"
        self.samtoolsDir = "/data3/zhaochen/project/colon_cancer/colon_chip/samtools"
        self.H3K27acDir = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K27ac"
        self.H3K27me3Dir = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K27me3"
        self.mm10_genebody_10kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K27me3/mm10_genebody_10000.txt"
        #self.merged_peak =

    # draw plotPCA to explain inflammation and cancer
    def multiBamSummary_bins(self, bamlist, labels):
        results = os.path.join(rc.poolDir,"25_results.npz")
        readCounts = os.path.join(rc.poolDir,"readCounts.tab")
        cmd = "multiBamSummary bins -b %s -l %s -e 147 --numberOfProcessors 21 -o %s --outRawCounts %s" % (bamlist, labels, results,readCounts)
        os.system(cmd)
        return cmd


    def plotCorrelation(self):
        readCounts = os.path.join(rc.poolDir,"25_results.npz")
        out1 = os.path.join(rc.poolDir,"heatmap_SpearmanCorr_readCounts.pdf")
        outFile1 = os.path.join(rc.poolDir,"SpearmanCorr_readCounts.tab")
        out2 = os.path.join(rc.poolDir, "heatmap_PearsonCorr_readCounts.pdf")
        outFile2 = os.path.join(rc.poolDir, "PearsonCorr_readCounts.tab")
        cmd1 = "plotCorrelation -in %s --skipZeros --corMethod spearman --whatToPlot heatmap --plotTitle 'Spearman Correlation of ChIP-seq Read Counts' " \
              "--colorMap YlOrRd -o %s --outFileCorMatrix %s" % (readCounts,out1,outFile1)
        cmd2 = "plotCorrelation -in %s --skipZeros --corMethod pearson --whatToPlot heatmap --plotTitle 'Pearson Correlation of ChIP-seq Read Counts' " \
               "--colorMap YlOrRd -o %s --outFileCorMatrix %s" % (readCounts, out2, outFile2)
        os.system(cmd1)
        os.system(cmd2)
        return cmd1,cmd2

    def plotPCA(self,marker):
        ### gobal
        #readCounts = os.path.join(rc.poolDir,"35_results.npz")
        #out  = os.path.join(rc.poolDir,"PCA_readCounts.pdf")
        ### each marker
        readCounts = os.path.join(gr.path,marker,"%s_results.npz" % marker)
        out = os.path.join(gr.path,marker,"PCA_%s_readCounts.pdf" % marker)
        cmd = "plotPCA -in %s --log2 -o %s -T 'PCA of %s ChIP-seq read counts'" % (readCounts,out,marker)
        os.system(cmd)
        return cmd

    # find each marker special region reads count for chip-seq maSigPro
    def multiBamSummary_BEDfile(self, marker, bamlist, BEDfile, labels):
        out = os.path.join(gr.path,marker,"%s_results.npz" % marker)
        RawCount = os.path.join(gr.path,marker,"%s_outRawCounts.tab" % marker)
        cmd = "multiBamSummary BED-file --BED %s --bamfiles %s --numberOfProcessors 21 --minMappingQuality 30 --labels %s -out %s " \
              "--outRawCounts %s" % (BEDfile, bamlist, labels, out, RawCount)
        os.system(cmd)
        return cmd



if __name__ == "__main__":
    step = 8
    gr = get_region()
    rc = calculate_reads_count()
    Times = ["ctrl", "2weeks", "4weeks", "7weeks", "10weeks"]
    Markers = ["Input", "H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"]

    if step < 1:
        gr.get_TSS()
        #gr.get_2kb()
        #gr.get_genebody(marker="H3K27me3",distance=10000)
        #gr.get_downstream()


    # all chip-seq data correlation and PCA
    if step < 2:
        bamlist = []
        labels = []
        for root,dirs,files in os.walk(rc.bamDir):
            for file in files:
                if "Input" not in file:
                    if "H3K9me2" not in file:
                        if os.path.splitext(file)[1] == ".bam":
                            bamlist.append(os.path.join(root,file))
                            label = file.split("_")[0]
                            labels.append(label)

        bamlist.sort()
        labels.sort()
        #print(bamlist)
        #print(labels)
        bamlist = " ".join(bamlist)
        labels = " ".join(labels)
        #print(labels)
#        rc.multiBamSummary_bins(bamlist=bamlist,labels=labels)
#        rc.plotPCA()
        rc.plotCorrelation()

    # all Markers read counts file perpare for masigpro
    if step < 3:
        labels = []
        bamlist = []
        for root,dirs,files in os.walk(rc.samtoolsDir):
            for file in files:
                if os.path.splitext(file)[1] == ".bam":
                    bamlist.append(os.path.join(root, file))
                    label = file.split("_")[0]
                    labels.append(label)

        bamlist.sort()
        labels.sort()
        bamlist = " ".join(bamlist)
        labels = " ".join(labels)

        rc.multiBamSummary_BEDfile(marker="all",BEDfile= rc.merged_peak,bamlist=bamlist,labels=labels)

    # signal marker perpare for masigpro input file
    if step < 4:
        labels = []
        bamlist = []
        for root,dirs,files in os.walk(rc.samtoolsDir):
            for file in files:
                if os.path.splitext(file)[1] == ".bam":
                    if Markers[1] in file or "Input" in file:
                        bamlist.append(os.path.join(root, file))
                        label = file.split("_")[0]
                        labels.append(label)


        bamlist.sort()
        labels.sort()
        print(bamlist)
        print(labels)
        bamlist = " ".join(bamlist)
        labels = " ".join(labels)

        # H3K27ac
        #rc.multiBamSummary_BEDfile(marker=Markers[1],BEDfile=rc.mm10_genebody_10kb, bamlist=bamlist, labels=labels)
        # H3K27me3
        #rc.multiBamSummary_BEDfile(marker=Markers[4], BEDfile=rc.mm10_genebody_10kb, bamlist=bamlist, labels=labels)
        # H3K4me1
        #rc.multiBamSummary_BEDfile(marker=Markers[2], BEDfile=rc.mm10_genebody_10kb, bamlist=bamlist, labels=labels)
        # H3K4me3
        #rc.multiBamSummary_BEDfile(marker=Markers[3], BEDfile=rc.mm10_genebody_10kb, bamlist=bamlist, labels=labels)
        # H3K9me3
        #rc.multiBamSummary_BEDfile(marker=Markers[5],BEDfile=gr.mm10_annotation_genebody,bamlist=bamlist,labels=labels)

        # signal marker plotPCA
        if step < 5:
            rc.plotPCA(marker="H3K27me3")

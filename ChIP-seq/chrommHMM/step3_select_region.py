#!/usr/bin/env python3
#######################################################################
# File Name: step2_select_region.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Sat 16 Jan 2021 11:31:06 AM CST
# Description: 
# History: 
#######################################################################
import os,sys
import pandas as pd
import numpy as np
import logging
#logging.basicConfig(level=logging.DEBUG,format=' %(asctime)s - %(levelname)s - %(message)s')
#logging.debug('Start of program')

chromDir = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM"
stateDir = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/bed/state"
genecountDir = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/genecount"
samtoolsDir = "/data3/zhaochen/project/colon_cancer/colon_chip/samtools"
chipDir = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/genecount/chip"

def get_genebody(length):
    file = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/mm10_annotation_genebody.bed"
    data = pd.read_csv(file,sep="\t",header=None,names=["chr","start","end","gene_name","score","strand","ensembl"])
    data["newstart"] = data["start"] - length
    data["newend"] = data["end"] + length
    df = data.loc[:,["chr","newstart","newend","gene_name","score","strand","ensembl"]]
    df.to_csv("/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/mm10_genebody_%sbp.bed" % length,sep="\t",index=False,header=None)

def get_tss():
    file = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/mm10_annotation_genebody.bed"
    df = pd.read_csv(file,sep="\t",header=None,names=["chr","start","end","gene_name","score","strand","ensembl"])

    def f(df,length=10000): # change this number
        if df["strand"] == "+":
            df["new_start"] = df["start"] - length
            df["new_end"] = df["start"] + length
        else:
            df["new_start"] = df["end"] - length
            df["new_end"] = df["end"] + length
        return df

    data = df.apply(f,axis=1)
    df = data.reindex(columns=["chr","new_start","new_end","gene_name","score","strand","ensembl"])
    df.to_csv(os.path.join(chromDir,"mm10_genetss_%sbp.bed" % 10000),sep="\t",index=False,header=None) # change this number


def chromatin_on_gene(genebody):
    segment_bin = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/merge_13_segments_200.txt"

    # try whether -wa and -wo have differences

    # -wa
    #output = genebody.replace(".bed","_segments.bed")
    #cmd = "bedtools intersect -a %s -b %s -wa -c > %s" % (genebody,segment_bin,output)

    # -wo
    output = genebody.replace(".bed", "_segments_wo.bed")
    cmd = "bedtools intersect -a %s -b %s -wo > %s" % (genebody,segment_bin,output)
    os.system(cmd)
    return cmd

# calculate state percentage
def calculate(segmentbed):
    df = pd.read_csv(segmentbed,sep="\t",
                     names=["chr","start","end","gene_name","score","strand","ensembl",
                            "chr","start","end","week0","week2","week4","week7","week10","overlap"])
    # we use week 0 state as whole genome gene stateannotation
    df_0 = df.groupby(['gene_name','ensembl', 'week0'])['overlap'].sum().reset_index(name="base")
    # df_0.to_csv("/home/zyang/Project/CRC/step47_geneExp_state/gene_base_summary.csv", index=False)
    idx = df_0.groupby(['gene_name','ensembl'])["base"].transform(max) == df_0['base']
    df_0 = df_0[idx]


    df_2 = df.groupby(['gene_name','ensembl', 'week2'])['overlap'].sum().reset_index(name="base")
    idx = df_2.groupby(['gene_name','ensembl'])["base"].transform(max) == df_2['base']
    df_2 = df_2[idx]
    df_merge = pd.merge(df_0, df_2, how="outer", on=["gene_name",'ensembl'])

    df_4 = df.groupby(['gene_name','ensembl', 'week4'])['overlap'].sum().reset_index(name="base")
    idx = df_4.groupby(['gene_name','ensembl'])["base"].transform(max) == df_4['base']
    df_4 = df_4[idx]
    df_merge = pd.merge(df_merge, df_4, how="outer", on=["gene_name",'ensembl'])

    df_7 = df.groupby(['gene_name','ensembl', 'week7'])['overlap'].sum().reset_index(name="base")
    idx = df_7.groupby(['gene_name','ensembl'])["base"].transform(max) == df_7['base']
    df_7 = df_7[idx]
    df_merge = pd.merge(df_merge, df_7, how="outer", on=["gene_name",'ensembl'])

    df_10 = df.groupby(['gene_name','ensembl','week10'])['overlap'].sum().reset_index(name="base")
    idx = df_10.groupby(['gene_name','ensembl'])["base"].transform(max) == df_10['base']
    df_10 = df_10[idx]
    df_merge = pd.merge(df_merge, df_10, how="outer", on=["gene_name",'ensembl'])

    df_merge = df_merge[["gene_name","ensembl","week0", "week2", "week4", "week7", "week10"]]
    basename = os.path.basename(segmentbed)
    filename = os.path.join(genecountDir,basename.replace(".bed","_geneCount_alltime.csv"))
    #filename = os.path.join("/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/test",basename.replace(".bed","_geneCount_alltime.csv"))
    df_merge.to_csv(filename, index=False)

def merge_rna_fpkm(genecount):
    fpkm = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/genes.fpkm_table.txt"
    file = pd.read_csv(fpkm,sep="\t")
    file['week0_mean'] = file[['control_0','control_1','control_2']].mean(axis=1)
    file['week2_mean'] = file[['2week_0', '2week_1', '2week_2']].mean(axis=1)
    file['week4_mean'] = file[['4week_0', '4week_1', '4week_2']].mean(axis=1)
    file['week7_mean'] = file[['7week_0', '7week_1', '7week_2']].mean(axis=1)
    file['week10_mean'] = file[['10week_0', '10week_1', '10week_2']].mean(axis=1)

    df = pd.read_csv(genecount)
    df_merge = file.merge(df, left_on="tracking_id", right_on="ensembl", how="inner")
    print(df_merge[0:10])
    df_merge.to_csv(genecount.replace("_segments_wo_geneCount_alltime.csv","_merge_rna.csv"), index=False)

def multiBamSummary_BEDfile(sample, bamlist, BEDfile, labels):
    out = os.path.join(genecountDir, "chip/%s_results.npz" % sample)
    RawCount = os.path.join(genecountDir,"chip/%s_outRawCounts.tab" % sample)
    cmd = "multiBamSummary BED-file --BED %s --bamfiles %s --numberOfProcessors 21 --minMappingQuality 30 --labels %s -out %s " \
          "--outRawCounts %s" % (BEDfile, bamlist, labels, out, RawCount)
    os.system(cmd)
    print(cmd)
    return cmd

def rpm(rawCount):
    dict = {
        'H3K27ac': pd.Series([20148549, 7537172, 6631827, 10447349, 14577061, 12625136, 14992457, 11270684, 10123389,
                              11656015, 12593325, 4441176, 17747602, 12298717, 3449216]),
        'H3K4me1': pd.Series([11729065, 8534148, 10775292, 8563526, 9182314, 11435226, 11238687, 11751678, 11434737,
                              13971382, 10284665, 7268299, 13824201, 3323698, 9741392]),
        'H3K4me3': pd.Series([14192207, 8556148, 10487448, 16700677, 11422956, 14070443, 10611359, 11320168, 12175416,
                              12193155, 15268777, 12149746, 5463889, 2921882, 10214387]),
        'H3K27me3': pd.Series([11635132, 8041169, 11855480, 10595682, 14524359, 10111478, 6878257, 14111538, 13076153,
                               11895055, 10714586, 7986872, 13541547, 3217629, 9287425]),
        'H3K9me2': pd.Series([21027858,7614592,13955063,13210419,10276422,15428638,11355885,10044634,10129604,
                              10078816,12072858,10351051,9979597,12811407]),
        'H3K9me3': pd.Series([9151767, 4726390, 6678253, 7936621, 9157008, 14095696, 8043469, 8884348, 11063354,
                              5976747, 6961629, 4203301, 13556435, 7952982, 13510368]),
        'Input': pd.Series([11623035, 7718653, 11267834, 13984409, 9180133, 11225237, 10965254, 12689517, 14403899,
                            10370311, 12236755, 15259467, 13896129, 10965790, 13972058])
    }

    Times = ["ctrl","2weeks","4weeks","7weeks","10weeks"]
    Markers = ["H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me3","Input"]

    columns = ['chr','start','end',
               '10weeks-1-H3K27ac','10weeks-1-H3K27me3','10weeks-1-H3K4me1','10weeks-1-H3K4me3','10weeks-1-H3K9me3','10weeks-1-Input',
               '10weeks-2-H3K27ac','10weeks-2-H3K27me3','10weeks-2-H3K4me1','10weeks-2-H3K4me3','10weeks-2-H3K9me2','10weeks-2-H3K9me3','10weeks-2-Input',
               '10weeks-3-H3K27ac','10weeks-3-H3K27me3','10weeks-3-H3K4me1','10weeks-3-H3K4me3','10weeks-3-H3K9me2','10weeks-3-H3K9me3','10weeks-3-Input',
               '2weeks-1-H3K27ac','2weeks-1-H3K27me3','2weeks-1-H3K4me1','2weeks-1-H3K4me3','2weeks-1-H3K9me2','2weeks-1-H3K9me3','2weeks-1-Input',
               '2weeks-2-H3K27ac','2weeks-2-H3K27me3','2weeks-2-H3K4me1','2weeks-2-H3K4me3','2weeks-2-H3K9me2','2weeks-2-H3K9me3','2weeks-2-Input',
               '2weeks-3-H3K27ac','2weeks-3-H3K27me3','2weeks-3-H3K4me1','2weeks-3-H3K4me3','2weeks-3-H3K9me2','2weeks-3-H3K9me3','2weeks-3-Input',
               '4weeks-1-H3K27ac','4weeks-1-H3K27me3','4weeks-1-H3K4me1','4weeks-1-H3K4me3','4weeks-1-H3K9me2','4weeks-1-H3K9me3','4weeks-1-Input',
               '4weeks-2-H3K27ac','4weeks-2-H3K27me3','4weeks-2-H3K4me1','4weeks-2-H3K4me3','4weeks-2-H3K9me2','4weeks-2-H3K9me3','4weeks-2-Input',
               '4weeks-3-H3K27ac','4weeks-3-H3K27me3','4weeks-3-H3K4me1','4weeks-3-H3K4me3','4weeks-3-H3K9me2','4weeks-3-H3K9me3','4weeks-3-Input',
               '7weeks-1-H3K27ac','7weeks-1-H3K27me3','7weeks-1-H3K4me1','7weeks-1-H3K4me3','7weeks-1-H3K9me2','7weeks-1-H3K9me3','7weeks-1-Input',
               '7weeks-2-H3K27ac','7weeks-2-H3K27me3','7weeks-2-H3K4me1','7weeks-2-H3K4me3','7weeks-2-H3K9me2','7weeks-2-H3K9me3','7weeks-2-Input',
               '7weeks-3-H3K27ac','7weeks-3-H3K27me3','7weeks-3-H3K4me1','7weeks-3-H3K4me3','7weeks-3-H3K9me2','7weeks-3-H3K9me3','7weeks-3-Input',
               'ctrl-1-H3K27ac','ctrl-1-H3K27me3','ctrl-1-H3K4me1','ctrl-1-H3K4me3','ctrl-1-H3K9me2','ctrl-1-H3K9me3','ctrl-1-Input',
               'ctrl-2-H3K27ac','ctrl-2-H3K27me3','ctrl-2-H3K4me1','ctrl-2-H3K4me3','ctrl-2-H3K9me2','ctrl-2-H3K9me3','ctrl-2-Input',
               'ctrl-3-H3K27ac','ctrl-3-H3K27me3','ctrl-3-H3K4me1','ctrl-3-H3K4me3','ctrl-3-H3K9me2','ctrl-3-H3K9me3','ctrl-3-Input']
    file = pd.read_csv(rawCount,sep="\t",comment="#",header=None,names=columns,index_col=False)
    basename = os.path.basename(rawCount)
    #print(basename)
    regionfile = os.path.join(chromDir, basename.replace("_outRawCounts.tab",".bed"))
    #print(regionfile)
    genebody = pd.read_csv(regionfile,sep="\t",header=None,names=["chr","start","end","gene_name","score","strand","ensembl"])
    merge = pd.merge(genebody,file,on=["chr","start","end"])
    df = merge.set_index("gene_name")
    df = df.drop(["chr","start","end","score","strand","ensembl"],axis=1)
    df.to_csv(rawCount.replace("_outRawCounts.tab","_merge.txt"),sep="\t")
    #print(df[0:3])


    # for marker in Markers[4]:
    #     print(marker)
        #print(file.loc[:,file.columns.str.contains(marker)])
    marker = Markers[4]
    print(marker)
    col_name = ["ctrl-1-%s" % marker,"ctrl-2-%s" % marker,"ctrl-3-%s" % marker,"2weeks-1-%s" % marker,
                    "2weeks-2-%s" % marker,"2weeks-3-%s" % marker,"4weeks-1-%s" % marker,"4weeks-2-%s" % marker,"4weeks-3-%s" % marker,
                    "7weeks-1-%s" % marker,"7weeks-2-%s" % marker,"7weeks-3-%s" % marker,"10weeks-1-%s" % marker,"10weeks-2-%s" % marker,
                    "10weeks-3-%s" % marker]
    data = df.reindex(columns=col_name)
    #print(data[0:3])
    bamCount = list(dict[marker] / 1000000)
    #print(data.loc[:, data.columns.str.contains(marker)])
    new_data = data.loc[:,data.columns.str.contains(marker)] / bamCount
    print(len(new_data))
    print(new_data.iloc[:, 0:3])
    name = os.path.join(chipDir,marker,basename.replace("_outRawCounts.tab","_%s_rpm.txt" % marker))
    new_data.to_csv(name, sep="\t", index=True,header=True)
    #d = pd.concat([data,new_data],axis=1)
    #print(d[0:2])

Times = ["ctrl","2weeks","4weeks","7weeks","10weeks"]
Markers = ["H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me3","Input"]

def merge_chromatin_rpm(segmentfile,sample):
    segment = pd.read_csv(segmentfile, sep=",")
    for marker in Markers[0:5]:
        filename = os.path.join(chipDir,marker,"%s_%s_rpm.txt" % (sample,marker))
        file = pd.read_csv(filename,sep="\t")
        data = pd.merge(file,segment,on="gene_name")
        output = filename.replace("_rpm.txt","_merge_chip.txt")
        data.to_csv(output,sep="\t",header=True,index=False)


if __name__ == "__main__":
    step = 10

    filelist = []
    segments = ["ctrl_13_segments.bed","2weeks_13_segments.bed","4weeks_13_segments.bed","7weeks_13_segments.bed","10weeks_13_segments.bed"]
    if step < 1:
        for segment in segments:
            file = os.path.join(stateDir,segment)
            filelist.append(file)

        merge_state()

    if step < 2:
        # get gene body upstream and downstream a special region
        lengths = [0,2000,4000,10000]
        #for length in lengths[2:4]:
            #get_genebody(length=length)

        # get gene tss upstream and downstream a special region
        get_tss()

    if step < 3:
        for root,dirs,files in os.walk(chromDir):
            for file in files:
                # step1 get gene region and segment bin overlap
                # if "bp.bed" in file:
                #     chromatin_on_gene(genebody=file)
                # step2 calculate
                if "bp_segments_wo.bed" in file:
                    #if "genetss" in file:
                    print(file)
                    calculate(segmentbed=file)

    if step < 4:
        for root,dirs,files in os.walk(genecountDir):
            for file in files:
                if "wo_geneCount_alltime" in file:
                    file = os.path.join(root,file)
                    merge_rna_fpkm(genecount=file)

    # from this step
    if step < 5:
        labels = []
        bamlist = []
        for root, dirs, files in os.walk(samtoolsDir):
            for file in files:
                if os.path.splitext(file)[1] == ".bam":
                    bamlist.append(os.path.join(root, file))
                    label = file.split("_")[0]
                    labels.append(label)

        bamlist.sort()
        labels.sort()
        bamlist = " ".join(bamlist)
        labels = " ".join(labels)

        # for root,dirs,files in os.walk(chromDir):
        #     for file in files:
        #         if "bp.bed" in file:
        #             sample = file.split(".")[0]
        #             file = os.path.join(root,file)
        #             print(file)
        #             multiBamSummary_BEDfile(sample=sample, BEDfile=file, bamlist=bamlist, labels=labels)

        #### this step for tss removed
        filelist = ["/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/mm10_genetss_2000bp.bed",
                    "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/mm10_genetss_4000bp.bed",
                    "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/mm10_genetss_10000bp.bed"]

        for file in filelist:
            filename = os.path.basename(file)
            sample = filename.split(".")[0]
            multiBamSummary_BEDfile(sample=sample, BEDfile=file, bamlist=bamlist, labels=labels)
            RawCount = os.path.join(genecountDir, "chip/%s_outRawCounts.tab" % sample)
            print(RawCount)
            rpm(rawCount=RawCount)

    if step < 6:
        for root,dirs,files in os.walk(chipDir):
            for file in files:
                if "_outRawCounts.tab" in file:
                    print(file)
                    file = os.path.join(root,file)
                    rpm(rawCount=file)

        #rpm(rawCount="/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/genecount/chip/mm10_genebody_0bp_outRawCounts.tab")

    if step > 7:
        for root,dirs,files in os.walk(genecountDir):
            for file in files:
                if "segments_wo_geneCount_alltime.csv" in file:
                    sample = file.split("_segments")[0]
                    file = os.path.join(root,file)
                    print(sample)
                    #print(file)
                    merge_chromatin_rpm(segmentfile=file,sample=sample)
        
    #logging.debug('End of program')

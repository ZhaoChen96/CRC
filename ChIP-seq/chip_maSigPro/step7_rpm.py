#!/usr/bin/env python3
#######################################################################
# File Name: step7_rpm.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Fri 11 Sep 2020 08:03:10 PM CST
# Description: 
# History: 
#######################################################################
import os,sys
import numpy as np
import pandas as pd


samtoolsDir = "/data3/zhaochen/project/colon_cancer/colon_chip/samtools"
tfDir = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/tf_region"
mm10_genebody_10kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K27me3/mm10_genebody_10000.txt"
mm10_tss_2kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/H3K4me3/mm10_TSS_2kb.txt"
mm10_tss_1500bp = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/tf_region/mm10_TSS_1500bp.txt"
mm10_tss_up10kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/tf_region/mm10_TSS_up10kb.txt"
mm10_tss_down10kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/tf_region/mm10_TSS_down10kb.txt"
mm10_tss_up100kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/tf_region/mm10_TSS_up100kb.txt"
mm10_tss_down100kb = "/data3/zhaochen/project/colon_cancer/colon_chip/chip_maSigPro/tf_region/mm10_TSS_down100kb.txt"

dict = {
    'H3K27ac': pd.Series([20148549,7537172,6631827,10447349,14577061,12625136,14992457,11270684,10123389,
                11656015,12593325,4441176,17747602,12298717,3449216]),
    'H3K4me1': pd.Series([11729065,8534148,10775292,8563526,9182314,11435226,11238687,11751678,11434737,
                        13971382,10284665,7268299,13824201,3323698,9741392]),
    'H3K4me3': pd.Series([14192207,8556148,10487448,16700677,11422956,14070443,10611359,11320168,12175416,
                          12193155,15268777,12149746,5463889,2921882,10214387]),
    'H3K27me3': pd.Series([11635132,8041169,11855480,10595682,14524359,10111478,6878257,14111538,13076153,
                           11895055,10714586,7986872,13541547,3217629,9287425]),
    'H3K9me3': pd.Series([9151767,4726390,6678253,7936621,9157008,14095696,8043469,8884348,11063354,
                          5976747,6961629,4203301,13556435,7952982,13510368]),
    'Input': pd.Series([11623035,7718653,11267834,13984409,9180133,11225237,10965254,12689517,14403899,
                        10370311,12236755,15259467,13896129,10965790,13972058])
}
Markers = ["H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me3"]

# find each marker special region reads count for chip-seq maSigPro
def multiBamSummary_BEDfile(marker, bamlist, BEDfile, labels, state, region):
    out = os.path.join(tfDir,"%s_results_%s_%s.npz" % (marker,state,region))
    RawCount = os.path.join(tfDir,"%s_outRawCounts_%s_%s.txt" % (marker,state,region))
    cmd = "multiBamSummary BED-file --BED %s --bamfiles %s --numberOfProcessors 21 --minMappingQuality 30 --labels %s -out %s " \
          "--outRawCounts %s" % (BEDfile, bamlist, labels, out, RawCount)
    print(cmd)
    os.system(cmd)
    return cmd

def multiBamSummary_BEDfile_promoter(marker, bamlist, BEDfile, labels):
    out = os.path.join(tfDir,"%s_results_1500bp.npz" % (marker))
    RawCount = os.path.join(tfDir,"%s_outRawCounts_1500bp.txt" % (marker))
    cmd = "multiBamSummary BED-file --BED %s --bamfiles %s --numberOfProcessors 21 --minMappingQuality 30 --labels %s -out %s " \
          "--outRawCounts %s" % (BEDfile, bamlist, labels, out, RawCount)
    print(cmd)
    os.system(cmd)
    return cmd

def rpm(marker,regionfile,outCount,name):
    Times = ['10weeks','2weeks','4weeks','7weeks','ctrl']
    sample_list = ["chr","start","end"]
    for time in Times:
        for i in range(1,4):
            sample = "%s-%s-%s" % (time,i,marker)
            sample_list.append(sample)

    times = ['ctrl', '2weeks', '4weeks', '7weeks', '10weeks']
    column_list = []
    for time in times:
        for i in range(1,4):
            sample = "%s-%s-%s" % (time, i, marker)
            column_list.append(sample)

    file = pd.read_csv(regionfile, sep="\t",header=None)
    file.columns = ["chr","start","end","gene_name","ensembl"]
    data = pd.read_csv(outCount,sep="\t",header=None,comment="#")
    data.columns = sample_list
    df = pd.merge(file,data,on=["chr","start","end"],validate="1:1")
    d = df.iloc[:,5:20]
    print(len(d))
    d.index = df.iloc[:,3]
    ### calculate rpm
    bamCount = list(dict[marker]/1000000)
    new_data = d/bamCount
    print(len(new_data))
    new_data.to_csv(name,sep="\t",columns=column_list,index=True)

def addfile(marker,condition):
    addrpm = os.path.join(tfDir, marker + "_addrpm_%s.txt" % condition)
    upfile = os.path.join(tfDir, marker + "_rpm_up_%s.txt" % condition)
    up = pd.read_csv(upfile, sep="\t")
    downfile = os.path.join(tfDir, marker + "_rpm_down_%s.txt" % condition)
    down = pd.read_csv(downfile, sep="\t")
    # print(len(down))
    # print(down[1:3])
    data = pd.merge(up,down,on="gene_name")
    # print(data.iloc[1:3,1:15])
    # print(data.iloc[1:3,16:31])
    up_df = pd.DataFrame(data.iloc[:,0:16])
    up_df.columns = up.columns
    up_df.set_index(["gene_name"],inplace=True)
    down_df = data.iloc[:,16:32]
    down_df.columns = up.columns[1:16]
    down_df.index = data.iloc[:,0]
    df = up_df + down_df
    # print(len(df))
    # print(df[1:3])
    # print(len(up_df))
    # print(up_df[1:3])
    # print(len(down_df))
    # print(down_df[1:3])
    df.to_csv(addrpm, sep="\t",index=True)



#upload file to mac and select target gene in R

step = 5

if step < 1:
    for marker in Markers:
        labels = []
        bamlist = []
        for root,dirs,files in os.walk(samtoolsDir):
            for file in files:
                if os.path.splitext(file)[1] == ".bam":
                    if marker in file:
                        sample = file.split("_")[0]
                        file = os.path.join(root, file)
                        bamlist.append(file)
                        labels.append(sample)

        bamlist.sort()
        labels.sort()
        bamlist = " ".join(bamlist)
        labels = " ".join(labels)
        # print(labels)
        #multiBamSummary_BEDfile(marker=marker, BEDfile=mm10_tss_up10kb, bamlist=bamlist, labels=labels, state="up",region="10kb")
        multiBamSummary_BEDfile_promoter(marker=marker, BEDfile=mm10_tss_1500bp, bamlist=bamlist, labels=labels)
        #multiBamSummary_BEDfile(marker=marker, BEDfile=mm10_tss_down10kb, bamlist=bamlist, labels=labels, state="down", region="10kb")
        #multiBamSummary_BEDfile(marker=marker, BEDfile=mm10_tss_up100kb, bamlist=bamlist, labels=labels, state="up",region="100kb")
        #multiBamSummary_BEDfile(marker=marker, BEDfile=mm10_tss_down100kb, bamlist=bamlist, labels=labels, state="down",region="100kb")

    # H3K27ac
    #multiBamSummary_BEDfile(marker=Markers[0],BEDfile=mm10_tss_2kb, bamlist=bamlist, labels=labels)
    # H3K4me1
    #multiBamSummary_BEDfile(marker=Markers[1], BEDfile=mm10_tss_2kb, bamlist=bamlist, labels=labels)
    # H3K4me3
    #multiBamSummary_BEDfile(marker=Markers[2], BEDfile=mm10_tss_2kb, bamlist=bamlist, labels=labels)
    # H3K27me3
    #multiBamSummary_BEDfile(marker=Markers[3], BEDfile=mm10_tss_2kb, bamlist=bamlist, labels=labels)
    # H3K9me3
    #multiBamSummary_BEDfile(marker=Markers[4], BEDfile=mm10_tss_2kb, bamlist=bamlist, labels=labels)

states = ["up","down"]
conditions = ["10kb","100kb"]
if step < 2:
    for marker in Markers:
        # test up10kb
        # state = states[1]
        # condition = conditions[0]
        # print(state,condition)
        # outCount = os.path.join(tfDir, marker + "_outRawCounts_%s_%s.txt" % (state, condition))
        # regionfile = os.path.join(tfDir, "mm10_tss_%s%s.txt" % (state, condition))
        # name = os.path.join(tfDir, marker + "_rpm_%s_%s.txt" % (state, condition))
        # print(outCount,regionfile)
        # rpm(marker=marker, regionfile=regionfile, outCount=outCount,name=name)

        # for every peak
        for state in states:
            for condition in conditions:
                outCount = os.path.join(tfDir, marker + "_outRawCounts_%s_%s.txt" % (state, condition))
                regionfile = os.path.join(tfDir,"mm10_TSS_%s%s.txt" % (state,condition))
                name = os.path.join(tfDir, marker + "_rpm_%s_%s.txt" % (state, condition))
                rpm(marker=marker, regionfile= regionfile, outCount=outCount, name=name)

# for tss promoter up and down 1500bp
if step < 3:
    for marker in Markers:
        outCount = os.path.join(tfDir, marker + "_outRawCounts_1500bp.txt")
        name = os.path.join(tfDir, marker + "_rpm_1500bp.txt")
        rpm(marker=marker, regionfile=mm10_tss_1500bp, outCount=outCount, name=name)



if step > 4:
    for marker in Markers:
        for condition in conditions:
            addfile(marker,condition=condition)



if step < 5:
    for marker in Markers:
        outCount = os.path.join(tfDir,marker + "_outRawCounts_1500bp.txt")
        rpm(marker=marker, regionfile=mm10_tss_1500bp, outCount=outCount, condition="1500bp")

    # for root, dirs, files in os.walk(Dir):
    #     for file in files:
    #         if "outRawCounts_1500bp.txt" in file:
    #             marker = file.split("_")[0]
    #             outCount = os.path.join(root,file)






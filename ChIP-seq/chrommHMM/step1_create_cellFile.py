#!/usr/bin/env python3
#######################################################################
# File Name: chromhmm.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Wed 29 Apr 2020 12:03:31 PM CST
# Description: 
# History: 
#######################################################################
import os,sys

bamdir = "/data3/zhaochen/project/colon_cancer/chip/every_sample_bam"
beddir = "/data3/zhaochen/project/colon_cancer/colon_chip/samtools/bed"
statedir = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/bed/state"
coorddir = "/data3/zhaochen/software/ChromHMM/COORDS/mm10"

def create_cellFileTable_bam(outputFile):
    filelist = os.listdir(bamdir)
    output = open(outputFile,"w")
    for file in filelist:
        array = file.split("-")
        cellType = array[0]
        mark = array[2].split("_")[0]
        if array[2].split(".")[-1] == "bai":
            continue
        control = cellType + "-" + array[1] + "-" + "Input_combined_1_paired.merged.nodup.bam"
        if mark == "Input":
            continue
        elif mark == "H3K9me2":
            continue
        else:
            output.write("%s\t%s\t%s\t%s\n" % (cellType,mark,file,control))
    output.close()

def create_cellFileTable_bed(outputFile):
    filelist = os.listdir(beddir)
    output = open(outputFile, "w")
    for file in filelist:
        array = file.split("-")
        cellType = array[0]
        mark = array[2].split(".")[0]
        control = cellType + "-" + array[1] + "-" + "Input.bed"
        if mark == "Input":
            continue
        elif mark == "H3K9me2":
            continue
        else:
            output.write("%s\t%s\t%s\t%s\n" % (cellType, mark, file, control))
    output.close()

def binary(chrlength):
    cellmarkfiletable = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/cellMarkFile_all_sample.txt"
    outputdir = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/binary"
    cmd = "java -mx4000M -jar /data3/zhaochen/anaconda3/pkgs/chromhmm-1.21-0/share/chromhmm-1.21-0/ChromHMM.jar BinarizeBam -b 200 %s %s %s %s" % (chrlength,bamdir,cellmarkfiletable,outputdir)
    os.system(cmd)
    return cmd

def BinarizeBed(chrlength):
    bed_cellmarkfile = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/bed/bed_cellMarkFile_all_sample.txt"
    outputdir = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/bed/binary"
    cmd = "java -mx4000M -jar /data3/zhaochen/anaconda3/pkgs/chromhmm-1.21-0/share/chromhmm-1.21-0/ChromHMM.jar BinarizeBed " \
          "%s %s %s %s" % (chrlength,beddir,bed_cellmarkfile,outputdir)
    os.system(cmd)
    return cmd

def LearnModel(chrlength):
    inputdir = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/bed/binary"
    outputdir = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/bed/state"
    cmd = "java -mx4000M -jar /data3/zhaochen/anaconda3/pkgs/chromhmm-1.21-0/share/chromhmm-1.21-0/ChromHMM.jar LearnModel -b 200 " \
          "-p 0 %s %s 13 mm10" % (inputdir,outputdir)
    os.system(cmd)
    #print(cmd)
    return cmd

def CompareModels():
    cmd = "java -mx4000M -jar /data3/zhaochen/anaconda3/pkgs/chromhmm-1.21-0/share/chromhmm-1.21-0/ChromHMM.jar CompareModels " \
          "referencemodel %s "

def OverlapEnrichment(inputsegment,sample):
    outfileprefix = os.path.join(statedir,sample)
    cmd = "java -mx4000M -jar /data3/zhaochen/anaconda3/pkgs/chromhmm-1.21-0/share/chromhmm-1.21-0/ChromHMM.jar OverlapEnrichment " \
          "-b 200 -t %s %s %s %s" % (sample,inputsegment,coorddir,outfileprefix)
    os.system(cmd)
    #print(cmd)
    return cmd

if __name__ == "__main__":
    step = 0
    chrlength = "/data3/zhaochen/reference/mm10/mm10.chrom.sizes"
    if step > 1:
        create_cellFileTable_bam(outputFile="cellMarkFile_all_sample.txt")


    if step > 2:
        binary(chrlength=chrlength)
        LearnModel(chrlength=chrlength)
        #CompareModels()

    if step > 3:
        create_cellFileTable_bed(outputFile="./bed/bed_cellMarkFile_all_sample.txt")
        BinarizeBed(chrlength=chrlength)
        LearnModel(chrlength=chrlength)

    if step < 4:
        for root,dirs,files in os.walk(statedir):
            for file in files:
                if "segments.bed" in file:
                    sample = file.split(".")[0]
                    file = os.path.join(root,file)
                    OverlapEnrichment(inputsegment=file,sample=sample)

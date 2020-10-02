#!/usr/bin/env python3
#######################################################################
# File Name: homer_findPeaks.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Wed 25 Dec 2019 04:46:14 PM CST
# Description: 
# History: 
#######################################################################
import os, sys
import re
import pandas as pd
from subprocess import *

find_TF = "/data3/zhaochen/project/colon_cancer/colon_chip/enhancer/find_TF"
pool_bed = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bed"
macs2 = "/data3/zhaochen/project/colon_cancer/colon_chip/macs2"


def run(cmd):
    p = Popen(cmd, shell=True)
    p.wait()
    return p


class DEEhancer():

    def __init__(self):
        self.rpmDir = "/data3/zhaochen/project/colon_cancer/colon_chip/enhancer/find_TF/rpm"
        self.DEEhancer = "/data3/zhaochen/project/colon_cancer/colon_chip/enhancer/find_TF/DEEhancer"

    def select_rm2kb(self, file, select):
        columns = ['PeakID', 'Chr', 'Start', 'End', 'Strand', 'Peak Score', \
                   'Focus Ratio/Region Size', 'Annotation', 'Detailed Annotation', 'Distance to TSS', \
                   'Nearest PromoterID', 'Entrez ID', 'Nearest Unigene', 'Nearest Refseq', 'Nearest Ensembl', \
                   'Gene Name', 'Gene Alias', 'Gene Description', 'Gene Type']
        df = pd.read_csv(file, header=0, sep="\t")
        df.columns = columns
        data = df[(df["Distance to TSS"] < -2000) | (df["Distance to TSS"] > 2000)]
        data = data[["Chr", "Start", "End", "PeakID"]]
        data["Start"] = df["Start"] - 1
        data.to_csv(select, header=0, sep="\t", index=False)
        return data

    def sort(self, file, sortBed):
        cmd = "sort -k1,1 -k2,2n %s > %s" % (file, sortBed)
        return cmd

    def bedtools_merge(self, sortBed, mergeBed):
        cmd = "bedtools merge -i %s > %s" % (sortBed, mergeBed)
        return cmd

    def intersect(self, mergeBed, treatBed, treatcount):
        cmd = "bedtools intersect -c -wa -a %s -b %s > %s" % (mergeBed, treatBed, treatcount)
        return cmd

    def rpm(self, sample, control, treatment):
        dict = {"ctrl-H3K27ac": 106.890540, "2weeks-H3K27ac": 121.182308, "4weeks-H3K27ac": 113.026669,
                "7weeks-H3K27ac": 88.780726, "10weeks-H3K27ac": 110.733206}
        df1 = pd.read_csv(control, sep="\t", header=None, names=["chr", "start", "end", "control"])
        df2 = pd.read_csv(treatment, sep="\t", header=None, names=["chr", "start", "end", "treat"])
        df = pd.merge(df1, df2, on=["chr", "start", "end"])
        columns = ["chr", "start", "end", "PeakID", "Not used", "strand", "treat", "control", "RPM"]
        if sample in dict:
            treat = dict[sample]
            df["RPM"] = df["treat"] * dict["ctrl-H3K27ac"] / (df["control"] * treat)
            df["PeakID"] = [sample + "_Input_peak_" + str(i) for i in range(1, len(df["chr"]) +1)]
            df["Not used"] =  "."
            df["strand"] = "+"
            dfname = os.path.join(d.DEEhancer, "total_" + sample + "_DEEnhancers.bed")
            df.to_csv(dfname, sep="\t", header=0, index=False,columns=columns)

            data = df[df["RPM"] > 2]
            data["RPM"] = data["treat"] * dict["ctrl-H3K27ac"] / (data["control"] * treat)
            data.loc[:,"PeakID"] = [sample + "_Input_peak_" + str(i) for i in range(1, len(data["chr"]) + 1)]
            data["Not used"] = "."
            data["strand"] = "+"
            filename = os.path.join(d.DEEhancer, "Up_" + sample + "_DEEnhancers.bed")
            data.to_csv(filename, sep="\t", header=0, index=False, columns=columns)

            down = df[df["RPM"] < 0.5]
            down["PeakID"] = [sample + "_Input_peak_" + str(i) for i in range(1,len(down["chr"]) +1)]
            down["Not used"] = "."
            down["strand"] = "+"
            downname = os.path.join(d.DEEhancer, "Down_" + sample + "_DEEnhancers.bed")
            down.to_csv(downname, sep="\t", header=0, index=False, columns=columns)
        return data

    def annotatePeak(self, file, annotateFile):
        cmd = "annotatePeak.pl mm10 %s > %s" % (file, annotateFile)
        return cmd


class without_NFR():

    def __init__(self):
        self.motifs = "/data3/zhaochen/project/colon_cancer/colon_chip/enhancer/find_TF/motifs"
        self.total_motifs = "/data3/zhaochen/project/colon_cancer/colon_chip/enhancer/find_TF/total_motifs"

    def findMotifsGenome(self, file, outputDir, size):
        cmd = "findMotifsGenome.pl %s mm10 %s -size %s -p 21" % (file, outputDir, size)
        return cmd


class withNfr():

    def __init__(self):
        self.nfr_motifs = "/data3/zhaochen/project/colon_cancer/colon_chip/enhancer/find_TF/nfr_motifs"
        self.samtools = "/data3/zhaochen/project/colon_cancer/colon_chip/samtools"
        self.makeTagDir = "/data3/zhaochen/project/colon_cancer/colon_chip/enhancer/find_TF/nfr_motifs/step1_makeTagDirectory"
        self.outputNfr = "/data3/zhaochen/project/colon_cancer/colon_chip/enhancer/find_TF/nfr_motifs/step2_outputNfr"

    def makeTagDirectory(self, tagDir, alignmentBam):
        cmd = "makeTagDirectory %s %s" % (tagDir, alignmentBam)
        return cmd

    def findPeaks(self, tagDir, controlDir, output, size):
        cmd = "findPeaks %s -i %s -o %s -nfr -size %s" % (tagDir, controlDir, output, size)
        return cmd

    def add_100(self, file, output):
        cmd = """pos2bed.pl %s | sed '1,40d' | sort -k1,1 -k2,2n | awk '{print $1,$2-100,$3+100,$4,$5,$6}' OFS="\t" > %s""" % (
        file, output)
        os.system(cmd)
        return cmd


if __name__ == "__main__":
    step = 6
    d = DEEhancer()
    if step < 1:
        for root, dirs, files in os.walk(macs2):
            for file in files:
                if "annotation.txt" in file:
                    if "H3K27ac" in file:
                        sample = file.split("_")[0]
                        file = os.path.join(root, file)
                        select = os.path.join(d.rpm, sample + "/" + sample + "_rm2kb.bed")
                        cmd = d.select_rm2kb(file=file, select=select)
                        run(cmd)
    #                        os.system("sed -i '1d' %s" % select)

    if step < 2:
        control = "/data3/zhaochen/project/colon_cancer/colon_chip/enhancer/find_TF/rpm/ctrl-H3K27ac/ctrl-H3K27ac_rm2kb.bed"
        for root, dirs, files in os.walk(d.rpmDir):
            for file in files:
                if "_rm2kb.bed" in file:
                    sample = file.split("_")[0]
                    if "ctrl" not in file:
                        treatment = os.path.join(root, file)
                        poolBed = os.path.join(root, sample + "_ctrl_pool.broadPeak")
                        mergeBed = os.path.join(root, sample + "_ctrl_merge.broadPeak")

                        #                        os.system("cat %s %s > %s" % (treatment,control,poolBed))
                        sortBed = os.path.join(root, sample + "_ctrl_sort.broadPeak")
                        #                        cmd = d.sort(file=poolBed,sortBed=sortBed)
                        #                        run(cmd)
                        #                        cmd1 = d.bedtools_merge(sortBed=sortBed,mergeBed=mergeBed)
                        #                        run(cmd1)
                        treatcount = os.path.join(root, sample + "_count.txt")
                        treatBed = os.path.join(pool_bed, sample + "_pool.bed")
                        #                        cmd2 = d.intersect(treatBed=treatBed, mergeBed=mergeBed, treatcount=treatcount)
                        #                        run(cmd2)
                        ctrlBed = os.path.join(pool_bed, "ctrl-H3K27ac_pool.bed")
                        ctrlcount = os.path.join(root, sample + "_control_count.txt")
    #                        cmd3 = d.intersect(mergeBed=mergeBed, treatBed=ctrlBed, treatcount=ctrlcount)
    #                        run(cmd3)
    #                        os.remove(poolBed)
    #                        os.remove(mergeBed)

    if step < 3:
        for root,dirs,files in os.walk(d.rpmDir):
            for file in files:
                if "_control_count.txt" in file:
                    ctrl = os.path.join(root,file)
                    sample = file.split("_")[0]
                    treatment = os.path.join(root, sample + "_count.txt")
                    cmd = d.rpm(sample=sample,control=ctrl,treatment=treatment)
                    run(cmd)

    w = without_NFR()
    if step < 4:
        samples = ["2weeks-H3K27ac", "4weeks-H3K27ac", "7weeks-H3K27ac", "10weeks-H3K27ac"]
        models = ["Up_", "total_"]
        for sample in samples:
            for model in models:
                file = os.path.join(d.DEEhancer, model + sample + "_DEEnhancers.bed")
                outputDir = os.path.join(w.motifs, model.split("_")[0] + "-" + sample)
                if not os.path.exists(outputDir):
                    os.makedirs(outputDir)
                cmd = w.findMotifsGenome(file=file, outputDir=outputDir, size=600)
                run(cmd)

    n = withNfr()
    if step < 6:
        times = ["ctrl", "2weeks", "4weeks", "7weeks", "10weeks"]
        for time in times:
            samples = []
            for i in range(1, 4):
                sample = "%s/%s-%s-H3K27ac_rmdup.bam" % (n.samtools, time, i)
                samples.append(sample)
            tagDir = os.path.join(n.makeTagDir, time + "-H3K27ac")
            alignmentBam = " ".join(samples)
            cmd = n.makeTagDirectory(tagDir=tagDir, alignmentBam=alignmentBam)
            run(cmd)

    if step < 6:
        samples = ["2weeks-H3K27ac", "4weeks-H3K27ac", "7weeks-H3K27ac", "10weeks-H3K27ac"]
        controlDir = os.path.join(n.makeTagDir, "ctrl-H3K27ac")
        for sample in samples:
            tagDir = os.path.join(n.makeTagDir, sample)
            output = os.path.join(tagDir, sample + "_peak_200.txt")
            cmd = n.findPeaks(tagDir=tagDir, controlDir=controlDir, output=output, size=200)
            run(cmd)

    if step < 7:
        for root, dirs, files in os.walk(n.makeTagDir):
            for file in files:
                if "_peak_200.txt" in file:
                    sample = file.split("_")[0]
                    file = os.path.join(root, file)
                    output = file.replace("_200.txt", "_400.bed")
#                    cmd1 = n.add_100(file=file,output=output)
#                    run(cmd1)
                    DEEhancer_Bed = os.path.join(d.DEEhancer, "Up_" + sample + "_DEEnhancers.bed")
                    treatcount = os.path.join(d.DEEhancer, sample + "_nfr_400.bed")
#                    cmd2 = d.intersect(mergeBed=output, treatBed=DEEhancer_Bed, treatcount=treatcount)
#                    run(cmd2)
#                    outputDir = os.path.join(n.outputNfr,sample)
                    outputDir = os.path.join("/data3/zhaochen/project/colon_cancer/colon_chip/enhancer/find_TF/nfr_motifs/step3_outputNfr_200",sample)
                    if not os.path.exists(outputDir):
                        os.makedirs(outputDir)
                    cmd3 = w.findMotifsGenome(file=treatcount,outputDir=outputDir,size=200)
                    run(cmd3)

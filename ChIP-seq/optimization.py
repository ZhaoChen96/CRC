#!/usr/bin/env python3
#######################################################################
# File Name: optimization.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Sun 08 Dec 2019 09:08:51 PM CST
# Description: 
# History: 
#######################################################################
import os,re
import sys
from subprocess import *

class quality_control():

    def __init__(self):
        self.raw_data = "/backup_raw/AOM-DSS/ChIP-seq"
        self.fastqc_dir = "/data3/zhaochen/project/colon_cancer/colon_chip/fastqc"
        self.clean_dir = "/data3/zhaochen/project/colon_cancer/colon_chip/clean"
        self.reference_genome = "/data4/genome_index/mm10/tophat2/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa"
        self.aligned_dir = "/data3/zhaochen/project/colon_cancer/colon_chip/aligned"
        self.samtools_dir = "/data3/zhaochen/project/colon_cancer/colon_chip/samtools"

    def run(self,cmd,wkdir=None):
        sys.stderr.write("Running: %s \n" % cmd)
        p = Popen(cmd,shell=True,cwd=wkdir)
        p.wait()
        return p

    def fastqc(self,sample):
        cmd = "fastqc -t 7 -o %s %s" % (qc.fastqc_dir,sample)
        return cmd

    def multiqc(self):
        multiqc_dir = "/data3/zhaochen/project/colon_cancer/colon_chip/fastqc/multiqc"
        cmd = "multiqc %s -o %s" % (qc.fastqc_dir, multiqc_dir)
        return cmd

    def cutadapt_PE(self,sample_R1,sample_R2):
        """
        PE mains pair end
        SE mains singal end
        """
        adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        name_R1 = os.path.basename(sample_R1)
        name_R2 = os.path.basename(sample_R2)
        output_R1 = os.path.join(qc.clean_dir,name_R1.replace("_R1.fastq.gz","_R1.fq.gz"))
        output_R2 = os.path.join(qc.clean_dir,name_R2.replace("_R2.fastq.gz","_R2.fq.gz"))
        cmd = "cutadapt -j 14 -a %s -A %s -u 3 -u -87 -U 3 -U -87 -m 30 -o %s -p %s %s %s" % \
              (adapter_1,adapter_2,output_R1,output_R2,sample_R1,sample_R2)
        return cmd

    def bwa(self,ref,rg,clean_R1,clean_R2,sample,threads=21):
        sam = os.path.join(qc.aligned_dir,"10weeks/" + sample + ".sam")
        cmd = "bwa mem -t %s %s -R '@RG\\tID:1\\tPL:ILLUMINA\\tSM:%s' %s %s > %s" %\
              (threads,ref,rg,clean_R1,clean_R2,sam)
        print(cmd)
        os.system(cmd)
        return cmd

    def samToview(self,sam_dir,sample):
        sam = os.path.join(sam_dir, sample + ".sam")
        view_bam = os.path.join(qc.samtools_dir,sample + "_view.bam")
        cmd = "samtools view -@ 7 -bS %s > %s" % (sam,view_bam)
        print(cmd)
        os.system(cmd)
        return cmd

    def viewTosort(self,sample):
        view_bam = os.path.join(qc.samtools_dir,sample + "_view.bam")
        sort_bam = os.path.join(qc.samtools_dir,sample + "_sort.bam")
        cmd = "samtools sort -@ 7 -O BAM -o %s %s" % (sort_bam,view_bam)
        print(cmd)
        os.system(cmd)
        return cmd

    def sortTormdup(self,sample):
        sort_bam = os.path.join(qc.samtools_dir, sample + "_sort.bam")
        rmdup_bam = os.path.join(qc.samtools_dir, sample + "_rmdup.bam")
        cmd = "samtools rmdup %s %s" % (sort_bam,rmdup_bam)
        print(cmd)
        os.system(cmd)
        return cmd

    def index(self,sample):
        rmdup_bam = os.path.join(qc.samtools_dir, sample + "_rmdup.bam")
        cmd = "samtools index -@ 7 -b %s" % (rmdup_bam)
        print(cmd)
        os.system(cmd)
        return cmd

    def bamTobed(self,sample):
        rmdup_bam = os.path.join(qc.samtools_dir, sample + "_rmdup.bam")
        bed = os.path.join(qc.samtools_dir, sample + ".bed")
        cmd = "bedtools bamtobed -i %s > %s" % (rmdup_bam,bed)
        print(cmd)
        os.system(cmd)
        return cmd

    def flagstat(self,sample):
        view_bam = os.path.join(qc.samtools_dir, sample + "_view.bam")
        flagstat = os.path.join(qc.samtools_dir, sample + "_flagstat.txt")
        cmd = "samtools flagstat -@ 7 %s > %s" % (view_bam, flagstat)
        print(cmd)
        os.system(cmd)
        return cmd
        

if __name__=="__main__":
    qc = quality_control()
    step = 1
    if step < 1:
        for root,dirs,files in os.walk(qc.raw_data):
            for file in files:
                sample = os.path.join(root, file)
                cmd = qc.fastqc(sample=sample)
                qc.run(cmd=cmd)
                qc.run(cmd=qc.multiqc())

    if step < 2:
        raw_data_10 = "/backup_raw/AOM-DSS/ChIP-seq/10_weeks"

        samples = ["10weeks-3-H3K27ac","10weeks-1-H3K4me3","10weeks-1-H3K9me2",
                   "10weeks-2-H3K4me1","10weeks-2-H3K27me3","10weeks-2-H3K4me3"]
        for sample in samples:
            sample_R1 = "%s_combined_R1.fastq.gz" % (sample)
            sample_R1 = os.path.join(raw_data_10,sample_R1)
            sample_R2 = "%s_combined_R2.fastq.gz" % (sample)
            sample_R2 = os.path.join(raw_data_10, sample_R2)
            cmd = qc.cutadapt_PE(sample_R1=sample_R1,sample_R2=sample_R2)
            qc.run(cmd=cmd)
            clean_R1 = os.path.join(qc.clean_dir,sample + "_combined_R1.fq.gz")
            clean_R2 = os.path.join(qc.clean_dir,sample + "_combined_R2.fq.gz")
            qc.bwa(ref=qc.reference_genome,rg=sample,clean_R1=clean_R1,clean_R2=clean_R2,sample=sample)
            sam_dir = "/data3/zhaochen/project/colon_cancer/colon_chip/aligned/10weeks"
#            qc.samToview(sam_dir=sam_dir,sample=sample)
            qc.viewTosort(sample=sample)
            qc.sortTormdup(sample=sample)
            qc.index(sample=sample)
            qc.bamTobed(sample=sample)
#            qc.flagstat(sample=sample)



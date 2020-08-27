#!/usr/bin/env python3
#######################################################################
# File Name: chip-seq_pipeline.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Sat 30 Nov 2019 05:00:09 PM CST
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
        self.bw_dir = "/data3/zhaochen/project/colon_cancer/colon_chip/bigwig"
        self.markduplicate_dir = "/data3/zhaochen/project/colon_cancer/colon_chip/step11_markduplicates"

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
        cmd = "cutadapt -j 21 -a %s -A %s -u 4 -u -35 -U 4 -U -35 -m 30 -o %s -p %s %s %s" % \
              (adapter_1,adapter_2,output_R1,output_R2,sample_R1,sample_R2)
        return cmd

    def bwa(self,ref,rg,clean_R1,clean_R2,sample,threads=21):
        sam = os.path.join(qc.aligned_dir,"10weeks/" + sample + ".sam")
        cmd = "bwa mem -t %s %s -R '@RG\\tID:1\\tPL:ILLUMINA\\tSM:%s' %s %s > %s" %\
              (threads,ref,rg,clean_R1,clean_R2,sam)
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

    def bamTobw(self,sample):
        """
        compare with input
        :param sample:
        :return:
        """
#        os.makedirs(qc.bw_dir)
        rmdup_bam = os.path.join(qc.samtools_dir, sample + "_rmdup.bam")
        output = os.path.join(qc.bw_dir,sample + ".bw")
        cmd = "bamCoverage -b %s --normalizeUsing RPKM -bs 100 --operation ratio -p 21 -o %s" % (rmdup_bam,output)
        print(cmd)
        os.system(cmd)
        return cmd

    def flagstat(self,sample):
        rmdup_bam = os.path.join(qc.samtools_dir,sample + "_rmdup.bam")
        flagstat = os.path.join(qc.samtools_dir,sample + "_flagstat.txt")
        cmd = "samtools flagstat -@ 7 %s > %s" % (rmdup_bam,flagstat)
        return cmd

    def markduplicate(self,sample,input_bam,marked_bam):
        metrics_txt = os.path.join(qc.markduplicate_dir,sample + "_matrix.txt")
        cmd = "picard MarkDuplicates I=%s O=%s M=%s REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=queryname " % \
              (input_bam,marked_bam,metrics_txt)
        return cmd

    def unique_mapping(self,input_bam,sort_bam):
        cmd = "samtools view -h -q 1 -F 4 -F 256 %s |grep -v XA:Z | grep -v SA:Z | samtools view -@ 10 -Sb - | " \
              "samtools sort -@ 10 > %s" %(input_bam, sort_bam)
        return cmd

class peak_calling():

    def __init__(self):
        self.macs2_dir = "/data3/zhaochen/project/colon_cancer/colon_chip/macs2"
        self.pool_bam = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bam"
        self.pool_bed = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bed"
        self.bw_InputDir = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input"
        self.bw_Dir = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_withoutInput"

    def merge(self,merge_bam,input_bam):
        cmd = "samtools merge -@ 14 %s %s" % (merge_bam,input_bam)
        return cmd

    def mergeTosort(self,merge_bam,sort_bam):
        cmd = "samtools sort -@ 14 -o %s %s" % (sort_bam,merge_bam)
        return cmd

    def bamTobed(self,sort_bam,poolbed):
        cmd = "bedtools bamtobed -i %s > %s" % (sort_bam, poolbed)
        print(cmd)
        os.system(cmd)
        return cmd

    def bamTobw(self,pool_bam,pool_bw):
        """
        not compare with Input
        """
        cmd = "bamCoverage -b %s --normalizeUsing RPKM -bs 100 -e 147 -p 21 -o %s" % (pool_bam,pool_bw)
        print(cmd)
        return cmd

    def bamcompare(self,treat_bam,sample,Input_bam):
        """
        compare with Input -b1 treatment BAM file -b2 control BAM file
        """
        #pool_bw = os.path.join(pc.bw_InputDir,sample + "_ratio.bw")
        #cmd = "bamCompare -b1 %s -b2 %s --operation ratio -bs 100 -p 21 --normalizeUsing RPKM -e 147 -o %s " % \
        #      (treat_bam,Input_bam,pool_bw)

        pool_bw = os.path.join(pc.bw_InputDir,sample + "_bs1bp_ratio.bw")
        cmd = "bamCompare -b1 %s -b2 %s --operation ratio -bs 1 -p 21 -e 147 -o %s " % \
              (treat_bam,Input_bam,pool_bw)
        return cmd

    def callpeak_broad(self,sample,treat1_bam,treat2_bam,treat3_bam,ctrl1_bam,ctrl2_bam,ctrl3_bam):
        cmd = "macs2 callpeak -t %s %s %s -c %s %s %s -f BAM --outdir %s -n %s_Input -g mm --nomodel --keep-dup all -p 1E-9 \
        --broad --broad-cutoff 1E-9 --extsize 147" % (treat1_bam,treat2_bam,treat3_bam,ctrl1_bam,ctrl2_bam,ctrl3_bam,pc.macs2_dir,sample)
        return cmd

    def callpeak_narrow(self,sample):
        cmd = "macs2 callpeak -t %s -c %s -f BAM --outdir %s -n %s_Input -g mm --nomodel --keep-dup all -p 1E-5 " \
              "--extsize 147" % (treatment,control,pc.macs2_dir,sample)
        return cmd

if __name__=="__main__":
    qc = quality_control()
    step = 18
    if step < 1:
        for root,dirs,files in os.walk(qc.raw_data):
            for file in files:
                sample = os.path.join(root, file)
#                cmd = qc.fastqc(sample=sample)
#                qc.run(cmd=cmd)
                qc.run(cmd=qc.multiqc())

    if step < 2:
        samples = []
        raw_data_10 = "/backup_raw/AOM-DSS/ChIP-seq/10_weeks"
        raw_data_ctrl = "/backup_raw/AOM-DSS/ChIP-seq/control"
        raw_data_2 = "/backup_raw/AOM-DSS/ChIP-seq/2_weeks"
        raw_data_4 = "/backup_raw/AOM-DSS/ChIP-seq/4_weeks"
        raw_data_7 = "/backup_raw/AOM-DSS/ChIP-seq/7_weeks"
        for root,dirs,files in os.walk(raw_data_10):
            for file in files:
                fi = file.split("_")[0]
                if fi not in samples:
                    samples.append(fi)

        for sample in samples:
            sample_R1 = "%s_combined_R1.fastq.gz" % (sample)
            sample_R1 = os.path.join(raw_data_10,sample_R1)
            sample_R2 = "%s_combined_R2.fastq.gz" % (sample)
            sample_R2 = os.path.join(raw_data_10, sample_R2)
            cmd = qc.cutadapt_PE(sample_R1=sample_R1,sample_R2=sample_R2)
            qc.run(cmd=cmd)

    if step < 3:
        for root, dirs, files in os.walk(qc.clean_dir):
            for file in files:
                sample = os.path.join(root, file)
                cmd = qc.fastqc(sample=sample)
                qc.run(cmd=cmd)
                qc.run(cmd=qc.multiqc())

    if step < 4:
        samples_ctrl = []
        samples_2weeks = []
        samples_4weeks = []
        samples_7weeks = []
        samples_10weeks = []
        for file in os.listdir(qc.clean_dir):
            fi = file.split("_")[0]
            if "ctrl" in file:
                if fi not in samples_ctrl:
                    samples_ctrl.append(fi)
            elif "2weeks" in file:
                if fi not in samples_2weeks:
                    samples_2weeks.append(fi)
            elif "4weeks" in file:
                if fi not in samples_4weeks:
                    samples_4weeks.append(fi)
            elif "7weeks" in file:
                if fi not in samples_7weeks:
                    samples_7weeks.append(fi)
            elif "10weeks" in file:
                if fi not in samples_10weeks:
                    samples_10weeks.append(fi)

        for sample in samples_10weeks:
            clean_R1 = "%s_combined_R1.fq.gz" % (sample)
            clean_R1 = os.path.join(qc.clean_dir,clean_R1)
            clean_R2 = "%s_combined_R2.fq.gz" % (sample)
            clean_R2 = os.path.join(qc.clean_dir, clean_R2)

            for file in os.listdir(qc.aligned_dir):
                if not re.search(".sam",file):
                    cmd = qc.bwa(ref=qc.reference_genome,rg=sample,clean_R1=clean_R1,clean_R2=clean_R2,sample=sample)
                    qc.run(cmd=cmd)

    if step < 5:
#        samFiles = os.listdir(qc.aligned_dir)
#        if len(samFiles) == 0:
#            exit(1)

        aligned = "/data3/zhaochen/project/colon_cancer/colon_chip/aligned"
        for root,dirs,files in os.walk(aligned):
            for file in files:
                sample = file.split(".")[0]
                sam = os.path.join(aligned,sample + ".sam")
                view = os.path.join(qc.samtools_dir,sample + "_view.bam")
                sort = os.path.join(qc.samtools_dir,sample + "_sort.bam")
                if not re.search(".sam",file):
                    continue
                else:
                    cmd_view = qc.samToview(sample=sample,sam_dir=aligned_ctrl)
                    qc.run(cmd_view)
                    cmd_sort = qc.viewTosort(sample=sample)
                    qc.run(cmd_sort)
                    cmd_rmdup = qc.sortTormdup(sample=sample)
                    qc.run(cmd_rmdup)
                    cmd_index = qc.index(sample=sample)
                    os.system("rm %s" % (sam))
                    os.system("rm %s" % (view))
                    os.system("rm %s" % (sort))
                    qc.run(cmd_index)
                    cmd_bed = qc.bamTobed(sample=sample)
                    qc.run(cmd_bed)
                    cmd_bw = qc.bamTobw(sample=sample)
                    qc.run(cmd_bw)

    if step < 10:
        for root,dirs,files in os.walk(qc.samtools_dir):
            for file in files:
                if os.path.splitext(file)[1] == ".bam":
                    sample = file.split("_")[0]
                    cmd = qc.flagstat(sample=sample)
                    qc.run(cmd)

    if step < 11:
        for root,dirs,files in os.walk(qc.samtools_dir):
            for file in files:
                if os.path.splitext(file)[1] == ".bam":
                    sample = file.split("_")[0]
                    input_bam = os.path.join(qc.samtools_dir,file)
                    marked_bam = os.path.join(qc.markduplicate_dir,sample + "_mkdup.bam")
                    cmd = qc.markduplicate(input_bam=input_bam,marked_bam=marked_bam,sample=sample)
                    qc.run(cmd)

    if step < 12:
        for root,dirs,files in os.walk(qc.samtools_dir):
            for file in files:
                if os.path.splitext(file)[1] == ".bam":
                    sample = file.split("_")[0]
                    input_bam = os.path.join(root,file)
                    sort_bam = os.path.join(root,sample + "_mkdup_unique.bam")
                    cmd = qc.unique_mapping(input_bam=input_bam,sort_bam=sort_bam)
                    qc.run(cmd)

    pc = peak_calling()

    if step < 13:
        bam_dict = {}
        for root,dirs,files in os.walk(qc.samtools_dir):
            for file in files:
                if os.path.splitext(file)[1] == ".bam":
                    sample = file.split("_")[0]
                    type_feature = sample.split("-")[0] + "-" + sample.split("-")[2]
                    bam_file = os.path.join(root,file)
                    if type_feature in bam_dict:
                        bam_dict[type_feature].append(bam_file)
                    else:
                        bam_dict[type_feature] = [bam_file]

        for key in bam_dict.keys():
            merge_bam = os.path.join(pc.pool_bam,key + "_merge.bam")
            input_bam = " ".join(bam_dict[key])
#            cmd = pc.merge(input_bam=input_bam,merge_bam=merge_bam)
#            qc.run(cmd)
            sort_bam = merge_bam.replace("_merge.bam","_sort.bam")
#            cmd_sort = pc.mergeTosort(sort_bam=sort_bam,merge_bam=merge_bam)
#            qc.run(cmd_sort)
            os.system("samtools index -@ 14 -b %s" % merge_bam)
            poolbed = os.path.join(pc.pool_bed,key + "_pool.bed")
#            cmd_bamTobed=pc.bamTobed(sort_bam=sort_bam,poolbed=poolbed)
#            qc.run(cmd_bamTobed)

    if step < 14:
        for root,dirs,files in os.walk(pc.pool_bam):
            for file in files:
                if file.split("_")[1] == "sort.bam":
                    pool_bam = os.path.join(root,file)
                    pool_bw = os.path.join(pc.bw_Dir,file.replace("_sort.bam",".bw"))
                    cmd = pc.bamTobw(pool_bam=pool_bam,pool_bw=pool_bw)
                    qc.run(cmd)

    if step > 15:
        modifications = ["H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K9me2","H3K9me3"]
        times = ["ctrl","2weeks","4weeks","7weeks","10weeks"]
        for time in times:
            for modification in modifications:
                sample = "%s-%s" % (time,modification)
                Input_bam = os.path.join(pc.pool_bam,time + "-Input_merge.bam")
                treat_bam = os.path.join(pc.pool_bam,sample + "_merge.bam")
                cmd = pc.bamcompare(treat_bam=treat_bam,sample=sample,Input_bam=Input_bam)
                qc.run(cmd)

    if step < 16:
        modifications = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "H3K9me2", "H3K9me3"]
        times = ["ctrl", "2weeks", "4weeks", "7weeks", "10weeks"]
        for time in times:
            for modification in modifications:
                sample = "%s-%s" % (time,modification)
                sample1 = "%s-1-%s" % (time,modification)
                sample2 = "%s-2-%s" % (time, modification)
                sample3 = "%s-3-%s" % (time, modification)
                t1 = os.path.join(qc.samtools_dir,sample1 + "_rmdup.bam")
                t2 = os.path.join(qc.samtools_dir, sample2 + "_rmdup.bam")
                t3 = os.path.join(qc.samtools_dir, sample3 + "_rmdup.bam")
                c1 = os.path.join(qc.samtools_dir,time + "-1-Input_rmdup.bam")
                c2 = os.path.join(qc.samtools_dir,time + "-2-Input_rmdup.bam")
                c3 = os.path.join(qc.samtools_dir, time + "-3-Input_rmdup.bam")
                cmd = pc.callpeak_broad(sample=sample,treat1_bam=t1,treat2_bam=t2,treat3_bam=t3,ctrl1_bam=c1,ctrl2_bam=c2,ctrl3_bam=c3)
                qc.run(cmd)







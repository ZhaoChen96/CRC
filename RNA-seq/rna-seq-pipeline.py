#!/usr/bin/env python3
#######################################################################
# File Name: cutadapter.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Mon 07 Oct 2019 10:52:41 AM CST
# Description: 
# History: 
#######################################################################
import os,sys
import numpy as np
from subprocess import *
import re


class work_flow():

    def __init__(self):
        self.rawdataDir = "/data3/zhaochen/Raw_data/zhangbo/rna-seq"
        self.outputdir = os.path.abspath("/data3/zhaochen/project/zhangbo/rna-seq/fastqc")
        self.cleandir = "/data3/zhaochen/project/zhangbo/rna-seq/clean"
        self.samtools_dir = "/data3/zhaochen/project/zhangbo/rna-seq/samtools"
        self.cufflinks_dir = "/data3/zhaochen/project/zhangbo/rna-seq/cufflinks"
        self.align_dir = "/data3/zhaochen/project/zhangbo/rna-seq/align"
        self.featureCount_dir = "/data3/zhaochen/project/zhangbo/rna-seq/samtools/featureCount"
        self.multiqc_dir = "/data3/zhaochen/project/zhangbo/rna-seq/fastqc/multiqc"

    def run(self,cmd):
        sys.stderr.write("Running: %s \n" % cmd)
        p = Popen(cmd, shell=True)
        p.wait()
        return p

    def fastqc(self,sample):
        cmd = "fastqc -t 7 -o %s %s" % (wf.outputdir,sample)
        os.system(cmd)
        return cmd

    def multiqc(self):
        cmd = "multiqc %s -o %s" % (wf.cleandir, wf.multiqc_dir)
        os.system(cmd)
        return cmd

#   def cutadapt_SE(self,sample):
#        adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
#       adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
#       cleandir = "/data3/zhaochen/NCBI/breast_cell_lines/rna-seq/clean"
#       output = os.path.join(cleandir,sample + "_fastq.gz")
#       files = os.path.join(fastq_dir,sample + ".fastq.gz")
#       cmd = "cutadapt -j 10 -a %s -u 10 -m 20 -o %s %s" % (adapter_1,output,files)
#       os.system(cmd)
#       return cmd 

    def cutadapt_PE(self,sample_R1,sample_R2):
        """
        PE mains pair end
        SE mains singal end
        """
        adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        name_R1 = os.path.basename(sample_R1)
        name_R2 = os.path.basename(sample_R2)
        output_R1 = os.path.join(wf.cleandir,name_R1.replace("_R1.fastq.gz","_R1.fq.gz"))
        output_R2 = os.path.join(wf.cleandir,name_R2.replace("_R2.fastq.gz","_R2.fq.gz"))
        cmd = "cutadapt -j 21 -a %s -A %s -u 10 -u -10 -U 10 -U -10 -m 30 -o %s -p %s %s %s" % \
              (adapter_1,adapter_2,output_R1,output_R2,sample_R1,sample_R2)
        return cmd

    def tophat(self,sample,clean_R1,clean_R2,species):
        if species == "mouse":
            gtf = "/data3/zhaochen/reference/mm10/mm10_annotation.gtf"
            tophat_index = "/data4/genome_index/mm10/tophat2/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
        elif species == "human":
            gtf = "/data3/zhaochen/reference/hg19/gencode.v32lift37.annotation.gtf"
            tophat_index = "/data4/genome_index/hg19/tophat2/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
        cmd = "/data3/zhaochen/anaconda3/envs/py2/bin/tophat -p 21 -o %s/%s_tophat -G %s %s %s %s" \
              % (self.align_dir,sample,gtf,tophat_index,clean_R1,clean_R2)
        os.system(cmd)
        return cmd

    def hisat2(self,sample,clean_R1,clean_R2):
        #hisat2_index = "/data3/zhaochen/reference/mm10/mm10_hisat_index/genome
        output = os.path.join(self.align_dir,sample + "-RNA.sam")
        cmd = "hisat2 -q --phred33 -p 21 -x %s -1 %s %s -S %s" % (hisat2_index,clean_R1,clean_R2,output)
        os.system(cmd)
        return cmd

    def STAR_build_index_human(self):
        cmd = "STAR --runMode genomeGenerate --runThreadN 21 --genomeDir /data3/zhaochen/reference/hg19/STAR_genome_index/star_genome " \
              "--genomeFastaFiles /data3/zhaochen/reference/hg19/STAR_genome_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa " \
              "--sjdbGTFfile /data3/zhaochen/reference/hg19/STAR_genome_index/Homo_sapiens.GRCh38.99.gtf --sjdbOverhang 149 --limitGenomeGenerateRAM 35000000000"
        os.system(cmd)
        return cmd

    def STAR_bulid_index_mouse(self):
        cmd = "STAR --runMode genomeGenerate --runThreadN 21 --genomeDir /data3/zhaochen/reference/mm10/STAR_genome_index/star_genome " \
              "--genomeFastaFiles /data3/zhaochen/reference/mm10/STAR_genome_index/Mus_musculus.GRCm38.dna.toplevel.fa " \
              "--sjdbGTFfile /data3/zhaochen/reference/mm10/STAR_genome_index/Mus_musculus.GRCm38.94.gtf --sjdbOverhang 149 --limitGenomeGenerateRAM 35000000000"
        os.system(cmd)
        return cmd

    def STAR(self,sample,clean_R1,clean_R2):
        prefix = os.path.join(wf.align_dir,sample)
        cmd = "STAR --genomeDir /data3/zhaochen/reference/mm10/STAR_genome_index/star_genome " \
              "--readFilesIn %s %s --readFilesCommand zcat --sjdbGTFfile /data3/zhaochen/reference/mm10/STAR_genome_index/Mus_musculus.GRCm38.94.gtf " \
              "--runThreadN 21 --outFileNamePrefix %s --outSAMtype BAM Unsorted" % (clean_R1,clean_R2,prefix)
        os.system(cmd)
        return cmd

    def samtools(self,sample):
        #bamdir = os.path.join(self.align_dir,sample + "_tophat")
        view_bam = os.path.join(wf.align_dir,'%sAligned.out.bam' % sample)
        sort_bam = os.path.join(wf.samtools_dir,sample + '_sorted.bam')
        bed = os.path.join(wf.samtools_dir,sample + ".bed")
        tmp = "/data3/zhaochen/project/zhangbo/rna-seq/samtools"
        os.system("samtools sort -@ 7 -O BAM -T %s -o %s %s" % (tmp,sort_bam,view_bam))
        os.system("samtools index -@ 7 -b %s" % (sort_bam))
        #os.system("bedtools bamtobed -i %s > %s" % (sort_bam,bed))
        os.system("rm %s" % view_bam)


    def cufflinks(self,prefix,sort_bam):
        cmd_cufflinks = "cufflinks -p 21 -o %s/%s %s" % (self.cufflinks_dir,prefix,sort_bam)
        os.system(cmd_cufflinks)
        return cmd_cufflinks
    
    def assemblies(self,prefix):
        cmd_assemblies = "ls -R %s/SRR*/transcripts.gtf > %s/assemblies.txt" % (self.cufflinks_dir,self.cufflinks_dir)
        os.system(cmd_assemblies)
        return cmd_assemblies 

    def cuffmerge(self,prefix):
        gtf = "/data4/genome_index/hg19/tophat2/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
        fasta = "/data4/genome_index/hg19/tophat2/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
        cmd_cuffmerge = "cuffmerge -p 21 -o %s -g %s -s %s %s/assemblies.txt" % (self.cufflinks_dir,gtf,fasta,self.cufflinks_dir)
        os.system(cmd_cuffmerge)
        return cmd_cuffmerge 

    def cuffdiff(self,sample_list,bam_list):
        fasta = "/data4/genome_index/hg19/tophat2/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
        cmd_cuffdiff = "cuffdiff -p 21 -o %s/diff_out -b %s -L %s -u %s/merged.gtf %s" % (self.cufflinks_dir,fasta,sample_list,wf.cufflinks_dir,bam_list)
        print(cmd_cuffdiff)
        os.system(cmd_cuffdiff)
        return cmd_cuffdiff 

    def htseqCount(self,sample,sort_bam):
        gtf = "/data3/zhaochen/reference/mm10/mm10_annotation.gtf"
        output = os.path.join(self.htseqCount_dir,sample + "_htseq_count.txt")
        cmd = "htseq-count -f bam -r name -s no -a 20 -m union %s %s > %s" % (sort_bam,gtf,output)
        os.system(cmd)
        return cmd

    def featureCounts(self,bamlist):
        gtf = "/data3/zhaochen/reference/mm10/mm10_annotation.gtf"
        output = os.path.join(wf.featureCount_dir, "KO-KDM2B_featurecounts.txt")
        cmd = "featureCounts -T 14 -a %s -o %s -p -t gene -g gene_name %s" % (gtf, output,bamlist)
        os.system(cmd)
        return cmd

if __name__=="__main__":
    import re
    wf = work_flow()
    step = 10
    sample_list = ["Ctrl-7","Ctrl-12","Ctrl-20","cKO-10","cKO-14","cKO-21"]

    if step < 1:
        for root,dirs,files in os.walk(wf.cleandir):
            for file in files:
                file = os.path.join(root,file)
                wf.fastqc(file)

        wf.multiqc()

###    step 1,cutadapter
    if step < 2:
        for sample in sample_list:
            sample_R1 = os.path.join(wf.rawdataDir,"%s_1.fq.gz" % (sample))
            sample_R2 = os.path.join(wf.rawdataDir,"%s_2.fq.gz" % (sample))
            cmd = wf.cutadapt_PE(sample_R1=sample_R1,sample_R2=sample_R2)
            wf.run(cmd=cmd)


    if step < 3:
        wf.STAR_build_index_human()

    if step < 4:
        sample_num = []
        for filename in os.listdir(wf.cleandir):
            if os.path.splitext(filename)[1] == '.gz':
                sample = filename.split('_')[0]
                if sample not in sample_num:
                    sample_num.append(sample)

        samples = sorted(sample_num)
        for sample in samples:
            clean_R1 = os.path.join(wf.cleandir, sample + "_1.fq.gz")
            clean_R2 = os.path.join(wf.cleandir, sample + "_2.fq.gz")
            # print(clean_R1)
            # print(clean_R2)
            #wf.STAR(sample=sample,clean_R1=clean_R1,clean_R2=clean_R2)
#            wf.tophat(sample=sample,clean_R1=clean_R1,clean_R2=clean_R2,species="human")
#            cutadapt_SE(sample=sample)
#            wf.hisat2(sample=sample,clean_R1=clean_R1,clean_R2=clean_R2)
            wf.samtools(sample=sample)
#            sort_bam = os.path.join(wf.samtools_dir, sample + '_sorted.bam')
#            wf.htseqCount(sample=sample,sort_bam=sort_bam)

    if step > 5:
        bamlist = []
        sample_list = []
        for root,dirs,files in os.walk(wf.samtools_dir):
            for file in files:
#                if not re.search("_sorted.bam",file):
#                    continue
                if os.path.splitext(file)[1] == ".bam":
                    sort_bam = os.path.join(wf.samtools_dir,file)
                    bamlist.append(sort_bam)
                    # prefix = file.replace("_sorted.bam","")
                    # sample_list.append(prefix)
#                    wf.cufflinks(prefix,sort_bam)
#            wf.assemblies(prefix)
#             sample_list = sorted(sample_list)
#             sample_list = ",".join(sample_list)
#            print(sample_list)
        bamlist = sorted(bamlist)
        bamlist = " ".join(bamlist)
        wf.featureCounts(bamlist=bamlist)
#            wf.cuffmerge(prefix)
#            wf.cuffdiff(sample_list,bam_list)


                





   



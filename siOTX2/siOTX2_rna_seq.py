import os,sys
import pandas as pd
import numpy as np

rawDir = "/backup_raw/AOM-DSS/siOtx2"
fastqcDir = "/data3/zhaochen/project/colon_cancer/colon_chip/siOTX2/fastqc"
multiqcDir = "/data3/zhaochen/project/colon_cancer/colon_chip/siOTX2/fastqc/multiqc"
fastpDir = "/data3/zhaochen/project/colon_cancer/colon_chip/siOTX2/fastp"
cleanDir = "/data3/zhaochen/project/colon_cancer/colon_chip/siOTX2/clean"
alignDir = "/data3/zhaochen/project/colon_cancer/colon_chip/siOTX2/align"
tophatDir = "/data3/zhaochen/project/colon_cancer/colon_chip/siOTX2/tophat"
samtoolsDir = "/data3/zhaochen/project/colon_cancer/colon_chip/siOTX2/samtools"
tophat_samtoolsDir = "/data3/zhaochen/project/colon_cancer/colon_chip/siOTX2/tophat_samtools"
cufflinksDir = "/data3/zhaochen/project/colon_cancer/colon_chip/siOTX2/cufflinks"
featureCountsDir = "/data3/zhaochen/project/colon_cancer/colon_chip/siOTX2/featureCounts"

def fastqc(file):
    cmd = "fastqc -t 7 -o %s %s" % (fastqcDir,file)
    os.system(cmd)
    return cmd

def multiqc():
    cmd = "multiqc %s -o %s" % (fastqcDir,multiqcDir)
    os.system(cmd)
    return cmd

def cutadapt_PE(sample_R1, sample_R2):
    """
    PE mains pair end
    SE mains singal end
    """
    adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    name_R1 = os.path.basename(sample_R1)
    name_R2 = os.path.basename(sample_R2)
    output_R1 = os.path.join(cleanDir, name_R1.replace("_R1.fastq.gz", "_R1.fq.gz"))
    output_R2 = os.path.join(cleanDir, name_R2.replace("_R2.fastq.gz", "_R2.fq.gz"))
    cmd = "cutadapt -j 8 -a %s -A %s -u 5 -u -15 -U 5 -U -15 -m 30 -o %s -p %s %s %s" % \
          (adapter_1, adapter_2, output_R1, output_R2, sample_R1, sample_R2)
    os.system(cmd)
    return cmd

def fastp(sample_R1,sample_R2):
    name_R1 = os.path.basename(sample_R1)
    name_R2 = os.path.basename(sample_R2)
    output_R1 = os.path.join(fastpDir, name_R1.replace("_R1.fastq.gz", "_R1.fq.gz"))
    output_R2 = os.path.join(fastpDir, name_R2.replace("_R2.fastq.gz", "_R2.fq.gz"))
    cmd = "fastp -w 14 -i %s -o %s -I %s -O %s" % (sample_R1,output_R1,sample_R2,output_R2)
    os.system(cmd)
    return cmd

def hisat2(sample,clean_R1,clean_R2):
    hisat2_index = "/data3/zhaochen/reference/hg19/hg19_hisat2_index/genome"
    output = os.path.join(alignDir,sample + "_RNA.sam")
    log = os.path.join(alignDir,sample + "_hisat2.txt")
    cmd = "hisat2 -q --phred33 -p 8 -x %s -1 %s -2 %s -S %s > %s" % (hisat2_index,clean_R1,clean_R2,output,log)
    print(cmd)
    os.system(cmd)
    view_bam = os.path.join(samtoolDir,sample + "_view.bam")
    cmd_view = "samtools view -@ 7 -bS -o %s %s" % (view_bam,output)
    os.system(cmd_view)
    os.system("rm %s" % output)
    return cmd

def tophat2(sample,cleanR1,cleanR2):
    tophat_index = "/data4/genome_index/hg19/tophat2/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
    gtf = "/data3/zhaochen/reference/hg19/hg19_annotation.gtf"
    cmd = "/data3/zhaochen/anaconda3/envs/py2/bin/tophat -p 21 -o %s/%s_tophat -G %s %s %s %s" \
              % (tophatDir,sample,gtf,tophat_index,cleanR1,cleanR2)
    os.system(cmd)
    return cmd

def samtools(view_bam,sample):
    sort_bam = os.path.join(tophat_samtoolsDir,sample + "_sort.bam")
    cmd_sort = "samtools sort -@ 7 -O BAM -o %s %s" % (sort_bam,view_bam)
    os.system(cmd_sort)
    cmd_index = "samtools index -@ 7 -b %s" % sort_bam
    os.system(cmd_index)
    return cmd_sort

def featureCounts(bam_list,output):
    #gtf = "/data3/zhaochen/reference/hg19/gencode.v32lift37.annotation.gtf"
    #gtf = "/data3/zhaochen/reference/hg19/STAR_genome_index/Homo_sapiens.GRCh38.99.gtf"
    gtf = "/data3/zhaochen/reference/hg19/hg19_annotation.gtf"
    cmd = "featureCounts -T 7 -a %s -o %s -p -t gene -g gene_name %s" % (gtf,output,bam_list)
    os.system(cmd)
    return cmd

def STAR(sample, clean_R1, clean_R2):
    prefix = os.path.join(alignDir, sample)
    cmd = "STAR --genomeDir /data3/zhaochen/reference/hg19/STAR_genome_index/star_genome " \
          "--readFilesIn %s %s --readFilesCommand zcat --sjdbGTFfile /data3/zhaochen/reference/hg19/STAR_genome_index/Homo_sapiens.GRCh38.99.gtf " \
          "--runThreadN 21 --outFileNamePrefix %s --outSAMtype BAM Unsorted" % (clean_R1, clean_R2, prefix)
    os.system(cmd)
    return cmd

def cufflinks(sample,sortbam):
    cmd = "cufflinks -p 21 -o %s/%s %s" % (cufflinksDir,sample,sortbam)
    os.system(cmd)
    return cmd

# should split file from assemblies, cuffdiff will report error
def assemblies():
    cmd ="ls -R %s/si*/transcripts.gtf > %s/assemblies.txt" % (cufflinksDir,cufflinksDir)
    os.system(cmd)
    return cmd

def cuffmerge():
    #gtf = "/data4/genome_index/hg19/tophat2/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
    #mask_gtf = "/data3/zhaochen/reference/hg19/hg19_rmsk.gtf"
    #fasta = "/data4/genome_index/hg19/tophat2/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"

    gtf = "/data3/zhaochen/reference/hg19/STAR_genome_index/Homo_sapiens.GRCh38.99.gtf"
    fasta = "/data3/zhaochen/reference/hg19/STAR_genome_index/star_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    #fasta = "/data3/zhaochen/reference/hg19/STAR_genome_index/star_genome/Genome"
    mergetxt = os.path.join(cufflinksDir,"assemblies.txt")
    cmd = "cuffmerge -p 14 -o %s/merged_asm -g %s -s %s %s" % (cufflinksDir,gtf,fasta,mergetxt)
    os.system(cmd)
    return cmd

def cuffdiff(sample_list,bam_list):
    fasta = "/data4/genome_index/hg19/tophat2/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
    #fasta = "/data3/zhaochen/reference/hg19/STAR_genome_index/star_genome/Genome"
    #fasta = "/data3/zhaochen/reference/hg19/STAR_genome_index/star_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    cmd = "cuffdiff -p 10 -o %s/diff_out -b %s -L %s -u %s/merged_asm/merged.gtf %s" % (cufflinksDir,fasta,sample_list,cufflinksDir,bam_list)
    print(cmd)
    os.system(cmd)
    return cmd

if __name__=="__main__":
    step = 8
    filelist = []
    if step < 1:
        for root,dirs,files in os.walk(cleanDir):
            for file in files:
                if ".gz" in file:
                    sample = file.split("_")[0]
                    filelist.append(sample)
                    file = os.path.join(root,file)
                    #print(file)
                    fastqc(file=file)
                    multiqc()

    if step < 2:
        for root, dirs, files in os.walk(rawDir):
            for file in files:
                if "fastq.gz" in file:
                    sample = file.split("_")[0]
                    if sample not in filelist:
                        filelist.append(sample)

        for sample in filelist:
            sample_R1 = os.path.join(rawDir, sample + "_combined_R1.fastq.gz")
            sample_R2 = os.path.join(rawDir, sample + "_combined_R2.fastq.gz")
            print(sample)
            #cutadapt_PE(sample_R1=sample_R1,sample_R2=sample_R2)
            #fastp(sample_R1=sample_R1,sample_R2=sample_R2)

    if step < 3:
        for root,dirs,files in os.walk(cleanDir):
            for file in files:
                sample = file.split("_")[0]
                if sample not in filelist:
                    filelist.append(sample)

        filelist = list(reversed(filelist))
        for sample in filelist:
            clean_R1 = os.path.join(cleanDir,sample + "_combined_R1.fq.gz")
            clean_R2 = os.path.join(cleanDir,sample + "_combined_R2.fq.gz")
            print(sample)
            #hisat2(sample=sample,clean_R1=clean_R1,clean_R2=clean_R2)
            #STAR(sample=sample,clean_R1=clean_R1,clean_R2=clean_R2)
            #tophat2(sample=sample,cleanR1=clean_R1,cleanR2=clean_R2)
            bamdir = os.path.join(tophatDir, sample + "_tophat")
            view_bam = os.path.join(bamdir, 'accepted_hits.bam')
            samtools(view_bam=view_bam,sample=sample)

    if step < 4:
        for root,dirs,files in os.walk(alignDir):
            for file in files:
                if ".out.bam" in file:
                    sample = file.split(".")[0]
                    file = os.path.join(root,file)
                    samtools(view_bam=file,sample=sample)

    if step < 5:
        bamlist = ['/data3/zhaochen/project/colon_cancer/colon_chip/siOTX2/samtool/siNC-1_sort.bam',
                   '/data3/zhaochen/project/colon_cancer/colon_chip/siOTX2/samtool/siNC-2_sort.bam']
        klf3list = ['/data3/zhaochen/project/colon_cancer/colon_chip/siOTX2/samtool/siNC-1_sort.bam',
                    '/data3/zhaochen/project/colon_cancer/colon_chip/siOTX2/samtool/siNC-2_sort.bam']
        sample_list = ["siNC-1","siNC-1"]
        klf3_sample_list = ["siNC-1","siNC-1"]
        for root,dirs,files in os.walk(samtoolDir):
            for file in files:
                if os.path.splitext(file)[1] == ".bam":
                    if "OTX2" in file:
                        sample =file.split("_")[0]
                        sample_list.append(sample)
                        file = os.path.join(root,file)
                        bamlist.append(file)
                    elif "KLF3" in file:
                        sample = file.split("_")[0]
                        klf3_sample_list.append(sample)
                        file = os.path.join(root,file)
                        klf3list.append(file)

        bamlist = sorted(bamlist)
        bamlist = " ".join(bamlist)
        output = os.path.join(featureCountsDir, "siOTX2_featureCounts_genecode.txt")
        featureCounts(bamlist=bamlist,output=output)
        klf3list = sorted(klf3list)
        klf3list = " ".join(klf3list)
        out = os.path.join(featureCountsDir,"siKLF3_featureCounts_genecode.txt")
        featureCounts(bamlist=klf3list,output=out)

        # sample_list = sorted(sample_list)
        # sample_list = " ".join(sample_list)
        # klf3_sample_list = sorted(klf3_sample_list)
        # klf3_sample_list = " ".join(klf3_sample_list)
        #
        # cuffdiff(tf="siOTX2",sample_list=sample_list,bam_list=bamlist)
        # cuffdiff(tf="siKLF3",sample_list=klf3_sample_list,bam_list=klf3list)

    if step < 6:
        sample_list = []
        bam_list = []
        for root,dirs,files in os.walk(samtoolsDir):
            for file in files:
                if os.path.splitext(file)[1] == ".bam":
                    sample = file.split("A")[0]
                    sample_list.append(sample)
                    file = os.path.join(root,file)
                    bam_list.append(file)
                    #print(sample)
                    #print(file)
                    #cufflinks(sample=sample,sortbam=file)

        #assemblies()
        #cuffmerge()

        sample_list = sorted(sample_list)
        sample_list = ",".join(sample_list)
        bam_list = sorted(bam_list)
        bam_list = " ".join(bam_list)
        output = os.path.join(featureCountsDir,"star_featureCounts.txt")
        #featureCounts(bam_list=bam_list,output=output)
        #print(sample_list)
        #print(bam_list)
        cuffdiff(sample_list=sample_list,bam_list=bam_list)

    if step > 7:
        sample_list = []
        bam_list = []
        for root, dirs, files in os.walk(tophat_samtoolsDir):
            for file in files:
                if os.path.splitext(file)[1] == ".bam":
                    sample = file.split("_")[0]
                    sample_list.append(sample)
                    file = os.path.join(root, file)
                    bam_list.append(file)
                    # print(sample)
                    # print(file)
                    # cufflinks(sample=sample,sortbam=file)

        # assemblies()
        # cuffmerge()

        sample_list = sorted(sample_list)
        sample_list = ",".join(sample_list)
        bam_list = sorted(bam_list)
        bam_list = " ".join(bam_list)
        output = os.path.join(featureCountsDir, "tophat_featureCounts.txt")
        featureCounts(bam_list=bam_list,output=output)

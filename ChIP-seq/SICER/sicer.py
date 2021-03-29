#!/usr/bin/env python3
#######################################################################
# File Name: bed.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Fri 01 Nov 2019 06:58:48 PM CST
# Description: 
# History: 
#######################################################################
import os,sys

class work_flow():
    
    def __init__(self):
        self.inputDir = "/data3/zhaochen/project/colon_cancer/colon_chip/SICER/Input"
        self.sicerDir = "/data3/zhaochen/project/colon_cancer/colon_chip/SICER"
        self.intensity_dir = "/data3/zhaochen/colon_chip/intensity"

    def bamTobed(self,bam_file,bed_file):
        cmd = "bedtools bamtobed -i %s > %s" % (bam_file,bed_file)
        os.system(cmd)
        return cmd 

    def SICER(self,treat_file,control_file,marker):
        treatDir = os.path.join(wf.sicerDir,marker)
        cmd = "sicer -t %s -c %s -s mm10 --redundancy_threshold 1 -f 150 -w 1000 -g 10000 -egf 0.77 -fdr 0.01 -cpu 14 " \
              "--significant_reads -o %s" % (treat_file,control_file,treatDir)
        #print(cmd)
        os.system(cmd)
        return cmd 

class signal_intensity():
    
    def __init__(self):
        self.intensity_dir = "/data3/zhaochen/colon_chip/intensity"
        
    def bamCompare(self,treat_bam,control_bam,output):
        cmd = "bamCompare -b1 %s -b2 %s --normalizeUsing RPKM -bs 100 -o %s --operation ratio -p 7" % (treat_bam,control_bam,output)
        os.system(cmd)
        return cmd

    def computeMatrix(self,peak_region,bigwig):
        cmd = "computeMatrix reference-point --referencePoint center -b 50000 -a 50000 --skipZeros -R %s -S %s -o %s -p 7" % (peak_region,bigwig,matrix)
        os.system(cmd)
        return cmd 

    def plotProfile(self,):
        cmd = "plotProfile -m %s -o %s.pdf" % (matrix,name)
        os.system(cmd)
        return cmd 

if __name__=="__main__":
    wf = work_flow()
    step = 8
###SICER call peaks               
    if step < 1:
        for filename in os.listdir(wf.treatDir):
            if os.path.splitext(filename)[1] == ".bed":
                treat_file = os.path.join("H3K9me3",filename)
                control_file = os.path.join("Input",filename)
                sample = filename.split('.')[0]
                dir1 = os.path.join(wf.sicerDir,"H3K9me3")
                dir2 = os.path.join(wf.sicerDir,"Input")
#                os.makedirs(dir1)
#                os.makedirs(dir2)
                wf.SICER(treat_file,control_file)

    poolbeddir = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bed"
    markers = ["H3K9me2","H3K9me3"]
    Times = ["ctrl","2weeks","4weeks","7weeks","10weeks"]
    if step > 5:
        for marker in markers:
            for time in Times:
                treat_file = os.path.join(poolbeddir,time + "-" + marker + "_pool.bed")
                control_file = os.path.join(poolbeddir,time + "-Input_pool.bed")
                wf.SICER(treat_file=treat_file,control_file=control_file,marker=marker)



###ChIP-seq intensity
    s = signal_intensity()

    if step < 2:
        for root,dirs,files in os.walk(s.bw_Input):
            for file in files:
                if "H3K9me3" in file:
                    broad = os.path.join(s.sicer, )
                    outName = os.path.join(s.intentityDir, "H3K9me3_center_matrix.mat.gz")
                    cmd = s.computeMatrix(bigwiglist=bigwiglist, broadPeaklist=broad, after=50000, brfore=50000,outName=outName)

    if step < 3:
        for filename in os.listdir(wf.treatDir):
            if os.path.splitext(filename)[1] == ".bam":
                treat_bam = os.path.join(wf.treatDir,filename)
                control_bam = os.path.join(wf.inputDir,filename)
                sample = filename.split('.')[0]
                output = "H3K9me3_Input_"+sample+"_ratio.bw"
#                si.bamCompare(treat_bam,control_bam,output)
                si.computeMatrix()


                




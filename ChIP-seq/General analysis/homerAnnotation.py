#!/usr/bin/env python3
#######################################################################
# File Name: homerAnnotation.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Sat 21 Dec 2019 05:03:03 PM CST
# Description: 
# History: 
#######################################################################
import os,sys
import re
import pandas as pd
import numpy as np
import subprocess

class homer():

    def __init__(self):
        self.macs2 = "/data3/zhaochen/project/colon_cancer/colon_chip/macs2"

    def run(self,cmd):
        p = subprocess.Popen(cmd,shell=True)
        p.wait()
        return p

    def annotatePeaks(self,broadPeak,txt):
        cmd = "annotatePeaks.pl %s mm10 > %s" % (broadPeak,txt)
        return cmd

if __name__=="__main__":
    h = homer()
    step = 0
    if step < 1:
        for root,dirs,files in os.walk(h.macs2):
            for file in files:
                if os.path.splitext(file)[1] == ".broadPeak":
                    sample = file.split("_")[0]
                    broadPeak = os.path.join(root,file)
                    txt = broadPeak.replace("_peaks.broadPeak","_annotation.txt")
#                    os.system("wc -l %s > peak_number.txt" % broadPeak)
                    cmd = h.annotatePeaks(broadPeak=broadPeak,txt=txt)
                    h.run(cmd)

    if step < 0:
        for root,dirs,files in os.walk(h.macs2):
            for file in files:
                if file.split("_")[1] == "Input_annotation.txt":
                    file = os.path.join(root,file)


                    





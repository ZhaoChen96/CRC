#!/usr/bin/env python3
#######################################################################
# File Name: matrix.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Wed 23 Oct 2019 08:42:02 PM CST
# Description: 
# History: 
#######################################################################
import pandas as pd
import numpy as np
import json
import re,os
import glob

dfs = []
filelist = []
labels = []
for root,dirs,files in os.walk('./'):
    for file in files:
        if os.path.splitext(file)[1] == '.counts':
            filelist.append(os.path.join(root,file))
            labels.append(file.split('.')[0])
            data = pd.read_csv(file,header=None,sep='\t',names=[str(file.split('.')[0])])
            dfs.append(data)

df = pd.concat(dfs,axis=1)
df.to_csv('normal_htseq_matrix.txt',sep='\t',header=True,index=True)

cancer_counts = "/data3/zhaochen/TCGA-count/raw_data/cancer/cancer_htseq_matrix.txt"
#frame = pd.read_csv(cancer_counts,nrows=3,sep='\t',)
#print(frame)


#!/usr/bin/env python3
#######################################################################
# File Name: step1_CRCSNPs.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Wed 13 Nov 2019 08:30:16 PM CST
# Description: 
# History: 
#######################################################################
import pandas as pd
import numpy as np

df = pd.read_csv("/data3/zhaochen/project/colon_cancer/GWAS/full",header=0,index_col=False,sep='\t',low_memory=False)

new_df = pd.DataFrame()
for index,row in df.iterrows():
    row["p_VALUE"] = float(row["P-VALUE"])
    if (row["DISEASE/TRAIT"].startswith("Colorectal") and row["p_VALUE"] <= 1E-8):
        new_df = new_df.append(row)

select_df = new_df.loc[:,["DISEASE/TRAIT","REGION","CHR_ID","CHR_POS","REPORTED GENE(S)","MAPPED_GENE","SNP_GENE_IDS","SNPS", "P-VALUE", "OR or BETA"]]

ensembl = []
for item in select_df["SNP_GENE_IDS"]:
    array = str(item).split(",")
    ensembl = ensembl + array
    for arr in array:
        print(arr)    
select_df.to_csv("CRC_GWAS_site.txt",header=True,sep='\t')

data = pd.read_csv("CRC_GWAS_site.txt",header=0,sep='\t',index_col=0)
data = data.loc[:,["CHR_ID","CHR_POS","SNPS"]]
data.to_csv("GWAS_site_3columns.txt",header=False,sep='\t',index=False)

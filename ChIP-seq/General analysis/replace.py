#!/usr/bin/env python3
#######################################################################
# File Name: replace.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Wed 14 Oct 2020 10:17:35 PM CST
# Description: 
# History: 
#######################################################################
import os,sys
import pandas as pd
import numpy as np

intensityDir = "/data3/zhaochen/project/colon_cancer/colon_chip/intensity"
Times = ["10weeks","4weeks","7weeks","ctrl"]
for root,dirs,files in os.walk(intensityDir):
    for file in files:
        if "noInput" in file:
            if "scale" in file:
                file = os.path.join(root,file)
                new_file = file.replace("scale","center")
                # print(file)
                # print(new_file)
                os.system("mv %s %s" % (file,new_file))


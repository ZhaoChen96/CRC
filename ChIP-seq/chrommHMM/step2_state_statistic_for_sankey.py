#!/usr/bin/env python3
#######################################################################
# File Name: step3_state_statistic.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: Sun 28 Feb 2021 08:55:17 PM CST
# Description: 
# History: 
#######################################################################
import os,sys
import pandas as pd
import numpy as np
import plotly.graph_objects as go

chromHMMDir = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/bed/state"
traceDir = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/trace"

def makewindows(bedfile):
    output = bedfile.replace("segments.bed","segments_200.bed")
    cmd = "bedtools makewindows -w 200 -b %s -i src > %s" % (bedfile,output)
    os.system(cmd)
    return cmd

def merge_state():
    fa = pd.read_csv(os.path.join(chromHMMDir,"ctrl_13_segments_200.bed"),header=None,sep="\t",names=["chr","start","end","week0"])
    fb = pd.read_csv(os.path.join(chromHMMDir,"2weeks_13_segments_200.bed"), header=None, sep="\t", names=["chr", "start", "end", "week2"])
    data = pd.merge(fa, fb, on=["chr", "start", "end"])
    fc = pd.read_csv(os.path.join(chromHMMDir,"4weeks_13_segments_200.bed"), header=None, sep="\t", names=["chr", "start", "end", "week4"])
    data = pd.merge(data, fc, on=["chr", "start", "end"])
    fd = pd.read_csv(os.path.join(chromHMMDir,"7weeks_13_segments_200.bed"), header=None, sep="\t", names=["chr", "start", "end", "week7"])
    data = pd.merge(data, fd, on=["chr", "start", "end"])
    fe = pd.read_csv(os.path.join(chromHMMDir,"10weeks_13_segments_200.bed"), header=None, sep="\t", names=["chr", "start", "end", "week10"])
    data = pd.merge(data, fe, on=["chr", "start", "end"])
    data.to_csv("/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/merge_13_segments_200.txt", sep="\t", header=True, index=False)

def sort(file,sortfile):
    cmd = "sort -k1,1 "
    return cmd

def statistic_states_count(binfile):
    df = pd.read_csv(binfile,sep="\t",names=["chr","start","end","week0","week2","week4","week7","week10"])

    change_state_list = []
    # ctrl_week2 = df.groupby(["week0","week2"]).size().reset_index(name="count")
    # ctrl_week2.to_csv(os.path.join(traceDir,"ctrl_week2.txt"),sep="\t")
    # change_state_list.append(ctrl_week2)
    #
    # week2_week4 = df.groupby(["week2","week4"]).size().reset_index(name="count")
    # week2_week4.to_csv(os.path.join(traceDir,"week2_week4.txt"),sep="\t")
    # change_state_list.append(week2_week4)
    #
    # week4_week7 = df.groupby(["week4", "week7"]).size().reset_index(name="count")
    # week4_week7.to_csv(os.path.join(traceDir, "week4_week7.txt"), sep="\t")
    # change_state_list.append(week4_week7)
    #
    # week7_week10 = df.groupby(["week7", "week10"]).size().reset_index(name="count")
    # week7_week10.to_csv(os.path.join(traceDir, "week7_week10.txt"), sep="\t")
    # change_state_list.append(week7_week10)

    all = df.groupby(["week0","week2","week4","week7",'week10']).size().reset_index(name="count")
    all.to_csv(os.path.join(traceDir,"all_timepoint.txt"),sep="\t")
    print(all)
    #print(change_state_list)
    return change_state_list

def node_trace():
    label = ["ctrl_E_" + str(i) for i in range(1,14)] + ["week2_E_" + str(i) for i in range(1,14)] + \
            ["week4_E_" + str(i) for i in range(1,14)] + ["week7_E_" + str(i) for i in range(1,14)] + ["week10_E_" + str(i) for i in range(1,14)]
    label_s = ["E" + str(i) for i in range(1,14)] + ["E" + str(i) for i in range(1,14)]

    print(label[0:13])
    source,target,value = [],[],[]

    ctrl_week2 = pd.read_csv("/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/trace/ctrl_week2.txt",sep="\t",index_col=0)
    for i in range(1,14):
        for j in range(1,14):
            if i == j and i != 8:
                continue
            source.append(i)
            target.append(j + 13)
            row = ctrl_week2[(ctrl_week2["week0"] == "E" + str(i)) & (ctrl_week2["week2"] == "E" + str(j))]
            value.append(int(row["count"]))

    week2_week4 = pd.read_csv("/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/trace/week2_week4.txt", sep="\t",
                             index_col=0)
    for i in range(1, 14):
        for j in range(1, 14):
            if i == j and i != 8:
                continue
            source.append(i + 13)
            target.append(j + 13 + 13)
            row = week2_week4[(week2_week4["week2"] == "E" + str(i)) & (week2_week4["week4"] == "E" + str(j))]
            if row.empty:
                continue
            #print(row["count"])
            value.append(int(row["count"]))

    week4_week7 = pd.read_csv("/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/trace/week4_week7.txt",sep="\t",index_col=0)
    for i in range(1, 14):
        for j in range(1, 14):
            if i == j and i != 8:
                continue
            source.append(i + 13 + 13)
            target.append(j + 13 + 13 + 13)
            row = week4_week7[(week4_week7["week4"] == "E" + str(i)) & (week4_week7["week7"] == "E" + str(j))]
            if row.empty:
                continue
            # print(row["count"])
            value.append(int(row["count"]))

    week7_week10 = pd.read_csv("/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/trace/week7_week10.txt",
                              sep="\t", index_col=0)
    for i in range(1, 14):
        for j in range(1, 14):
            if i == j and i != 8:
                continue
            source.append(i + 13 + 13 + 13)
            target.append(j + 13 + 13 + 13 + 13)
            row = week7_week10[(week7_week10["week7"] == "E" + str(i)) & (week7_week10["week10"] == "E" + str(j))]
            if row.empty:
                continue
            # print(row["count"])
            value.append(int(row["count"]))

    result = []
    result.append(source)
    result.append(target)
    result.append(value)
    print(result)
    return result

def sankey():
    result = node_trace()
    fig = go.Figure(data=[go.Sankey(
        node= dict(
            pad=15,
            thickness = 20,
            line = dict(color = "black", width=0.5),
            color = "blue"
        ),
        link = dict(
            source = result[0],
            target = result[1],
            value = result[2]
        ),
        textfont = dict(size=1)
    )])

    fig.update_layout(title_text="Basic Sankey Diagram", font_size=10)
    fig.write_image("./fig.png",width=800,height=800)





if __name__ == "__main__":
    step = 5

    if step < 1:
        for root,dirs,files in os.walk(chromHMMDir):
            for file in files:
                if "_13_segments.bed" in file:
                    file = os.path.join(root,file)
                    makewindows(bedfile=file)

    if step < 2:
        merge_state()

    if step < 3:
        binfile = "/data3/zhaochen/project/colon_cancer/colon_chip/chromHMM/merge_13_segments_200.txt"
        statistic_states_count(binfile=binfile)

    if step > 4:
        node_trace()
        sankey()
import sys,os
import pandas as pd
import numpy as np

def creat_config():
    #data_dir = "/home/zhluo/Project/CRC/data_nazhang/step28_checkbigWig/bigwig"
    data_dir = "/data3/zhaochen/project/colon_cancer/colon_chip/pool_bigwig/bw_Input"
    times = ["ctrl", "2weeks", "4weeks", "7weeks", "10weeks"]
    element = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"]
    
    sample_list = []
    for one_marker in element[0:3]:
        #enhancer
        marker = one_marker
        
        for one_p in times:
            #sample = "%s-1-%s_combined_1_paired.merged.nodup.pooled_x_%s-1-Input_combined_1_paired.merged.nodup.pooled.fc.signal.bigwig" %(one_p, marker, one_p)
            sample = "%s-%s_bs1bp_ratio.bw" % (one_p, marker)
            sample_list.append(os.path.join(data_dir, sample))
        
    sample_str = " ".join(sample_list)
    #anno = "/home/zyang/Project/CRC/step49_track/gencode.vM25.basic.annotation.sort.bed12"
    anno = "/data3/zhaochen/project/colon_cancer/colon_chip/genomeTrack/gencode.vM25.basic.annotation.sorted.bed12"
    #anno = "/data3/zhaochen/project/colon_cancer/colon_chip/genomeTrack/gencode.vM25.basic.annotation.gtf"
    cmd = "make_tracks_file --trackFiles %s %s -o zhao_all_tracks.ini"  %(sample_str, anno)
    os.system(cmd)

def pyGenomeTracks(region,name):
    cmd = "pyGenomeTracks --tracks zhao_all_tracks.ini --region chr%s -o %s.png --fontSize 20 --dpi 300 --trackLabelFraction 0.15" % \
          (region,name)
    os.system(cmd)
    return cmd

def change_symbol():
    sortbed = pd.read_csv("/data3/zhaochen/project/colon_cancer/colon_chip/genomeTrack/gencode.vM25.basic.annotation.sort.bed12",
                          sep="\t",header=None,names=["chr","start","end","transcript_id","score","strand","a","b","c","d","e","f"])
    file = pd.read_csv("/data3/zhaochen/project/colon_cancer/colon_chip/genomeTrack/gene_name_and_transcriptID.txt",sep="\t")
    data = pd.merge(sortbed,file,on="transcript_id")
    df = data.loc[:,["chr","start","end","gene_name","score","strand","a","b","c","d","e","f"]]
    df.to_csv("/data3/zhaochen/project/colon_cancer/colon_chip/genomeTrack/gencode.vM25.basic.annotation.sorted.bed12",sep="\t",
              index=False,header=False)

if __name__ == "__main__":
    #change_symbol()
    #creat_config()
    """
    make_tracks_file --trackFiles <bigwig file> <bed file> etc. -o tracks.ini
    
    pyGenomeTracks --tracks tracks.ini --region chr1:2,500,000-3,000,000 -o bigwig.png
    pyGenomeTracks --tracks enahcner_tracks.ini --region chr11:112780224-112789760 -o bigwig_bed.png
    pyGenomeTracks --tracks all_tracks.ini --region chr13:116296262-116312689 -o Isl1.png --fontSize 20 --trackLabelFraction 0.2

    """
    #pyGenomeTracks(region="15:36565057-36709190",name="Pabpc1")
    #pyGenomeTracks(region="5:91946482-91972773",name="Thap6")
    #pyGenomeTracks(region="12:21,375,358-21,458,565",name="Ywhaq")
    #pyGenomeTracks(region="9:45,034,619-45,078,073",name="Mpzl3")
    #pyGenomeTracks(region="7:140,966,314-140,972,609",name="Ifitm1")
    #pyGenomeTracks(region="17:83,687,398-83,818,985",name="Mta3")
    #pyGenomeTracks(region="8:60908429-60989846",name="Clcn3")
    #pyGenomeTracks(region="7:140,931,029-141,031,769",name="Ifitm1")
    pyGenomeTracks(region="5:90759360-90761624",name="Cxcl5")













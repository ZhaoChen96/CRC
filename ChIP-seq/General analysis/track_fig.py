import os,sys
import pandas as pd
from gtfparse import read_gtf


def bed12_process():
    #file_path = "/home/zyang/Project/CRC/step49_track/gencode.vM25.basic.annotation.gtf"
    file_path = "/data3/zhaochen/project/colon_cancer/colon_chip/genomeTrack/gencode.vM25.basic.annotation.gtf"
    df = read_gtf(file_path)
    print(df.columns)
    print(df[0:10])
    df_genes = df[df["feature"] == "transcript"]
    df_sub = df_genes[['gene_name', 'transcript_id']]
    df_sub.to_csv("gene_name_and_transcriptID.txt", sep="\t", index=False)
    
    
    #bed_file = "/home/zyang/Project/CRC/step49_track/gencode.vM25.basic.annotation.sort.bed12"
    bed_file = "/data3/zhaochen/project/colon_cancer/colon_chip/genomeTrack/gencode.vM25.basic.annotation.sort.bed12"
    bed_df = pd.read_csv(bed_file, sep="\t",\
                names=["chr", "start", "end", "name", "score", "strand", "thick_start", "thick_end",\
                 "rgb", "block_count", "block_size", "block_start"])
    new_bed = bed_df.merge(df_sub, left_on = "name", right_on = "transcript_id", how="left")
    new_bed =new_bed[["chr", "start", "end", "gene_name", "score", "strand", "thick_start", "thick_end",\
                 "rgb", "block_count", "block_size", "block_start"]]
    new_bed.to_csv("sorted.changeName.bed", sep="\t", index=False, header=False)



def draw_track_fluff():
    data_dir = "/home/zhluo/Project/CRC/data_nazhang/step28_checkbigWig/bigwig"
    times = ["ctrl", "2weeks", "4weeks", "7weeks", "10weeks"]
    element = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "H3K9me2", "H3K9me3"]
    
    
    #enhancer
    marker = "H3K27ac"
    sample_list = []
    for one_p in times:
        sample = "%s-1-%s_combined_1_paired.merged.nodup.pooled_x_%s-1-Input_combined_1_paired.merged.nodup.pooled.fc.signal.bigwig" %(one_p, marker, one_p)
        sample_list.append(os.path.join(data_dir, sample))
        
    sample_str = " ".join(sample_list)
    anno = "/home/zyang/Project/CRC/step49_track/gencode.vM25.basic.annotation.sort.bed12"
    cmd = "fluff profile -i chr6:41112015-41135714 -d %s -a %s -f 0 -s 1:4 -n -o test" %(sample_str, anno)
    print(cmd)
    
    
if __name__ == "__main__":
    bed12_process()
    #draw_track_fluff()
    # """
    # fluff profile -i chr6:41112015-41135714 -d 10weeks-1-H3K27ac_combined_1_paired.merged.nodup.pooled_x_10weeks-1-Input_combined_1_paired.merged.nodup.pooled.fc.signal.bw -a /home/zyang/Project/CRC/step49_track/gencode.vM25.basic.annotation.bed12 -f 0   -n -o test
    # """
rm(list = ls())
setwd("~/project/colon cancer/chip-seq/chromHMM/bed/")
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")

emission <- read.table("~/project/colon cancer/chip-seq/chromHMM/bed/emissions_13.txt",sep = "\t",header = TRUE)

display.brewer.all()
display.brewer.pal(n =11, name = 'RdBu')
color_list = brewer.pal(n = 11, name = "RdBu")
col_state = colorRamp2(c(0, 1), colors=color_list[c(6,11)])

element = c("Weakly active enhancer","Active enhancer","Poised enhancer", "Flanking active TSS","Flanking active TSS","Active TSS",
            "Transcription at gene 5' & 3'","Poised promoter","Repressed polycomb", "Repressed polycomb", "Bivalent enhancer", "Quiescent", "Heterochromatin"
             )

heterochromatin_color = brewer.pal(n = 11, name = "PRGn")[3]
promoter_color = brewer.pal(n = 9, name = "Greens")[8]
enhancer_color = brewer.pal(n = 9, name = "Reds")[8]
repressed_color = brewer.pal(n = 9, name = "Blues")[7]
quit_color = brewer.pal(n = 9, name = "Greys")[4]

state_color = c("1"=enhancer_color, "2"=enhancer_color, "3"=enhancer_color, "4"=promoter_color, "5"=promoter_color, 
                "6"=promoter_color, "7"=promoter_color, "8"=promoter_color, "9"=repressed_color, "10"=repressed_color, 
                "11"=enhancer_color, "12"=quit_color, "13"=heterochromatin_color)
text <- rep("",13)
text[2] = "Enhancer"
text[13] = "Heterochromatin"
text[9] = "Repressed"
text[12] = "Quiescent"
text[4] = "Promoter"
text_col = c("1"=enhancer_color, "2"=enhancer_color, "3"=enhancer_color, "4"=promoter_color, "5"=promoter_color, 
             "6"=promoter_color, "7"=promoter_color, "8"=promoter_color, "9"=repressed_color, "10"=repressed_color, 
             "11"=enhancer_color, "12"=quit_color, "13"=heterochromatin_color)
ha = rowAnnotation(foo = anno_text(element, location = 0.5, just = "center",
                                   gp = gpar(fill = state_color, col = "white", border = FALSE,fontsize = 8),
                                   width = max_text_width(element)*0.7))
right_annotation <- rowAnnotation(
  foo = anno_text(text,location = 0.9,just = "right",gp = gpar(col = text_col)),
  stat = 1:13,
  col= list(stat=state_color),
  simple_anno_size = unit(2,units = 'mm'),show_annotation_name = FALSE,show_legend = FALSE
)

#heatmap
ht <- Heatmap(as.matrix(emission[,2:6]), name = "Chromatin state", col=col_state,
              cluster_rows=FALSE, cluster_columns = FALSE, 
              show_row_names = FALSE,#column_labels = column_labels,
              left_annotation = right_annotation,
              right_annotation = ha,
              row_order = c(12,13,9,10,1,11,3,2,6,7,4,5,8),
              column_order = c("H3K4me3","H3K27ac","H3K4me1","H3K27me3","H3K9me3"),
              column_names_rot = 45,
              height = unit(10,'cm'),width = unit(4,'cm'),gap = unit(1,'mm'),
              heatmap_legend_param = list(title  = "State",title_position="topcenter",
                                          legend_width = unit(4,units = "cm"),legend_height = unit(0.3,units = "cm"),
                                          direction="horizontal")
)
pdf("emission state.pdf",width = 4.8,height = 5.5)
draw(ht,padding = unit(c(t = 0.1,r = 0.1,b = 0.1,l = 0.1),"cm"),heatmap_legend_side = "top")
dev.off()

# each overlapenrichment segment ------------------------------------------
file0 <- read.table("ctrl_13_segments.txt",sep = "\t",header = TRUE,row.names = 1)
file2 <- read.table("2weeks_13_segments.txt",sep = "\t",header = TRUE,row.names = 1)
file4 <- read.table("4weeks_13_segments.txt",sep = "\t",header = TRUE,row.names = 1)
file7 <- read.table("7weeks_13_segments.txt",sep = "\t",header = TRUE,row.names = 1)
file10 <- read.table("10weeks_13_segments.txt",sep = "\t",header = TRUE,row.names = 1)

plotHeatmap <- function(file,time) {
  names(file) <- c("Genome %","CpG","Exon","Gene","TES","TSS","TSS2kb")
  ht0 <- Heatmap(matrix = as.matrix(file[-14,]),
                 column_title = paste(time,"FE",sep = " "),
                 column_title_gp = gpar(face="sans",fontsize=14),
                 show_row_names = FALSE,
                 col = colorRamp2(c(0, 12), colors=color_list[c(6,11)]),
                 cluster_columns = FALSE,cluster_rows = FALSE,show_heatmap_legend = FALSE,
                 row_order = c(12,13,9,10,1,11,3,2,6,7,4,5,8),
                 column_order = c("Genome %","CpG","TSS2kb","TSS","Exon","TES","Gene"),
                 column_names_rot = 45,
                 height = unit(10,'cm'),width = unit(5.45,'cm'),gap = unit(1,'mm'),
                 heatmap_legend_param = list(title  = "2weeks FC",title_position="topcenter",
                                             legend_width = unit(4,units = "cm"),legend_height = unit(0.3,units = "cm"),
                                             direction="horizontal")
  )
  pdf(paste(time,"overlapenrichment.pdf"),width = 3.2,height = 5)
  draw(ht0,padding = unit(c(t = 0.1,r = 0.1,b = 0.1,l = 0.1),"cm"))#,ht_gap = unit(1,"mm"))
  dev.off()
}

plotHeatmap(file = file2,time = "2weeks")
plotHeatmap(file = file4,time = "4weeks")
plotHeatmap(file = file7,time = "7weeks")
plotHeatmap(file = file10,time = "10weeks")




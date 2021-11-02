

#setwd("~/projects/brain_aging_mouse/analysis/snapATAC/9mo+18mo/")
dat = list()
#setwd("../OLG/")
#files = list.files(".",pattern = ".txt",full.names = T)

files = list.files("age_diff_edgeR.snap",pattern = ".txt",full.names = T)
for (f in files){
  tmp = read.csv(f)
  tmp = tmp[which(tmp$PValue<0.01),]
  if(nrow(tmp) == 0) {
    next()
  }
  
  tmp = tmp[,-c(2,3,4,5)]
  if(ncol(tmp) != 6) {
    next()
  }
  outF =gsub(".edger.txt", "", f)
  outF =gsub("age_diff_edgeR.snap/", "", outF)
  tmp$id = outF 
  dat[[outF]] = tmp
}

diff = do.call(rbind, dat)
write.table(diff, "all_diff_peaks_concat.txt", sep ="\t", quote = F)
diff = read.table( "all_diff_peaks_concat.txt", sep ="\t")

pval = 0.001

t2 = table(diff$PValue < pval, diff$logFC>0, diff$id)
mt2 =melt(t2)

# # diff peak per cell type/tissue bar plot
setwd("age_diff_edgeR.snap/diff_peak_overlap/")
pdf(paste("all_tissues_number_of_diff_peaks_",pval,".pdf", sep = ""),height = 18,width = 10)
ggplot(subset(mt2,Var1==TRUE)) +
  #  geom_col(aes(x=factor(Var3, levels = rev.Vector(unique(diff$id ))),
  geom_col(aes(x=factor(Var3, levels = ),
               y=value,fill=Var2),position="dodge") +
  xlab("Cluster") +
  scale_fill_discrete(name="Down in Aging") +
  coord_flip()
dev.off()


diff = diff[which(diff$logFC<0),]
mat = matrix(0,nrow=length(unique(diff$id)),ncol=length(unique(diff$id)))
x = 1
for (i in unique(diff$id )){
  y=1
  set1 = diff[diff$id ==i,"X"] 
  #coors = set1
  # chr = sapply(strsplit(as.character(coors), ":"), "[[", 1)
  # start = sapply(strsplit(as.character(coors), ":"), "[[", 2)
  # end = as.numeric(paste(sapply(strsplit(as.character(start), "-"), "[[", 2)))
  # start = as.numeric(paste(sapply(strsplit(as.character(start), "-"), "[[", 1)))
  # gr1 <- GRanges(seqnames=chr,
  #                ranges=IRanges(start,end))
  for (j in unique(diff$id )){
    print(paste(i,j))
    set2 = diff[diff$id ==j,"X"]
    #coors = set2
    # chr = sapply(strsplit(as.character(coors), ":"), "[[", 1)
    # start = sapply(strsplit(as.character(coors), ":"), "[[", 2)
    # end = as.numeric(paste(sapply(strsplit(as.character(start), "-"), "[[", 2)))
    # start = as.numeric(paste(sapply(strsplit(as.character(start), "-"), "[[", 1)))
    # gr2 <- GRanges(seqnames=chr,
    #                ranges=IRanges(start,end))
    #mat[x,y] = length(GenomicRanges::intersect(gr1,gr2))/length(union(gr1,gr2))
    mat[x,y] = length(intersect(set1,set2))/length(union(set1,set2))
    y = y+1
  }
  x = x+1
}
pval = 0.01
#write.csv(mat,paste("All_tissues.diffpeak_share.table_p",pval,".csv", sep = ""))

#mat = read.csv("All_tissues.diffpeak_share.table_p.01.csv",row.names=1)
mat[which(mat==1, arr.ind = T)] = 0.03
mat[which(mat>0.03, arr.ind = T)] = 0.03
rownames(mat) = unique(diff$id)
colnames(mat) = unique(diff$id)
logmat = log10(mat+0.0001)

cluster = sapply(strsplit(as.character(rownames(mat)), "_"), "[[", 1)
annotation_col = data.frame(
   Cluster = paste0("c", cluster)
 )
rownames(annotation_col) = make.unique(as.character(rownames(mat)))

cluster_colors=c(
  "grey", "#E31A1C", "#FFD700", "#771122", "#777711", "#1F78B4", "#68228B", "#AAAA44",
  "#60CC52", "#771155", "#DDDD77", "#774411", "#AA7744", "#AA4455", "#117744", 
  "#000080", "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9", "#BC80BD", "#FFED6F",
  "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
  "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
  "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
  "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3",
  "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2",
  "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC",
  "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628",
  "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
  "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
  "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5"
)
my_colour = list(
  Cluster = c(c1 =cluster_colors[1],c2=cluster_colors[2] ,c3 = cluster_colors[3],
              c4 =cluster_colors[4],c5=cluster_colors[5] ,c6 = cluster_colors[6],
              c7 =cluster_colors[7],c8=cluster_colors[8] ,c9 = cluster_colors[9],
              c10 =cluster_colors[10],c11=cluster_colors[11] ,c12 = cluster_colors[12],
              c13 =cluster_colors[13],c14=cluster_colors[14] ,c15 = cluster_colors[15],
              c16 =cluster_colors[16],c17=cluster_colors[17] ,c18 = cluster_colors[18],
              c19 =cluster_colors[19],c20=cluster_colors[20] ,c21 = cluster_colors[21],
              c22 =cluster_colors[22],c23=cluster_colors[23] ,c24 = cluster_colors[24],
              c25 =cluster_colors[25],c26=cluster_colors[26] ,c27 = cluster_colors[27]))

# jaccard similarity for diff peaks bar plot
#setwd("/mnt/tscc/lamaral/projects/Aging/aging_share/figures/diffpeak_num_overlap/")
pdf(paste("all_tissues_jaccard_heat_ph_",pval,"_max0.01.pdf",sep = ""), height = 15, width = 15)
#pheatmap(logmat, main = paste( "All tissues Jaccard of differential peaks \n log10(jaccard+0.0001)"),legend_labels = "log10(jaccard +0.0001)",annotation_row =annotation_col,annotation_colors = my_colour)
pheatmap(mat, main = paste( "All tissues Jaccard of differential peaks"),annotation_row =annotation_col,annotation_colors = my_colour)#,annotation_colors = my_colour,show_colnames = F )
pheatmap(mat, main = paste( "All tissues Jaccard of differential peaks"), annotation_colors = my_colour,show_colnames = T,clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
         annotation_row =annotation_col)
dev.off()




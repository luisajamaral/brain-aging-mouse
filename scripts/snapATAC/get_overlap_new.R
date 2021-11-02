
library(data.table)
#setwd("~/projects/brain_aging_mouse/analysis/snapATAC/9mo+18mo/")
dat = list()
#setwd("../OLG/")
#files = list.files(".",pattern = ".txt",full.names = T)

files = list.files("age_diff_edgeR.snap",pattern = "8wk_vs_9mo.edger.txt",full.names = T)
for (f in files){
  tmp = fread(f,sep = ",")
  tmp = tmp[which(tmp$fdr<0.05),]
  if(nrow(tmp) == 0) {
    next
  }
  #tmp = tmp[,-c(2,3,4,5)]
  rn = tmp$V1
  tmp  = tmp[,(ncol(tmp)-4):(ncol(tmp))]
  tmp$peak = rn
  outF =gsub(".edger.txt", "", f)
  outF =gsub("age_diff_edgeR.snap/", "", outF)
  tmp$id = outF 
  dat[[outF]] = tmp
}

diff = do.call(rbind, dat)
write.table(diff, "all_diff_peaks_concat.txt", sep ="\t", quote = F)
diff = read.table( "all_diff_peaks_concat.txt", sep ="\t")


d = table(diff$id)
names(d) = gsub("_8wk_vs_9mo", "", names(d))
names(d)[which(names(d) %in% c(13))]="CA1GL"
names(d)[which(names(d) %in% c(31))]="CA3GL"
names(d)[which(names(d) %in% c(12))]="PVGA"
names(d)[which(names(d) %in% c(15))]="ITHGL"
names(d)[which(names(d) %in% c(14))]="L6bGL"
#names(d)[which(names(d) %in% c(25))]="OBGA2"
#names(d)[which(names(d) %in% c(17))]="OLFGL"
names(d)[which(names(d) %in% c(21))]="ITL5GL"
names(d)[which(names(d) %in% c(22))]="ITL5GL2"

names(d)[which(names(d) %in% c(26))]="PER"

names(d)[which(names(d) %in% c(24))]="VEC"
names(d)[which(names(d) %in% c(7))]="VLMC"
names(d)[which(names(d) %in% c(16))]="VPIA"
m = cbind(table(cluster)[names(d)],d)

colnames(m) = c("#cells", "# DAR")
plot(m,pch=19, cex=2,col="lightblue", main = "8wk vs 9mo")
text(m, labels=rownames(m),data=m, cex=0.9, font=2)

num_diff = c(0,0,2,56,480,1040, 1414)
num_cells = c(1000,2500,5000,25000,50000,75000,100333)

num_diff = c(5,214,761,13513, 27656, 37502, 43650, 44059)
num_cells = c(1000,2500,5000,25000,50000,75000,100000,100333)

num_diff = c(0,57,421,8563,15991,24571,31417,31956)
num_cells = c(1000,2500,5000,25000,50000,75000,100000,100333)


num_diff = c(3,19,92,995,5888,16175,37106)
num_cells = c(250,500,1000,2500,5000,10000,25000)

#8wk vs 9mo
num_diff = c(0,1,17,331,2421, 9613, 27771)

num_diff = c(0,1,2,10,24,116, 950)

num_cells = c(250,500,1000,2500,5000,10000,15000,20000,25000)

num_diff = c(3,19,110,1200,5500,16510,24291, 31779, 37000)
num_diff = c(0,1,15,400,2420,9700,16140,22589,27777)
num_diff = c(0,0,2,9,24,116,305,748,1053)

m = cbind(num_cells, num_diff)
colnames(m) = c("#cells", "# DAR")

plot(m,pch=19, cex=2,col="lightblue", main = "9mo vs 18mo OLG cluster # DAR / # of cells sampled")



fdr = c(0.503281866641127, 0.179234133003917, 0.251052000326504, 0.000112065166832834, 6.80309491922082e-26,6.08654336100654e-19,9.51098154834489e-32)


m = cbind(num_cells, fdr)
colnames(m) = c("#cells", "fdr")

plot(m,pch=19, cex=2,col="lightblue", main = "9mo vs 18mo OLG cluster fdr of Hist1H locus")



pval = 0.01

t2 = table(diff$fdr < pval, diff$logFC>0, diff$id)
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

bigdiff = diff
diff = bigdiff[which(bigdiff$logFC>0),]
mat = matrix(0,nrow=length(unique(diff$id)),ncol=length(unique(diff$id)))
x = 1
for (i in unique(diff$id )){
  y=1
  set1 = diff[diff$id ==i,"peak"] 
  #coors = set1
  # chr = sapply(strsplit(as.character(coors), ":"), "[[", 1)
  # start = sapply(strsplit(as.character(coors), ":"), "[[", 2)
  # end = as.numeric(paste(sapply(strsplit(as.character(start), "-"), "[[", 2)))
  # start = as.numeric(paste(sapply(strsplit(as.character(start), "-"), "[[", 1)))
  # gr1 <- GRanges(seqnames=chr,
  #                ranges=IRanges(start,end))
  for (j in unique(diff$id )){
    print(paste(i,j))
    set2 = diff[diff$id ==j,"peak"]
    #coors = set2
    # chr = sapply(strsplit(as.character(coors), ":"), "[[", 1)
    # start = sapply(strsplit(as.character(coors), ":"), "[[", 2)
    # end = as.numeric(paste(sapply(strsplit(as.character(start), "-"), "[[", 2)))
    # start = as.numeric(paste(sapply(strsplit(as.character(start), "-"), "[[", 1)))
    # gr2 <- GRanges(seqnames=chr,
    #                ranges=IRanges(start,end))
    #mat[x,y] = length(GenomicRanges::intersect(gr1,gr2))/length(union(gr1,gr2))
    mat[x,y] = length(intersect(set1$peak,set2$peak))/length(union(set1$peak,set2$peak))
    y = y+1
  }
  x = x+1
}
pval = 0.01
#write.csv(mat,paste("All_tissues.diffpeak_share.table_p",pval,".csv", sep = ""))

#mat = read.csv("All_tissues.diffpeak_share.table_p.01.csv",row.names=1)
mat[which(mat==1, arr.ind = T)] = 0.1
mat[which(mat>0.1, arr.ind = T)] = 0.1
rownames(mat) = unique(diff$id)
colnames(mat) = unique(diff$id)
logmat = log10(mat+0.0001)

cluster = sapply(strsplit(as.character(rownames(mat)), "_"), "[[", 1)
cluster[which(cluster %in% c(28,29))]="D2MSN"
cluster[which(cluster %in% c(13))]="CA1GL"
cluster[which(cluster %in% c(31))]="CA3GL"
cluster[which(cluster %in% c(12))]="PVGA"
cluster[which(cluster %in% c(15))]="ITHGL"
cluster[which(cluster %in% c(14))]="L6bGL"
#cluster[which(cluster %in% c(25))]="OBGA2"
#cluster[which(cluster %in% c(17))]="OLFGL"
cluster[which(cluster %in% c(21))]="ITL5GL"
cluster[which(cluster %in% c(22))]="ITL5GL2"

cluster[which(cluster %in% c(26))]="PER"

cluster[which(cluster %in% c(24))]="VEC"
cluster[which(cluster %in% c(7))]="VLMC"
cluster[which(cluster %in% c(16))]="VPIA"

annotation_col = data.frame(
  Cluster = make.names(cluster)
)

rownames(mat) = make.names(cluster)
colnames(mat) = make.names(cluster)
rownames(annotation_col) = make.unique(as.character(rownames(mat)))

down_mat = mat
up_mat = mat
# jaccard similarity for diff peaks bar plot
#setwd("/mnt/tscc/lamaral/projects/Aging/aging_share/figures/diffpeak_num_overlap/")
pdf(paste("all_tissues_jaccard_heat_ph_8wk_vs_9mo_",pval,"_max0.01.pdf",sep = ""), height = 8, width = 8)
#pheatmap(logmat, main = paste( "All tissues Jaccard of differential peaks \n log10(jaccard+0.0001)"),legend_labels = "log10(jaccard +0.0001)",annotation_row =annotation_col,annotation_colors = my_colour)
pheatmap(down_mat, main = paste( "All tissues Jaccard of differential peaks (Down peaks)"))
pheatmap(up_mat, main = paste( "All tissues Jaccard of differential peaks (Up peaks)"))

dev.off()




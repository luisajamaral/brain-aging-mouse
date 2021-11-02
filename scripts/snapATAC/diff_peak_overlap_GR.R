################## count the number of differential peaks using the same cut-off.
library(pheatmap)
library(reshape2)
library(ggplot2)
library(GenomicRanges)
library(data.table)
library(RColorBrewer)

pval = as.numeric(commandArgs(trailing=T)[1])
cat(pval,"\n")
key = read.csv("/mnt/silencer2/home/shz254/projects/mouse_aging/aging_share/figures/celltype_annotation.txt", sep = "")
key$id = paste0(key$Tissue,"_",key$Name)

# gets differential peak data for all cell types in tissue in one table
getDiffTable = function(key, pval) {
  dat = list()
  for ( i in 1:nrow(key)) {
    cat(i)
    tissue = key$Tissue[i]
    cluster = key$Name[i]
    clade = key$Clade[i]
    clustern = key$Cluster[i]
    id = key$id[i]
    curr = fread(paste0("/mnt/tscc/yanxiao/projects/mouse_aging/analysis/snapATAC/",tissue,"/age_diff_edgeR.snap/",clustern,".edger.txt"))
    curr$cluster = clustern
    curr$cell_type = cluster
    curr$clade = clade
    curr$tissue = tissue
    curr$name = key$id[i]
    # Select the top 1% of peaks as diff. 
    #    dat[[cl]] = dat[[cl]][which(dat[[cl]]$PValue < quantile(dat[[cl]]$PValue,0.01)),]
    curr = curr[which(curr$PValue<pval),]
    colnames(curr) = gsub(paste(tissue,"_",sep = ""), "", colnames(curr))
    dat[[key$id[i]]] = curr
  }
  big_dat = do.call(rbind, dat)
  return(big_dat)
}
#all tissues diff table
diff = getDiffTable(key, pval)


diff$tissue_celltype = paste0(diff$tissue,"_",diff$cell_type, sep="")
diff$tissue_cluster = paste0(diff$tissue,"_",diff$cluster, sep="")

peak_cnt = data.frame(table(diff$tissue_celltype))
peak_cnt = peak_cnt[order(-peak_cnt$Freq),]

t2 = table(diff$PValue < pval, diff$logFC>0, diff$tissue_celltype, diff$tissue)
mt2 =melt(t2)

# # diff peak per cell type/tissue bar plot
setwd("/mnt/tscc/lamaral/projects/Aging/aging_share/diff_peak_overlap/")

pdf(paste("all_tissues_number_of_diff_peaks_",pval,".pdf", sep = ""),height = 12,width = 10)
ggplot(subset(mt2,Var1==TRUE)) +
  #  geom_col(aes(x=factor(Var3, levels = rev.Vector(unique(diff$tissue_celltype))),
  geom_col(aes(x=factor(Var3, levels = rev(peak_cnt$Var1)),
               y=value,fill=Var2),position="dodge") +
  xlab("Cluster") +
  scale_fill_discrete(name="Up in Aging") +
  coord_flip()
dev.off()


mat = matrix(0,nrow=length(unique(diff$tissue_celltype)),ncol=length(unique(diff$tissue_celltype)))
x = 1
for (i in unique(diff$tissue_celltype)){
  y=1
  set1 = diff[diff$tissue_celltype==i,"V1"] 
  coors = set1$V1
  chr = sapply(strsplit(as.character(coors), ":"), "[[", 1)
  start = sapply(strsplit(as.character(coors), ":"), "[[", 2)
  end = as.numeric(paste(sapply(strsplit(as.character(start), "-"), "[[", 2)))
  start = as.numeric(paste(sapply(strsplit(as.character(start), "-"), "[[", 1)))
  gr1 <- GRanges(seqnames=chr,
                 ranges=IRanges(start,end))
  for (j in unique(diff$tissue_celltype)){
    print(paste(i,j))
    set2 = diff[diff$tissue_celltype==j,"V1"]
    coors = set2$V1
    chr = sapply(strsplit(as.character(coors), ":"), "[[", 1)
    start = sapply(strsplit(as.character(coors), ":"), "[[", 2)
    end = as.numeric(paste(sapply(strsplit(as.character(start), "-"), "[[", 2)))
    start = as.numeric(paste(sapply(strsplit(as.character(start), "-"), "[[", 1)))
    gr2 <- GRanges(seqnames=chr,
                   ranges=IRanges(start,end))
    mat[x,y] = length(GenomicRanges::intersect(gr1,gr2))/length(union(gr1,gr2))
    #mat[x,y] = length(intersect(set1,set2))/length(union(set1,set2))
    y = y+1
  }
  x = x+1
}

melted = melt(mat)

write.csv(mat,paste("All_tissues.diffpeak_share.table_p",pval,".csv", sep = ""))
mat = read.csv("All_tissues.diffpeak_share.table_p.001.csv",row.names=1)
mat[which(mat==1, arr.ind = T)] = 0.01
mat[which(mat>0.01, arr.ind = T)] = 0.01
rownames(mat) = unique(diff$tissue_celltype)
colnames(mat) = unique(diff$tissue_celltype)
logmat = log10(mat+0.0001)

annotation_col = data.frame(
  Tissue = key$Tissue,
  Clade = key$Clade
)
rownames(annotation_col) = make.unique(as.character(key$id))
cols = brewer.pal(6, "Set1")
cols2 = brewer.pal(4, "Set2")
my_colour = list(
  Tissue = c(DH = cols2[1], FC = cols2[2],LM = cols2[3],HT = cols2[4]),
  Clade = c(ExN = cols[1], Glia = cols[2], Immune = cols[3], InN = cols[4], Muscle = cols[5], Other = "gray58")
)

# jaccard similarity for diff peaks bar plot
#setwd("/mnt/tscc/lamaral/projects/Aging/aging_share/figures/diffpeak_num_overlap/")
pdf(paste("all_tissues_jaccard_heat_ph_",pval,"_max0.01.pdf",sep = ""), height = 10, width = 13)
#pheatmap(logmat, main = paste( "All tissues Jaccard of differential peaks \n log10(jaccard+0.0001)"),legend_labels = "log10(jaccard +0.0001)")
pheatmap(mat, main = paste( "All tissues Jaccard of differential peaks"),annotation_col =annotation_col,annotation_colors = my_colour,show_colnames = F )
pheatmap(mat, main = paste( "All tissues Jaccard of differential peaks"), annotation_colors = my_colour,show_colnames = F,
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",annotation_col =annotation_col)
dev.off()

# 
# reorder_cormat <- function(cormat){
#   # Use correlation between variables as distance
#   dd <- as.dist((1-cormat)/2)
#   hc <- hclust(dd)
#   cormat <-cormat[hc$order, hc$order]
# }
# get_upper_tri <- function(cormat){
#   cormat[lower.tri(cormat)]<- NA
#   return(cormat)
# }
# cormat <- reorder_cormat(mat)
# upper_tri <- get_upper_tri(cormat)
# # Melt the correlation matrix
# melted_cormat <- melt(upper_tri, na.rm = TRUE)
# # Create a ggheatmap
# ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
#   geom_tile(color = "white")+
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                        midpoint = 0, limit = c(0,0.02), space = "Lab", 
#                        name="Pearson\nCorrelation") +
#   theme_minimal()+ # minimal theme
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                    size = 12, hjust = 1))+
#   coord_fixed()
# # Print the heatmap
# print(ggheatmap)


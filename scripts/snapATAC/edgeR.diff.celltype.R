library(SnapATAC)
library(RColorBrewer)
library(ggplot2)
library(gplots)
library(ggrepel)
library(pheatmap)
library(dplyr)


c25 <- c(
  "dodgerblue2",
  "#E31A1C",
  "green4",
  "#6A3D9A",
  "#FF7F00",
  "gold1"
)
tissue = "9mo+18mo"
RData = "all.cluster.pmat.RData"
setwd(paste0("../../analysis/snapATAC/", tissue))

load(RData)

#system("mkdir -p age_diff_edgeR.snap_celltype_celltype/")

#x.sp =  addPmatToSnap(x.sp,do.par = T,num.cores = 6)
#save(x.sp,file = "all.cluster.pmat.RData")

#peaks_table = read.table(peaks_file)
# library(GenomicRanges)
# peaks_table = read.table(peaks_file)
# peaks.gr = GRanges(
#   peaks_table[,1],
#   IRanges(peaks_table[,2], peaks_table[,3]), name = peaks_table[,4], score = peaks_table[,5]
# );
# x.sp = createPmat(x.sp, peaks = peaks.gr, num.cores = 10)

pmat = x.sp@pmat
peak = x.sp@peak$name
cluster = x.sp@cluster
x.sp@metaData$sample = gsub("[+]", "_",x.sp@metaData$sample)
#sample = x.sp@sample
sample = x.sp@metaData$sample
stage = x.sp@metaData$stage
#tissue = gsub("18mo_", "", sample)
#tissue = gsub("9mo_", "", tissue)
#tissue = gsub("_1$", "", tissue)
#tissue = gsub("_2$", "", tissue)


cluster_sample = paste(cluster, sample, sep = ".")
#rownames(pmat) = cluster_sample
#cl = "6"
max_cluster = max(as.numeric(cluster))
#max_cluster = 15
library(edgeR)
library(doParallel)
registerDoParallel(cores = max_cluster)

#cluster = rep("",length(cluster))
#cluster[which(x.sp@cluster%in%c(3,20))] = "OLG"
#cluster_sample = paste(cluster, sample, sep = ".")

foreach(cl = 1:max_cluster) %dopar% {
#for (cl in 1:max_cluster) {
  print(cl)
  samples = unique(sample)
  dat = list()
  for (s in samples) {
    idx1 = which(cluster_sample == paste(cl, s, sep = "."))
    if (length(idx1) < 2) {
      cat(cl, s, " BAD\n")
      next
    } else {
      dat[[s]] = colSums(pmat[idx1, ])
    }
  }
  mat = do.call(cbind, dat)
  if(is.null(mat) || ncol(mat) < 4) {
    cat(cl, " BAD 2\n")
  } else {
    rownames(mat) = peak
    
    grps = substr(colnames(mat), 1, 2)
    tissue = gsub("18mo_", "", colnames(mat))
    tissue = gsub("9mo_", "", tissue)
    tissue = gsub("_1$", "", tissue)
    tissue = gsub("_2$", "", tissue)
    
    y = DGEList(mat)
    y = y[which(rowSums(cpm(y) > 1) >= 2), ]
    y = calcNormFactors(y)
    design = model.matrix( ~tissue+grps)
    y <- estimateCommonDisp(y)
    y <- estimateGLMTagwiseDisp(y, design)
    fit_tag = glmFit(y, design)
    lrt = glmLRT(fit_tag)#, coef = "grps9m")
    fdr = p.adjust(lrt$table$PValue, method = "BH")
    out = cbind(cpm(y), lrt$table, fdr)
    out = out[order(out$PValue), ]
    show(head(out))
    write.csv(out,
              paste0("age_diff_edgeR.snap_celltype/", cl, ".edger.txt"))
    pdf(paste(cl,"_heatmaps.pdf", sep = ""), width = 7, height = 7)
    
    o = out[which(out$fdr<0.05),1:(ncol(out)-5)]
    
    hr <- hclust(as.dist(1 - cor(t(o))), method = "ward.D2")
    
    col = grps
    col[which(col=="9m")] = "blue"
    col[which(col=="18")] = "red"
    
    my_palette = colorpanel(100, "darkblue", "white", "red")
    heatmap.2(
      data.matrix(o),
      cexCol = 1,
      margins = c(15, 5),
      Rowv = as.dendrogram(hr),
      Colv = T,
      col = my_palette,
      scale = "row",
      labRow = F,
      trace = "none",
      #RowSideColors = mycolhc,
      ColSideColors = col,
      main = paste("cluster ",cl," (", nrow(o), " fdr <0.05)", sep =
                     "")
    )
    
    
    o = out[1:50,1:(ncol(out)-5)]
    
    hr <- hclust(as.dist(1 - cor(t(o))), method = "ward.D2")
    
    ##Forcing 4 groups (can change)
    mycl <- cutree(hr, k = 4)
    mycolhc <- c25[1:length(unique(mycl))]
    mycolhc <- mycolhc[as.vector(mycl)]
    col = grps
    col[which(col=="9m")] = "blue"
    col[which(col=="18")] = "red"
    
    my_palette = colorpanel(100, "darkblue", "white", "red")
    heatmap.2(
      data.matrix(o),
      cexCol = 1,
      margins = c(15, 15),
      Rowv = as.dendrogram(hr),
      Colv = T,
      col = my_palette,
      scale = "row",
      trace = "none",
      #RowSideColors = mycolhc,
      ColSideColors = col,
      main = paste("cluster ",cl," (top 50)", sep =
                     "")
    )
    
    dev.off()
    
  }
}

## write differential peaks into bed files.
#files = list.files("age_diff_edgeR.snap_celltype",pattern = ".txt",full.names = T)
for (cl in 1:max_cluster){
  f = paste("age_diff_edgeR.snap_celltype/", cl, ".edger.txt",sep ="")
  tmp = read.csv(f)
  sig = tmp[which(tmp$PValue < 0.01), ]
  #  sig = tmp[which(tmp$PValue< quantile(tmp$PValue,0.01)),]
  chr = sub("(chr.*):(.*)-(.*)", "\\1", sig$X)
  start = sub("(chr.*):(.*)-(.*)", "\\2", sig$X)
  end = sub("(chr.*):(.*)-(.*)", "\\3", sig$X)
  out = data.frame(chr, start, end, sig$logFC, as.integer(-log10(sig$PValue)))
  #  write.table(out, paste0("../age_diff_edgeR.snap_celltype/",cl,".both.bed"),row.names=F,col.names=F,quote=F,sep="\t")
  outF =gsub(".txt", "", f)
  write.table(
    subset(out, sig.logFC > 0),
    paste0(outF,".up.bed"),
    row.names = F,
    col.names = F,
    quote = F,
    sep = "\t"
  )
  write.table(
    subset(out, sig.logFC < 0),
    paste0(outF, ".down.bed"),
    row.names = F,
    col.names = F,
    quote = F,
    sep = "\t"
  )
  
}

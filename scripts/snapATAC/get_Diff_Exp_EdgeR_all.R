library(SnapATAC)
library(RColorBrewer)
library(ggplot2)
library(gplots)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(tictoc)


c25 <- c(
  "dodgerblue2",
  "#E31A1C",
  "green4",
  "#6A3D9A",
  "#FF7F00",
  "gold1"
)
tissue = "all66"
RData = "all.cluster.pmat.RData"
setwd(paste0("../../analysis/snapATAC/", tissue))

load(RData)

#system("mkdir -p age_diff_edgeR.snap/")
#tic("add Pmat")
#x.sp =  addPmatToSnap(x.sp,do.par = T,num.cores = 8)
#save(x.sp,file = "all.cluster.pmat.RData")
#toc()

pmat = x.sp@pmat
peak = x.sp@peak$name
cluster = x.sp@metaData$cluster_names
cluster[which(cluster %in% c(28,29))]="D2MSN"
cluster[which(cluster %in% c(13))]="CA1GL"
cluster[which(cluster %in% c(31))]="CA3GL"
cluster[which(cluster %in% c(12))]="PVGA"
cluster[which(cluster %in% c(15))]="ITHGL"
cluster[which(cluster %in% c(14))]="L6bGL"
#cluster[which(cluster %in% c(25))]="OBGA2"
#cluster[which(cluster %in% c(17))]="OLFGL"
cluster[which(cluster %in% c(21,22))]="ITL5GL"
cluster[which(cluster %in% c(26))]="PER"

cluster[which(cluster %in% c(24))]="VEC"
cluster[which(cluster %in% c(7))]="VLMC"
cluster[which(cluster %in% c(16))]="VPIA"

cluster[which(cluster %in% c("CRC"))]="19"

cluster[which(cluster %in% c("LAMGA"))]="9"


plotViz(
  obj= x.sp,
  method="umap",
  main="Cluster",
  point.color=cluster,
  point.size=0.2,
  point.shape=19,
  text.add=T,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)







x.sp@metaData$cluster_names = cluster
#x.sp@metaData$sample = gsub("[+]", "_",x.sp@metaData$sample)
#sample = x.sp@sample
sample = x.sp@metaData$sample_name
stage = x.sp@metaData$stage
tissue = x.sp@metaData$tissue
tissue[which(tissue %in% c("3F","4E"))] = "3F+4E"
tissue[which(tissue %in% c("5E", "6E"))] = "5E+6E"
tissue[which(tissue %in% c("7H","8H", "9G"))] = "7H+8H+9G"
tissue[which(tissue %in% c("13D","14C"))] = "13D+14C"
tissue[which(tissue %in% c("11E","11F","12E"))] ="11E+11F+12E"
tissue[which(tissue %in% c("8E","9H","8J","9J"))] = "8E+9H+8J+9J"

cluster_sample = paste(cluster, sample, sep = ".")
age_cluster_tissue_rep = paste(cluster,stage,tissue, x.sp@metaData$replicate, sep = "_")
age_tissue_rep = paste(stage,tissue, x.sp@metaData$replicate, sep = "_")
cluster_tissue = paste(cluster,tissue,sep = "_")
rm(x.sp)

#rownames(pmat) = cluster_sample
#cl = "6"
max_cluster = length(unique(cluster))
#max_cluster = 15
library(edgeR)
library(doParallel)
get_diff <- function(st1, st2, bigmat, cl) {
  mat = bigmat[,c(grep(st1, colnames(bigmat)), grep(st2, colnames(bigmat)))]
  grps = substr(colnames(mat), 1, 3)
  tissue = gsub("18mo_", "", colnames(mat))
  tissue = gsub("9mo_", "", tissue)
  tissue = gsub("8wk_", "", tissue)
  tissue = gsub("_rep1$", "", tissue)
  tissue = gsub("_rep2$", "", tissue)
  cat(colnames(mat), "\n")
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
            paste0("age_diff_edgeR.snap/", cl, "_",st1,"_vs_",st2, ".edger.txt"))
  return(out)
}

registerDoParallel(cores = 8)


foreach(i = 1:max_cluster) %dopar% {
#foreach(i = c(25,31,"D2MSN")) %dopar% {
#for (i in c("D2MSN",31,25)) {
  print(i)
  cl = unique(cluster)[i]
  #cl = i
  print(cl)
  dat = list()
  for (s in unique(age_tissue_rep)) {
    idx1 = which(age_cluster_tissue_rep == paste(cl, s, sep = "_"))
    cl_tiss = paste(cl,strsplit(s,"_")[[1]][2], sep = "_")
    if (length(idx1) < 200 || as.numeric(table(cluster_tissue)[cl_tiss])<600) {
      cat(cl, s, " not enough cells\n")
      next
    } else {
      dat[[s]] = colSums(pmat[idx1, ])
    }
  }
  mat = do.call(cbind, dat)
  if(is.null(mat) || ncol(mat) < 6) {
    cat(cl, " not enough samples 2\n")
  } else {
    
    rownames(mat) = peak
    o_8wkvs18mo = get_diff("8wk","18mo",mat, cl)
    o_8wkvs9mo = get_diff("8wk","9mo",mat, cl)
    o_9movs18mo = get_diff("9mo","18mo",mat, cl)
    
    
    peaks = unique(c(rownames(o_8wkvs18mo)[which(o_8wkvs18mo$fdr<0.05)],rownames(o_8wkvs9mo)[which(o_8wkvs9mo$fdr<0.05)],rownames(o_9movs18mo)[which(o_9movs18mo$fdr<0.05)]))
    y = DGEList(mat)
    y = y[which(rowSums(cpm(y) > 1) >= 4), ]
    
    out = cpm(y)
    write.csv(out, file = paste(cl,"cpm.txt", sep = "_"), quote = F)
    pdf(paste(cl,"_heatmaps.pdf", sep = ""), width = 12, height = 7)
    main =paste("cluster ",cl," (", length(peaks), " fdr <0.05)", sep = "")
    
    if (length(peaks)>5000) {
      peaks = unique(c(rownames(o_8wkvs18mo)[which(o_8wkvs18mo$fdr<0.00001)],rownames(o_8wkvs9mo)[which(o_8wkvs9mo$fdr<0.00001)],rownames(o_9movs18mo)[which(o_9movs18mo$fdr<0.00001)]))
      main =paste("cluster ",cl," (", length(peaks), " fdr <0.00001)", sep = "")
    }
    o = out[peaks,]
    
    hr <- hclust(as.dist(1 - cor(t(o))), method = "ward.D2")
    
   
    o = o[,order(colnames(o))]
    oder = c(colnames(o)[grep("8wk", colnames(o))],colnames(o)[grep("9mo", colnames(o))],colnames(o)[grep("18mo", colnames(o))])
    o = o[,oder]
    col = substr(colnames(o), 1,3)
    col[which(col=="9mo")] = "blue"
    col[which(col=="18m")] = "red"
    col[which(col=="8wk")] = "yellow"
    my_palette = colorpanel(100, "darkblue", "white", "red")
    heatmap.2(
      data.matrix(o),
      cexCol = 1,
      margins = c(15, 5),
      Rowv = as.dendrogram(hr),
      Colv = F,dendrogram = "row",
      col = my_palette,
      scale = "row",
      labRow = F,
      trace = "none",
      #RowSideColors = mycolhc,
      ColSideColors = col,
      main = main
    )
    legend("topright", title = "age",legend=c("8wk","9mo", "18mo"), 
           fill=c("yellow","blue","red"), cex=0.8, box.lty=0)
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
      main = main
    )
    legend("topright", title = "age",legend=c("8wk","9mo", "18mo"), 
           fill=c("yellow","blue","red"), cex=0.8, box.lty=0)
    peaks = unique(c(rownames(o_8wkvs18mo)[1:50],rownames(o_8wkvs9mo)[1:50],rownames(o_9movs18mo)[1:50]))
    
    o = out[peaks,]
    
    hr <- hclust(as.dist(1 - cor(t(o))), method = "ward.D2")
    o = o[,order(colnames(o))]
    oder = c(colnames(o)[grep("8wk", colnames(o))],colnames(o)[grep("9mo", colnames(o))],colnames(o)[grep("18mo", colnames(o))])
    o = o[,oder]
    col = substr(colnames(o), 1,3)
    col[which(col=="9mo")] = "blue"
    col[which(col=="18m")] = "red"
    col[which(col=="8wk")] = "yellow"
    my_palette = colorpanel(100, "darkblue", "white", "red")
    heatmap.2(
      data.matrix(o),
      cexCol = 1,
      margins = c(15, 15),
      Rowv = as.dendrogram(hr),
      Colv = F,dendrogram = "row",
      col = my_palette,
      scale = "row",
      trace = "none",
      #RowSideColors = mycolhc,
      ColSideColors = col,
      main = paste("cluster ",cl," (top 50)", sep =
                     "")
    )
    legend("topright", title = "age",legend=c("8wk","9mo", "18mo"), 
           fill=c("yellow","blue","red"), cex=0.8, box.lty=0)
    
    dev.off()
    
  }
}

## write differential peaks into bed files.
#files = list.files("age_diff_edgeR.snap_celltype",pattern = ".txt",full.names = T)
for (cl in 1:max_cluster){
  f = paste("age_diff_edgeR.snap/", cl, ".edger.txt",sep ="")
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





# 
# tab = data.frame(sample=age_tissue_rep,cluster=cluster,tissue = tissue)
# summ = plyr::count(tab[,1])
# tab2 = plyr::count(tab)
# tab2$total = summ$freq[match(tab2$sample,summ$x)]
# tab2$frac = tab2$freq/tab2$total
# tab2$age = substr(tab2$sample,1,3)
# tab2$age = factor(tab2$age,levels = c("8wk","9mo","18m"), ordered = T)
# tab2 = tab2[order(tab2$age),]
# pdf("celltype_frac_per_tissue.pdf",height=5,width=10)
# for(i in unique(tissue)) {
#   p = ggplot(tab2[which(tab2$tissue==i),]) + 
#     geom_col(aes(x=cluster,y=frac,fill=age, color = sample),position="dodge") +
#     theme_bw() +ggtitle(i)
#   print(p)
# }
# dev.off()

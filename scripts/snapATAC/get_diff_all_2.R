


library(SnapATAC)
library(RColorBrewer)
library(ggplot2)
library(gplots)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(tictoc)
library(doParallel)


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
#x.sp@metaData$sample = gsub("[+]", "_",x.sp@metaData$sample)
#sample = x.sp@sample
big_sample = x.sp@metaData$sample_name
stage = x.sp@metaData$stage
tissue = x.sp@metaData$tissue
rm(x.sp)

cluster_sample = paste(cluster, big_sample, sep = ".")
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
  #show(head(out))
  write.csv(out,
            paste0("age_diff_edgeR.snap_2/", cl, "_",st1,"_vs_",st2, ".edger.txt"))
  return(out)
}
max_cluster = length(unique(cluster))
registerDoParallel(cores = 10)
big_pmat = pmat
#for (cl  in unique(cluster)) {
foreach(i = 1:max_cluster) %dopar% {
  cl = unique(cluster)[i]
  cat(cl, "\n")
  sample = big_sample[which(cluster==cl)]
  pmat = big_pmat[which(cluster==cl),]
  
  samples = unique(sample)
  dat = list()
  for (s in samples) {
    idx1 = which(sample == s)
    if (length(idx1) < 2) {
      cat(cl, s, " BAD\n")
      next
    } else {
      dat[[s]] = colSums(pmat[idx1, ])
    }
  }
  mat = do.call(cbind, dat)
  if(is.null(mat) || ncol(mat) < 6) {
    cat(cl, " BAD 2\n")
  } else {
    
    rownames(mat) = peak
    o_8wkvs18mo = get_diff("8wk","18mo",mat, cl)
    o_8wkvs9mo = get_diff("8wk","9mo",mat, cl)
    o_9movs18mo = get_diff("9mo","18mo",mat, cl)
    
    
    cat("8wk vs 18 mo ",length(rownames(o_8wkvs18mo)[which(o_8wkvs18mo$fdr<0.05)]))
    cat("\n 8wk vs 9 mo ",length(rownames(o_8wkvs9mo)[which(o_8wkvs9mo$fdr<0.05)]))
    cat("\n 9mo vs 18 mo ",length(rownames(o_9movs18mo)[which(o_9movs18mo$fdr<0.05)]))
    cat("\n_____________________\n")
    
  }
}


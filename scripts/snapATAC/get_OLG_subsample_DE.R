


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
#x.sp@metaData$sample = gsub("[+]", "_",x.sp@metaData$sample)
#sample = x.sp@sample
sample = x.sp@metaData$sample_name
stage = x.sp@metaData$stage
tissue = x.sp@metaData$tissue
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

cluster_sample = paste(cluster, sample, sep = ".")
#rownames(pmat) = cluster_sample
#cl = "6"
max_cluster = length(unique(cluster))
#max_cluster = 15
library(edgeR)
library(doParallel)
get_diff <- function(st1, st2, bigmat, cl, cell_num) {
  mat = bigmat[,c(grep(st1, colnames(bigmat)), grep(st2, colnames(bigmat)))]
  grps = substr(colnames(mat), 1, 3)
  tissue = gsub("18mo_", "", colnames(mat))
  tissue = gsub("9mo_", "", tissue)
  tissue = gsub("8wk_", "", tissue)
  tissue = gsub("_rep1$", "", tissue)
  tissue = gsub("_rep2$", "", tissue)
  
  y = DGEList(mat)
  y = y[which(rowSums(cpm(y) > 1) >= 4), ]
  y = calcNormFactors(y)
  design = model.matrix( ~tissue+grps)
  
  #design = model.matrix( ~0+grps)
  y <- estimateCommonDisp(y)
  y <- estimateGLMTagwiseDisp(y, design)
  fit_tag = glmFit(y, design)
  lrt = glmLRT(fit_tag)#, coef = "grps9m")
  
  fdr = p.adjust(lrt$table$PValue, method = "BH")
  out = cbind(cpm(y), lrt$table, fdr)
  out = out[order(out$PValue), ]
  #show(head(out))
  write.csv(out,
            paste0("age_diff_edgeR.snap/",cell_num,"_", cl, "_",st1,"_vs_",st2, ".edger.txt"))
  return(out)
}

registerDoParallel(cores = 3)
i = "OLG"
cl = "OLG"

sample = sample[which(cluster=="OLG")]
pmat = pmat[which(cluster=="OLG"),]
age_cluster_tissue_rep_big = age_cluster_tissue_rep[which(cluster=="OLG")]
age_tissue_rep_big = age_tissue_rep[which(cluster=="OLG")]
big_pmat = pmat
cl = "OLG"

for (num_cells in c(1000,2500,5000,10000,25000,50000,75000,100000, length(age_cluster_tissue_rep_big))) {
  cat(num_cells, "\n")
  samp = sample.int(nrow(big_pmat), num_cells)
  pmat = big_pmat[samp,]
  sample_samp = sample[samp]
  age_cluster_tissue_rep = age_cluster_tissue_rep_big[samp]
  age_tissue_rep_big = age_tissue_rep_big[samp]
  dat = list()
  for (s in unique(age_tissue_rep)) {
    #idx1 = which(sample_samp == s)
    idx1 = which(age_cluster_tissue_rep == paste(cl, s, sep = "_"))
    
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
    o_8wkvs18mo = get_diff("8wk","18mo",mat, cl, num_cells)
    o_8wkvs9mo = get_diff("8wk","9mo",mat, cl,num_cells)
    o_9movs18mo = get_diff("9mo","18mo",mat, cl,num_cells)
    
    
    cat("8wk vs 18 mo ",length(rownames(o_8wkvs18mo)[which(o_8wkvs18mo$fdr<0.05)]))
    cat("\n 8wk vs 9 mo ",length(rownames(o_8wkvs9mo)[which(o_8wkvs9mo$fdr<0.05)]))
    cat("\n 9mo vs 18 mo ",length(rownames(o_9movs18mo)[which(o_9movs18mo$fdr<0.05)]))
    cat(" \n_____________________\n")
    
    
  }
}

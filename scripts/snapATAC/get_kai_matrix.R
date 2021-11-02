


library(SnapATAC)
tissue = "all66"
RData = "all.cluster.pmat.RData"
setwd(paste0("../../analysis/snapATAC/", tissue))

load(RData)


pmat = x.sp@pmat 


peak = x.sp@peak$name
cluster_stage = paste(x.sp@metaData$cluster_names,"_",x.sp@metaData$stage, sep = "")
olg_cpm = read.csv("cpm_tables/OLG_cpm.txt")
olg_peaks = olg_cpm$X


gt = function(x) {
  return(paste(which(x != 0)-1, x[which(x!=0)], sep = ",", collapse = "\t"))
}

pmat_OLG_9mo = pmat[which(cluster_stage=="OLG_9mo"),which(peak %in% olg_peaks)]
pmat_OLG_18mo = pmat[which(cluster_stage=="OLG_18mo"),which(peak %in% olg_peaks)]

bi = which(colSums(pmat_OLG_9mo)==0)
bi18 = which(colSums(pmat_OLG_18mo)==0)

rm_bi = unique(c(bi,bi18))
pmat_OLG_9mo = pmat_OLG_9mo[,-rm_bi]
pmat_OLG_18mo = pmat_OLG_18mo[,-rm_bi]

olg_peaks = olg_peaks[-rm_bi]
dim(pmat_OLG_9mo)
dim(pmat_OLG_18mo)

rm(x.sp)

all_res = c()
for (i in seq(1,ncol(pmat_OLG_9mo),500)){
  mini = min((i+499),ncol(pmat_OLG_9mo))
  cat(i,"-",mini," pmat_OLG_9mo \n")  
  a1 = apply(pmat_OLG_9mo[,i:mini],2,gt)
  all_res = c(all_res,a1)
  
}
cat(length(all_res), length(olg_peaks), "\n")

dat = data.frame(all_res,row.names=olg_peaks[1:length(all_res)])
#head(dat)
colnames(dat) = paste("sparse matrix:", dim(pmat_OLG_9mo)[2],"x",dim(pmat_OLG_9mo)[1])
write.table(dat, file = "kai_matrix_DE/pmat_OLG_9mo.txt", quote = F, col.names = T, sep = "\t")


all_res = c()
for (i in seq(1,ncol(pmat_OLG_18mo),500)){
  mini = min((i+499),ncol(pmat_OLG_18mo))
  cat(i,"-",mini," pmat_OLG_18mo \n")
  a1 = apply(pmat_OLG_18mo[,i:mini],2,gt)
  all_res = c(all_res,a1)
  
}

dat = data.frame(all_res,row.names=olg_peaks[1:length(all_res)])
colnames(dat) = paste("sparse matrix:", dim(pmat_OLG_18mo)[2],"x",dim(pmat_OLG_18mo)[1])
write.table(dat, file = "kai_matrix_DE/pmat_OLG_18mo.txt", quote = F, col.names = T, sep = "\t")







library(doParallel)
library(SnapATAC)
tissue = "all66"
RData = "all.cluster.pmat.RData"
setwd(paste0("../../analysis/snapATAC/", tissue))

load(RData)


pmat = x.sp@pmat 


peak = x.sp@peak$name
x.sp@metaData$cluster_names[which(x.sp@metaData$cluster_names %in% c(28,29))]="D2MSN"
cluster_stage_region = paste(x.sp@metaData$cluster_names,"_",x.sp@metaData$stage, "_", x.sp@metaData$tissue, sep = "")

clusters = x.sp@metaData$cluster_names
rm(x.sp)


gt = function(x) {
  return(paste(which(x != 0)-1, x[which(x!=0)], sep = ",", collapse = "\t"))
}

get_table = function(mat) {
  all_res = c()
  for (i in seq(1,ncol(mat),500)){
    mini = min((i+499),ncol(mat))
    a1 = apply(mat[,i:mini],2,gt)
    all_res = c(all_res,a1)
  }
  return(all_res)
}

#get peaks for each cluster
all_peaks = list()
files = list.files(path = "cpm_tables/",pattern = ".txt")
for(f in files) {
  cl = strsplit(f, "_")[[1]][1]
  cpm = read.csv(paste("cpm_tables/",f, sep = ""))
  all_peaks[[cl]] = cpm$X
  write.table(cpm$X, paste("kai_matrix_DE/cluster_all_peaks/",cl,"_peaks.txt",sep = ""),row.names = F,quote=F,col.names = F)
}

registerDoParallel(cores = 10)
cluster_stage_region_options = cluster_stage_region[-grep("8wk",cluster_stage_region)]

foreach(i = 1:length(unique(cluster_stage_region_options))) %dopar% {
#for(i in 59:length(unique(cluster_stage_region_options))) {
  c = unique(cluster_stage_region_options)[i]
  cluster = strsplit(c,"_")[[1]][1]
  peaks = all_peaks[[cluster]]
  if(length(peaks)>1) {
    cat(c,"\n")
    if (length(which(cluster_stage_region==c))>10) {
    mat =  pmat[which(cluster_stage_region==c),which(peak %in% peaks)]
    cat(dim(mat))
    mat = t(mat)
    res = get_table(mat)
    
    #dat = data.frame(res,row.names=peak[which(peak %in% peaks)])
    dat = data.frame(res,row.names = paste(names(res)))
    colnames(dat) = paste("sparse matrix:", dim(mat)[2],"x",dim(mat)[1])
    
    #show(head(dat, 20))
    write.table(dat, file = paste("kai_matrix_DE/pmat_",c,".txt",sep = ""), quote = F, col.names = T, sep = "\t")
    }
  }
}

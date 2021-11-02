
library(SnapATAC)
library(parallel)
setwd("~/projects/brain_aging_mouse/analysis/snapATAC/all66/")        
load("all21.cluster.RData")
plotViz(
  obj= x.sp,
  method="umap", 
  main="Cluster",
  point.color=x.sp@cluster, 
  point.size=0.2, 
  point.shape=19, 
  text.add=TRUE,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
);

x.sp@metaData$cluster_names = as.character(x.sp@cluster)

x.sp@metaData$cluster_names[which(x.sp@metaData$cluster_names %in% as.character(c(2,5,11)))] = "DG"
x.sp@metaData$cluster_names[which(x.sp@metaData$cluster_names %in% as.character(c(10,18,23)))] = "OLG"
x.sp@metaData$cluster_names[which(x.sp@metaData$cluster_names %in% as.character(c(4,20,30)))] = "AST"
x.sp@metaData$cluster_names[which(x.sp@metaData$cluster_names %in% as.character(c(1)))] = "OPC"
x.sp@metaData$cluster_names[which(x.sp@metaData$cluster_names %in% as.character(c(8)))] = "MCG"

plotViz(
  obj= x.sp,
  method="umap", 
  main="Cluster",
  point.color=x.sp@metaData$cluster_names, 
  point.size=0.2, 
  point.shape=19, 
  text.add=TRUE,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
);



clusters.sel = names(table(x.sp@metaData$cluster_names))[which(table(x.sp@metaData$cluster_names) > 200)];
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  runMACS(
    obj=x.sp[which(x.sp@metaData$cluster_names==clusters.sel[i]),], 
    output.prefix=paste0("atac_v1_adult_brain_fresh_5k.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/projects/ps-renlab/lamaral/software/miniconda3/envs/py27/bin/snaptools",
    path.to.macs="/projects/ps-renlab/lamaral/software/miniconda3/envs/py27/bin/macs2",
    gsize="mm", # mm, hs, etc
    buffer.size=500, 
    num.cores=1,
    macs.options="--shift 75 --extsize 150 --nomodel --call-summits --SPMR --keep-dup all -q 0.01",
    tmp.folder=tempdir()
  );
}, mc.cores=10)


peaks.names = system("ls | grep narrowPeak", intern=TRUE);
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = reduce(Reduce(c, peak.gr.ls));
peak.gr
peaks.df = as.data.frame(peak.gr)[,1:3];
write.table(peaks.df,file = "peaks.combined.bed",append=FALSE,
              quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
              fileEncoding = "")

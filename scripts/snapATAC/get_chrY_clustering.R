library(tictoc)
library(SnapATAC)
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("doParallel"))
library("umap")

setwd("~/projects/brain_aging_mouse/analysis/snapATAC/all66/")
RData = "all.cluster.meta.final.RData"
load(RData)
black_list <- "/projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed.gz"
cpus = 5
pc_num = 50
dims = 40
sample_list = read.csv("~/projects/brain_aging_mouse/scripts/snapATAC/doublet_thresholds_manual.csv")
sample.list <- sample_list$Sample

plotViz(
  obj= x.sp,
  method="umap",
  main="Type",
  point.color=x.sp@metaData$final_clusters,
  point.size=0.2,
  point.shape=19,
  text.add=T,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)
#peaks = x.sp@peak$name
#chr = sapply(strsplit(as.character(peaks), ":"), "[[", 1)
#x.sp <- x.sp[,which(chr == 'chrY'), mat = "pmat"]

dims = 15
pc_num = 35
x.sp <- x.sp[which(x.sp@metaData$final_clusters%in%c("CTGL", "ITL23GL","ITHGL", "ITL6GL")), ]
x.sp <- x.sp[-which(x.sp@sample %in% names(which(table(x.sp@sample)<5))),]

bin_size = 5000
tic("addBmatToSnap")
x.sp = addBmatToSnap(x.sp, bin.size=bin_size);
toc()

tic("make Binary")
x.sp = makeBinary(x.sp, mat="bmat");
toc()
#x.sp.sub <- x.sp.sub[-which(x.sp.sub@sample %in% names(which(table(x.sp.sub@sample)<40))),]
outF = paste("chrY_exN_Y_", sep = "")

x.sp

tic("filterBins")
black_list = read.table(black_list);
black_list.gr = GRanges(
  black_list[,1], 
  IRanges(black_list[,2], black_list[,3])
);
idy = queryHits(
  findOverlaps(x.sp@feature, black_list.gr)
);
if(length(idy) > 0){
  x.sp = x.sp[,-idy, mat="bmat"];
};
x.sp

chr.include = seqlevels(x.sp@feature)[grep("chrY", seqlevels(x.sp@feature))];
idy = grep(paste(chr.include, collapse="|"), x.sp@feature);
if(length(idy) > 0){
  x.sp = x.sp[,idy, mat="bmat"]
};
x.sp

bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
#pdf(paste(outF, ".BinCoverage.pdf",sep=""))
hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 5)
)
#dev.off()

bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
x.sp = x.sp[, idy, mat="bmat"]
x.sp
bin.cov = log10(Matrix::rowSums(x.sp@bmat)+1);
idy = which(bin.cov > 0)
x.sp = x.sp[idy, ,mat="bmat"]
x.sp



x.sp = runDiffusionMaps(
  obj=x.sp,
  input.mat="bmat", 
  num.eigs=50
);

save(x.sp, file = paste(outF, "x.sp.RData", sep = ""))
load(paste(outF, "x.sp.RData", sep = ""))

tic("runKNN")
x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:dims,
  k=15
);
toc()


## Clustering


## R-igraph
tic("runCluster_R-igraph")
x.sp = runCluster(obj=x.sp, 
                      tmp.folder=tempdir(), 
                      louvain.lib="R-igraph",
                      seed.use=10
);
toc()

## Visulization

tic("runViz_umap")
x.sp = runViz(
  obj=x.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:dims, 
  method="umap",
  seed.use=10
);
toc()

outfname = paste(outF, ".cluster.RData",sep="")
save(x.sp, file=outfname)

pdf(paste(outF, ".cluster.pdf",sep=""))
plotViz(
  obj= x.sp,
  method="umap",
  main="Cluster",
  point.color=x.sp@metaData$final_clusters,
  point.size=0.2,
  point.shape=19,
  text.add=T,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)

plotViz(
  obj= x.sp,
  method="umap", 
  main="Cluster",
  point.color=x.sp@metaData$stage, 
  point.size=0.2, 
  point.shape=19, 
  text.add=F,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=T
)
plotViz(
  obj= x.sp,
  method="umap", 
  main="Cluster orig",
  point.color=x.sp@metaData$cluster_names, 
  point.size=0.2, 
  point.shape=19, 
  text.add=TRUE,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)
plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData[,"log10UQ"],
  method="umap",
  main="Read Depth",
  point.size=0.2,
  point.shape=19,
  down.sample=10000,
  quantiles=c(0.01, 0.99)
);

plotViz(
  obj= x.sp,
  method="umap", 
  main="Landmark",
  point.size=0.2, 
  point.shape=19, 
  point.color=x.sp@metaData[,"landmark"], 
  text.add=FALSE,
  text.size=1.5,
  text.color="black",
  down.sample=10000,
  legend.add=TRUE
);

plotViz(
  obj= x.sp,
  method="umap",
  main="stage",
  point.size=0.2,
  point.shape=19,
  point.color=x.sp@metaData[,"stage"],
  text.add=FALSE,
  text.size=1.5,
  text.color="black",
  down.sample=10000,
  legend.add=TRUE
);

plotViz(
  obj= x.sp,
  method="umap",
  main="Replicates",
  point.size=0.2,
  point.shape=19,
  point.color=x.sp@metaData[,"replicate"],
  text.add=FALSE,
  text.size=1.5,
  text.color="black",
  down.sample=10000,
  legend.add=TRUE
);

plotViz(
  obj= x.sp,
  method="umap",
  main="region",
  point.size=0.2,
  point.shape=19,
  point.color=x.sp@metaData[,"region"],
  text.add=F,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=F
);

plotViz(
  obj= x.sp,
  method="umap",
  main="region",
  point.size=0,
  point.shape=19,
  point.color=x.sp@metaData[,"region"],
  text.add=FALSE,
  text.size=1.5,
  text.color="black",
  down.sample=10000,
  legend.add=TRUE
);
dev.off()

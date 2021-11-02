
load("~/projects/Aging/all_tissues/all_merged_40kL.cluster.RData")

plotViz(
  obj= x.sp,
  method="umap",
  main="Cluster",
  point.color=x.sp@cluster,
  point.size=0.2,
  point.shape=19,
  text.add=F,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)
#Endo:
# DH.13
# FC.15
# LM.6
# HT.7
# HT.3
# BM.15

library(tictoc)
library(SnapATAC)
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("doParallel"))
library("umap")

setwd("~/projects/brain_aging_mouse/scripts/snapATAC/")
outF = "~/projects/Aging/all_tissues/endo/Endo_"
black_list <- "/projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed.gz"
cpus = 5
pc_num = 50
dims = 15
x.sp@metaData$cluster = x.sp@cluster
#sample_list = read.csv("doublet_thresholds_manual.csv")
#sample.list <- sample_list$Sample


x.sp@metaData$endo = ""
x.sp@metaData$endo[which(x.sp@metaData$tissue_cluster %in% c("DH_13","FC_15","LM_6","HT_7","HT_3","BM_15"))] = "Endo"

plotViz(
  obj= x.sp,
  method="umap",
  main="Type",
  point.color=x.sp@metaData$endo,
  point.size=0.2,
  point.shape=19,
  text.add=T,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)


get_subclustering <- function(type) {
  
  dims = 20
  pc_num = 35
  x.sp.sub <- x.sp[which(x.sp@metaData$endo == "Endo"), ]
  #x.sp.sub <- x.sp.sub[-which(x.sp.sub@sample %in% names(which(table(x.sp.sub@sample)<40))),]
  outF = paste("~/projects/Aging/all_tissues/endo/Endo_", sep = "")
  set.seed(10)
  

  ## 1. identify usable features
  bin_size = 5000
  tic("addBmatToSnap")
  x.sp.sub = addBmatToSnap(x.sp.sub, bin.size=bin_size,do.par = T,num.cores = 5);
  toc()
  
  tic("makeBinary")
  x.sp.sub = makeBinary(x.sp.sub, mat="bmat");
  toc()
  
  tic("filterBins")
  black_list = read.table(black_list);
  black_list.gr = GRanges(
    black_list[,1], 
    IRanges(black_list[,2], black_list[,3])
  );
  idy = queryHits(
    findOverlaps(x.sp.sub@feature, black_list.gr)
  );
  if(length(idy) > 0){
    x.sp.sub = x.sp.sub[,-idy, mat="bmat"];
  };
  x.sp.sub
  
  chr.exclude = seqlevels(x.sp.sub@feature)[grep("random|chrM", seqlevels(x.sp.sub@feature))];
  idy = grep(paste(chr.exclude, collapse="|"), x.sp.sub@feature);
  if(length(idy) > 0){
    x.sp.sub = x.sp.sub[,-idy, mat="bmat"]
  };
  x.sp.sub
  
  bin.cov = log10(Matrix::colSums(x.sp.sub@bmat)+1);
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
  x.sp.sub = x.sp.sub[, idy, mat="bmat"]
  x.sp.sub
  
  #idx = which(Matrix::rowSums(x.landmark.sp@bmat) > 500);
  #x.landmark.sp = x.landmark.sp[idx,];
  #x.landmark.sp
  toc()
  
  #summarySnap(x.landmark.sp)
  x.sp.sub

  # 2. embedding
  tic("runDiffusionMaps")
  x.sp.sub = runDiffusionMaps(
    obj= x.sp.sub,
    input.mat="bmat",
    num.eigs=50
  );
  toc()
  
  
 
  

  pdf(paste(outF, ".pre.plotDimResuce.pdf",sep=""))
  plotDimReductElbow(
    obj=x.sp.sub,
    point.size=1.5,
    point.shape=19,
    point.color="red",
    point.alpha=1,
    pdf.file.name=NULL,
    pdf.height=7,
    pdf.width=7,
    labs.title="PCA Elbow plot",
    labs.subtitle=NULL
  )
  
  plotDimReductPW(
    obj=x.sp.sub,
    eigs.dims=1:pc_num
  )
  dev.off()
  
  
  
  
  
  # KNN
  tic("runKNN")
  x.sp.sub = runKNN(
    obj=x.sp.sub,
    eigs.dims=1:dims,
    k=15
  );
  toc()
  
  
  ## Clustering
  
  
  ## R-igraph
  tic("runCluster_R-igraph")
  x.sp.sub = runCluster(obj=x.sp.sub, 
                        tmp.folder=tempdir(), 
                        louvain.lib="R-igraph",
                        seed.use=10
  );
  toc()
  
  ## Visulization
  
  tic("runViz_umap")
  x.sp.sub = runViz(
    obj=x.sp.sub, 
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:dims, 
    method="umap",
    seed.use=10
  );
  toc()
  
  outfname = paste(outF, ".cluster.RData",sep="")
  save(x.sp.sub, file=outfname)
  
  pdf(paste(outF, ".cluster.pdf",sep=""))
  plotViz(
    obj= x.sp.sub,
    method="umap",
    main="Cluster",
    point.color=x.sp.sub@cluster,
    point.size=0.2,
    point.shape=19,
    text.add=T,
    text.size=1,
    text.color="black",
    down.sample=10000,
    legend.add=FALSE
  )
  
  plotViz(
    obj= x.sp.sub,
    method="umap", 
    main="tissue combined Cluster",
    point.color=x.sp.sub@metaData$cluster, 
    point.size=0.2, 
    point.shape=19, 
    text.add=TRUE,
    text.size=1,
    text.color="black",
    down.sample=10000,
    legend.add=FALSE
  )
  plotViz(
    obj= x.sp.sub,
    method="umap", 
    main="Tissue Cluster",
    point.color=x.sp.sub@metaData$tissue_cluster, 
    point.size=0.2, 
    point.shape=19, 
    text.add=TRUE,
    text.size=1,
    text.color="black",
    down.sample=10000,
    legend.add=FALSE
  )
  plotViz(
    obj= x.sp.sub,
    method="umap", 
    main="Tissue Cluster",
    point.color=x.sp.sub@metaData$tissue_cluster, 
    point.size=0.2, 
    point.shape=19, 
    text.add=F,
    text.size=1,
    text.color="black",
    down.sample=10000,
    legend.add=T
  )
  # plotFeatureSingle(
  #   obj=x.sp.sub,
  #   feature.value=x.sp.sub@metaData[,"log10UQ"],
  #   method="umap", 
  #   main="Read Depth",
  #   point.size=0.2, 
  #   point.shape=19, 
  #   down.sample=10000,
  #   quantiles=c(0.01, 0.99)
  # );
  
  
  
  plotViz(
    obj= x.sp.sub,
    method="umap",
    main="stage",
    point.size=0.2,
    point.shape=19,
    point.color=x.sp.sub@metaData[,"stage"],
    text.add=FALSE,
    text.size=1.5,
    text.color="black",
    down.sample=10000,
    legend.add=TRUE
  );
  
  plotViz(
    obj= x.sp.sub,
    method="umap",
    main="Replicates",
    point.size=0.2,
    point.shape=19,
    point.color=x.sp.sub@metaData[,"replicate"],
    text.add=FALSE,
    text.size=1.5,
    text.color="black",
    down.sample=10000,
    legend.add=TRUE
  );
 
  dev.off()
  
  
}
#get_subclustering("InhN")
#get_subclustering("Glial")
#get_subclustering("ExN")
#get_subclustering("D2MSN")
#get_subclustering("AST")
get_subclustering("OLG")
get_subclustering("MCG")
get_subclustering("END")





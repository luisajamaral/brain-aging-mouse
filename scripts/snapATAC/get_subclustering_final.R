library(tictoc)
library(SnapATAC)
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("doParallel"))
library("umap")

setwd("~/projects/brain_aging_mouse/scripts/snapATAC/")
outF = "../../analysis/snapATAC/all66/"
load(paste(outF, "all.cluster.meta.final.RData",sep=""))
black_list <- "/projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed.gz"
cpus = 5
pc_num = 50
dims = 15
x.sp@metaData$cluster = x.sp@metaData$cluster_names
sample_list = read.csv("doublet_thresholds_manual.csv")
sample.list <- sample_list$Sample


x.sp@metaData$macro_type = ""
x.sp@metaData$macro_type[which(x.sp@metaData$cluster_names %in% c("L6bGL","DG","CA1GL","ITHGL","ITL5GL","U4","CA3GL"))] = "ExN"
x.sp@metaData$macro_type[which(x.sp@metaData$cluster_names %in% c("D2MSN","U1","U2","PVGA","U3"))] = "InhN"
#x.sp@metaData$macro_type[which(x.sp@metaData$cluster_names %in% c("OLG","PER","VEC","VLMC","VPIA","OPC","MCG","AST"))] = "Glial"
x.sp@metaData$macro_type[which(x.sp@metaData$cluster_names %in% c("D2MSN"))] = "D2MSN"
x.sp@metaData$macro_type[which(x.sp@metaData$cluster_names %in% c("AST"))] = "AST"
x.sp@metaData$macro_type[which(x.sp@metaData$cluster_names %in% c("OLG","OPC"))] = "OLG"
x.sp@metaData$macro_type[which(x.sp@metaData$cluster_names %in% c("MCG"))] = "MCG"
x.sp@metaData$macro_type[which(x.sp@metaData$cluster_names %in% c("PER","VEC","VLMC","VPIA"))] = "VASC"

plotViz(
  obj= x.sp,
  method="umap",
  main="Type",
  point.color=x.sp@metaData$macro_type,
  point.size=0.2,
  point.shape=19,
  text.add=T,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)


get_subclustering <- function(type) {
  
  dims = 15
  pc_num = 35
  x.sp.sub <- x.sp[which(x.sp@metaData$macro_type == type), ]
  #x.sp.sub <- x.sp.sub[-which(x.sp.sub@sample %in% names(which(table(x.sp.sub@sample)<40))),]
  outF = paste("../../analysis/snapATAC/all66/sub_",type, sep = "")
  set.seed(10)
  
  ## 0. sampling
  sampleSize <- 10000
  
  tic("sampling landmark")
  row.covs.dens <- density(
    x = x.sp.sub@metaData[,"log10UQ"],
    bw = 'nrd', adjust = 1
  );
  
  sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = x.sp.sub@metaData[,"log10UQ"])$y + .Machine$double.eps);
  set.seed(2022)
  idx.landmark.ds <- sort(sample(x = seq(nrow(x.sp.sub)), size = sampleSize, prob = sampling_prob));
  x.landmark.sp = x.sp.sub[idx.landmark.ds,];
 # #x.landmark.sp <- x.landmark.sp[-which(x.landmark.sp@sample %in% names(which(table(x.landmark.sp@sample)<5))),]
  
  x.query.sp = x.sp.sub[-idx.landmark.ds,];
  toc()
  
  
  ## 1. identify usable features
  bin_size = 5000
  tic("addBmatToSnap")
  x.landmark.sp = addBmatToSnap(x.landmark.sp, bin.size=bin_size,do.par = T,num.cores = 5);
  toc()
  
  tic("makeBinary")
  x.landmark.sp = makeBinary(x.landmark.sp, mat="bmat");
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
    x.landmark.sp = x.landmark.sp[,-idy, mat="bmat"];
  };
  x.landmark.sp
  
  chr.exclude = seqlevels(x.landmark.sp@feature)[grep("random|chrM", seqlevels(x.landmark.sp@feature))];
  idy = grep(paste(chr.exclude, collapse="|"), x.landmark.sp@feature);
  if(length(idy) > 0){
    x.landmark.sp = x.landmark.sp[,-idy, mat="bmat"]
  };
  x.landmark.sp
  
  bin.cov = log10(Matrix::colSums(x.landmark.sp@bmat)+1);
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
  x.landmark.sp = x.landmark.sp[, idy, mat="bmat"]
  x.landmark.sp
  
  #idx = which(Matrix::rowSums(x.landmark.sp@bmat) > 500);
  #x.landmark.sp = x.landmark.sp[idx,];
  #x.landmark.sp
  toc()
  
  #summarySnap(x.landmark.sp)
  x.landmark.sp
  outfname = paste(outF, ".landmark.RData",sep="")
  #save(x.landmark.sp, file=outfname)
  
  #summarySnap(x.query.sp)
  x.query.sp
  outfname = paste(outF, ".query.RData",sep="")
  #save(x.query.sp, file=outfname)
  
  
  # 2. embedding
  tic("runDiffusionMaps")
  x.landmark.sp = runDiffusionMaps(
    obj= x.landmark.sp,
    input.mat="bmat",
    num.eigs=50
  );
  x.landmark.sp@metaData$landmark = 1;
  toc()
  
  
  # 3.extension
  tic("extension")
  num.files = length(sample.list);
  
  registerDoParallel(cpus)
  sample.list = unique(x.landmark.sp@sample)
  num.files = length(sample.list);
  
  #x.query.ls <- foreach (i=1:num.files, .combine=snapRbind, .inorder=TRUE) %dopar% {
  x.query.ls <- foreach (i=1:num.files, .inorder=TRUE) %dopar% {
    #for (i in 1:num.files) {
    print(sample.list[i])
    x.query.sub <- x.query.sp[which(x.query.sp@sample == as.character(sample.list[i])), ]
    x.query.sub = addBmatToSnap(x.query.sub, do.par=T, num.cores=cpus, bin.size=bin_size);
    x.query.sub = makeBinary(x.query.sub);
    
    idy = unique(queryHits(findOverlaps(x.query.sub@feature, x.landmark.sp@feature)));
    x.query.sub = x.query.sub[,idy, mat="bmat"];
    x.query.sub = runDiffusionMapsExtension(
      obj1=x.landmark.sp,
      obj2=x.query.sub,
      input.mat="bmat"
    );
    x.query.sub@metaData$landmark = 0;
    x.query.sub = rmBmatFromSnap(x.query.sub)
    return(x.query.sub)
  }
  closeAllConnections()
  
  x.query.sp = Reduce(snapRbind, x.query.ls);
  x.landmark.sp = rmBmatFromSnap(x.landmark.sp)
  x.sp.sub = snapRbind(x.landmark.sp, x.query.sp);
  x.sp.sub = x.sp.sub[order(x.sp.sub@metaData[,"sample"])];
  toc()
  
  rm(x.landmark.sp)
  rm(x.query.sp)
  
  
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
    main="Cluster",
    point.color=x.sp.sub@metaData$final_clusters, 
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
    main="Cluster orig",
    point.color=x.sp.sub@metaData$cluster_names, 
    point.size=0.2, 
    point.shape=19, 
    text.add=TRUE,
    text.size=1,
    text.color="black",
    down.sample=10000,
    legend.add=FALSE
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
    main="Landmark",
    point.size=0.2, 
    point.shape=19, 
    point.color=x.sp.sub@metaData[,"landmark"], 
    text.add=FALSE,
    text.size=1.5,
    text.color="black",
    down.sample=10000,
    legend.add=TRUE
  );
  
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
  
  plotViz(
    obj= x.sp.sub,
    method="umap",
    main="region",
    point.size=0.2,
    point.shape=19,
    point.color=x.sp.sub@metaData[,"region"],
    text.add=F,
    text.size=1,
    text.color="black",
    down.sample=10000,
    legend.add=F
  );
  
  plotViz(
    obj= x.sp.sub,
    method="umap",
    main="region",
    point.size=0,
    point.shape=19,
    point.color=x.sp.sub@metaData[,"region"],
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



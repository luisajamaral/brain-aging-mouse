library(tictoc)
library(SnapATAC)
setwd("~/projects/brain_aging_mouse/scripts/snapATAC/")
outF = "../../analysis/snapATAC/all66/all21"
load(paste(outF, ".cluster.RData",sep=""))
x.sp@metaData$cluster = x.sp@cluster

x.sp@metaData$macro_type = ""
x.sp@metaData$macro_type[which(x.sp@metaData$cluster_names %in% c("L6bGL","DG","CA1GL","ITHGL","ITL5GL","U4","CA3GL"))] = "ExN"
x.sp@metaData$macro_type[which(x.sp@metaData$cluster_names %in% c("D2MSN","U1","U2","PVGA","U3"))] = "InhN"
x.sp@metaData$macro_type[which(x.sp@metaData$cluster_names %in% c("OLG","PER","VEC","VLMC","VPIA","OPC","MCG","AST"))] = "Glial"

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




#x.sp.sub <- x.sp[which(x.sp@cluster %in% c(2,5,11,13,14,15,17,21,22,27,31)), ]
#x.sp.sub <- x.sp[which(x.sp@cluster %in% c(3,9,12,19,25,28,29)), ]


get_subclustering <- function(type) {
  

x.sp.sub <- x.sp[which(x.sp@metaData$macro_type == type), ]
outF = paste("../../analysis/snapATAC/all66/sub_",type)
#outF = "../../analysis/snapATAC/all66/Inh"
set.seed(2021)
x.sp.sub <- x.sp.sub[sample(1:nrow(x.sp.sub),10000), ]
x.sp.sub = addBmatToSnap(x.sp.sub, bin.size=5000,do.par = T,num.cores = 5)
x.sp.sub = makeBinary(x.sp.sub, mat="bmat")
x.sp.sub = runDiffusionMaps(
  obj= x.sp.sub,
  input.mat="bmat",
  num.eigs=50
);
pc_num = 50
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




dims = 12
# KNN
tic("runKNN")
x.sp.sub = runKNN(
  obj=x.sp.sub,
  eigs.dims=1:dims,
  k=15
);
toc()


## Clustering

#library("leiden")
#tic("runCluster_leiden")
#x.sp.sub = runCluster(
#    obj=x.sp.sub,
#    tmp.folder=tempdir(),
#    louvain.lib="leiden",
#    resolution = resolution,
#    seed.use=10
#    );
#toc()

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


pdf(paste(outF, ".cluster.pdf",sep=""))
plotViz(
  obj= x.sp.sub,
  method="umap",
  main="Cluster",
  point.color=x.sp.sub@cluster,
  point.size=0.2,
  point.shape=19,
  text.add=FALSE,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)

plotViz(
  obj= x.sp.sub,
  method="umap", 
  main="Cluster",
  point.color=x.sp.sub@cluster, 
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
  main="Cluster",
  point.color=x.sp.sub@metaData$cluster_names, 
  point.size=0.2, 
  point.shape=19, 
  text.add=TRUE,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
);


plotFeatureSingle(
  obj=x.sp.sub,
  feature.value=x.sp.sub@metaData[,"log10UQ"],
  method="umap", 
  main="Read Depth",
  point.size=0.2, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99)
);

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
dev.off()
outfname = paste(outF, ".cluster.RData",sep="")
save(x.sp.sub, file=outfname)

}
get_subclustering("InhN")
get_subclustering("Glial")
get_subclustering("ExN")

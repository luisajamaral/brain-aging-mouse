
suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("doParallel"))
library("umap")
library("leiden")
library("tictoc")




#inputF <- "heroin.DNA.snap.list"
bin_size <- 5000
#datf <- "heroin.dataset.list"
#sumf <- "heroin.barcode.sum.txt"
black_list <- "/projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed.gz"
pc_num <- 50
cpus <- 5
#outF <- "heroin.DNA"
seed.use <- 2020
dims <- 20
geneF <- "gene.markers.txt"
resolution <- 0.5
outF = "../../analysis/snapATAC/all21"
#sample.list = c("8E-9H-8J-9J_09_rep1", "8E-9H-8J-9J_09_rep2", "8E-9H-8J-9J_18_rep1", "8E-9H-8J-9J_18_rep2", "8E_03_rep1", "8E_03_rep2", "8J_03_rep1", "8J_03_rep2", "9H_03_rep1", "9H_03_rep2", "9J_03_rep1",  "9J_03_rep2")
sumF = data.frame()




#inputf <- read.table(inputF,sep="\t",header=F)
#file.list <- as.character(inputf[,2])
#sample.list <- as.character(inputf[,1])

sample_list = read.csv("doublet_thresholds_manual.csv")
rownames(sample_list) = sample_list$Sample
sample_list$Sample = as.character(sample_list$Sample)
sample.list <- sample_list$Sample
file.list = c()
for (i in sample.list) {
  file.list = c(file.list, paste("../../analysis/by_sample.snapATAC/",i,"/snapFiles/",i,".snap", sep = ""))
}


tic("createSnap")
x.sp.ls = lapply(seq(file.list), function(i){
  x.sp = createSnap(file=file.list[i], sample=sample.list[i], do.par = TRUE, num.cores=cpus);
  x.sp
})
names(x.sp.ls) = sample.list;
x.sp = Reduce(snapRbind, x.sp.ls);
x.sp@metaData["sample"] = x.sp@sample;
toc()
rm(x.sp.ls)

tmp <- do.call(rbind, strsplit(x.sp@metaData$sample, "_"))
x.sp@metaData$replicate = sample_list[paste(x.sp@metaData$sample),"rep"]
x.sp@metaData$stage <- sample_list[paste(x.sp@metaData$sample),"stage"]

for (i in sample.list) {
  cur =  read.table(paste("../../analysis/by_sample.snapATAC/",i,"/snapFiles/",i,".fitDoublets.txt", sep = "") , header = T)
  newthresh = sample_list[paste(i),"manual.threshold"]
  cur$manual.threshold = newthresh
  cur$new.pred.doub = FALSE
  cur$new.pred.doub[which(cur$doublet_scores>newthresh)]=TRUE
  sumF = rbind(sumF, cur)
}
x.sp@metaData$quality <- "lowQ"

idx <- which(paste(x.sp@sample, x.sp@barcode, sep=".") %in% paste(sumF$sample, sumF$barcode, sep="."))
x.sp@metaData[idx, "quality"] <- "passQC"

selF <- subset(sumF, sumF$new.pred.doub=="TRUE")
idx <- which(paste(x.sp@sample, x.sp@barcode, sep=".") %in% paste(selF$sample, selF$barcode, sep="."))
x.sp@metaData[idx, "quality"] <- "doublet"

x.sp@metaData <- join(x.sp@metaData, sumF[, c("sample","barcode","TSS_enrich","doublet_scores")], by=c("sample","barcode"))
x.sp@metaData$log10UQ <- log10(x.sp@metaData$UQ)

outfname = paste(outF, ".raw.RData",sep="")
save(x.sp, file=outfname)

outmetaf <- paste(outF, ".raw.meta.txt", sep="")
outmetamx <- x.sp@metaData
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)

pdf(paste(outF, ".raw.sta.pdf",sep=""))
plotBarcode(x.sp,
            pdf.file.name=NULL,
            pdf.width=7,
            pdf.height=7,
            col="grey",
            border="grey",
            breaks=50
)
dev.off()


# filter used cells
idx <- which(x.sp@metaData$quality == "passQC")
x.sp <- x.sp[idx, ]

outfname = paste(outF, ".sel.RData",sep="")
save(x.sp, file=outfname)

outmetaf <- paste(outF, ".sel.meta.txt", sep="")
outmetamx <- x.sp@metaData
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)

pdf(paste(outF, ".sel.sta.pdf",sep=""))
plotBarcode(x.sp,
            pdf.file.name=NULL,
            pdf.width=7,
            pdf.height=7,
            col="grey",
            border="grey",
            breaks=50
)
dev.off()





x.sp <- get(load(paste(outF, ".sel.RData",sep="")))

## 0. sampling
sampleSize <- 10000

tic("sampling landmark")
row.covs.dens <- density(
  x = x.sp@metaData[,"log10UQ"],
  bw = 'nrd', adjust = 1
);

sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = x.sp@metaData[,"log10UQ"])$y + .Machine$double.eps);
set.seed(2020);
idx.landmark.ds <- sort(sample(x = seq(nrow(x.sp)), size = sampleSize, prob = sampling_prob));
x.landmark.sp = x.sp[idx.landmark.ds,];
x.query.sp = x.sp[-idx.landmark.ds,];
toc()


## 1. identify usable features
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
  findOverlaps(x.sp@feature, black_list.gr)
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
pdf(paste(outF, ".BinCoverage.pdf",sep=""))
hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 5)
)
dev.off()

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
save(x.landmark.sp, file=outfname)

#summarySnap(x.query.sp)
x.query.sp
outfname = paste(outF, ".query.RData",sep="")
save(x.query.sp, file=outfname)


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
#x.query.ls <- foreach (i=1:num.files, .combine=snapRbind, .inorder=TRUE) %dopar% {
x.query.ls <- foreach (i=1:num.files, .inorder=TRUE) %dopar% {
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
x.sp = snapRbind(x.landmark.sp, x.query.sp);
x.sp = x.sp[order(x.sp@metaData[,"sample"])];
toc()

rm(x.landmark.sp)
rm(x.query.sp)

outfname = paste(outF, ".pre.RData",sep="")
save(x.sp, file=outfname)

outmeta = paste(outF, ".pre.meta.txt",sep="")
write.table(x.sp@metaData, file=outmeta, sep="\t", quote=F, col.names=T, row.names=F)


pdf(paste(outF, ".pre.plotDimResuce.pdf",sep=""))
plotDimReductElbow(
  obj=x.sp,
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
  obj=x.sp,
  eigs.dims=1:pc_num
)
dev.off()





# KNN
tic("runKNN")
x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:dims,
  k=15
);
toc()


## Clustering

#library("leiden")
#tic("runCluster_leiden")
#x.sp = runCluster(
#    obj=x.sp,
#    tmp.folder=tempdir(),
#    louvain.lib="leiden",
#    resolution = resolution,
#    seed.use=10
#    );
#toc()

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


pdf(paste(outF, ".cluster.pdf",sep=""))
plotViz(
  obj= x.sp,
  method="umap",
  main="Cluster",
  point.color=x.sp@cluster,
  point.size=0.2,
  point.shape=19,
  text.add=FALSE,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)

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
dev.off()

## Heretical clustering of the clusters
# calculate the ensemble signals for each cluster
#ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
#    Matrix::colMeans(x.sp@bmat[x,])
#    })

# cluster using 1-cor as distance
#hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");
#par(mfrow=c(1,1))
#pdf(paste(outF, ".hc.pdf",sep=""))
#plot(hc, hang=-1, xlab="");
#dev.off()

## Create chromatin lanscape and identify cis-elements for each cluster seperately.

outmetaf <- paste(outF, ".cluster.meta.txt", sep="")
outmetamx <- cbind(x.sp@sample, x.sp@metaData, x.sp@cluster, x.sp@umap)
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)

outfname = paste(outF, ".cluster.RData",sep="")
save(x.sp, file=outfname)

graph.knn <- x.sp@graph@mat

outfname = paste(outF, ".knn.mmtx", sep="")
writeMM(graph.knn,file=outfname)



# plot markers

set.seed(2020)
x.sp.sub <- x.sp[sample(1:nrow(x.sp),10000), ]

x.sp.sub <- addGmatToSnap(x.sp.sub, do.par = T, num.cores = 4)


# normalize the cell-by-gene matrix
x.sp.sub = scaleCountMatrix(
  obj=x.sp.sub, 
  cov=x.sp.sub@metaData$UQ + 1,
  mat="gmat",
  method = "RPM"
);

# smooth the cell-by-gene matrix
x.sp.sub = runMagic(
  obj=x.sp.sub,
  input.mat="gmat",
  step.size=3
);


geneF <- read.table(geneF, sep="\t", header=F)
marker.genes <- geneF[, 1]
marker.genes <- marker.genes[marker.genes %in% colnames(x.sp.sub@gmat)]

pdf(paste(outF, ".umap2marker.pdf",sep=""))
par(mfrow = c(2, 2));
for(i in 1:length(marker.genes)){
  plotFeatureSingle(
    obj=x.sp.sub,
    feature.value=x.sp.sub@gmat[, marker.genes[i]],
    method="umap", 
    main=marker.genes[i],
    point.size=0.1, 
    point.shape=19, 
    down.sample=10000,
    quantiles=c(0.01, 0.98)
  )}
dev.off()




# get sub cluster

RDataF <- "./heroin.DNA.refineCluster.RData"
x.sp <- get(load(RDataF))

typeF <- read.table("heroin.DNA.subCluster.meta.txt", sep="\t", header=T)
annoF <- read.table("heroin.DNA.subCluster.anno.txt", sep="\t", header=T)

typeF$MajorType <- NULL
colnames(typeF) <- c("sample","barcode", "cluster")
typeF <- join(typeF, annoF, by="cluster")
typeF$cluster <- NULL

x.sp@metaData$snapATAC_clusters <- x.sp@cluster

x.sp@metaData <- join(x.sp@metaData, typeF, by=c("sample","barcode"))

outfname = paste(outF, ".cluster2anno.RData",sep="")
save(x.sp, file=outfname)

outmetaf <- paste(outF, ".cluster2anno.meta.txt", sep="")
outmetamx <- cbind(x.sp@metaData, x.sp@umap)
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)

load(paste(outF, ".cluster2anno.RData",sep=""))
pdf(paste(outF, ".cluster2anno.attr.pdf",sep=""))
plotViz(
  obj= x.sp,
  method="umap",
  main="Cluster",
  point.color=x.sp@metaData$MajorType,
  point.size=0.2,
  point.shape=19,
  text.add=FALSE,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
);

plotViz(
  obj= x.sp,
  method="umap", 
  main="Cluster",
  point.color=x.sp@metaData$MajorType, 
  point.size=0.2, 
  point.shape=19, 
  text.add=TRUE,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
);

plotViz(
  obj= x.sp,
  method="umap",
  main="Cluster",
  point.color=x.sp@metaData$SubType,
  point.size=0.2,
  point.shape=19,
  text.add=FALSE,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
);

plotViz(
  obj= x.sp,
  method="umap", 
  main="Cluster",
  point.color=x.sp@metaData$SubType, 
  point.size=0.2, 
  point.shape=19, 
  text.add=TRUE,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
);


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
  main="Treatment",
  point.size=0.2,
  point.shape=19,
  point.color=x.sp@metaData[,"treat"],
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
dev.off()




setwd("~/projects/brain_aging_mouse/scripts/snapATAC/")
library(SnapATAC)
library("GenomicRanges")
geneF <- "gene.markers.txt"
geneF <- "CEMBA_annotation.csv"

outF = "../../analysis/snapATAC/all66/all21"
load(paste(outF, ".cluster.RData",sep=""))
set.seed(2020)
x.sp.sub <- x.sp[sample(1:nrow(x.sp),10000), ]

x.sp.sub <- addGmatToSnap(x.sp.sub, do.par = T, num.cores = 5)


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


geneF <- read.table(geneF, sep=",", header=F)
marker.genes <- geneF[, 3]
marker.genes <- marker.genes[marker.genes %in% colnames(x.sp.sub@gmat)]
marker.genes <- c("Prox1","Gpr161","Rnf182","Garnl3","Hs3st4","Tshz2","Olfm3","Rorb","Satb2","Cux2","Foxp2","Erbb4")

marker.genes <- c("Tshz1","Ndnf","Vip","Rgs12","Npas1","Lamp5","Npy","Pnoc","Sst","Elfn1","Pvalb","Erbb4","Drd2","Foxp2","Dgkg","Tshz2","Fam19a2","Prox1")

pdf(paste(outF, ".umap2marker_cemba.pdf",sep=""))
par(mfrow = c(2,2));
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
  )
  }
dev.off()




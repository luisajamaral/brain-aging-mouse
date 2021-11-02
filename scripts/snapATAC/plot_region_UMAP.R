library(ggplot2)
setwd("../../analysis/snapATAC/all66")
sample_list = read.csv("../../../scripts/snapATAC/doublet_thresholds_manual.csv")
rownames(sample_list) = sample_list$Sample
outF = "all21"
load(paste(outF, ".cluster.RData",sep=""))


x.sp@metaData$tissue = ""
x.sp@metaData$sample_name = ""
for(n in sample_list$Sample) {
  n = paste(n)
  x.sp@metaData[which(paste(x.sp@metaData$sample) == n),"tissue"] = paste(sample_list[n,"tissue"])
  samp_name = paste(sample_list[n,"stage"], sample_list[n,"tissue"], sample_list[n,"rep"], sep="_")
  x.sp@metaData[which(paste(x.sp@metaData$sample) == n),"sample_name"] = samp_name
}

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
pdf(paste(outF, ".tissuelegend.pdf",sep=""), height = 7, width = 7)

plotViz(
  obj= x.sp,
  method="umap",
  main="tissue",
  point.color=x.sp@metaData$tissue,
  point.size=0,
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
  main="sample_name",
  point.color=x.sp@metaData$sample_name,
  point.size=0,
  point.shape=19,
  text.add=F,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=T
)

dev.off()

pdf(paste(outF, ".region.pdf",sep=""), height = 7, width = 7)

plotViz(
  obj= x.sp,
  method="umap",
  main="region",
  point.color=x.sp@metaData$tissue,
  point.size=.3,
  point.shape=19,
  text.add=F,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=F
)

plotViz(
  obj= x.sp,
  method="umap",
  main="sample_name",
  point.color=x.sp@metaData$sample_name,
  point.size=.3,
  point.shape=19,
  text.add=F,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=F
)

dev.off()
####

x.sp@metaData$region = x.sp@metaData$tissue

tab = data.frame(sample=x.sp@metaData$sample_name,cluster=x.sp@cluster,region = x.sp@metaData$tissue)
summ = plyr::count(tab[,1])
tab2 = plyr::count(tab)
tab2$total = summ$freq[match(tab2$sample,summ$x)]
tab2$frac = tab2$freq/tab2$total

pdf("celltype_frac_per_region.pdf",height=5,width=10)
for(i in unique(x.sp@metaData$tissue)) {
  p = ggplot(tab2[which(tab2$region==i),]) + 
    geom_col(aes(x=cluster,y=frac,fill=sample),position="dodge") +
    theme_bw() +ggtitle(i)
  print(p)
}
dev.off()




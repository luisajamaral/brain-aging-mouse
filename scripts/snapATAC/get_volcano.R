################## volcano plot


tissue =  "9mo+18mo"
setwd(paste("~/projects/brain_aging_mouse/analysis/snapATAC/",tissue,"/age_diff_edgeR.snap_celltype",sep = ""))
## read the meta info. 
meta = read.delim("../all.cluster.meta.txt")

files = list.files(pattern=".edger.txt")
max_cluster = length(files)
#dat = list()

for (f in files) {
 # dat[[f]] = read.csv(f)
 # cl = strsplit(f, "_")[[1]][1]
  cl = gsub(".edger.txt", "", f)
  res = read.csv(f)
  cols = rep("black", nrow(res))
  cols[which(res$fdr<0.05)] = "red"
  png(paste("~/projects/brain_aging_mouse/analysis/snapATAC/",tissue,"/age_diff_edgeR.snap_celltype/volcano/",cl,"_Volcanoplot.png",sep = ""))
  print(plot(res$logFC, -log10(res$fdr), col=cols, panel.first=grid(),
             main=paste("Volcano plot", cl), xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
             pch=20, cex=0.6))
  print(abline(v=0))
  print(abline(v=c(-1,1), col="brown"))
  alpha = 0.05
  print(abline(h=-log10(alpha), col="brown"))
  dev.off()
}

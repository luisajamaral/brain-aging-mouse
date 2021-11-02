
library(SnapATAC)
library(GenomicRanges)
library(tictoc)

setwd("~/projects/brain_aging_mouse/analysis/snapATAC/all66/")
tic("load things")
load("Ast.cluster.RData")
#setwd("~/projects/brain_aging_mouse/scripts/snapATAC/")        
peak.df = read.table("all21.Ast_peaks.narrowPeak")
peak.df = GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
toc()
tic("add Pmat")
x.sp.sub = createPmat(x.sp.sub,peaks =peak.df,do.par = T,num.cores = 5 )
toc()

save(x.sp.sub, "Ast.cluster.Pmat.RData")
library(SnapATAC)

tissue = commandArgs(trailing = T)[1] #"9mo+18mo"
RData = commandArgs(trailing = T)[2] # "all.cluster.pmat.RData"
setwd(paste0("../../analysis/snapATAC/", tissue))

load(RData)

#system("mkdir -p age_diff_edgeR.snap/")

#x.sp =  addPmatToSnap(x.sp,do.par = T,num.cores = 6)
#save(x.sp,file = "all.cluster.pmat.RData")

#peaks_table = read.table(peaks_file)
# library(GenomicRanges)
# peaks_table = read.table(peaks_file)
# peaks.gr = GRanges(
#   peaks_table[,1],
#   IRanges(peaks_table[,2], peaks_table[,3]), name = peaks_table[,4], score = peaks_table[,5]
# );
# x.sp = createPmat(x.sp, peaks = peaks.gr, num.cores = 10)

pmat = x.sp@pmat
peak = x.sp@peak$name
cluster = x.sp@cluster
x.sp@metaData$sample = gsub("[+]", "_",x.sp@metaData$sample)
#sample = x.sp@sample
sample = x.sp@metaData$sample
stage = x.sp@metaData$stage
tissue = gsub("18mo_", "", sample)
tissue = gsub("9mo_", "", tissue)
tissue = gsub("_1$", "", tissue)
tissue = gsub("_2$", "", tissue)


cluster_sample = paste(cluster, sample, sep = ".")
#rownames(pmat) = cluster_sample
#cl = "6"
max_cluster = max(as.numeric(cluster))
library(edgeR)
library(doParallel)
registerDoParallel(cores = max_cluster)

#cluster = rep("",length(cluster))
#cluster[which(x.sp@cluster%in%c(3,20))] = "OLG"
#cluster_sample = paste(cluster, sample, sep = ".")

#foreach(cl = 1:max_cluster) %dopar% {
for (cl in 1:max_cluster) {
  print(cl)
  samples = unique(sample)
 # tissues = unique(tissue)
  tissues = c("3F_4E" ,      "5E_6E"  )
  for (t in tissues) {
    dat = list()
    tsamples = samples[grep(t, samples, fixed = T)]
    for (ss in tsamples) {
      idx1 = which(cluster_sample == paste(cl, ss, sep = "."))
      if (length(idx1) < 2) {
        next
        #dat[[ss]] = rep(0, length(peak))
      } else {
        dat[[ss]] = colSums(pmat[idx1, ])
      }
    }
    mat = do.call(cbind, dat)
    if(is.null(mat) || ncol(mat) != 4) {
      next
    }
    rownames(mat) = peak
    
    grps = substr(colnames(mat), 1, 2)
    
    y = DGEList(mat)
    y = y[which(rowSums(cpm(y) > 1) >= 2), ]
    y = calcNormFactors(y)
    design = model.matrix( ~ grps)
    y <- estimateCommonDisp(y)
    y <- estimateGLMTagwiseDisp(y, design)
    fit_tag = glmFit(y, design)
    
    #contrast.matrix = matrix(c(c(-1,1,0),c(0,-1,1),c(-1,0,1)),nrow=3)
    #lrt = glmLRT(fit_tag, contrast =contrast.matrix)
    #lrt = glmLRT(fit_tag, contrast =c(-1,0,1))
    lrt = glmLRT(fit_tag)#, coef = "grps9m")
    fdr = p.adjust(lrt$table$PValue, method = "BH")
    out = cbind(cpm(y), lrt$table, fdr)
    out = out[order(out$PValue), ]
    show(head(out))
    write.csv(out,
              paste0("age_diff_edgeR.snap/", cl, "_", t, ".edger.txt"))
  }
}

## write differential peaks into bed files.
files = list.files("age_diff_edgeR.snap",pattern = ".txt",full.names = T)
for (f in files){
  tmp = read.csv(f)
  sig = tmp[which(tmp$PValue < 0.01), ]
  #  sig = tmp[which(tmp$PValue< quantile(tmp$PValue,0.01)),]
  chr = sub("(chr.*):(.*)-(.*)", "\\1", sig$X)
  start = sub("(chr.*):(.*)-(.*)", "\\2", sig$X)
  end = sub("(chr.*):(.*)-(.*)", "\\3", sig$X)
  out = data.frame(chr, start, end, sig$logFC, as.integer(-log10(sig$PValue)))
  #  write.table(out, paste0("../age_diff_edgeR.snap/",cl,".both.bed"),row.names=F,col.names=F,quote=F,sep="\t")
  outF =gsub(".txt", "", f)
  write.table(
    subset(out, sig.logFC > 0),
    paste0(outF,".up.bed"),
    row.names = F,
    col.names = F,
    quote = F,
    sep = "\t"
  )
  write.table(
    subset(out, sig.logFC < 0),
    paste0(outF, ".down.bed"),
    row.names = F,
    col.names = F,
    quote = F,
    sep = "\t"
  )
  
}

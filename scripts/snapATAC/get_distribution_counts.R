setwd("~/projects/brain_aging_mouse/analysis/snapATAC/all66/")
RData = "all.cluster.meta.final.pmat.RData"
load(RData)
peak = x.sp@peak$name
x.sp@metaData$final_clusters[which(x.sp@metaData$final_clusters %in% c("D1MSN","D2MSN"))] = "MSN"
cluster = x.sp@metaData$final_clusters
cluster = gsub("/","-",cluster)
sample = x.sp@metaData$sample_name
stage = x.sp@metaData$stage
tissue = x.sp@metaData$region
tissue = gsub(" ", ".", tissue)
cluster_sample = paste(cluster, sample, sep = ".")
age_cluster_tissue_rep = paste(cluster,stage,tissue, x.sp@metaData$replicate, sep = "_")
age_tissue_rep = paste(stage,tissue, x.sp@metaData$replicate, sep = "_")
cluster_tissue = paste(cluster,tissue,sep = "_")
pmat = x.sp@pmat 
peaks = x.sp@peak$name
pmat = pmat/rowSums(pmat)
pmat = pmat[which(cluster == 'OLG'), ]
stage = x.sp@metaData$stage
stage = stage[which(cluster == "OLG")]
d = pmat[,which(peaks == "chr13:21809701-21810979")]
d = data.frame(stage, d)
d$stage = factor(d$stage, levels = c("18mo","9mo","8wk"))
ggplot(d, aes(x=d, color=stage)) +
  geom_histogram(fill="white", alpha=0.1, position="identity") #+
  ylim(c(0,12000))

ggplot(d, aes(x=d, color=stage)) +
  geom_histogram(fill="white") +
  ylim(c(0,2000))

length(which(d$d!=0))/length(d$d)
length(which(d$d[which(d$stage=="18mo")]!=0))/length(d$d[which(d$stage=="18mo")])
length(which(d$d[which(d$stage=="9mo")]!=0))/length(d$d[which(d$stage=="9mo")])
length(which(d$d[which(d$stage=="8wk")]!=0))/length(d$d[which(d$stage=="8wk")])





hist(pmat[,which(peaks == "chr13:21809701-21810979")])
hist(pmat[which(stage == "8wk"),which(peaks == "chr13:21809701-21810979")], main = "8wk")
hist(pmat[which(stage == "9mo"),which(peaks == "chr13:21809701-21810979")], main = "9mo")
hist(pmat[which(stage == "18mo"),which(peaks == "chr13:21809701-21810979")], main = "18mo")
